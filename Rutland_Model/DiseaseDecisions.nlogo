;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Disease Decisions, version 3.1.
;; This is a specially truncated version for reproducing
;; the LSHTM compartmental model with POLYMOD contact matrices only.
;; All other methods for simulating contacts (Activity-Locations
;; using Time-Diaries survey) have been removed from menus for now.
;; (C) Christopher J Watts, 2021.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

extensions [array matrix table rnd csv profiler]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; General structure (2020/05/30):

;; Globals
;; Breeds
;; Agent attributes
;; Setup, RNG
;; Debugging aids: Trace-Turtles and Error messages
;; Setup Disease state and transitions, and how to access them.
;; Setup households and their inhabitant people.
;; City setup , including locations, and people.
;; City of Activity-Locations to determine contacts, instead of input contacts matrices.
;; Create new locations and people, and give them their initial attributes.
;; Code related to Contacts Matrices.
;; Setup initial diseased population.
;; Useful reporters for converting times
;; Go
;; DES engine: Management of Time-Bound events list (b-events)
;; Main scheduling procedure
;; Decisions about locations
;; Process infections and disease-state transitions.
;; Recording changes visually
;; Plots and statistics.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Globals

globals [
  prev-seed-setup
  prev-seed-go
  sorted-patches
  cur-person ; Used to inform code run by transitions and disease-states, who the current person is (and hence access their location, age, sex etc.)
  deaths ; Cumulative deaths per day from sim-time = 0

  current-view

  sim-time ; Starts at 0. In minutes.
  seed-infected ; Most recent person to be infected by extrinsic event, i.e. either during setup-disease-people, or by button press.
;  first-seed-day ; Actual first day of seed infections (as opposed to Seed-Start-Day, which might not be in use).

  b-events ; Time-Bound events.
  sim-stopping? ; Used for control.
  sim-stopping-reason ; E.g. "Infections Impossible", "Error..."
  start-date-as-num ; start-daye parsed

  ; LSHTM disease states
  susceptible
  exposed
  compartments-sorted
  new-case-transition

  i-preclinical
  i-clinical
  i-subclinical
  hospitalize
  h-icu
  h-non-icu
  recovered
  dead

  parameters-prop ; Various proportion parameters contained in a table.
  symptomatic-rates ; Sampled from Markov Chains giving Joint Posterior distribution for Symptomatic rates. Used with the LSHTM model.

  ; Now using turtle-sets for these. Breeds merged to locations.
  workplaces
  schools
  households
  hospitals
  loc-type-table

  school-holiday-dates
  school-holiday?

  age-groups ; Nobody ages in this model, so might as well store them at start.

  ; Contacts-Matrices
;  contacts-matrix-all
  contacts-matrix-home
  contacts-matrix-school
  contacts-matrix-work
  contacts-matrix-other
  contacts-matrix-weekday-school
  contacts-matrix-weekday-no-school
  contacts-matrix-weekend
  susceptibility

  ; Used for weighted sampling in schedule-next-random-event.
  total-infection-rate
;  total-dst-rate
  loc-type-infection-rate ; Arrays. 1 item for each location type.
;  loc-type-dst-rate
  loc-type-num-infections-here ; Array. 1 item for each location type.

  ; See calc-stats for definitions.
  max-new-exposed ; i.e. peak infection rate
  max-new-cases ; i.e. peak flow Ip -> Ic
  max-new-hospitalize
  max-new-deaths
  num-new-infected
  num-new-cases
  num-new-icu-beds
  num-new-deaths
  total-cases
  total-deaths
  num-peak-week-cases
  num-peak-week-deaths
  num-peak-icu-beds
  num-peak-non-icu-beds
  num-weeks-to-peak-week-cases
  total-infected

  superspreaders
  max-infected-by-one
  post-disease-people

  num-sus-net-components
  largest-sus-net-component

  num-days-in-lockdown
  locked-down?

  ; Useful debigging aids.
  trace-turtles

  time-taken-setup
  time-taken-go
]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Breeds

breed [locations location]
breed [people person]


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

directed-link-breed [habitants habitant] ; from people to households
directed-link-breed [staffings staffing] ; from people to workplaces
directed-link-breed [studyings studying] ; from people to schools
directed-link-breed [visitings visiting] ; from people to schools
undirected-link-breed [friendships friendship] ; with people

breed [disease-states disease-state]
directed-link-breed [transitions transition] ; from state to state

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Agent attributes

people-own [
  p-location
  p-home ; A household
  p-workplace ; A location
  p-state ; Disease state
  p-time-of-last-change ; When did I last change disease-state or location?
  p-infection-time
  p-ds-history ; [[time new-disease-state]]
  p-age ; Fixed at start. (Though we could schedule birthdays...)
  p-sex ; Fixed at start.
  p-health ; Controls how quickly I get better.
  p-contacts ; Controls how often I interact.
  p-susceptibility ; Controls how easily I catch disease. E.g. do I wash my hands before touching my face?
  p-infectiousness ; Controls how easily I transmit disease. E.g. do I cough into a handkerchief?
  p-off-school?
  p-off-work?
  p-num-infected ; How many I have infected so far.
  p-component ; Network component. E.g. potential transmission network between susceptibles.
  p-options ; List of options for choosing next activity.
  p-data ; Copy of data row used to create person.
]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

locations-own [
  l-type ; ID for Type of location (school, workplace, household etc.)
  l-type-name ; Recording this saves time when inspecting locations.
  l-source ; If location based on data file, give reference so it can be traced back. E.g. Row ID
  l-contents ; Array of disease-states' people at this location.
;  l-health-boost ; Multiplier to increase rate of health recovery.
  l-infection-boost ; Multiplier to increase rate of infections.
  l-sc-rate ; Total weighted contacts rates for susceptibles here.
  l-ic-rate ; Total weighted contacts rates for infectious here.
  l-nc-rate ; Total contacts rates for all non-dead here.
  l-infection-rate ; Total infection rates for all  non-dead here.
  l-dst-rate ; Total disease-state transitions' rates for all non-dead here.
  l-num-infections-here ; Number of infections occurring here since start.

  l-min-age ; Schools only.
  l-max-age ; Schools only.

;  l-people-here ; If we weren't moving people to locations, we wouldn't have people-here .

]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

disease-states-own [
  ds-name
  ds-color
  ds-cur-num
  ds-max-num
  ds-max-at-time
  ds-rel-inf ; Relative infectiousness
]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

transitions-own [
  ts-perc-weight
;  ts-condition
  ts-ie-time
  ts-time-param
  ts-rate ; = 1.0 / ts-time-param
  ts-event-commands
  ts-num-events
  ts-num-new-events
]



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Setup and RNG

to setup
  clear-all
  setup-folders
  if Using-Profiler? [profiler:reset profiler:start]
  set time-taken-setup timer
  set time-taken-go false

  if print-scheduling? or print-processing-b-events? [print "" print "New Setup:" print ""]
  let new-world-width max (list 2 (ceiling sqrt (population-size * 100 / (2.0 * perc-patches-with-household))))
  if world-width != new-world-width [
    resize-world 0 (new-world-width - 1) 0 (new-world-width - 1)
    set-patch-size 13 * 30 / new-world-width
  ]
  set sorted-patches sort patches
  foreach sorted-patches [p -> ask p [set pcolor white]]

  set current-view "City"
  setup-disease ; Setup the disease-states and transitions.

  setup-rng "Seed-Setup" ; To fix the contents of the city.

  if generate-random-r0? [set r0 random-normal 2.675739 0.5719293] ; Numbers from LSHTM's UK.R script.

  setup-population ; Setup households and their inhabitant people.
  setup-city ; Setup city structure, including other locations.

  ; Initialise the discrete-event simulation engine.
  reset-ticks ; Will also initialise view, plots, monitors. Any plotted or monitored information must be ready for this.
  set sim-time 0 ; We enter at midnight.
  set sim-stopping? false
  set sim-stopping-reason ""

  b-events-setup

;  setup-rng "Seed-Go". Would mean, given fixed seed-setup, different people could be initial infected person.

  setup-disease-people ; Initialise the population's disease-states.
  if Calculate-Susceptibles-Net? [calc-susceptibles-network] ; Calculates uninteresting(?) statistic. Delete?

  setup-school-holidays ; Defines school holidays, and schedules next closure/reopening.

  set num-days-in-lockdown 0
  set locked-down? false

  if social-interactions = "Activity-Locations" [
    setup-personal-schedules ; Give each person their personal list of options for choosing activities.
    schedule-initial-decisions
    setup-locations-rates
    ; TO DO: Schedule interventions
  ]

  if social-interactions = "Contacts-Matrices" [
    setup-contacts-matrices
    setup-Intervention "Base"
    set susceptibility susceptibility-given-R0 r0 ; NB: assuming susceptibility is not a function of age group.
    schedule (sim-time + 0 * 60) nobody [-> schedule-contacts-infections] nobody "Infections via contacts matrices."

    ; Intervention start and end. Scheduled to occur 1-2 minutes after next day.
    schedule (sim-time + 1 + 1440 * (max (list 0 intervention-day))) nobody [-> setup-Intervention Intervention] nobody (word "Intervention: " Intervention)
    schedule (sim-time + 2 + 1440 * (max (list 0 (intervention-day + intervention-duration)))) nobody [-> setup-Intervention "Base"] nobody "Return to Base scenario."

  ]
;  schedule sim-time nobody [-> schedule-next-day] nobody "Next Day"

  ; Initialize stats
  set max-new-exposed 0
  set max-new-cases 0
  set max-new-hospitalize 0
  set max-new-deaths 0

  run (word "setup-seed-infectors-" seed-infectors) ; Schedule the seed infection events. Should be last use of RN stream, since could vary in number of seed infectors..

  setup-trace-turtles

  setup-initial-plots

  schedule-next-day ; which includes calculating stats and doing plots at time 0.
  set sim-time -1 ; Having done plots at time 0, go back 1 minute.
  setup-rng "Seed-Go"
  ; So if seed-setup is fixed,
  ; go processes same city structure, locations, people, social networks,
  ; and initial infected person.
;  set trace-turtles sort runresult trace-turtles-reporter ; Which seed do you want this based on?
  print-trace-turtles "Initialized."

  set time-taken-setup timer - time-taken-setup

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-rng [given-variable-name]
  ifelse 0 = runresult given-variable-name [
    run (word "set prev-" given-variable-name " " new-seed)
  ]
  [
    run (word "set prev-" given-variable-name " " given-variable-name)
  ]
  random-seed runresult (word "prev-" given-variable-name)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-folders
  if input-data-folder = "" [stop]
  if not file-exists? input-data-folder [
    user-message "Please update input-data-folder"
    setup-input-folder
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-input-folder
  let tmp-dir user-directory
  if tmp-dir = false [stop]
  set input-data-folder replace-string tmp-dir "\\" "\\\\"
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report replace-string [str-to-be-searched str-to-be-replaced str-replacement]
  let cur-pos 0
  let str-to-return ""
  foreach n-values (length str-to-be-searched) [i -> item i str-to-be-searched] [x ->
    set str-to-return (word str-to-return (ifelse-value (x = str-to-be-replaced) [str-replacement] [x]))
  ]
  report str-to-return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to change-view
  ; NB: Use of lists and sort. We don't want to use any random numbers here.
  ifelse current-view = view-of-city [
    set current-view view-of-disease-states
  ]
  [
    if current-view = view-of-disease-states [
      set current-view view-of-city
    ]
  ]
  if any? turtles [
    foreach sort turtles [p ->
      ask p [
        if breed != disease-states [
          set hidden? (current-view != view-of-city)
          foreach sort (my-links with [end1 = myself]) [s ->
            ask s [
              set hidden? Hide-Peoples-Links? or (current-view != view-of-city)
            ]
          ]
        ]
      ]
    ]
  ]
  foreach compartments-sorted [c ->
    ask c [
      set hidden? (Current-View != view-of-disease-states)
      foreach sort (my-out-transitions) [t ->
        ask t [
          set hidden? (Current-View != view-of-disease-states)
        ]
      ]
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report view-of-city
  report "City"
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report view-of-disease-states
  report "Disease States"
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-contacts-matrices-from-files
  set contacts-matrix-home contacts-from-file "home"
  set contacts-matrix-work contacts-from-file "work"
  set contacts-matrix-school contacts-from-file "school"
  set contacts-matrix-other contacts-from-file "other"
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report contacts-from-file [mat-name]
  if not file-exists? (word input-data-folder "cm_" mat-name ".csv") [
    user-message (word "Could not find file: \n" input-data-folder "cm_" mat-name ".csv")
    report false
  ]
  let tmp-mat map [r -> but-first r] but-first csv:from-file (word input-data-folder "cm_" mat-name ".csv")
  if 16 != length tmp-mat [user-message (word "ERROR in contacts-matrix-" mat-name "!\n\nMatrix is " (length tmp-mat) " rows long.")]
  if 16 != length filter [x -> 16 = length x] tmp-mat [user-message (word "ERROR in contacts-matrix-" mat-name "!\n\nSome rows do not have 16 columns.")]
  report matrix:from-row-list tmp-mat
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Debugging aids: Trace-Turtles and Error messages

to setup-trace-turtles
  ; Trace-turtles will show their progress through the simulation.
  ; Useful debugging aid (speaks from experience).
  set trace-turtles runresult trace-turtles-reporter ; Which seed do you want this based on?
  ifelse not is-list? trace-turtles and not is-turtle-set? trace-turtles [
    ifelse is-turtle? trace-turtles [
      set trace-turtles (list trace-turtles) ; Single turtle?
    ]
    [
      set trace-turtles (list )
    ]
  ]
  [
    set trace-turtles sort trace-turtles
  ]
  if trace-turtles? [
    print ""
    print "Tracing Turtles through Simulation Run:"
    print ""
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to print-trace-turtles [given-text]
  if trace-turtles? [
    foreach trace-turtles [t ->
      ask t [
        if is-person? self [
          show (word "Ticks=" ticks " : Time=" sim-time " : Day=" (sim-days) " " sim-time-as-hh-mm " at "([item l-type location-type-names] of p-location) " " p-location " : Feeling=" ([ds-name] of p-state) " : " given-text)
        ]
        if is-location? self [
          show (word "(" (item l-type location-type-names) ") Ticks=" ticks " : Time=" sim-time " : Day=" (sim-days) " " sim-time-as-hh-mm " : Count=" (count people-here) " : sc=" l-sc-rate ", ic=" l-ic-rate ", nc=" l-nc-rate ", IR=" l-infection-rate ", DSTR=" l-dst-rate " : " given-text)
        ]
      ]
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to print-error [given-message]
  set sim-stopping-reason (word given-message " at Sim-Time=" sim-time ", Seed-Setup=" prev-seed-setup ", Seed-Go=" prev-seed-go)
  print sim-stopping-reason
  set sim-stopping? true
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to show-error [given-message]
  set sim-stopping-reason (word given-message " from " self " at Sim-Time=" sim-time ", Seed-Setup=" prev-seed-setup ", Seed-Go=" prev-seed-go)
  show sim-stopping-reason
  set sim-stopping? true
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to user-message-error [given-message]
  set sim-stopping-reason (word given-message " at Sim-Time=" sim-time ", Seed-Setup=" prev-seed-setup ", Seed-Go=" prev-seed-go)
  set sim-stopping? true
  user-message sim-stopping-reason
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Setup disease-states, their transitions.

; Other choices of states could include:
; SIR
; SEIR
; MSEIRS
; Susceptible Exposed Infectious Symptoms Incapacited Hospitalised Ventilated Dead Recovered

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-disease
  run (word "setup-disease-" disease-model)
  setup-parameters-prop
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-disease-LSHTM-Covid-19
  ; LSHTM use different labels in Lancet paper to those in their code corona.cpp :)
  ; Here I'm giving states the labels from Lancet paper (where the state is mentioned in the paper, esp. fig. 1),
  ; However, I'm merging the transmission, health burden, and deaths models,
  ; because, while cloning flows and feeding them into parallel processes may be tolerable in an SD model,
  ; cloning agents in an ABM is weird and wasteful,
  ; and I like my agents to leave infectious, before they enter hospital, and leave hospital before they enter death. ;-)

  let mid-x mean (list max-pxcor min-pxcor)
  let x-step world-width / 8
  let y-step world-height / 8

  ; new-compartment name relative-infectiousness xcor ycor color
  ; new-compartment name relative-infectiousness xcor ycor
  let ds-S new-disease-state "Susceptible" 0.0 mid-x (max-pycor - 1 * y-step) green
  let ds-E new-disease-state "Exposed" 0.0 mid-x (max-pycor - 2 * y-step) yellow
  let ds-Ip new-disease-state "I-Preclinical" 1.00 (mid-x + 1 * x-step) (max-pycor - 3 * y-step) orange
  let ds-Ic new-disease-state "I-Clinical" 1.00 (mid-x + 1 * x-step) (max-pycor - 4 * y-step) red
  let ds-Is new-disease-state "I-Subclinical" 0.50 (mid-x - 1 * x-step) (max-pycor - 4 * y-step) pink
  let ds-H new-disease-state "Hospitalize?" 0.0 (mid-x + 1 * x-step) (max-pycor - 5 * y-step) (violet + 1)
  let ds-ICU new-disease-state "ICU" 0.0 (mid-x + 2 * x-step) (max-pycor - 6 * y-step) (sky)
  let ds-Non-ICU new-disease-state "Non-ICU" 0.0 (mid-x + 1 * x-step) (max-pycor - 6 * y-step) (blue - 2)
  let ds-R new-disease-state "Recovered" 0.0 (mid-x - 1 * x-step) (max-pycor - 7 * y-step) brown
  let ds-D new-disease-state "Dead" 0.0 (mid-x + 2 * x-step) (max-pycor - 7 * y-step) (gray - 4)

  ; Setup globals.
  set susceptible ds-S
  set exposed ds-E
  set compartments-sorted (list ds-S ds-E ds-Ip ds-Ic ds-Is ds-H ds-ICU ds-Non-ICU ds-R ds-D)

  ; make-transition-to destination-compartment weight maturity-time-reporter mean-maturity-time on-entry-procedure
  ask ds-S [make-transition-to ds-E [-> 1.0] [-> ] "inf" [-> ] ]
  ask ds-E [make-transition-to ds-Ip [-> Perc-Show-Symptoms] [-> random-event-time-gamma 4.0 4.0] 4.0 [-> ] ]
  ask ds-E [make-transition-to ds-Is [-> 1.0 - Perc-Show-Symptoms] [-> random-event-time-gamma 4.0 4.0] 4.0 [-> ] ]
  ask ds-Ip [make-transition-to ds-Ic [-> 1.0] [-> random-event-time-gamma 1.5 4.0] 1.5 [-> ] ]
  ask ds-Ic [make-transition-to ds-H [-> 1.0] [-> random-event-time-gamma 3.5 3.5] 3.5 [-> ] ]
  ask ds-Is [make-transition-to ds-R [-> 1.0] [-> random-event-time-gamma 5.0 4.0] 5.0 [-> ] ]

  ask ds-H [make-transition-to ds-R [-> 1.0 - Perc-Need-Hospitalizing] [-> 0] 0 [-> add-to-post-disease-people] ]
  ask ds-H [make-transition-to ds-ICU [-> Perc-Need-Hospitalizing * Perc-Need-ICU] [-> random-event-time-gamma 3.5 3.5] 3.5 [-> goto-hospital] ]
  ask ds-H [make-transition-to ds-Non-ICU [-> Perc-Need-Hospitalizing * (1.0 - Perc-Need-ICU)] [-> random-event-time-gamma 3.5 3.5] 3.5 [-> goto-hospital] ]
  ask ds-ICU [make-transition-to ds-R [-> 1.0 - Perc-Die-From-ICU] [-> random-event-time-gamma 10.0 10.0] 10.0 [-> add-to-post-disease-people] ]
  ask ds-ICU [make-transition-to ds-D [-> Perc-Die-From-ICU] [-> random-event-time-gamma 10.0 10.0] 10.0 [-> add-to-post-disease-people] ]
  ask ds-Non-ICU [make-transition-to ds-R [-> 1.0 - Perc-Die-From-Non-ICU] [-> random-event-time-gamma 8.0 8.0] 8.0 [-> add-to-post-disease-people] ]
  ask ds-Non-ICU [make-transition-to ds-D [-> Perc-Die-From-Non-ICU] [-> random-event-time-gamma 8.0 8.0] 8.0 [-> add-to-post-disease-people] ]

  ;ask ds-R [make-transition-to ds-S 1.0 [-> random-event-time-gamma 365 365] 365 [-> remove-from-post-disease-people]] ; If immunity wears off.

  ask ds-Ip [ask out-transition-to ds-Ic [set new-case-transition self]] ; Define a new (recorded?) case.

  ; Setup globals for convenience(?).
  ; Comments contain labels from corona.cpp.
  set i-preclinical ds-Ip ; Ip: Presymptomatic
  set i-clinical ds-Ic ; Is: Symptomatic
  set i-subclinical ds-Is ; Ia: Asymptomatic
  ; NB: The Lancet paper does not distinguish between the different types of "Removed" from population.
  ; (Whether Covid-19 hospitalized people were really "removed" from the population is another question... ;) )
  set hospitalize ds-H ; Assuming hospital staff protected from infectious patient.
  ; The Lancet paper just has a Removed state.
  ; The Supplementary Information describes additional holding and beds states HG, BG, HI, BI for (General ward and ICU)
  ; to represent health burden processes.
  ; The numbers added to these compartments are cloned from Ip-Ic, not flows from a stock.
  ; corona.cpp has one additional state H (for hospitalized). This is a stock.
  ; Deaths are also based on duplication from Ip-Ic, not flow between stocks.
  ; For an ABM, consistency between states and time order is important.
  set h-non-icu ds-Non-ICU
  set h-icu ds-ICU
  set dead ds-D
  set recovered ds-R

  relabel-all-transitions

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-disease-Susceptible-Infectious
  ; The simplest model: Susceptible-Infectious
  let mid-x mean (list max-pxcor min-pxcor)
  let x-step world-width / 8
  let y-step world-height / 4

  ; new-compartment name relative-infectiousness xcor ycor color
  let ds-S new-disease-state "Susceptible" 0.0 mid-x (max-pycor - 1 * y-step) green
  let ds-I new-disease-state "Infectious" 1.0 mid-x (max-pycor - 2 * y-step) red

  set susceptible ds-S
  set exposed ds-I
  set compartments-sorted (list ds-S ds-I)

  ; make-transition-to destination-compartment weight maturity-time-reporter
  ask ds-S [make-transition-to ds-I [-> 1.0] [-> ] "inf" [-> ]]

  ask ds-S [ask out-transition-to ds-I [set new-case-transition self]]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-disease-SIR
  ; The classic model introduced by Kermack & McKendrick (1927): Susceptible-Infectious-Removed
  let mid-x mean (list max-pxcor min-pxcor)
  let x-step world-width / 8
  let y-step world-height / 4

  ; new-compartment name relative-infectiousness xcor ycor color
  let ds-S new-disease-state "Susceptible" 0.0 mid-x (max-pycor - 1 * y-step) green
  let ds-I new-disease-state "Infectious" 1.0 mid-x (max-pycor - 2 * y-step) red
  let ds-R new-disease-state "Removed" 0.0 (mid-x - 0 * x-step) (max-pycor - 3 * y-step) brown

  set susceptible ds-S
  set exposed ds-I
  set compartments-sorted (list ds-S ds-I ds-R)

  ; make-transition-to destination-compartment weight maturity-time-reporter
  ask ds-S [make-transition-to ds-I [-> 1.0] [-> ] "inf" [-> ]]
  ask ds-I [make-transition-to ds-R [-> 1.0] [-> random-event-time-gamma 8.0 8.0] 8.0 [-> add-to-post-disease-people]]

  ask ds-S [ask out-transition-to ds-I [set new-case-transition self]]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report rescalar

  ; Rescale coordinates
  report world-width / 30.0
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report scenario-params-LSHTM
  ; See Table 1 in Lancet paper.
  report table:from-list (list
    ; Scenario-Name Home-Contacts Work-Contacts School-Contacts Other-Contacts Relative-Infectiousness-Of-I-Clinical
    ; Work contacts may be of form [[a p] q] where reports p if age < a, q otherwise.
    ; School contacts may be of [p q] where reports p if schools open, q if schools closed.
    (list "Base" [100 100 100 100 100])
    (list "School Closures" [100 100 0 100 100])
    (list "Social Distancing" [100 50 100 50 100])
    (list "Elderly Shielding" [100 [100 25] 100 [100 25] 100])
    (list "Self-Isolation" [100 100 100 100 65])
    (list "Combination" [100 [50 25] 0 [50 25] 65])
    (list "Intensive-Interventions" [100 [65 25] [100 0] [59 16] 65])
    (list "Lockdown" [100 10 [10 0] 10 65])
  )
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-intervention [given-Intervention]
  let cur-Intervention-params map [p -> ifelse-value (is-number? p) [p / 100] [map [p2 -> p2 / 100] p]] table:get scenario-params-LSHTM given-Intervention
  ; Home
  set contacts-matrix-weekday-school matrix:times contacts-matrix-home (item 0 cur-Intervention-params)
  set contacts-matrix-weekday-no-school matrix:times contacts-matrix-weekday-school 1.0 ; Shouldn't we have higher home contacts when there's no school?
  set contacts-matrix-weekend matrix:times contacts-matrix-weekday-no-school 1.0

  let split-item int (70 / 5) ; Split off the elderly (70+)
  let num-items 16 ; 16 age groups in contacts matrix

  ; Work
  let cur-entry (item 1 cur-Intervention-params)
  ifelse is-number? cur-entry [
    set contacts-matrix-weekday-school matrix:plus contacts-matrix-weekday-school (matrix:times contacts-matrix-work cur-entry)
    set contacts-matrix-weekday-no-school matrix:plus contacts-matrix-weekday-no-school (matrix:times contacts-matrix-work cur-entry)
    set contacts-matrix-weekend matrix:plus contacts-matrix-weekend (matrix:times contacts-matrix-work (Weekend-Work * cur-entry / 100))
  ]
  [
    ; Weekday, School
    set contacts-matrix-weekday-school matrix:plus contacts-matrix-weekday-school
    (matrix:times-element-wise contacts-matrix-work matrix:from-row-list
      sentence
      n-values split-item [->
        sentence
        n-values split-item [-> first cur-entry]
        n-values (num-items - split-item) [-> last cur-entry]
      ]
      n-values (num-items - split-item) [->
        n-values num-items [-> last cur-entry]
      ]
    )

    ; Weekday, No School
    set contacts-matrix-weekday-no-school matrix:plus contacts-matrix-weekday-no-school
    (matrix:times-element-wise contacts-matrix-work matrix:from-row-list
      sentence
      n-values split-item [->
        sentence
        n-values split-item [-> first cur-entry]
        n-values (num-items - split-item) [-> last cur-entry]
      ]
      n-values (num-items - split-item) [->
        n-values num-items [-> last cur-entry]
      ]
    )

    ; Weekend
    set contacts-matrix-weekend matrix:plus contacts-matrix-weekend
    (matrix:times-element-wise contacts-matrix-work matrix:from-row-list
      sentence
      n-values split-item [->
        sentence
        n-values split-item [-> (Weekend-Work * first cur-entry / 100)]
        n-values (num-items - split-item) [-> (Weekend-Work * last cur-entry / 100)]
      ]
      n-values (num-items - split-item) [->
        n-values num-items [-> (Weekend-Work * last cur-entry / 100)]
      ]
    )
  ]

  ; School
  set cur-entry (item 2 cur-Intervention-params)
  ifelse is-number? cur-entry [
    set contacts-matrix-weekday-school matrix:plus contacts-matrix-weekday-school (matrix:times contacts-matrix-school cur-entry)
  ]
  [
    set contacts-matrix-weekday-school matrix:plus contacts-matrix-weekday-school (matrix:times contacts-matrix-school first cur-entry)
  ]

  ; Other
  set cur-entry (item 3 cur-Intervention-params)
  ifelse is-number? cur-entry [
    set contacts-matrix-weekday-school matrix:plus contacts-matrix-weekday-school (matrix:times contacts-matrix-other cur-entry)
    set contacts-matrix-weekday-no-school matrix:plus contacts-matrix-weekday-no-school (matrix:times contacts-matrix-other cur-entry)
    set contacts-matrix-weekend matrix:plus contacts-matrix-weekend (matrix:times contacts-matrix-other cur-entry) ; Transport should change if no school or weekend!
  ]
  [
    ; Weekday, School
    set contacts-matrix-weekday-school matrix:plus contacts-matrix-weekday-school
    (matrix:times-element-wise contacts-matrix-other matrix:from-row-list
      sentence
      n-values split-item [->
        sentence
        n-values split-item [-> first cur-entry]
        n-values (num-items - split-item) [-> last cur-entry]
      ]
      n-values (num-items - split-item) [->
        n-values num-items [-> last cur-entry]
      ]
    )

    ; Weekday, No School
    set contacts-matrix-weekday-no-school matrix:plus contacts-matrix-weekday-no-school
    (matrix:times-element-wise contacts-matrix-other matrix:from-row-list
      sentence
      n-values split-item [->
        sentence
        n-values split-item [-> first cur-entry]
        n-values (num-items - split-item) [-> last cur-entry]
      ]
      n-values (num-items - split-item) [->
        n-values num-items [-> last cur-entry]
      ]
    )

    ; Weekend
    set contacts-matrix-weekend matrix:plus contacts-matrix-weekend
    (matrix:times-element-wise contacts-matrix-other matrix:from-row-list
      sentence
      n-values split-item [->
        sentence
        n-values split-item [-> first cur-entry]
        n-values (num-items - split-item) [-> last cur-entry]
      ]
      n-values (num-items - split-item) [->
        n-values num-items [-> last cur-entry]
      ]
    )


  ]

  ; Relative Infectiousness of I-Clinical people.
  ask i-clinical [set ds-rel-inf (item 4 cur-Intervention-params)]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-parameters-prop
  if Disease-Model = "LSHTM-Covid-19" [setup-parameters-Prop-LSHTM stop]
  setup-parameters-Prop-LSHTM
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-parameters-Prop-LSHTM
  ; Taken from file UK.R, line 79+
  let info-by-rows csv:from-string (word
    "Age,Prop_symptomatic,IFR,Prop_inf_hosp,Prop_inf_critical,Prop_critical_fatal,Prop_noncritical_fatal,Prop_symp_hospitalised,Prop_hospitalised_critical\n"
    "10,0.66,8.59E-05,0.002361009,6.44E-05,0.5,0,0,0.3\n"
    "20,0.66,0.000122561,0.003370421,9.19E-05,0.5,9.47E-04,0.007615301,0.3\n"
    "30,0.66,0.000382331,0.010514103,0.000286748,0.5,0.001005803,0.008086654,0.3\n"
    "40,0.66,0.000851765,0.023423527,0.000638823,0.5,0.001231579,0.009901895,0.3\n"
    "50,0.66,0.001489873,0.0394717,0.001117404,0.5,0.002305449,0.018535807,0.3\n"
    "60,0.66,0.006933589,0.098113786,0.005200192,0.5,0.006754596,0.054306954,0.3\n"
    "70,0.66,0.022120421,0.224965092,0.016590316,0.5,0.018720727,0.150514645,0.3\n"
    "80,0.66,0.059223786,0.362002579,0.04441784,0.5,0.041408882,0.332927412,0.3\n"
    "100,0.66,0.087585558,0.437927788,0.065689168,0.5,0.076818182,0.617618182,0.3"
  )
  set parameters-Prop table:from-list (n-values (length first info-by-rows) [i -> (list (item i first info-by-rows) (reformat-Prop-LSHTM map [y -> item i y] but-first info-by-rows))])

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report reformat-Prop-LSHTM [given-list]
  ; Line 91 in UK.R
  ; 70-74,3388.488  75-79,2442.147  80-84,1736.567  85-89,1077.555  90-94,490.577  95-99,130.083  100+,15.834
  ; x = c(P[1:7], weighted.mean(c(P[8], P[9]), c(3388.488 + 2442.147, 1736.567 + 1077.555 + 490.577 + 130.083 + 15.834)));
  let age-weights (list (3388.488 + 2442.147) (1736.567 + 1077.555 + 490.577 + 130.083 + 15.834))
  report (reduce sentence map [x -> (list x x)] (lput (sum (map [[w x] -> w * x / sum age-weights] age-weights (sublist given-list 7 9))) (sublist given-list 0 7) ))
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report Prop_symptomatic
  ; The LSHTM model samples a set from a joint posterior estimated using MCMC.
  ; How best to represent this in NetLogo?

  ; Use this if you have collated the actual covy values used for your LSHTM runs.
;  report array:from-list sublist (item LSHTM-Run csv:from-file "Runs_covy.csv") 2 10

  ; Use this to sample values from the Markov Chains, i.e. to sample from the estimate of the joint posterior generated by the LSHTM model.
  if symptomatic-rates = 0 [
    set symptomatic-rates array:from-list sublist (one-of but-first csv:from-file "2-linelist_symp_fit_fIa0.5.csv") 5 13 ; NB: The columns for f_00:f_70
  ]
  report symptomatic-rates

  ; These are the means from the (Markov Chains) in 2-linelist_symp_fit_fIa0.5.qs
  ; NB: a combination of means might be very unrepresentative of the joint posterior distribution.
;  report array:from-list [0.147371843	0.073028676	0.296122577	0.419103564	0.444586716	0.563571982	0.816944346	0.750559872]


  ; These are the max likelihood estimates from the Markov Chains in 2-linelist_symp_fit_fIa0.5.qs
  ; NB: Just because they had the highest likelihood, does not mean they are good representations of the PDF as a whole.
  ; There may be plenty of combinations with important likelihood far from this combination.
;  report array:from-list [0.132377732	0.06060549	0.267664196	0.37323943	0.406284315	0.521538669	0.759989905	0.664960788]

  ; Just returns 0.66 for every double age group
;  report array:from-list table:get parameters-Prop "Prop_symptomatic"
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report Perc-Show-Symptoms
  ; This version for use with the age-varying symptomatic rate in 2-linelist_symp_fit_fIa0.5.qs
  report [array:item Prop_symptomatic ifelse-value (p-age >= 70) [7] [int (p-age / 10)]] of cur-person

  ; This version for use with parameters-Prop-LSHTM
;  report [array:item Prop_symptomatic p-age-group] of cur-person
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report Prop_symp_hospitalised
  report array:from-list table:get parameters-Prop "Prop_symp_hospitalised"
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report Perc-Need-Hospitalizing
  report [array:item Prop_symp_hospitalised p-age-group] of cur-person
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report Prop_hospitalised_critical
  report array:from-list table:get parameters-Prop "Prop_hospitalised_critical"
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report Perc-Need-ICU
  report [array:item Prop_hospitalised_critical p-age-group] of cur-person
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report Prop-NonCritical-Dead
  report array:from-list table:get parameters-Prop "Prop_noncritical_fatal"
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report Prop-Critical-Dead
  ; NB: Lines 101, 119 in UK.R: LSHTM model uses Prop_noncritical_fatal for both ICU and Non-ICU patients.
  ; (Presumably this is because they didn't manage to get data that distinguished between ICU and Non-ICU deaths.
  ; Though it might just be that they simplified, used one probability only, but retained the now misleading label.)
  report array:from-list table:get parameters-Prop "Prop_noncritical_fatal"
;  report array:from-list table:get parameters-Prop "Prop_critical_fatal"
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report Perc-Die-From-ICU
  ; NB: The LSHTM model works with flows straight from Ip to Deaths.
  ; If we are going to branch to Dead after Hospitalization and ICU/Non-ICU,
  ; then we have to adjust the probabilities a bit.
  report [(array:item Prop-Critical-Dead p-age-group) / (Perc-Need-Hospitalizing)] of cur-person
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report Perc-Die-From-Non-ICU
  ; NB: The LSHTM model works with flows straight from Ip to Deaths.
  ; If we are going to branch to Dead after Hospitalization and ICU/Non-ICU,
  ; then we have to adjust the probabilities a bit.
  report [(array:item Prop-NonCritical-Dead p-age-group) / (Perc-Need-Hospitalizing)] of cur-person
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report new-disease-state [given-name given-relative-infectiousness given-xcor given-ycor given-color]
  let retobj nobody
  create-disease-states 1 [
    set hidden? (current-view = "City")
    set shape "circle"
    set size 2 * rescalar
    set ds-color given-color
    set color ds-color
    set ds-name given-name
    set label (word given-name "     ")
    set label-color black
    setxy (given-xcor) (given-ycor)
    set ds-cur-num 0
    set ds-max-num 0
    set ds-max-at-time 0
    set ds-rel-inf given-relative-infectiousness
    relabel-disease-state
    set retobj self
  ]
  report retobj

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to make-transition-to [target-state percentage-weight time-reporter time-param event-commands]
  create-transition-to target-state [
    set hidden? (current-view = "City")
    set thickness 0.2 * rescalar
    set shape "transition"
    set color grey
    set label-color grey
    set label (word "   " time-param)
    set ts-perc-weight percentage-weight
    set ts-ie-time time-reporter
    ifelse time-param = "inf" [
      set ts-time-param infinity
      set ts-rate 0
    ]
    [
      ifelse time-param = 0 [
        set ts-time-param time-param ; i.e. 0
        set ts-rate infinity ; Not entirely convinced this will work. No crashes so far.
      ]
      [
        set ts-time-param time-param ; In days.
        set ts-rate round (100.0 / time-param) ; NB: Events-Per-100-Days and integer.
      ]
    ]
    set ts-event-commands event-commands
    set ts-num-events 0
    set ts-num-new-events (list 0)
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report transition-between [from-state-name to-state-name]
  ; Assumes there is such a link.
  report [out-transition-to to-state-name] of from-state-name
  report [out-transition-to min-one-of (disease-states with [ds-name = to-state-name]) [who]] of min-one-of (disease-states with [ds-name = from-state-name]) [who]
;  if not any? disease-states with [ds-name = from-state-name] [report nobody] ; Shouldn't be needed, once bug-free. :)
;  if not any? disease-states with [ds-name = to-state-name] [report nobody]
;  report transition ([who] of min-one-of (disease-states with [ds-name = from-state-name]) [who]) ([who] of min-one-of (disease-states with [ds-name = to-state-name]) [who])
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report susceptible?
  report p-state = susceptible
;  report [ds-name = "Susceptible"] of p-state
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report ds-allows-infection?
  report ds-rel-inf > 0
;  report member? ds-name ["Infectious" "Symptoms" "Hospitalize?"] ; Will need to update if disease states change.
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report exposed?
  ; Have been exposed to the disease, and am incubating it, but not yet infectious.
  report p-state = exposed
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report i-preclinical?
  ; Now infectious, and will soon show symptoms, but not yet.
  report p-state = i-preclinical
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report i-clinical?
  ; Now infectious, and showing symptoms.
  report p-state = i-clinical
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report i-subclinical?
  ; Now infectious, but will remain asymptomatic.
  report p-state = i-subclinical
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report infectious?
  ; Now infectious.
  report [0 < ds-rel-inf] of p-state
  ;report [member? ds-name ["Infectious" "I-Preclinical" "I-Subclinical"]] of p-state
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report symptoms?
  ; Now infectious, and showing symptoms - feeling ill.
  report [member? ds-name ["Symptoms" "I-Clinical"]] of p-state
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report hospitalize?
  ; Infectious, and ill enough to need to be in hospital.
  report p-state = hospitalize
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report h-non-icu?
  ; Need hospitalizing, but in a General Ward.
  report p-state = h-non-icu
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report h-icu?
  ; Need hospitalizing in an Intensive Care Unit.
  ; (Acute respiratory distress syndrome (ARDS))
  report p-state = h-icu
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report infector?
;  report [member? ds-name ["Infectious" "Symptoms" "Hospitalize?"]] of p-state
  report [ds-rel-inf > 0] of p-state
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report dead?
  ; Dead.
  report p-state = dead
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report recovered?
  ; Recovered, and currently immune.
  report p-state = recovered
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to add-to-post-disease-people
  set post-disease-people (turtle-set post-disease-people self)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to remove-from-post-disease-people
  set post-disease-people post-disease-people with [self != myself]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report superspreader?
  report superspreader-threshold >= p-num-infected
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to add-to-superspreaders?
  ; 7 was chosen because it picked out the outliers
  ; when testing a city of 250 people.
  ; May need to scale it.
  if superspreader-threshold = p-num-infected [set superspreaders (turtle-set superspreaders self)]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report infectors
  report people with [infector?]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report num-susceptible
  report [ds-cur-num] of susceptible
;  report count people with [susceptible?]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report num-exposed
  report [ds-cur-num] of exposed
;  report count people with [exposed?]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report total-exposed
  report [ts-num-events] of transition-between Susceptible Exposed
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report num-infectious
  report sum map [ds -> [ifelse-value (ds-rel-inf > 0) [ds-cur-num] [0]] of ds] compartments-sorted
;  report count people with [infectious?]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report num-i-preclinical
  report [ds-cur-num] of i-preclinical
;  report count people with [i-preclinical?]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report num-i-clinical
  report [ds-cur-num] of i-clinical
;  report count people with [i-clinical?]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report num-i-subclinical
  report [ds-cur-num] of i-subclinical
;  report count people with [i-subclinical?]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report num-symptoms
  report num-i-clinical
;  report count people with [symptoms?]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report num-hospitalize
  report [ds-cur-num] of hospitalize + [ds-cur-num] of h-icu + [ds-cur-num] of h-non-icu
;  report count people with [hospitalize?]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report total-icu
  report ([ts-num-events] of transition-between Hospitalize H-ICU)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report total-hospitalize
  report ([ts-num-events] of transition-between Hospitalize H-ICU) +
   [ts-num-events] of transition-between Hospitalize H-Non-ICU
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report num-h-non-icu
  report [ds-cur-num] of h-non-icu
;  report count people with [h-non-icu?]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report num-h-icu
  report [ds-cur-num] of h-icu
;  report count people with [h-icu?]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report num-recovered
  report [ds-cur-num] of recovered
;  report count people with [recovered?]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report num-dead
  report [ds-cur-num] of dead
;  report count people with [dead?]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Setup households and their inhabitant people.

to setup-population
  set age-groups array:from-list n-values 16 [-> array:from-list n-values (count disease-states) [-> (list )]]
  if Population-Generator = "Nuclear-Household City" [setup-households-of-fixed-size stop]
  if Population-Generator = "Survey Data" [setup-households-from-survey-data-file stop]
  if Population-Generator = "Demographic Data : UK" [setup-households-with-demog-data demographic-data-UK stop]
  if Population-Generator = "Demographic Data" [setup-households-with-demog-data demographic-data stop]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-households-of-fixed-size
  let cur-household nobody
  let cur-hh-id 0
  set cur-person nobody
  set households (turtle-set )
;  foreach sort n-of round (Perc-Patches-With-Household * count patches / 100) patches [cur-patch ->
  let cur-patch nobody
  while [population-size > count people] [
    set cur-patch vacant-patch
    set cur-household new-household
    set cur-hh-id cur-hh-id + 1
    ask cur-household [
      set l-source cur-hh-id
      move-to cur-patch
    ]
    foreach sublist (reduce sentence n-values (ceiling (People-To-1-Household / 4)) [? -> shuffle (list
      (list (40 + round random-normal 0 2) 1)
      (list (35 + round random-normal 0 2) 0)
      (list (11 + random 7) (random 2))
      (list (10 - random 6) (random 2))
    )]) 0 People-To-1-Household [? ->
      set cur-person new-person (item 0 ?) (item 1 ?)
      ask cur-person [create-habitant-to cur-household [set hidden? true]]
    ]
    ask cur-household [setup-cur-household]
  ]

  ; Safer, but slower, to do this in new-household.
  ;set households locations with [l-type-name = "household"]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-households-from-survey-data-file
  ; Loads the household survey data file model-pop2.csv
  ; Creates desired number of households, sampling from the data in the file.
  ; Creates people for the households, using the data in the file.
  let cur-household nobody
  set cur-person nobody

  set households (turtle-set )

  ; Load household survey data.
  if not file-exists? "model-pop2.csv" [user-message "FATAL! Tried to create people from file model-pop2.csv, but cannot find it in current folder." stop]
  let cur-data but-first (csv:from-file "model-pop2.csv" ",")
  let num-hhs length remove-duplicates map [r -> item 5 r] cur-data
;  print (word "Number of distinct Household IDs in file = " num-hhs)
  ; Assuming model-pop2.csv is sorted by household id.
  ; if not: set cur-data sort-by [[x y] item 5 x < item 5 y ] cur-data

  ; Build list of row numbers for finding each household's people in cur-data.
  let hh-first-rows (list 0)
  let cur-item 0
  let cur-hh-id item 5 first cur-data
  foreach cur-data [r ->
    if cur-hh-id != item 5 r [
      set hh-first-rows fput cur-item hh-first-rows
      set cur-hh-id item 5 r
    ]
    set cur-item cur-item + 1
  ]
  set hh-first-rows reverse hh-first-rows
  ;print hh-first-rows ; Debugging.
  set cur-data array:from-list cur-data

  ; Create desired number of households, sampling from cur-data, create their people.
  let cur-hh-item 0
  let cur-ditem 0
  let cur-drow []
;  foreach sort n-of round (Perc-Patches-With-Household * count patches / 100) patches [cur-patch ->
  let cur-patch nobody
  while [population-size > count people] [
    set cur-patch vacant-patch
    set cur-household new-household
    ask cur-household [
      move-to cur-patch
    ]

    ; Sample household from cur-data
    if cur-hh-item = 0 [set hh-first-rows shuffle hh-first-rows] ; (Re)shuffle the pointers to households in cur-data
    set cur-ditem item cur-hh-item hh-first-rows
    set cur-drow array:item cur-data cur-ditem
    set cur-hh-id item 5 cur-drow
    ask cur-household [set l-source cur-hh-id]
    while [cur-hh-id = item 5 cur-drow] [ ; Do all people with that household ID
      ;print (word "cur-hh-item = " cur-hh-item " : cur-ditem = " cur-ditem " : cur-hh-id = " cur-hh-id) ; For debugging.
      ; new-person like cur-drow
      set cur-person new-person (item 1 cur-drow) (ifelse-value ("Female" = item 2 cur-drow) [1] [0])
      ask cur-person [
        set p-location cur-household
        create-habitant-to cur-household [set hidden? true]
        set p-data map [x -> x] cur-drow ; Store personal data somewhere. We'll process the rest of it later.
;        show (word "cur-ditem = " cur-ditem " : cur-hh-id = " cur-hh-id " : cur-hh-item = " cur-hh-item)
      ]
      set cur-ditem cur-ditem + 1 ; Look for next person in cur-data
      ifelse cur-ditem < array:length cur-data [ ; Not at end of cur-data yet?
        set cur-drow array:item cur-data cur-ditem
      ]
      [
        set cur-hh-id ""
      ]
    ]
    ask cur-household [setup-cur-household]
    set cur-hh-item (cur-hh-item + 1) mod length hh-first-rows
  ]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-households-with-demog-data [given-demog-data]
  ; Uses demographic data, broken down by age group and sex.
  ; Allocates people to arbitrary households.
  ; Use this if you don't care about households as locations for transmission.
  ; E.g. if attempting to replicate behaviour of a compartmental model.

  ; If desired population-size >= sum of numbers in demog-data
  ; then use the numbers to create copy of demog-data population.
  ; Make up any shortfall by sampling from demog-data, weighted by those numbers.
  ; Remember: Who is in which household and with whom is completely random.

  let cur-household nobody
  set cur-person nobody
  set households (turtle-set )

  let demog-data reverse given-demog-data ; NB: reverse, so that we do old people first.
;  let total-demog sum map [cur-person-type -> int first cur-person-type] demog-data
  let total-demog sum map [cur-person-type -> first cur-person-type] demog-data

  repeat int (population-size / total-demog) [
    foreach demog-data [cur-person-type ->
      repeat int first cur-person-type [
        if population-size > count people [ ; NB: Slight penalty to last age group & sex: If person created for lots of fractions. Is it better to be unbiased or to generate requested pop size?
          if cur-household = nobody [set cur-household new-household-for-demog]
          if People-To-1-Household <= [count in-habitant-neighbors] of cur-household [
            ask cur-household [setup-cur-household]
            set cur-household new-household-for-demog
          ]
          set cur-person new-person ((random 5) + item 1 cur-person-type) (item 2 cur-person-type)
          ask cur-person [create-habitant-to cur-household [set hidden? true]]
        ]
      ]

      ; One extra for fraction?
      if (first cur-person-type) mod 1 != 0 [
        if population-size > count people [
          if (first cur-person-type) mod 1 > random-float 1 [
            if cur-household = nobody [set cur-household new-household-for-demog]
            if People-To-1-Household <= [count in-habitant-neighbors] of cur-household [
              ask cur-household [setup-cur-household]
              set cur-household new-household-for-demog
            ]
            set cur-person new-person ((random 5) + item 1 cur-person-type) (item 2 cur-person-type)
            ask cur-person [create-habitant-to cur-household [set hidden? true]]
          ]
        ]
      ]

    ] ; Next item of demog-data
  ] ; Next pass through demog-data

  let cur-person-type []
  while [population-size > count people] [
    if cur-household = nobody [set cur-household new-household-for-demog]
    if People-To-1-Household <= [count in-habitant-neighbors] of cur-household [
      ask cur-household [setup-cur-household]
      set cur-household new-household-for-demog
    ]
    set cur-person-type (rnd:weighted-one-of-list demog-data [x -> first x])
    set cur-person new-person ((random 5) + item 1 cur-person-type) (item 2 cur-person-type)
    ask cur-person [create-habitant-to cur-household [set hidden? true]]
  ]

  if cur-household != nobody [
    ask cur-household [setup-cur-household]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report new-household-for-demog
  let cur-patch vacant-patch
  let cur-household new-household
  let cur-hh-id count households
  ask cur-household [
    set l-source cur-hh-id
    move-to cur-patch
  ]
  report cur-household
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-cur-household
;  let next-heading nobody
  (foreach
    (sort in-habitant-neighbors)
    (n-values count in-habitant-neighbors [? -> (? - 0.5) / (count in-habitant-neighbors)]) [[?1 ?2] ->
      ask ?1 [
        set p-location myself
        move-to myself
        set heading (360 * ?2)
        fd 0.25
        ;set heading 180
        ;fd 1
      ]
  ])
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report demographic-data
  ; Source: wpp2019_pop2020.rds as used by LSHTM
  if not file-exists? (word input-data-folder "Demog.csv") [user-message "FATAL: Cannot find file Demog.csv in Input-Data-Folder!" report false]
  let demog-by-age-group-and-sex but-first (csv:from-file (word input-data-folder "Demog.csv") ",")
  report (sentence
    (map [[dat ag] -> (list (1000 * item 4 dat) (ag * 5) 0)] demog-by-age-group-and-sex (n-values length demog-by-age-group-and-sex [? -> ?]))
    (map [[dat ag] -> (list (1000 * item 5 dat) (ag * 5) 1)] demog-by-age-group-and-sex (n-values length demog-by-age-group-and-sex [? -> ?]))
    )
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


to-report demographic-data-UK
  ; Constructed in Excel from filtered wpp2019_pop2020.csv

  let demog-by-age-group-and-sex (list
    ["age group: 0-4" 1915.127 2009.363]
    ["age group: 5-9" 2011.016 2108.55]
    ["age group: 10-14" 1933.97 2022.37]
    ["age group: 15-19" 1805.522 1880.611]
    ["age group: 20-24" 2001.966 2072.674]
    ["age group: 25-29" 2208.929 2275.138]
    ["age group: 30-34" 2345.774 2361.054]
    ["age group: 35-39" 2308.36 2279.836]
    ["age group: 40-44" 2159.877 2148.253]
    ["age group: 45-49" 2167.778 2128.343]
    ["age group: 50-54" 2353.119 2281.421]
    ["age group: 55-59" 2306.537 2232.388]
    ["age group: 60-64" 1985.177 1919.839]
    ["age group: 65-69" 1734.37 1647.391]
    ["age group: 70-74" 1763.853 1624.635]
    ["age group: 75-79" 1304.709 1137.438]
    ["age group: 80-84" 969.611 766.956]
    ["age group: 85-89" 638.892 438.663]
    ["age group: 90-94" 320.625 169.952]
    ["age group: 95-99" 95.559 34.524]
    ["age group: 100+" 12.818 3.016]
  )
  report (sentence
    (map [[dat ag] -> (list (item 1 dat) (ag * 5) 0)] demog-by-age-group-and-sex (n-values length demog-by-age-group-and-sex [? -> ?]))
    (map [[dat ag] -> (list (item 2 dat) (ag * 5) 1)] demog-by-age-group-and-sex (n-values length demog-by-age-group-and-sex [? -> ?]))
    )
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Setup city, its locations, and people.

to setup-city
  if social-interactions = "Activity-Locations" [setup-city-activity-locations]
  if social-interactions = "Contacts-Matrices" [setup-city-contact-matrices]

  set loc-type-num-infections-here array:from-list n-values (length location-types) [-> 0]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; City of Activity-Locations to determine contacts, instead of input matrices.

to setup-city-activity-locations
  setup-hospitals
  setup-schools
  setup-workplaces
  setup-friendships
  if any? patches with [1 < count locations-here] [
    user-message-error "FATAL! Some patches have multiple locations on them."
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-friendships
  ask people [
    let ego self
    ask other people [
      if who < [who] of myself [
        if (random-float 1.0) < exp (0 - (
          (friendship-distance * distance myself) +
          (friendship-not-same-sex * ifelse-value (p-sex = [p-sex] of myself) [0] [1]) +
          (Friendship-Age-Gap * abs (p-age - [p-age] of myself)) +
          (friendship-not-same-school * ifelse-value (any? out-studying-neighbors) [ifelse-value (same-school-as? ego) [0] [1]] [0]) +
          (friendship-not-same-workplace * ifelse-value (any? out-staffing-neighbors) [ifelse-value (same-workplace-as? ego) [0] [1]] [0]) +
          friendship-base-param
          ))
        [
          create-friendship-with myself [
            set hidden? true
            set color magenta
          ]
        ]
      ]

    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report same-school-as? [given-alter]
  ; run by person. Reports if self has a school in common with given-alter.
  if not any? out-studying-neighbors [report false]
  report any? out-studying-neighbors with [member? self [out-studying-neighbors] of given-alter]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report same-workplace-as? [given-alter]
  ; run by person. Reports if self has a workplace in common with given-alter.
  if not any? out-staffing-neighbors [report false]
  report any? out-staffing-neighbors with [member? self [out-staffing-neighbors] of given-alter]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-workplaces
  set workplaces (turtle-set )
  repeat min (list (max (list 1 round ((count people) / People-To-1-Workplace))) (count patches with [not any? locations-here])) [
;  repeat max (list 1 round ((count people) / People-To-1-Workplace)) [
    ask new-workplace "workplace" "factory" blue [move-to vacant-patch]
  ]

  ask people with [p-age >= 18 and p-age < 67] [
    if 90 > random-float 100 [ ; Some people might be unemployed.
      create-staffing-to one-of workplaces [
        set hidden? true
        set color blue + 1
      ]
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-schools
  ; Primary schools
  set schools (turtle-set )
  repeat max (list 1 round ((count people) / People-To-1-Primary-School)) [
    ask new-school "house ranch" [
      move-to vacant-patch
      set l-min-age 5
      set l-max-age 10
    ]
  ]
  ; Secondary School
  repeat max (list 1 round ((count people) / People-To-1-Secondary-School)) [
    ask new-school "house ranch" [
      move-to vacant-patch
      set l-min-age 11
      set l-max-age 18
    ]
  ]

  ask people with [p-age >= 5 and p-age <= 18] [
    if 95 > random-float 100 [ ; Some children might not be at school here.
      create-studying-to one-of schools with [l-min-age <= [p-age] of myself and l-max-age >= [p-age] of myself] [
        set hidden? true
        set color yellow + 1
      ]
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-hospitals
  set hospitals (turtle-set )
  let cur-patch nobody
  repeat max (list 1 round ((count people) / People-To-1-Hospital)) [
    set cur-patch vacant-patch
    ask new-hospital [
      move-to cur-patch
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report vacant-patch
  if any? patches with [not any? (locations-here with [not hidden?])] [
    report one-of (patches with [not any? locations-here with [not hidden?]])
  ]
  report false
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Create new locations and people, and give them their initial attributes.

to-report location-type-names
  report (list "household" "school" "workplace" "hospital")
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report location-types
  ; Obviously, the order should match that of names above.
  ; Every location's l-type will index this.
  report (list households schools workplaces hospitals)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to initialize-new-location
  set l-sc-rate 0
  set l-ic-rate 0
  set l-nc-rate 0
  set l-infection-rate 0
  set l-dst-rate 0
  set l-num-infections-here 0
  set l-contents array:from-list n-values (count disease-states) [(list )]
  set label-color black
  relabel-venue
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report new-household
  let retobj nobody
  create-locations 1 [
    set hidden? not (current-view = "City")
    set l-type position "household" location-type-names
    set l-type-name item l-type location-type-names
    set households (turtle-set households self)
    set color grey - 2
    set shape "house"
    set size 1 * base-size
    initialize-new-location
;    set l-health-boost round (1.0 * 10) ; To be used to alter ds transition rates. NB: Multiplied by 10 and converted to integer.
    ifelse No-Isolation-In-Household? [
      set l-infection-boost round (1.0 * 10) ; To be used to alter infection rates. NB: Multiplied by 10 and converted to integer.
    ]
    [
      set l-infection-boost round (0.0 * 10) ; To be used to alter infection rates. NB: Multiplied by 10 and converted to integer.
    ]
    set retobj self
  ]
  report retobj
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report new-person [age sex]
  let retobj nobody
  create-people 1 [
    set hidden? not (current-view = "City")
    set p-state susceptible
    set p-time-of-last-change -1 ; Be careful if you analyze this. In reality, people have been susceptible all their lives...
    set p-ds-history []
    set p-infection-time infinity
    set p-num-infected 0
    set p-age age
;    array:set age-groups p-age-group (fput self array:item age-groups p-age-group)
    set p-sex sex ; let 1=female, 0=male.
    set p-susceptibility round (100 * sqrt (base-transmission-chance / 100)) ; NB: Converted back to integer.
    set p-infectiousness round (100 * sqrt (base-transmission-chance / 100)) ; NB: Converted back to integer.
    set p-contacts round (Contacts-Per-Hour * 24) ; Converted to per-day in order to match other rates.
    set color ifelse-value (p-sex = 1) [sky] [pink]
    set size ifelse-value (p-age < 20) [0.25 * base-size] [0.5 * base-size]
    set shape "person"
    set p-off-school? false
    set p-off-work? false
    set retobj self
  ]
  report retobj
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report new-workplace [given-type given-shape given-color]
  let retobj nobody
  create-locations 1 [
    set hidden? not (current-view = "City")
    set l-type position "workplace" location-type-names
    set l-type-name item l-type location-type-names
    set workplaces (turtle-set workplaces self)
    set color given-color
    set size 1 * base-size
    set shape given-shape
    initialize-new-location
;    set l-health-boost round (1.0 * 10) ; To be used to alter ds transition rates. NB: Multiplied by 10 and converted to integer.
    set l-infection-boost round (1.0 * 10) ; To be used to alter infection rates. NB: Multiplied by 10 and converted to integer.
    set retobj self
  ]
  report retobj
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report new-hospital
  let retobj nobody
  create-locations 1 [
    set hidden? not (current-view = "City")
    set l-type position "hospital" location-type-names
    set l-type-name item l-type location-type-names
    set hospitals (turtle-set hospitals self)
    set color cyan
    set size 1 * base-size
    set shape "hospital"
    initialize-new-location
;    set l-health-boost round (1.0 * 10) ; To be used to alter ds transition rates. NB: Multiplied by 10 and converted to integer.
    set l-infection-boost round (1.0 * 10) ; To be used to alter infection rates. NB: Multiplied by 10 and converted to integer.
    set retobj self
  ]
  report retobj
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report new-school [given-shape]
  let retobj nobody
  create-locations 1 [
    set hidden? not (current-view = "City")
    set l-type position "school" location-type-names
    set l-type-name item l-type location-type-names
    set schools (turtle-set schools self)
    set color yellow
    set size 1 * base-size
    set shape given-shape
    initialize-new-location
;    set l-health-boost round (1.0 * 10) ; To be used to alter ds transition rates. NB: Multiplied by 10 and converted to integer.
    set l-infection-boost round (1.0 * 10) ; To be used to alter infection rates. NB: Multiplied by 10 and converted to integer.
    set retobj self
  ]
  report retobj
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to relabel-venue
  set label (word l-num-infections-here "     ")
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Code related to Contacts Matrices.

to setup-contacts-matrices
;  set contacts-matrix-home matrix:from-row-list contacts-matrix-home-reporter
;  set contacts-matrix-school matrix:from-row-list contacts-matrix-school-reporter
;  set contacts-matrix-work matrix:from-row-list contacts-matrix-work-reporter
;  set contacts-matrix-other matrix:from-row-list contacts-matrix-other-reporter

  setup-contacts-matrices-from-files

  ; These are now replaced by setup-intervention "":
;  set contacts-matrix-weekday-school contacts-matrix-all-reporter 100 100 100 100 ; home school work other
;  set contacts-matrix-weekday-no-school contacts-matrix-all-reporter 100 0 100 100
;  set contacts-matrix-weekend contacts-matrix-all-reporter 100 0 Weekend-Work 100 ; How much to adjust Other?

  ;  set contacts-matrix-all contacts-matrix-weekday-school
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report contacts-matrix-all-reporter [h s w o]
  report (matrix:plus
    (matrix:times contacts-matrix-home (h / 100))
    (matrix:times contacts-matrix-school (s / 100))
    (matrix:times contacts-matrix-work (w / 100))
    (matrix:times contacts-matrix-other (o / 100))
    )
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report p-age-group
  ; Convert p-age into age-group for use with contacts matrices.
  if p-age >= 75 [report 15]
  report int (p-age / 5)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report contacts-matrix-home-reporter
  ; Constructed in Excel from 'United Kingdom of Great Britain-home.csv'
  let tmp-mat map [s -> but-last but-first (csv:from-row s "\t")] but-first (list
    "	group: 0-4	group: 5-9	group: 10-14	group: 15-19	group: 20-24	group: 25-29	group: 30-34	group: 35-39	group: 40-44	group: 45-49	group: 50-54	group: 55-59	group: 60-64	group: 65-69	group: 70-74	group: 75+	"	 ; 	0
    "	0.4788128	0.55185414	0.334323605	0.132361228	0.138531588	0.281604887	0.406440259	0.493947983	0.113301081	0.074682641	0.041964034	0.017983199	0.005536943	0.001421873	0	0.000505582	"	 ; 	group: 0-4
    "	0.263264243	0.918274813	0.524179768	0.116285331	0.025055685	0.16985875	0.448042262	0.577627602	0.32478196	0.072276845	0.026533231	0.005856645	0	0.002420402	0.000285809	0	"	 ; 	group: 5-9
    "	0.168812075	0.535714614	1.080219323	0.388591513	0.039314522	0.011057129	0.21680897	0.591908717	0.485585642	0.133959469	0.044333387	0.017595164	0.004676779	0.00751153	0.001198292	0.000239288	"	 ; 	group: 10-14
    "	0.093901281	0.153999088	0.417215826	0.979076097	0.128063144	0.033498119	0.061002703	0.253606638	0.421886057	0.206528142	0.074921369	0.048022516	0.007194711	0.024556583	0.002397444	0	"	 ; 	group: 15-19
    "	0.167946317	0.08083609	0.07284187	0.356303405	0.80479534	0.205901995	0.074874642	0.041384182	0.164106139	0.283182758	0.080998683	0.085310811	0.01643175	0.007560766	0.001051554	0.001579202	"	 ; 	group: 20-24
    "	0.489661848	0.296565218	0.039897703	0.066153309	0.12500062	0.65954203	0.210168617	0.025338351	0.007441191	0.046801553	0.167113552	0.11084633	0.033629319	0.000839367	0	2.94523E-06	"	 ; 	group: 25-29
    "	0.319984866	0.472630659	0.269616873	0.075588594	0.046645734	0.088492899	0.642245625	0.148784529	0.0325759	0.004944669	0.012322182	0.039788976	0.019188409	0.000965681	0	7.67613E-13	"	 ; 	group: 30-34
    "	0.37824281	0.70078354	0.564193426	0.196423246	0.023690027	0.008989378	0.087075144	0.590987777	0.153469321	0.003244555	0.023579649	0.008028595	0.011731299	0.006121734	0	0	"	 ; 	group: 35-39
    "	0.166028735	0.516161356	0.730664425	0.415001753	0.067417032	0.004119461	0.090163489	0.191874369	0.457516113	0.104913169	0.017196255	0.001054882	0.025936005	0.012079583	0.008403247	0.001241388	"	 ; 	group: 40-44
    "	0.127256603	0.146586243	0.334933702	0.706502619	0.352848641	0.083114422	0.010816395	0.076359468	0.071951891	0.533831753	0.112210304	0.028418428	0.012267441	3.04005E-06	0.001471558	0.007002297	"	 ; 	group: 45-49
    "	0.121649534	0.099598995	0.175445499	0.257908401	0.309365024	0.07016612	0.07047664	0.0364506	0.07931335	0.055175027	0.371421901	0.127825451	0.01491615	6.01014E-06	0.004000386	0.00243485	"	 ; 	group: 50-54
    "	0.018613581	0.005513579	0.10908009	0.278496861	0.229404069	0.18324847	0.063586133	0.015694582	0.006869126	0.077729635	0.094560049	0.389562252	0.097670212	0.00308734	0	0	"	 ; 	group: 55-59
    "	0.021279381	0	0.047773559	0.034614117	0.085021554	0.1007503	0.11050783	0.069982275	0.083090798	0.017293876	0.044752237	0.113862629	0.47874739	0.057291891	0.014530004	0	"	 ; 	group: 60-64
    "	0.045960895	0.070106321	0.11455988	0.150158357	0.015328847	0.007949634	0.021608141	0.096096571	0.216056745	0.027157126	0.017570266	0.047951496	0.082105837	0.510463324	0.054842564	0	"	 ; 	group: 65-69
    "	0	0.02235594	0.163429174	0.245762052	0.008405939	0	0	0	0.296636658	0.084490549	0.034410222	0	0.058216007	0.142592716	0.156981294	0.089754718	"	 ; 	group: 70-74
    "	0.020607383	0	0.036166193	0	0.018581732	1.01979E-06	0.041371276	0	0.070873007	0.054252436	0.076604655	0	0	0	0.110469936	0.270621588	"	 ; 	group: 75+
  )
  if 16 != length tmp-mat [user-message (word "ERROR in contacts-matrix-home!\n\nMatrix is " (length tmp-mat) " rows long.")]
  if 16 != length filter [x -> 16 = length x] tmp-mat [user-message (word "ERROR in contacts-matrix-home!\n\nSome rows do not have 16 columns.")]
  report tmp-mat
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report contacts-matrix-school-reporter
  let tmp-mat map [s -> but-last but-first (csv:from-row s "\t")] but-first (list
    "	group: 0-4	group: 5-9	group: 10-14	group: 15-19	group: 20-24	group: 25-29	group: 30-34	group: 35-39	group: 40-44	group: 45-49	group: 50-54	group: 55-59	group: 60-64	group: 65-69	group: 70-74	group: 75+	"	 ; 	0
    "	0.974577996	0.151369805	0.008748809	0.026279091	0.011128161	0.089104305	0.125477587	0.088318278	0.03718242	0.02940927	5.10911E-38	1.13982E-32	0.007584287	0.001516368	1.23262E-50	5.97486E-64	"	 ; 	group: 0-4
    "	0.240133743	2.345158876	0.063993981	0.003888028	0.018033721	0.036688664	0.063607893	0.060756045	0.062134577	0.047317782	0.042025891	0.00436578	0.003438226	5.95795E-35	0.000493943	2.0766E-118	"	 ; 	group: 5-9
    "	0.002820885	1.120104417	2.561841358	0.118940767	0.014715788	0.053985706	0.041360786	0.122566825	0.085524651	0.073356995	0.032159515	0.024845813	0.012029306	2.65862E-39	6.67045E-80	1.18947E-80	"	 ; 	group: 10-14
    "	0.041037894	0.081500583	1.168804006	4.14022281	0.062569848	0.116615502	0.07523341	0.076614828	0.073089363	0.054652401	0.044335663	0.026557116	0.004125758	1.38647E-79	1.26699E-64	5.2748E-157	"	 ; 	group: 15-19
    "	4.19673E-10	0.126450247	2.15658E-12	0.273508417	0.231775202	0.013723373	0.022558539	0.031110981	0.014231203	0.005304263	0.013296088	0.004265201	2.41843E-38	0.002491266	1.8147E-104	9.56827E-94	"	 ; 	group: 20-24
    "	0.10559814	0.066700269	0.014972679	9.32704E-18	4.45743E-16	0.055416506	0.038981135	0.055495541	0.042656861	0.025310017	5.22686E-74	7.65054E-50	4.3106E-63	3.20011E-58	1.1043E-112	4.704E-138	"	 ; 	group: 25-29
    "	9.08419E-13	0.101860702	0.010341629	0.011748017	0.005535669	2.27561E-16	0.037791249	1.02626E-21	0.017300561	2.19613E-48	5.26097E-54	4.26686E-73	0.00440714	1.31285E-38	9.04985E-75	5.9746E-166	"	 ; 	group: 30-34
    "	0.035954028	0.170694149	0.073485368	0.030045224	2.9199E-19	0.024287958	0.052545181	0.072197497	0.066340666	0.045916143	2.12417E-39	0.012848423	0.00384942	0.011717019	4.88904E-97	2.56518E-88	"	 ; 	group: 35-39
    "	0.070026978	0.101401068	0.034330893	1.13431E-28	0.043343021	0.058747861	0.03242495	0.058476294	0.020853638	0.007773407	0.006905305	0.00617465	4.38008E-46	3.90102E-68	4.6908E-106	2.7621E-130	"	 ; 	group: 40-44
    "	0.015334051	2.79942E-63	0.017465773	0.207384215	4.23521E-50	3.34923E-36	0.045138288	0.009703647	0.067049402	0.048510635	0.062647584	0.039001299	0.007406986	9.9152E-119	5.1507E-112	5.2085E-116	"	 ; 	group: 45-49
    "	1.66847E-64	6.53088E-32	1.33939E-24	4.12769E-65	5.7437E-78	1.74433E-90	2.13522E-93	2.49943E-19	0.018665995	3.93009E-28	3.02016E-31	0.040803405	1.37379E-46	5.46886E-68	4.1042E-139	1.5455E-131	"	 ; 	group: 50-54
    "	0.047243691	0.139679917	0.050231035	4.67713E-29	2.06363E-39	4.0382E-105	0.031512485	1.27298E-35	0.150283699	0.033482205	2.15611E-42	0.076583168	1.43175E-28	1.469E-106	9.60833E-84	1.3197E-132	"	 ; 	group: 55-59
    "	1.20464E-55	0.128480957	3.64534E-28	3.6267E-61	0.019213852	0.007529818	0.046279882	0.112103189	0.020537965	0.034763688	0.007303726	0.00666334	1.06836E-56	0.004734946	6.15573E-85	5.0531E-99	"	 ; 	group: 60-64
    "	6.77119E-75	0.064070384	0.038699614	1.38552E-41	9.87473E-58	3.59626E-62	0.035009878	0.034592457	4.58619E-66	3.70697E-37	2.54826E-55	0.034236651	6.92907E-79	3.38216E-97	2.81377E-76	4.1457E-142	"	 ; 	group: 65-69
    "	1.0777E-128	1.83356E-67	3.22342E-58	2.5405E-102	1.86466E-74	8.7473E-179	2.11147E-55	3.80704E-77	3.3382E-103	5.8212E-122	2.7422E-101	5.36543E-88	5.8783E-114	2.81586E-98	3.1149E-135	2.0247E-131	"	 ; 	group: 70-74
    "	5.1675E-145	3.9448E-99	2.86063E-60	0.058663049	8.06704E-60	4.7091E-107	1.7098E-133	2.3279E-102	7.4213E-124	4.6866E-136	1.2654E-131	1.9791E-91	1.023E-133	8.9234E-124	6.5736E-126	9.049E-118	"	 ; 	group: 75+
  )
  if 16 != length tmp-mat [user-message (word "ERROR in contacts-matrix-school!\n\nMatrix is " (length tmp-mat) " rows long.")]
  if 16 != length filter [x -> 16 = length x] tmp-mat [user-message (word "ERROR in contacts-matrix-school!\n\nSome rows do not have 16 columns.")]
  report tmp-mat
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report contacts-matrix-work-reporter
  let tmp-mat map [s -> but-last but-first (csv:from-row s "\t")] but-first (list
    "	group: 0-4	group: 5-9	group: 10-14	group: 15-19	group: 20-24	group: 25-29	group: 30-34	group: 35-39	group: 40-44	group: 45-49	group: 50-54	group: 55-59	group: 60-64	group: 65-69	group: 70-74	group: 75+	"	 ; 	0
    "	0	0	0	0	0	0	0	0	0	0	0	0	0	1.849E-124	8.0627E-118	1.5103E-125	"	 ; 	group: 0-4
    "	0	0	0	0	0	0	0	0	0	0	0	0	0	4.9902E-109	3.4459E-132	1.2686E-134	"	 ; 	group: 5-9
    "	0	0	0.017070652	3.39774E-58	5.21109E-55	1.77192E-62	0.103290434	1.96122E-41	0.039028736	1.97799E-44	1.59468E-41	3.90847E-47	1.67091E-52	3.05322E-67	4.5715E-105	1.71365E-89	"	 ; 	group: 10-14
    "	0	0	0.011126026	1.082806495	0.866967626	0.277814213	0.107287816	0.281420074	0.151004085	0.482846088	0.322855848	0.061746494	0.017070012	9.9437E-06	3.67139E-06	2.32484E-60	"	 ; 	group: 15-19
    "	0	0	0.020101138	0.43798595	0.640587498	0.709850066	0.557489171	0.712702933	0.350716589	0.371190606	0.173785224	0.162704445	0.018825423	1.03983E-05	3.77158E-25	6.28789E-06	"	 ; 	group: 20-24
    "	0	0	0.008744471	0.487143168	0.613381662	0.836808081	0.766883222	0.623447246	0.947190892	0.539282641	0.549818808	0.090587144	0.022345947	3.10747E-05	2.16644E-05	3.05975E-46	"	 ; 	group: 25-29
    "	0	0	0.048424685	0.098423131	0.592950691	0.626772329	0.646893403	1.045224552	0.581682804	0.687889559	0.270550441	0.214812973	0.018043267	1.16739E-05	1.78221E-06	2.15436E-17	"	 ; 	group: 30-34
    "	0	0	0.01874034	0.327483491	0.366842833	0.591763122	0.536313765	0.702885834	0.624741019	0.660862381	0.523040167	0.193609278	0.007908794	5.88072E-06	9.26733E-06	9.30126E-06	"	 ; 	group: 35-39
    "	0	0	0.00558367	0.229767122	0.392219793	0.533086777	0.551139706	0.830081835	0.877503953	0.768856044	0.349715761	0.249994838	0.024571673	1.36056E-05	4.58282E-06	4.59779E-06	"	 ; 	group: 40-44
    "	0	0	0.263051547	0.28006664	0.39335553	0.757348378	0.756714423	0.679089053	0.84906079	1.044428958	0.424243495	0.194960632	0.029805995	1.32082E-05	1.71616E-05	5.85948E-06	"	 ; 	group: 45-49
    "	0	0	1.74585E-18	2.05216E-31	0.121657066	0.425861003	0.47269345	0.475417257	0.881593841	0.517705523	0.415224535	0.156205347	0.053766869	2.68951E-05	2.64826E-05	3.43401E-05	"	 ; 	group: 50-54
    "	0	0	0.005371844	0.085803538	0.329750962	0.357451427	0.312520128	0.249804092	0.292846051	0.327971964	0.295271255	0.140112625	0.013452652	2.68206E-06	1.79208E-05	1.11046E-19	"	 ; 	group: 55-59
    "	0	0	1.31286E-63	6.70807E-30	0.031718461	0.067430268	0.031800233	0.088491662	0.074881205	0.052168593	0.059257257	0.041984247	0.002035294	6.73993E-18	4.66687E-06	4.52861E-06	"	 ; 	group: 60-64
    "	2.5841E-111	7.392E-130	5.33141E-87	3.78778E-05	6.22554E-05	8.71081E-05	6.24939E-05	3.79228E-05	8.74564E-05	3.81468E-05	6.26558E-05	8.79221E-05	6.20689E-05	1.36475E-05	3.77999E-05	1.98536E-47	"	 ; 	group: 65-69
    "	2.5064E-126	2.9992E-135	1.21231E-89	9.79383E-70	5.30157E-66	1.31128E-65	8.71796E-63	2.46221E-70	4.22266E-70	1.58532E-82	1.40542E-86	1.17749E-85	2.02032E-67	5.50643E-66	2.02546E-87	1.10077E-89	"	 ; 	group: 70-74
    "	4.2709E-122	8.4473E-119	7.9353E-138	2.4183E-104	1.4525E-126	1.8901E-117	9.0217E-119	4.4564E-114	2.0025E-100	5.0168E-115	1.7303E-131	2.388E-101	2.9045E-115	5.8158E-106	1.2348E-102	6.0273E-149	"	 ; 	group: 75+
  )
  if 16 != length tmp-mat [user-message (word "ERROR in contacts-matrix-work!\n\nMatrix is " (length tmp-mat) " rows long.")]
  if 16 != length filter [x -> 16 = length x] tmp-mat [user-message (word "ERROR in contacts-matrix-work!\n\nSome rows do not have 16 columns.")]
  report tmp-mat
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report contacts-matrix-other-reporter
  let tmp-mat map [s -> but-last but-first (csv:from-row s "\t")] but-first (list
    "	group: 0-4	group: 5-9	group: 10-14	group: 15-19	group: 20-24	group: 25-29	group: 30-34	group: 35-39	group: 40-44	group: 45-49	group: 50-54	group: 55-59	group: 60-64	group: 65-69	group: 70-74	group: 75+	"	 ; 	0
    "	0.257847576	0.100135168	0.045803677	0.127084549	0.187303683	0.257979215	0.193228849	0.336594917	0.30922329	0.070538523	0.152218422	0.113554852	0.061577148	0.040429874	0.037356499	0.006697816	"	 ; 	group: 0-4
    "	0.178249137	0.770818807	0.126353007	0.089393832	0.057506323	0.219815579	0.227399212	0.177127045	0.263450303	0.170585267	0.158471291	0.062805024	0.088913253	0.055617614	0.038755915	0.013488697	"	 ; 	group: 5-9
    "	0.145914285	0.352969647	0.882339649	0.310008842	0.087803316	0.211294343	0.09111423	0.165645507	0.25510578	0.133087803	0.109222892	0.051837944	0.061615692	0.024815363	0.0392821	0.052295804	"	 ; 	group: 10-14
    "	0.041505375	0.260686913	0.681929499	1.671008065	0.266344449	0.204639942	0.157765096	0.372781663	0.254490869	0.208705381	0.158881713	0.082891945	0.058742615	0.036233733	0.012924995	0.000626015	"	 ; 	group: 15-19
    "	0.152449194	0.040352731	0.17734556	0.956584354	0.747473605	0.365269979	0.242523037	0.297250794	0.162561984	0.318278622	0.09373735	0.135322516	0.08498937	0.027664918	0.021409597	0.049716507	"	 ; 	group: 20-24
    "	0.246624662	0.093251822	0.079041606	0.204867788	0.992753107	0.853320503	0.457489561	0.265745555	0.287944489	0.321539788	0.199559888	0.095429894	0.093161241	0.054518163	0.017075012	3.48407E-06	"	 ; 	group: 25-29
    "	0.16330337	0.16928499	0.090773094	0.180129515	0.317586122	0.439556514	0.668057369	0.444733217	0.231203718	0.228177542	0.220039677	0.246935289	0.101243142	0.057276834	0.009740622	0.024821296	"	 ; 	group: 30-34
    "	0.12896771	0.219607305	0.095860135	0.099670201	0.239546405	0.298953187	0.306342303	0.521367017	0.510712253	0.251296074	0.147519291	0.196720545	0.207242065	0.10056861	0.111409484	0.038388513	"	 ; 	group: 35-39
    "	3.75239E-11	0.252235248	0.150369063	0.1459222	0.266348797	0.177376873	0.363119379	0.348387652	0.3938907	0.389781423	0.238353809	0.197666911	0.09668136	0.023403792	0.053226524	7.30242E-05	"	 ; 	group: 40-44
    "	0.000138622	0.007329721	0.02827579	0.129584707	0.18685834	0.131553248	0.169202815	0.393752639	0.464534939	0.556651247	0.318427668	0.122253182	0.162523947	0.099046429	0.061370074	0.059666354	"	 ; 	group: 45-49
    "	0.028355129	0.030461309	0.101989787	0.386949156	0.298423992	0.477952248	0.236001646	0.464053305	0.406314308	0.603299664	0.231973686	0.291053827	0.197333552	0.14377182	0.124694792	0.009826119	"	 ; 	group: 50-54
    "	0.106247959	0.062183545	0.058521372	0.108802565	0.356437338	0.574277189	0.546699177	0.462661146	0.486305996	0.220824496	0.462143721	0.639590054	0.397758991	0.204284068	0.123780529	0.109684868	"	 ; 	group: 55-59
    "	6.06193E-08	0.094685263	0.065270552	0.025313227	0.286419963	0.216646288	0.283354109	0.312674071	0.32230567	0.270592125	0.277529477	0.398381793	0.25876186	0.094789369	0.113889277	0.081778403	"	 ; 	group: 60-64
    "	0.018954296	0.110966265	0.062998679	0.063055707	0.107816765	0.493773233	0.288114222	0.462696678	0.417082281	0.430372098	0.387893725	0.516293634	0.289415248	0.215854635	0.136359067	0.056493948	"	 ; 	group: 65-69
    "	0.072816898	0.026073871	0.090931519	0.242115267	0.506491395	0.439531838	0.24289862	0.138245804	0.488467243	0.356689765	0.337426065	0.467690646	0.482265927	0.551656694	0.241565408	0.24924717	"	 ; 	group: 70-74
    "	1.35425E-09	0.000688316	0.04139986	0.111795121	0.068446769	0.038603864	0.129207583	0.215447768	0.198027155	0.406778809	0.352954451	0.000198892	0.466582938	0.000338157	0.192071526	0.467710373	"	 ; 	group: 75+
  )
  if 16 != length tmp-mat [user-message (word "ERROR in contacts-matrix-other!\n\nMatrix is " (length tmp-mat) " rows long.")]
  if 16 != length filter [x -> 16 = length x] tmp-mat [user-message (word "ERROR in contacts-matrix-other!\n\nSome rows do not have 16 columns.")]
  report tmp-mat
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-school-holidays
  ; Dates taken from UK.R
  set start-date-as-num parsed-date Start-Date
  set school-holiday-dates (map [[c o] -> (list c o)]
    map [x -> (parsed-date x) - start-date-as-num] ["2020-2-16" "2020-4-05" "2020-5-24" "2020-7-22" "2020-10-25" "2020-12-20" "2021-02-14" "2021-04-01" "2021-05-30" "2021-07-25"]
    map [x -> (parsed-date x) - start-date-as-num] ["2020-2-22" "2020-4-18" "2020-5-30" "2020-9-01" "2020-10-31" "2021-01-02" "2021-02-20" "2021-04-17" "2021-06-05" "2021-09-01"]
  )
  schedule-school-holiday
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report next-school-holiday
  ; school-holidays is a sorted list of 2-item lists [school-closed-date re-opened-date]
  let cur-day sim-day
  foreach school-holiday-dates [h ->
    if cur-day <= last h [report h]
  ]
  report (list infinity infinity)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to schedule-school-holiday
  ; Works out dates of current/next school holiday.
  ; Set flag school-holiday?
  ; If currently a holiday, schedules reopening.
  ; Else schedules closure.
  let next-hols next-school-holiday
  set school-holiday? sim-day >= first next-hols
  ifelse school-holiday? [
    schedule ((1 + last next-hols) * 1440) nobody [-> schedule-school-holiday] nobody "School reopens after holiday."
  ]
  [
    schedule (first next-hols * 1440) nobody [-> schedule-school-holiday] nobody "School closes for holiday."
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report parsed-date [given-string]
  ; Convert date string in format yyyy-m-d
  ; to a number. Used with start-date for parsing school holiday dates.
  let sep1 position "-" given-string
  let sep2 sep1 + 1 + position "-" substring given-string (sep1 + 1) (length given-string)
  let y cnum substring given-string 0 sep1
  let m cnum substring given-string (sep1 + 1) sep2
  let d cnum substring given-string (sep2 + 1) (length given-string)
  ; Treat 2020 as year 0. 2020 was a leapyear.
  let ret sum n-values (y - 2020) [x -> ifelse-value (0 = x mod 4) [366] [365]]
  set ret ret + ifelse-value (0 = y mod 4) [sum sublist [31 29 31 30 31 30 31 31 30 31 30 31] 0 (m - 1)] [sum sublist [31 28 31 30 31 30 31 31 30 31 30 31] 0 (m - 1)]
  ; 2020-1-1 is day 0
  set ret ret + d - 1
  report ret
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report date-string [given-sim-time]
  let tmp-sum start-date-as-num + int (given-sim-time / 1440)
  ; Day0 is 2020-1-1
  let cur-year 2020
  set tmp-sum tmp-sum - ifelse-value (0 = cur-year mod 4) [366] [365]
  while [tmp-sum > 0] [
    set cur-year cur-year + 1
    set tmp-sum tmp-sum - ifelse-value (0 = cur-year mod 4) [366] [365]
  ]
  set tmp-sum tmp-sum + ifelse-value (0 = cur-year mod 4) [366] [365]
  let cur-month 0
  set tmp-sum tmp-sum - item cur-month [31 29 31 30 31 30 31 31 30 31 30 31]
  while [tmp-sum > 0] [
    set cur-month cur-month + 1
    set tmp-sum tmp-sum - item cur-month [31 29 31 30 31 30 31 31 30 31 30 31]
  ]
  set tmp-sum tmp-sum + item cur-month [31 29 31 30 31 30 31 31 30 31 30 31]
  report (word cur-year "-" (cur-month + 1) "-" (1 + int tmp-sum)) ; yyyy-m-d
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report cnum [given-string]
  let ret 0
  let cur-item 0
  repeat length given-string [
    set ret 10 * ret + position (item cur-item given-string) ["0" "1" "2" "3" "4" "5" "6" "7" "8" "9"]
    set cur-item cur-item + 1
  ]
  report ret
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-city-contact-matrices
  ; If using contacts matrices then no need for locations other than hospital.

  setup-hospitals

  ;setup-contacts-matrices ; Using hard-coded procedures instead.

  if any? patches with [1 < count locations-here] [
    user-message-error "FATAL! Some patches have multiple locations on them."
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-city-contacts-generator
  ; (Not yet ready to use.)
  setup-hospitals

  foreach (list
    (list 50 "busroutes" "bus" green "transport")
    (list 10 "tramroutes" "tram" green "transport")
    (list 20 "pub-cafes" "food" violet "non-work")
    (list 10 "coachroutes" "bus" green "transport")
    (list 10 "trainroutes" "train" green "transport")
    (list 10 "sports-centres" "ball football" violet "non-work")
    (list 20 "rec-centres" "ball tennis" violet "non-work")
    (list 10 "hotels" "building colonial" violet "non-work")
    (list 5 "parks" "tree" violet "parks")
    (list 50 "other-locations" "triangle" violet "parks")
    (list 50 "shops" "building store" violet "non-work")

    (list 10 "schools-pre" "house ranch" yellow "non-work")
    (list 10 "schools-primary" "house ranch" yellow "non-work")
    (list 10 "schools-secondary" "house ranch" yellow "non-work")
    (list 10 "schools-FE" "building institution" yellow "non-work")
    (list 10 "schools-HE" "building institution" yellow "non-work")

    (list 200 "MT-workplaces" "person business" blue "work")
    (list 0 "NA" "wheel" blue "work")
    (list 200 "S-workplaces" "person service" blue "work")
    (list 200 "P-workplaces" "person doctor" blue "work")
    (list 200 "SM-workplaces" "person construction" blue "work")
    (list 200 "SnM-workplaces" "person graduate" blue "work")
    (list 200 "U-workplaces" "person farmer" blue "work")
    ) [cur-type ->

    create-locations (item 0 cur-type) [
      move-to vacant-patch
      set l-type item 1 cur-type
      set shape item 2 cur-type
      set color item 3 cur-type
    ]
  ]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Setup initial diseased population.

to setup-disease-people
  ; Base case scenario:
  ; When simulated time starts, there is just 1 infectious person in the population.

  set deaths []

  ask people [
    set p-state susceptible
    array:set (array:item age-groups p-age-group) ([who] of p-state) (fput self array:item (array:item age-groups p-age-group) ([who] of p-state))
    ask p-location [array:set l-contents ([[who] of p-state] of myself) (fput myself array:item l-contents ([[who] of p-state] of myself))]

    ask p-state [adjust-cur-num 1]
    set p-ds-history fput (list sim-time p-state) p-ds-history
    recolor-by-state
  ]
  ;infect-one-susceptible
  ; Now create seed infection events later in Setup.
  set post-disease-people (turtle-set )
  set superspreaders (turtle-set )
  set max-infected-by-one 0

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-seed-infectors-n-per-day-for-t-days-from-d [n t d]
  let initial-day 0
  ifelse is-number? d [
    set initial-day d
  ]
  [
    set initial-day runresult d  ; LSHTM used 7 for London Boroughs and 21 for other counties.
  ]
  let a int n
  let b n - a
  set seed-start-day initial-day
  foreach n-values t [x -> x] [y ->
;    print (word "y = " y)
    foreach n-values a [z -> z + 1] [i ->
;      print (word "i = " i)
      schedule (last-midnight + ((initial-day + y) * 1440) + random 1440) nobody [-> infect-one-susceptible] nobody (word "Seed infector day " y " addition " i)
    ]
    if b > random-float 1 [
      schedule (last-midnight + ((initial-day + y) * 1440) + random 1440) nobody [-> infect-one-susceptible] nobody (word "Seed infector day " y ", chance addition.")
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to infect-one-susceptible
  ; Make one randomly chosen susceptible infectious.
  ; Infected by someone outside the system.
  ; Called during Setup. Maybe called from a button.

  if not any? seed-age-susceptibles [stop] ; Could happen if called from user button, or if population very small.
  ask one-of seed-age-susceptibles [
    set seed-infected self ; To study relation between epidemic outcome and seed person's attributes.
    schedule sim-time self [-> become-infected-by nobody p-location] p-location "One-off infection."
;    schedule (sim-time + random (1440 * 7)) self [-> become-infected-by nobody p-location] p-location "One-off infection."
    ; NB: This event is only scheduled. It still needs processing.
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report seed-age-susceptibles
  ; See use of dist_seed_ages in UK.R, section "# 4a. Set parameters" for age range.
  ; Bit arbitrary, no?
;  report (reduce sentence (map [ag-id -> array:item (array:item age-groups ag-id) 0] (range 5 10)))
  report (turtle-set (map [ag-id -> array:item (array:item age-groups ag-id) 0] (range 5 10)))
;  report people with [susceptible? and p-age >= 25 and p-age < 50]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to adjust-cur-num [given-adjustment]
  ; Run by disease-state
  set ds-cur-num ds-cur-num + given-adjustment
  if ds-cur-num > ds-max-num [
    set ds-max-num ds-cur-num
    set ds-max-at-time sim-time
  ]
  relabel-disease-state
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to recolor-by-state
  set color [ds-color] of p-state
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Useful reporters for converting times

to-report sim-time-as-hh-mm
  ; Reports sim-time as a time of day in neat format.
  ; For display purposes, not experiment output.
  let cur-hour (word "0" sim-hour)
  let cur-min (word "0" sim-minutes)
  report (word
    (substring cur-hour (-2 + length cur-hour) (length cur-hour))
    ":"
    (substring cur-min (-2 + length cur-min) (length cur-min))
    )
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report sim-day
  ; Assuming sim-time is in minutes.
  report int (sim-time / (60 * 24))
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report sim-hour
  ; Used by clock-time.
  ; Assuming sim-time is in minutes.
  report int ((sim-time / 60) mod 24)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report sim-minutes
  ; Used by clock-time.
  ; Assuming sim-time is in minutes.
  report int (sim-time mod 60)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report time-as-days [given-time]
  ; Converts minutes into days.
  ; Reports a float.
  ; To compute a particular day, use
  ; int time-as-days given-time
  ; To convert a duration to days, use
  ; round time-as-days given-duration
  report given-time / (24 * 60)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report sim-days
  ; How many complete days have passed since sim start.
  report int time-as-days sim-time
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report last-midnight
  ; Midnight before the current time.
  ; Add 24 * 60 minutes to this to get the next midnight,
  ; or hh * 60 + mm to get a particular time of day hh:mm.
  report (24 * 60) * int (sim-time  / (24 * 60))
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report sim-time-of-day
  ; Minutes past midnight.
  report sim-time mod (24 * 60)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report day
  ; Reports the current day of the week as a text abbreviation.
  ; NB: Monday is day 0, Saturday day 5, Sunday day 6
  report item (( sim-time / (24 * 60)) mod 7) ["Mon" "Tue" "Wed" "Thu" "Fri" "Sat" "Sun"]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report workday?
  ; NB: Monday is day 0, Saturday day 5, Sunday day 6
  report (((sim-time / (24 * 60)) mod 7) < 5)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report weekend?
  ; NB: Monday is day 0, Saturday day 5, Sunday day 6
  report (((sim-time / (24 * 60)) mod 7) >= 5)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Implements the 3-phase (ABC) approach to discrete-event simulation (DES) of K.D. Tocher.

to go [step-time]
  if time-taken-go = false [set time-taken-go timer] ; Start timing.
  if print-scheduling? or Print-Processing-B-Events? [print "" tick]

  while [sim-time < step-time and not sim-stopping?] [
    ; Step A: Advance clock.
    set sim-time b-events-next-time
    ;  b-events-shuffle-next-events ; Shuffle order of events at this time. NB: N simultaneous events will need N random numbers for the permutation.

    ; Step B: Do Time-Bound Events.
    ; People relocating.
    ; Disease progression (other than transmission).
    ; Events were scheduled using
    ; b-events-put event-time (list given-agent given-task given-object given-reason sim-time)
    let cur-event []
    while [sim-time = b-events-next-time] [
      set cur-event last b-events-get
      if Print-Processing-B-Events? [
        print (word ticks ": " sim-time ": Event=" cur-event ", IR=" total-infection-rate ".")
      ]
      ifelse is-turtle? item 0 cur-event [
        ; Run as turtle (person, household, etc.)
        ask item 0 cur-event [
          ;        ifelse p-time-of-last-change <= item 4 cur-event [ ; Else event now obsolete so ignore it? Not debugged yet.
          run item 1 cur-event
          ;        ]
          ;        [
          ;          show-error (word "WARNING: Obsolete event: " cur-event)
          ;        ]
        ]
      ]
      [
        ; Run as observer.
        run item 1 cur-event
        if sim-stopping? [
          set time-taken-go timer - time-taken-go
          if Using-Profiler? [profiler:stop print profiler:report]
          stop
        ] ; In case sim-stopping? set by code just run.
          ;      if sim-stopping? [set b-events but-first b-events stop]
      ]
      ;    set b-events but-first b-events
    ]

    if sim-stopping? [stop]

    ; Step C: Do Conditional Events.
    if halt-when = "Infections Impossible" [if (0 = [ds-cur-num] of susceptible) [set sim-stopping-reason (word "Infections Impossible") set sim-stopping? true]]
  if halt-when = "All Post-Infectious" [if (count people) = [ds-cur-num] of recovered + [ds-cur-num] of dead [set sim-stopping-reason (word "All Post-Infectious") set sim-stopping? true]]
    if halt-when = "No Susceptible" [if (0 = [ds-cur-num] of susceptible) [set sim-stopping-reason (word "No Susceptible") set sim-stopping? true]]
    if halt-when = "6 Months" [if sim-time >= (1440 * 366 / 2) [set sim-stopping-reason (word "End of 6 months") set sim-stopping? true]]
    if halt-when = "9 Months" [if sim-time >= (1440 * 366 * 3 / 4) [set sim-stopping-reason (word "End of 9 months") set sim-stopping? true]]
    if halt-when = "2 Years" [if sim-time >= (1440 * (366 + 365)) [set sim-stopping-reason (word "End of 9 months") set sim-stopping? true]]
    if sim-stopping? [
      calc-stats
      set time-taken-go timer - time-taken-go
      if Using-Profiler? [profiler:stop print profiler:report]
      stop
    ]

    ;tick ; Remember: Plots etc. update on tick.
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to go-1-day
  set sim-stopping? false
;  schedule-stop-at sim-time + (24 * 60)
  go last-midnight + 1440
  if sim-stopping? [stop]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to go-1-week
  set sim-stopping? false
;  schedule-stop-at sim-time + (7 * 24 * 60)
  go last-midnight + (7 * 1440)
  if sim-stopping? [stop]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Discrete-Event Simulation
;; Code for processing list of time-Bound events (b-events).
;; b-events represents a time-priority queue.
;; Currently, b-events is structured as a tree:
;; 1 branch for each day on which events will occur,
;; 1 subbranch for each time at which events will occur,
;; 1 leaf for each event, at which the details of that event are stored.

; Methods: b-events-setup, b-events-shuffle-next-events, b-events-put
; Reporters: b-events-get, b-events-on, b-events-at, b-events-next-time
; Other b-events- procedures can be considered "internal". You shouldn't access them directly.

to b-events-setup
  set b-events []
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report b-events-next-time
  ; Reports time of first event.
  if empty? b-events [report false]

  report first first last first b-events
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report b-events-next-but-one-time
  ; Reports time of first event later than the time of first event.
  if empty? b-events [report false]
  ; if first day has 1 time, then need a second day.
  if 1 = length last first b-events [
    if 1 = length b-events [report false] ; There is only 1 day.
    report first first last item 1 b-events
  ]
  ; else report second time on first day.
  report first item 1 last first b-events
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to b-events-shuffle-next-events
  ; Shuffle the events scheduled to run at next time.
  ; (Running simultaneous in the order they were scheduled in
  ; might cause artifacts. If simultaneous events then scheduled simultaneous events for some future time,
  ; the order would repeated.
  ; If events schedule instantaneous events, these latter will wait until last and not be re-shuffled.
  ; - But then what sense does it make to have instantaneous events?
  ; Preserving the order may be useful when debugging. So maybe comment the next three lines out.
  let cur-day first first b-events
  let cur-time first first last first b-events
  set b-events fput (list cur-day fput (list cur-time shuffle last first last first b-events) but-first last first b-events) but-first b-events
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report b-events-on [given-day]
  ; Reports times list for given-day.
  ; b-events is (list day (times list)), sorted ascending by day.
  foreach b-events [day-times-pair ->
    if given-day = first day-times-pair [ report last day-times-pair ] ; Found it, so report it.
    if given-day < first day-times-pair [ report [] ] ; Missed it, so it can't be there.
  ]
  report [] ; All days earlier than given day.
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report b-events-at [given-time]
  ; Reports events list for given-time.
  ; b-events is (list day (list time (list event-details))), sorted ascending by day, ascending by time, reverse order inserted.
  let cur-day int (given-time / (60  * 24))
  foreach b-events-on cur-day [time-events-pair ->
    if given-time = first time-events-pair [ report last time-events-pair ] ; Found it, so report it.
    if given-time < first time-events-pair [ report [] ] ; Missed it, so it can't be there.
  ]
  report [] ; Either nothing on day, or failed to find time.
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report b-events-before [given-time]
  ; Reports events list before given-time but on same day of given-time.
  ; b-events is (list day (list (list time (list event-details)))), sorted ascending by day, ascending by time, reverse order inserted.
  let cur-day int (given-time / (60  * 24))
;  report filter [given-time > first ?] b-events-on cur-day ; Assume this is slower for days with lots of times (but is it?)

  let found-pos b-events-found-position-of-time given-time
  report sublist (b-events-on cur-day) 0 last found-pos
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report b-events-after [given-time]
  ; Reports events list after given-time but on same day of given-time.
  ; b-events is (list day (list (list time (list event-details)))), sorted ascending by day, ascending by time, reverse order inserted.
  let cur-day int (given-time / (60  * 24))
;  report filter [given-time < first ?] b-events-on cur-day ; Assume this is slower for days with lots of times (but is it?)

  let found-pos b-events-found-position-of-time given-time
  if first found-pos [
    report sublist (b-events-on cur-day) (1 + last found-pos) length (b-events-on cur-day)
  ]
  report sublist (b-events-on cur-day) (last found-pos) length (b-events-on cur-day)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report b-events-get
  ; Removes first event from b-events.
  ; Reports that event, or [] if no events.
  if empty? b-events [report []] ; No events yet.
  let event-to-return (list
    (first first last first b-events) ; Time of first time on first day.
    (first last first last first b-events) ; Details of first event at first time on first day.
    )

  let cur-day first first b-events
  let cur-time first first last first b-events
  let cur-events but-first last first last first b-events
  ifelse empty? but-first last first last first b-events [
    let cur-times but-first last first b-events
    ifelse empty? cur-times [ ; No other events on this day.
      set b-events but-first b-events
    ]
    [ ; No more events at this time, but other times on this day.
      set b-events fput (list cur-day but-first last first b-events) but-first b-events
    ]
  ]
  [ ; More events at this time.
    set b-events fput (list cur-day fput (list cur-time but-first last first last first b-events) but-first last first b-events) but-first b-events
  ]

  report event-to-return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to b-events-put [given-time given-event-details]
  ; Insert given event in b-events in the correct place.
  ; May already be events on same day.
  ; May already be events at same time.
  let cur-day int (given-time / (60  * 24))
  let found-pos b-events-found-position-of-day cur-day
  ifelse first found-pos [ ; Day was found at position.
    set b-events sentence
    (sublist b-events 0 last found-pos)
    ; NB: lput for lifo, fput for fifo. fput much faster if many events have same time.
    fput (list cur-day (sentence (b-events-before given-time) fput (list given-time (fput given-event-details b-events-at given-time)) (b-events-after given-time)))
    (sublist b-events (1 + last found-pos) (length b-events))
  ]
  [ ; Day was not found, but should be inserted at position.
    set b-events sentence
    (sublist b-events 0 last found-pos)
    ; NB: lput for lifo, fput for fifo. fput much faster if many events have same time.
    fput (list cur-day (sentence (b-events-before given-time) fput (list given-time (fput given-event-details b-events-at given-time)) (b-events-after given-time)))
    (sublist b-events (last found-pos) (length b-events))
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report b-events-test
  ; Used for debugging.
  report (list
    (list 0 ; Day 0
      (list ; List of times and events on that day.
        (list 60 ; Time 60 minutes
          (list ; List of events at that time.
            (list "A" "a")
            (list "B" "b")
            )
        )
        (list 90 ; Time 90 minutea
          (list ; List of events at that time.
            (list "C" "c")
            )
        )
      )
    )

    (list 1 ; Day 1
      (list ; List of times and events on that day.
        (list 1500 ; Time 1500 minutes
          (list ; List of events at that time.
            (list "D" "d")
            )
        )
        (list 1530 ; Time 1500 minutes
          (list ; List of events at that time.
            (list "E" "e")
            (list "F" "f")
            (list "G" "g")
            )
        )
      )
    )

  )
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to b-events-print
  ; Prints the b-events in readable format.
  ; Useful for debugging.
  let cur-day 0
  let cur-time 0
  foreach b-events [day-times-pair ->
    set cur-day first day-times-pair
    foreach last day-times-pair [time-events-pair ->
      set cur-time first time-events-pair
      foreach last time-events-pair [cur-event ->
        print (word "Day: " cur-day ", Time:" cur-time ", Details:" cur-event)
      ]
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report b-events-found-position-of-day [given-day]
  ; Reports (list day-found? position-to-insert-at)
  (foreach b-events (n-values (length b-events) [? -> ?]) [[day-times-pair cur-pos] ->
    if given-day = first day-times-pair [ report (list true cur-pos) ] ; Found it. Report position.
    if given-day < first day-times-pair [ report (list false cur-pos) ] ; Missed it. Report how many come before.
  ])
  report (list false (length b-events))
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report b-events-found-position-of-time [given-time]
  ; Reports (list time-found? position-to-insert-at)
  ; Not currently used and not tested.
  let cur-day int (given-time / (60  * 24))
  (foreach (b-events-on cur-day) (n-values (length (b-events-on cur-day)) [? -> ?]) [[time-events-pair cur-pos] ->
    if given-time = first time-events-pair [ report (list true cur-pos) ] ; Found it. Report position.
    if given-time < first time-events-pair [ report (list false cur-pos) ] ; Missed it. Report how many come before.
  ])
  report (list false (length (b-events-on cur-day)))
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Main scheduling procedure

to schedule [event-time given-agent given-task given-object given-reason]
  ; Insert the event in the list of time-bound events.
  if event-time < sim-time [
    user-message-error (word "ERROR! Trying to schedule an event that should be past.\n Event-time=" event-time ", Agent=" given-agent ", Object=" given-object ", Reason=" given-reason)
  ]

  b-events-put event-time (list given-agent given-task given-object given-reason sim-time)

  if print-scheduling? [
    ; Next line for debugging purposes:
    print (word ticks ": " sim-time ": Scheduled " given-agent " for time " event-time " to do " given-task " with " given-object " because of \"" given-reason "\"." " IR=" total-infection-rate ".")
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to schedule-next-day
  ; Run at midnight.
  ; Determines the day's main events.
  ; Also schedules itself for the next midnight.

  schedule (last-midnight + 24 * 60) nobody [-> schedule-next-day] nobody "Next Day"

  calc-stats

  do-time-series-plots ; Do these once a day, not on every tick.
  ; NB: Plots with their own code, and monitors will update on tick.
  if not print-scheduling? and not print-processing-b-events? [tick]

  ; Reset counters to 0.
  ask transitions [set ts-num-new-events fput 0 ts-num-new-events]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to schedule-stop-at [given-time]
  ; How to setup sim to halt after a fixed length of simulated time (e.g. 6 months).
  schedule given-time nobody [-> set sim-stopping? true set sim-stopping-reason (word "Scheduled stop made at " given-time)] nobody "STOP"
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to schedule-contacts-infections
  ; Uses random-binomial to sample how many susceptibles infected,
  ; then uses n-of to sample which susceptibles they are,
  ; and therefore uses fewer random numbers (1 Binom + p * N_s uniform distrib integers)
  let new-infecteds []
  let AGI age-groups-total-infectiousness ; Calc here, just to make sure we only do it once, not for every person.
  let AGN age-groups-total-people
  let lambda 0.0 ; Force of infection = Rate at which infection events occur to each susceptible during 1 unit of time (i.e. 1 day)
  let prob 0.0
;  ask people with [susceptible?] [ ; If we needed to do them in random order. Can always shuffle final new-infecteds instead.
  let ag-id 0
  foreach array:to-list age-groups [ag ->
    set lambda force-of-infection ag-id AGI AGN
    set prob ifelse-value (lambda > 0) [1.0 - exp (0 - lambda * 0.25)] [0]
    set ag-id ag-id + 1

    set new-infecteds sentence new-infecteds n-of (random-binomial (length (array:item ag 0)) prob) (array:item ag 0)
  ]

;  foreach shuffle new-infecteds [p -> ; Order shouldn't matter, no?
  foreach new-infecteds [p ->
    if member? p trace-turtles [print-trace-turtles (word "Became infected at " sim-time)]
    ask p [become-infected]
  ]

  schedule (sim-time + 6 * 60) nobody [-> schedule-contacts-infections] nobody "Infections via contacts matrices."
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to schedule-contacts-infections-bern
  ; Samples whether infected for each individual susceptible.
  ; Hence, if there are N_s susceptibles, need N_s Bernoulli-distributed random numbers.
  let new-infecteds []
  let AGI age-groups-total-infectiousness ; Calc here, just to make sure we only do it once, not for every person.
  let AGN age-groups-total-people
  let lambda 0.0 ; Force of infection = Rate at which infection events occur to each susceptible during 1 unit of time (i.e. 1 day)
  let prob 0.0
;  ask people with [susceptible?] [ ; If we needed to do them in random order. Can always shuffle final new-infecteds instead.
  let ag-id 0
  foreach array:to-list age-groups [ag ->
    set lambda force-of-infection ag-id AGI AGN
    set prob ifelse-value (lambda > 0) [1.0 - exp (0 - lambda * 0.25)] [0]
    set ag-id ag-id + 1

    foreach (array:item ag 0) [p -> ; Susceptibles only (disease-state 0)
      ask p [
        if member? self trace-turtles [print-trace-turtles (word "Force of Infection = " lambda)]

        ; Is the time to next infection event < time step? If so, then infection occurs in that time step.
        ; (Markovian process = Memoryless, so if event doesn't occur in this time step, forget the time and resample in next time step.)
        ; Recall: contacts matrix is per day, so force of infection is events per day, hence time step of 6 hours means t = 0.25 days.
        ; Sample a time to next infection event from exponential distribution with rate parameter lambda = force of infection.
        ; Alternatively: calculate critical probability for time step, and compare with sampled probability.
        ; Which is faster?
        ; if lambda > 0 [if 0.25 > 0 - (1 / lambda) * ln (1 - random-float 1) [set new-infecteds fput self new-infecteds]]
        ;if lambda > 0 [if 0.25 > random-exponential (1.0 / lambda) [set new-infecteds fput self new-infecteds]]
        if prob > random-float 1 [set new-infecteds fput self new-infecteds]
      ]
    ]
  ]

;  foreach shuffle new-infecteds [p -> ; Order shouldn't matter, no?
  foreach new-infecteds [p ->
    if member? p trace-turtles [print-trace-turtles (word "Became infected at " sim-time)]
    ask p [become-infected]
  ]

  schedule (sim-time + 6 * 60) nobody [-> schedule-contacts-infections] nobody "Infections via contacts matrices."
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report random-binomial [n p]
  ; Not currently used.
  ; Could use to sample which Susceptibles become infected,
  ; and which of the Exposed become symptomatic.
  ; Will be slow if n large and p = 0.5.
  ; NB: Knowing that x of n people undergo some process
  ; does not tell which x they are.
  ; A random permutation will still require x random numbers.
  let q 1.0 - p
  if p > q [report n - random-binomial n q]
  let u random-float 1.0
  let x q ^ n
  let pq p / q
  let r 0
  set u u - x
  while [u >= 0] [
    set r r + 1
;    set stored-r r
    set x x * (pq * (n - r + 1) / r)
    set u u - x
  ]
  report r
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report age-groups-total-people
  ; Returns list of age groups' total number of people in any disease-state.
  ; NB: This includes the Dead and the Hositalized!
  report map [ag-ds-array -> sum map [dspeople -> length dspeople] array:to-list ag-ds-array] array:to-list age-groups

;  report map [ag -> length ag] array:to-list age-groups

;  let non-dead-people (people with [not dead?]) ; Cuts down calls to dead?
;  report array:from-list n-values 16 [ag -> count non-dead-people with [p-age-group = ag]] ; Beats all the arithmetic in the next method.

  ;  report array:from-list reduce [[a b] -> (map [[x y] -> x + y] a b)] [replace-item p-age-group (n-values 16 [0]) (1)] of non-dead-people
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report age-groups-total-infectiousness
  ; Returns list of age groups' total relative infectiousness
  let ds-list filter [ds -> [ds-rel-inf > 0] of ds] compartments-sorted
  report map [ag-ds-array -> sum map [ds -> ([ds-rel-inf] of ds) * length array:item ag-ds-array [who] of ds] ds-list] array:to-list age-groups

;  report (map [ag -> (sum map [pp -> [[ds-rel-inf] of p-state] of pp] ag)] array:to-list age-groups)

  ;  let non-dead-people (people with [not dead?]) ; Cuts down calls to dead?
;  report n-values 16 [ag -> sum [[ds-rel-inf] of p-state] of non-dead-people with [p-age-group = ag]] ; Beats all the arithmetic in the next method.

  ;  report reduce [[a b] -> (map [[x y] -> x + y] a b)] [replace-item p-age-group (n-values 16 [0]) ([ds-rel-inf] of p-state)] of non-dead-people
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report susceptibility-given-r0 [given-r0]
  ; Hacked from cmmid:
;  report 0.08 * given-R0 / 3.211871916
  ; i.e. u = base-u * R0s[run] / cm_calc_R0(parametersUK1, 1)
  ; see lines 364, 374 in UK.R

  ; See supplementary information from LSHTM model.
  ; R0 = the absolute value of the dominant eigenvalue of the next generation matrix NGM.
  ; NGM_ij = u * c_ij * (y_j * E(d_P + d_C) + (1 - y_j) * f * E(d_S))
  ; where f = relative infectiousness of Subclinical,
  ; c_ij is the element i,j of the contacts matrix,
  ; E(d_X) is the expected duration in disease state X.

  let dp ([ds-rel-inf] of i-preclinical) * [ts-time-param] of transition-between I-Preclinical I-Clinical
  let dc ([ds-rel-inf] of i-clinical) * [ts-time-param] of transition-between I-Clinical Hospitalize ; What if you go to Dead instead? Or Recovered??
  let ds ([ds-rel-inf] of i-subclinical) * [ts-time-param] of transition-between I-Subclinical Recovered

  let A map [y -> y * (dp + dc) + (1 - y) * ds] n-values (first matrix:dimensions contacts-matrix-weekday-school) [i -> array:item Prop_symptomatic int (i / 2)]
;  print (word "Eigenvalues = " matrix:real-eigenvalues matrix:from-row-list map [r -> (map [[x y] -> x * y] A r)] matrix:to-row-list contacts-matrix-all)
  let dominant-ev max map abs matrix:real-eigenvalues matrix:from-row-list map [r -> (map [[x y] -> x * y] A r)] matrix:to-row-list contacts-matrix-weekday-school
  report given-R0 / dominant-ev
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report r0-given-susceptibility [given-susceptibility]
  ; See supplementary information from LSHTM model.
  ; R0 = the absolute value of the dominant eigenvalue of the next generation matrix NGM.
  ; NGM_ij = u * c_ij * (y_j * E(d_P + d_C) + (1 - y_j) * f * E(d_S))
  ; where f = relative infectiousness of Subclinical,
  ; c_ij is the element i,j of the contacts matrix,
  ; E(d_X) is the expected duration in disease state X.

  let dp ([ds-rel-inf] of i-preclinical) * [ts-time-param] of transition-between I-Preclinical I-Clinical
  let dc ([ds-rel-inf] of i-clinical) * [ts-time-param] of transition-between I-Clinical Hospitalize ; What if you go to Dead instead? Or Recovered??
  let ds ([ds-rel-inf] of i-subclinical) * [ts-time-param] of transition-between I-Subclinical Recovered

  let contacts-matrix-all contacts-matrix-weekday-school
  let A map [y -> y * (dp + dc) + (1 - y) * ds] n-values (first matrix:dimensions contacts-matrix-all) [i -> array:item Prop_symptomatic int (i / 2)]
  let dominant-ev max map abs matrix:real-eigenvalues matrix:from-row-list map [r -> (map [[x y] -> x * y] A r)] matrix:to-row-list contacts-matrix-all
  report given-susceptibility * dominant-ev
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-LSHTM-R0-and-Intervention
  ; LSHTM-Run is the run number from the LSHTM model, and it ranges from 1 to 200.
  let run-row []
  if Population-Generator = "Demographic Data" [
    set run-row (item (LSHTM-Run - 1) Runs-R0-SeedStart-Peakt)
    set r0 item 1 run-row
    set intervention-day (intervention-shift - (7 * 12 / 2) + last run-row)
    set seed-start-day item 2 run-row
    stop
  ]
  if Population-Generator = "Demographic Data : UK" [
    set r0 item 1 (item (LSHTM-Run - 1) r0-list-uk)
    set intervention-day (intervention-shift - (7 * 12 / 2) + last (item (LSHTM-Run - 1) r0-list-uk))
    set seed-start-day 0
    stop
  ]
  ; Default option:
  set r0 item 1 (item (LSHTM-Run - 1) r0-list-uk)
  set intervention-day (intervention-shift - (7 * 12 / 2) + last (item (LSHTM-Run - 1) r0-list-uk))
  set seed-start-day 0
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report Runs-R0-SeedStart-Peakt
  if not file-exists? (word input-data-folder "Runs_R0_SeedStart_Peakt.csv") [user-message "FATAL: Cannot find file Runs_R0_SeedStart_Peakt.csv in Input-Data-Folder." report false]
  report map [r -> but-first r] but-first csv:from-file (word input-data-folder "Runs_R0_SeedStart_Peakt.csv")
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report r0-list-UK
  ; These are the first 200 R0 values generated by the LSHTM model.
  ; Also included are the times of peak cases in their "Base" scenario.
  ; This is used by them to calculate the launch date of interventions
  ; (Peak Day - 7 * 12 / 2 + Intervention-Shift).
  ; Use these numbers to dock with their model.
  report (list
    [	1	3.260446601	69	]		[	2	2.014948188	110	]		[	3	2.566780636	89	]		[	4	2.622673171	85	]		[	5	2.686515811	81	]
    [	6	2.018551289	109	]		[	7	3.016344781	75	]		[	8	1.637423395	145	]		[	9	3.0412625	75	]		[	10	2.661695084	86	]
    [	11	3.880042323	63	]		[	12	1.829586253	118	]		[	13	2.614126358	81	]		[	14	1.992594982	109	]		[	15	2.473145606	90	]
    [	16	2.720384223	79	]		[	17	3.77726208	67	]		[	18	2.525646317	89	]		[	19	4.184524672	61	]		[	20	2.638017327	86	]
    [	21	3.024906879	73	]		[	22	3.330937782	70	]		[	23	2.624308704	87	]		[	24	2.533031307	88	]		[	25	2.137467441	102	]
    [	26	2.465403907	87	]		[	27	2.955696535	75	]		[	28	2.843714137	77	]		[	29	4.957751528	57	]		[	30	2.426514886	90	]
    [	31	3.129246088	74	]		[	32	2.763787931	78	]		[	33	2.337549531	94	]		[	34	2.647235151	86	]		[	35	2.184763885	101	]
    [	36	2.760946674	79	]		[	37	1.811946828	118	]		[	38	2.937667545	75	]		[	39	2.162668089	104	]		[	40	2.845462207	79	]
    [	41	3.482179008	69	]		[	42	2.342551866	96	]		[	43	3.271409068	70	]		[	44	3.073886042	72	]		[	45	1.909989271	116	]
    [	46	2.955558142	74	]		[	47	2.105063024	104	]		[	48	2.937414178	77	]		[	49	2.83795786	78	]		[	50	2.489721663	91	]
    [	51	3.248007593	71	]		[	52	3.143735522	71	]		[	53	2.525334751	88	]		[	54	2.866452903	79	]		[	55	2.644319547	86	]
    [	56	2.241790722	100	]		[	57	2.833067799	77	]		[	58	1.805248666	118	]		[	59	2.577458214	87	]		[	60	1.964922024	111	]
    [	61	2.35253024	96	]		[	62	3.865898834	64	]		[	63	1.913539678	115	]		[	64	2.227856539	98	]		[	65	2.923445211	76	]
    [	66	3.389194526	69	]		[	67	2.167449764	104	]		[	68	2.387825497	95	]		[	69	2.721739083	86	]		[	70	3.641653697	67	]
    [	71	2.773508537	79	]		[	72	2.573169672	85	]		[	73	1.926157223	114	]		[	74	3.202884592	72	]		[	75	3.685383305	67	]
    [	76	2.789539353	77	]		[	77	2.515312377	88	]		[	78	2.641421709	86	]		[	79	2.756179284	80	]		[	80	3.259565991	70	]
    [	81	2.517050426	89	]		[	82	3.434097023	69	]		[	83	3.536806498	68	]		[	84	3.122301782	72	]		[	85	2.124820548	106	]
    [	86	2.524938238	90	]		[	87	2.603834559	86	]		[	88	2.195221244	102	]		[	89	2.783977205	77	]		[	90	3.754825843	66	]
    [	91	1.502411823	166	]		[	92	2.232364949	99	]		[	93	2.265763348	98	]		[	94	1.824078589	117	]		[	95	2.881853021	78	]
    [	96	1.417284696	176	]		[	97	2.791898437	79	]		[	98	3.634280145	65	]		[	99	2.644299338	81	]		[	100	2.584165259	86	]
    [	101	2.237672714	100	]		[	102	2.220789354	99	]		[	103	2.735258598	82	]		[	104	2.401270767	93	]		[	105	2.392507023	91	]
    [	106	2.585074627	89	]		[	107	2.951745222	75	]		[	108	3.511121818	68	]		[	109	4.054639581	63	]		[	110	1.928980539	114	]
    [	111	2.754748169	80	]		[	112	2.469844227	91	]		[	113	4.062332673	61	]		[	114	2.898067238	76	]		[	115	2.657707969	85	]
    [	116	1.320411861	184	]		[	117	2.508230966	88	]		[	118	1.626192607	145	]		[	119	2.279911367	97	]		[	120	3.074343568	74	]
    [	121	2.638784726	87	]		[	122	3.353004377	70	]		[	123	2.223401856	101	]		[	124	2.603188154	87	]		[	125	2.267257768	98	]
    [	126	2.983930852	75	]		[	127	1.73892926	128	]		[	128	2.266154873	96	]		[	129	2.495622126	88	]		[	130	2.430715457	90	]
    [	131	3.273911495	71	]		[	132	2.129245821	104	]		[	133	2.907519472	75	]		[	134	2.557461398	88	]		[	135	1.973113443	114	]
    [	136	2.581965838	86	]		[	137	1.821365754	119	]		[	138	3.274205429	70	]		[	139	1.807873209	118	]		[	140	2.118014367	104	]
    [	141	2.615536606	87	]		[	142	3.056669705	74	]		[	143	3.811935076	66	]		[	144	2.721100063	79	]		[	145	1.662575237	138	]
    [	146	2.051446689	107	]		[	147	2.277404042	98	]		[	148	2.217994544	101	]		[	149	3.672599483	68	]		[	150	3.124256405	73	]
    [	151	2.361836307	95	]		[	152	3.33156975	69	]		[	153	3.705112361	68	]		[	154	4.145382921	62	]		[	155	3.517971555	69	]
    [	156	2.935184159	76	]		[	157	3.215770782	71	]		[	158	1.713655307	130	]		[	159	2.447642405	90	]		[	160	3.229922709	71	]
    [	161	2.515623885	89	]		[	162	3.10328301	72	]		[	163	3.369963969	69	]		[	164	3.33911909	70	]		[	165	2.855339475	77	]
    [	166	2.713280183	82	]		[	167	2.133034363	103	]		[	168	1.7272594	127	]		[	169	2.401783899	91	]		[	170	3.422691653	70	]
    [	171	3.21521907	72	]		[	172	3.017227294	74	]		[	173	2.829045362	79	]		[	174	2.373970457	93	]		[	175	3.063607151	75	]
    [	176	3.546606095	68	]		[	177	2.219268232	100	]		[	178	2.453929361	90	]		[	179	2.655052756	82	]		[	180	2.778664156	79	]
    [	181	2.207988586	100	]		[	182	3.525786281	68	]		[	183	1.991189693	112	]		[	184	2.442106572	90	]		[	185	2.538184964	87	]
    [	186	2.607408608	86	]		[	187	1.626256624	145	]		[	188	2.743506448	83	]		[	189	3.22788306	70	]		[	190	1.869097927	118	]
    [	191	1.986994467	111	]		[	192	3.587960782	68	]		[	193	2.592383937	85	]		[	194	2.581159843	86	]		[	195	1.500931348	163	]
    [	196	2.627195047	85	]		[	197	1.990289946	112	]		[	198	3.624266659	67	]		[	199	2.702071976	80	]		[	200	2.33547892	93	]
  )
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report force-of-infection [given-age-group-id given-AGI given-AGN]
  if Simulating-Weekends? [
    ifelse not weekend? [
      ifelse school-holiday? [
        report force-of-infection-weekday-no-school given-age-group-id  given-AGI given-AGN
      ]
      [
        report force-of-infection-weekday-school given-age-group-id  given-AGI given-AGN
      ]
    ]
    [
      report force-of-infection-weekend given-age-group-id given-AGI given-AGN
    ]
  ]

  ; Not recognizing weekends as different.
  if school-holiday? [
    report force-of-infection-weekday-no-school given-age-group-id given-AGI given-AGN
  ]
  report force-of-infection-weekday-school given-age-group-id given-AGI given-AGN

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report force-of-infection-weekday-school [given-age-group-id given-AGI given-AGN]
  report
   susceptibility *
  (sum (map [[c i n] -> ifelse-value (n = 0) [0] [c * i / n]] (matrix:get-row contacts-matrix-weekday-school given-age-group-id) given-AGI given-AGN))
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report force-of-infection-weekday-no-school [given-age-group-id given-AGI given-AGN]
  report
   susceptibility *
  (sum (map [[c i n] -> ifelse-value (n = 0) [0] [c * i / n]] (matrix:get-row contacts-matrix-weekday-no-school given-age-group-id) given-AGI given-AGN))
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report force-of-infection-weekend [given-age-group-id given-AGI given-AGN]
  report
   susceptibility *
  (sum (map [[c i n] -> ifelse-value (n = 0) [0] [c * i / n]] (matrix:get-row contacts-matrix-weekend given-age-group-id) given-AGI given-AGN))
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Decisions about locations

to schedule-initial-decisions
  ask people [
    let chosen-option option-with-max-penalty
    schedule 0 self ([-> do-decision-in p-location]) p-location "Initial Deciding."
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-personal-schedules
  ask people [
    if any? my-out-staffings [
      set p-options schedule-of-ordinary-adult (min-one-of out-habitant-neighbors [who]) (min-one-of out-staffing-neighbors [who]) (min-one-of hospitals [distance myself])
    ]
    if any? my-out-studyings [
      set p-options schedule-of-ordinary-child (min-one-of out-habitant-neighbors [who]) (min-one-of out-studying-neighbors [who]) (min-one-of hospitals [distance myself])
    ]
    if p-options = 0 [
      ifelse p-age < 16 [
        set p-options schedule-of-ordinary-child (min-one-of out-habitant-neighbors [who]) (min-one-of out-habitant-neighbors [who]) (min-one-of hospitals [distance myself])
      ]
      [
        set p-options schedule-of-ordinary-adult (min-one-of out-habitant-neighbors [who]) (min-one-of out-habitant-neighbors [who]) (min-one-of hospitals [distance myself])
      ]
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report schedule-of-ordinary-adult [given-household given-workplace given-hospital]
  report (list
    ; location condition start-time end-time penalty
    ; Between start-time and end-time, if I am not at location, and condition is true, then I incur the penalty.
    ; Choose location that minimizes penalty.
    (list given-household [-> true] 0 (24 * 60) 1) ; If nowhere better to go, go home.
    (list given-household [-> true] (23 * 60) ((24 + 7) * 60) 10) ; Sleep at home until next day.
    (list given-workplace [-> Go-To-Work? and not weekend?] (8 * 60 + 30) (17 * 60 + 30) 100) ; At work.
    (list given-household [-> ifelse-value symptoms? [Stay-At-Home-With-Symptoms? > random-float 100] [false]] 0 (24 * 60) 200) ; Stay at home, you're ill.
    (list [-> friends-house] [-> Go-To-Friend? and weekend?] (19 * 60 + 30) (23 * 60 + 30) 100) ; Visiting friend.
    (list given-hospital [-> hospitalize?] 0 (24 * 60) 1000) ; Go to hospital!
    (list given-hospital [-> dead?] 0 (24 * 60) 10000) ; Go to hospital (morgue)!
    )
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report schedule-of-ordinary-child [given-household given-school given-hospital]
  report (list
    ; location condition start-time end-time penalty
    ; Between start-time and end-time, if I am not at location, and condition is true, then I incur the penalty.
    ; Choose location that minimizes penalty.
    (list given-household [-> true] 0 (24 * 60) 1) ; If nowhere better to go, go home.
    (list given-household [-> true] (23 * 60) ((24 + 7) * 60) 10) ; Sleep at home until next day.
    (list given-school [-> Go-To-School? and not weekend?] (8 * 60 + 30) (15 * 60 + 30) 100) ; At school.
    (list given-household [-> ifelse-value symptoms? [Stay-At-Home-With-Symptoms? > random-float 100] [false]] 0 (24 * 60) 200) ; Stay at home, you're ill.
    (list [-> friends-house] [-> Go-To-Friend? and weekend?] (10 * 60 + 00) (12 * 60 + 30) 100) ; Visiting friend.
    (list [-> friends-house] [-> Go-To-Friend? and weekend?] (14 * 60 + 30) (17 * 60 + 30) 100) ; Visiting friend.
    (list given-hospital [-> hospitalize?] 0 (24 * 60) 1000) ; Go to hospital!
    (list given-hospital [-> dead?] 0 (24 * 60) 10000) ; Go to hospital (morgue)!
    )
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report friends-house
  ; Run by person.
  ; NB: The friend might not be at home!
  if symptoms? [report p-location] ; You feel ill. Stay home.
  if hospitalize? or dead? [report p-location] ; Forget it!

  if not any? friendship-neighbors with [not dead? and not hospitalize? and not symptoms?] [
    report p-location ; All your friends are sick. Stay where you are. (Presumably at home.)
  ]

  let destinations [p-location] of friendship-neighbors with [not dead? and not hospitalize? and not symptoms?]
  set destinations filter [loc -> [not any? people-here with [symptoms? or hospitalize?]] of loc] destinations
  if empty? destinations [report p-location] ; Nowhere to go. Stay put.
  report one-of destinations ; Found somewhere to go.
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to do-decision-in [given-location]
  ; Run by person.
  ; Option has form (list location condition start-time end-time penalty).
  ; Between start-time and end-time, if I am not at location, and condition is true, then I incur the penalty.
  ; Choose location that minimizes penalty.
  if given-location != p-location [stop] ; Decision obsolete, because I'm not where I was expecting to be.
  let chosen-option option-with-max-penalty
  if chosen-option = [] [stop]
  let next-loc item 0 chosen-option
  if is-anonymous-reporter? next-loc [
    set next-loc runresult next-loc
  ]
  relocate-to next-loc
;  schedule sim-time self ([-> relocate-to next-loc]) next-loc "Relocating."
  schedule min (list
    (last-midnight + item 3 chosen-option)
    next-option-start-time
    )
    self ([-> do-decision-in next-loc]) next-loc "Deciding."

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report next-option-start-time
  ; Work out what the next option would be (if its condition holds),
  ; and report its start time.
  ; Option has form (list location condition start-time end-time penalty).
  let earliest-options filter [cur-option -> sim-time-of-day < item 2 cur-option] p-options
  if empty? earliest-options [
    report last-midnight + (24 * 60) + min map [cur-option -> item 2 cur-option] p-options
  ]
  report last-midnight + min map [cur-option -> item 2 cur-option] earliest-options
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report option-with-max-penalty
  ; The option that is valid and currently incurs the largest penalty
  ; is the best option to carry out.
  let options []
  let cur-time-of-day sim-time-of-day
  foreach p-options [cur-option ->
;    if p-location != item 0 cur-option [
      if runresult item 1 cur-option [
        if cur-time-of-day >= item 2 cur-option [
          if cur-time-of-day < item 3 cur-option [
            set options fput cur-option options
          ]
        ]
      ]
;    ]
  ]
  if empty? options [
    show (word "WARNING: I've got nothing to do at time " sim-time "!")
    report []
  ]
  let max-penalty max map [x -> last x] options
  report one-of filter [x -> max-penalty = last x] options
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to goto-hospital
  let old-location p-location
  set p-location one-of hospitals
  move-to p-location
  if trace-turtles? [
    if member? self trace-turtles [print-trace-turtles (word "Went to hospital from " ([item l-type location-type-names] of old-location) " " old-location)]
    if member? old-location trace-turtles [print-trace-turtles (word self " left for " ([item l-type location-type-names] of p-location) " " p-location )]
    if member? p-location trace-turtles [print-trace-turtles (word self " arrived from " ([item l-type location-type-names] of old-location) " " old-location )]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to goto-home
  let old-location p-location
  set p-location min-one-of out-habitant-neighbors [who]
  move-to p-location
  if trace-turtles? [
    if member? self trace-turtles [print-trace-turtles (word "Went home from hospital " ([item l-type location-type-names] of old-location) " " old-location)]
    if member? old-location trace-turtles [print-trace-turtles (word self " left for " ([item l-type location-type-names] of p-location) " " p-location )]
    if member? p-location trace-turtles [print-trace-turtles (word self " arrived from " ([item l-type location-type-names] of old-location) " " old-location )]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to relocate-to [given-location]
  if given-location = p-location [schedule-next-random-event stop] ; No change.

  let old-location p-location
  move-to given-location
  set p-location given-location
  recalc-rates-after-relocation-from old-location
  set p-time-of-last-change sim-time

  if trace-turtles? [
    if member? self trace-turtles [print-trace-turtles (word "Relocated from " ([item l-type location-type-names] of old-location) " " old-location)]
    if member? old-location trace-turtles [print-trace-turtles (word self " left for " ([item l-type location-type-names] of p-location) " " p-location )]
    if member? p-location trace-turtles [print-trace-turtles (word self " arrived from " ([item l-type location-type-names] of old-location) " " old-location )]
  ]

  ; Anything else needed here?
  schedule-next-random-event
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Process infections and disease-state transitions.
;; Calculates rates of occurrence for random events.
;;
; DOTS: Duration infectious Opportunity Transmission probability Susceptibility
; Wuhan: 15 contacts per day => R = 2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report random-event-time-expon [inv-lambda]
  ; Where mean=1/lambda, in days.
  ; Returns in minutes.
  report round (1440 * random-exponential inv-lambda)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report random-event-time-gamma [mu k]
  ; Where mean=mu, shape=k, in days.
  ; Returns in minutes.
  report round (1440 * random-gamma k (k / mu))
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to schedule-next-random-event
  ; It is possible to do the whole thing with just two random numbers.
  ; Maybe later...

  if (total-infection-rate) = 0 [stop] ; Nothing happens.

  ; Sample time of the next random event.
  let sampled-time random-event-time
  if sampled-time >= b-events-next-time [stop] ; World might change at the next event, whereupon we will need to resample.

  ; Ok. This will be the next event. Now sample what type of event it will be, where it will be, who it will involve...
  ; Infection
  let chosen-loc-type last rnd:weighted-one-of-list (map [[r t] -> (list r t)] (array:to-list loc-type-infection-rate) location-types) [p -> first p]
  let chosen-location rnd:weighted-one-of chosen-loc-type [l-infection-rate]
  ; NB: Using people-here, so mustn't have multiple locations on same patch.
  let chosen-susceptible [rnd:weighted-one-of-list (array:item l-contents 0) [p -> [p-susceptibility * p-contacts] of p]] of chosen-location
;  let chosen-susceptible [rnd:weighted-one-of (people-here with [susceptible?]) [p-susceptibility * p-contacts]] of chosen-location
  if chosen-susceptible = nobody [ask chosen-location [show-error (word "ERROR! Managed to sample nobody to become infected.") ]]
  ; NB: Using people-here, so mustn't have multiple locations on same patch.
  let chosen-infector [rnd:weighted-one-of-list (sentence (map [ds -> array:item l-contents [who] of ds] filter [ds -> [ds-rel-inf > 0] of ds] compartments-sorted)) [p -> [p-infectiousness * p-contacts] of p]] of chosen-location
;  let chosen-infector [rnd:weighted-one-of (people-here with [infector?]) [p-infectiousness * p-contacts]] of chosen-location
  schedule sampled-time chosen-susceptible [-> become-infected-by chosen-infector chosen-location] chosen-location "Infection."

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report ds-transition-rate
  if not any? my-out-transitions [report 0]
  report sum [(runresult ts-perc-weight) * ts-rate] of my-out-transitions
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report rate-conversion
  ; MUST BE REVIEWED! All percentages now ~[0, 1]

  ; Infection rate (in events per day) has been adjusted by:
  ; 100 for p-susceptibility (in %)
  ; 100 for p-infectiousness (in %)
  ; 24 for p-contacts (per hour) of susceptibles
  ; 24 for p-contacts of infectors
  ; 1/24 for p-contacts of people-here
  ; 10 for infection-bonus
  ; =100*100*10 = 100000. So /100000 to get back to events per day.

  ; dst rate (in events per day) has been adjusted by:
  ; 1 for ts-perc-weight (no longer in %)
  ; 100 for ts-rate to be in events per 100 day.
  ; So /(100 * 1) to get back to events per day.
  ; So to bring dst-rate into line with infection-rate, multiple l-dst-rate by
  report 1000 ; =100000 / 100*1
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report random-event-time
  ; Need time in minutes.
  ; 1440 = 24 * 60
  ; Given rate in events per 100000 days.
  ; So multiply by 1 / (24 * 60 * 100000) to get events per minute.

;  print (word sim-time ": IR=" total-infection-rate ", DST=" total-dst-rate)
  report sim-time + round random-exponential (144000000 /  total-infection-rate )
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-locations-rates
  ask locations [
;  ask locations with [any? people with [p-location = myself]] [ ; Misses the point.
    set l-sc-rate sum map [p -> [p-susceptibility * p-contacts] of p] array:item l-contents 0
    set l-ic-rate sum map [p -> [p-infectiousness * p-contacts] of p] (sentence map [ds -> array:item l-contents [who] of ds] filter [ds -> [ds-rel-inf > 0] of ds] compartments-sorted)
;    set l-ic-rate sum [p-infectiousness * p-contacts] of people-here with [infector?]
    set l-nc-rate sum map [p -> [p-contacts] of p] (sentence map [ds -> array:item l-contents [who] of ds] filter [ds -> ds != dead] compartments-sorted)
;    set l-nc-rate sum [p-contacts] of people-here with [not dead?]
    calc-infection-rate-here
  ]

  set loc-type-infection-rate array:from-list n-values (length location-types) [-> 0]

  (foreach location-types (n-values (length location-types) [? -> ?]) [[cur-type cur-pos] ->
    array:set loc-type-infection-rate cur-pos sum [l-infection-rate] of cur-type
  ])
  set total-infection-rate sum array:to-list loc-type-infection-rate
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to check-locations-rates
  print ""
  print "Performing check on location rates."
  print ""
  foreach sort locations [? -> ask ? [
    ; NB: Using people-here, so mustn't have multiple locations on same patch.
    if l-sc-rate != sum [p-susceptibility * p-contacts] of people-here with [susceptible?] [show (word "WARNING: Failed SC-Rate: " l-sc-rate)]
    if l-ic-rate != sum [p-infectiousness * p-contacts] of people-here with [infector?] [show (word "WARNING: Failed IC-Rate: " l-ic-rate)]
    if l-nc-rate != sum [p-contacts] of people-here with [not dead?] [show (word "WARNING: Failed NC-Rate: " l-nc-rate)]
    calc-infection-rate-here
  ]]

  (foreach location-types (n-values (length location-types) [? -> ?]) [[cur-type cur-pos] ->
    if array:item loc-type-infection-rate cur-pos != sum [l-infection-rate] of cur-type [print (word "WARNING: Failed loc-type-infection-rate item " cur-pos ": " array:to-list loc-type-infection-rate)]
  ])
  if total-infection-rate != sum array:to-list loc-type-infection-rate  [print (word "WARNING: Failed total-infection-rate" total-infection-rate)]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to recalc-rates-after-transition-from [old-state]
  ; Run by person, after they have updated p-state.
  let cur-susc-rate p-susceptibility * p-contacts
  let cur-inf-rate p-infectiousness * p-contacts

  ask p-location [
    ; Remove old weights
    array:set loc-type-infection-rate l-type (array:item loc-type-infection-rate l-type) - l-infection-rate
    set total-infection-rate total-infection-rate - l-infection-rate
    ; Add new weights
    ifelse old-state = susceptible [
      set l-sc-rate l-sc-rate - cur-susc-rate
    ]
    [ ; Old-state not susceptible. One of the states with random transitions out.
    ]
    if [susceptible?] of myself [set l-sc-rate l-sc-rate + cur-susc-rate]
    if [ds-allows-infection?] of old-state [set l-ic-rate l-ic-rate - cur-inf-rate]
    if [infector?] of myself [set l-ic-rate l-ic-rate + cur-inf-rate]
    if [dead?] of myself [set l-nc-rate l-nc-rate - [p-contacts] of myself]

    calc-infection-rate-here

    ; Update weights for location types.
    array:set loc-type-infection-rate l-type (array:item loc-type-infection-rate l-type) + l-infection-rate
    set total-infection-rate total-infection-rate + l-infection-rate
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to recalc-rates-after-relocation-from [old-location]
  recalc-rates-at old-location -1
  recalc-rates-at p-location 1
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to recalc-rates-at [given-location adjustment]
  ; Run by person, after they have updated p-state.
  ; Will either remove person's rates from those at given-location (adjustment = -1)
  ; or add them to those at given-location (adjustment = 1).
  let cur-susc-rate p-susceptibility * p-contacts
  let cur-inf-rate p-infectiousness * p-contacts

  ask given-location [
    ; Remove old rates
    array:set loc-type-infection-rate l-type (array:item loc-type-infection-rate l-type) - l-infection-rate
    set total-infection-rate total-infection-rate - l-infection-rate

;    ; Debugging purposes.
;    if self = location 273 [
;      show (word myself ": " ([ds-name] of [p-state] of myself) ", sc=" l-sc-rate ", ic=" l-ic-rate ", nc=" l-nc-rate ", IR=" l-infection-rate)
;    ]

    if [susceptible?] of myself [set l-sc-rate l-sc-rate + adjustment * cur-susc-rate]
    if [infector?] of myself [set l-ic-rate l-ic-rate + adjustment * cur-inf-rate]
    if [not dead?] of myself [set l-nc-rate l-nc-rate + adjustment * [p-contacts] of myself]
    if [not susceptible? and not dead?] of myself [set l-dst-rate l-dst-rate + adjustment * rate-conversion * ([[ds-transition-rate] of p-state] of myself) ]

    calc-infection-rate-here

;    if self = location 273 [
;      show (word myself ": " ([ds-name] of [p-state] of myself) ", sc=" l-sc-rate ", ic=" l-ic-rate ", nc=" l-nc-rate ", IR=" l-infection-rate )
;    ]

    ; Add recomputed rates
    array:set loc-type-infection-rate l-type (array:item loc-type-infection-rate l-type) + l-infection-rate
    set total-infection-rate total-infection-rate + l-infection-rate
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to calc-infection-rate-here
  ; Recompute infection rate at this location.
  ; SI model: dI/dt = contact-rate * susceptibility * infectiousness * S * I / N
  ; Heterogeneous agents:
  ; lambda =
  ; location_boost
  ; * sum_all_S(My-Susceptibility * My-Contact-Rate)
  ; * sum_all_I(My-Infectiousness * My-Contact-Rate)
  ; / sum_all_people-here(My-Contact-Rate)
  ifelse l-nc-rate = 0 [
    set l-infection-rate 0
  ]
  [
;    set l-infection-rate l-infection-boost * l-sc-rate * l-ic-rate / l-nc-rate
    set l-infection-rate round (l-infection-boost * l-sc-rate * l-ic-rate / l-nc-rate)
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to become-infected-by [given-infector given-location]
  ; Run as person. Should be susceptible.
  ; Update infection and ds transition rates for person, their location, and the location's type.
  if not susceptible? [user-message-error (word "ERROR! Someone not susceptible is trying to become infected.")]
  if given-location != p-location [show-error "WARNING! Susceptible has moved before infection event." stop]

  if is-person? given-infector [
    if given-location != [p-location] of given-infector [show-error "WARNING! Infector has moved before infection event." stop]
    if [dead?] of given-infector [show-error "WARNING! Infector is dead before infection event." stop]
    ask given-infector [
      set p-num-infected p-num-infected + 1
      if p-num-infected > max-infected-by-one [set max-infected-by-one p-num-infected]
      add-to-superspreaders?
    ]
  ]

  if trace-turtles? [
    if member? self trace-turtles [print-trace-turtles (word "Infected by " given-infector)]
    if member? given-infector trace-turtles [print-trace-turtles (word "Infected  " self)]
    if member? p-location trace-turtles [print-trace-turtles (word self " was infected here by " given-infector)]
  ]

  become-infected

  ask p-location [
    set l-num-infections-here l-num-infections-here + 1
    array:set loc-type-num-infections-here l-type 1 + array:item loc-type-num-infections-here l-type
    relabel-venue
  ]

;  recalc-rates-after-transition-from susceptible

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to become-infected
  let route-out one-of [my-out-transitions] of p-state ; Presumably to Exposed.

  do-ds-transition route-out
  set p-infection-time sim-time
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to become-dead
  goto-hospital ; For want of a better place to send them. No morgue in this city?
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report infinity
  ; Arbitrarily large number.
  ; Could represent time of events that are far in the future,
  ; beyond the end of simulated time.
  report (2 ^ 31) - 1
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to do-ds-transition [given-transition]
  ; Run by person.
  set cur-person self

  array:set (array:item age-groups p-age-group) ([who] of p-state) (remove self array:item (array:item age-groups p-age-group) ([who] of p-state))

  ; Count the occurrence of this type of ds transition.
  ask given-transition [
    set ts-num-events ts-num-events + 1
    set ts-num-new-events fput (1 + first ts-num-new-events) but-first ts-num-new-events
    relabel-transition
    ask end1 [adjust-cur-num -1] ; Old state has 1 fewer people in it.
    ask end2 [adjust-cur-num 1 ] ; New state has 1 extra person in it.
  ]
  set p-state [end2] of given-transition
  array:set (array:item age-groups p-age-group) ([who] of p-state) (fput self array:item (array:item age-groups p-age-group) ([who] of p-state))
  set p-time-of-last-change sim-time
  set p-ds-history fput (list sim-time p-state) p-ds-history ; Each person knows their ds history.

  if social-interactions != "Contacts-Matrices" [
    recalc-rates-after-transition-from [end1] of given-transition ; Location recomputes its rates.
  ]

  recolor-by-state
  if trace-turtles? [
    if member? self trace-turtles [print-trace-turtles (word "Transitioned from " ([ds-name] of [end1] of given-transition))]
    if member? p-location trace-turtles [print-trace-turtles (word self " transitioned from " ([ds-name] of [end1] of given-transition) " to " ([ds-name] of p-state))]
  ]
  run [ts-event-commands] of given-transition ; Run any special code for this transition.

  if [any? my-out-transitions] of p-state [ ; If it's not a terminating state,
    let sampled-transition [rnd:weighted-one-of my-out-transitions [runresult ts-perc-weight]] of p-state
    let sampled-time sim-time + [runresult ts-ie-time] of sampled-transition
    schedule sampled-time self [-> do-ds-transition sampled-transition] sampled-transition "DS Transition."
  ]
  schedule-next-random-event ; If next event is a random one, schedule it.
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Recording changes visually

to relabel-transition
  if Transition-Labels ="Count" [set label (word ts-num-events "     ") stop]
  if Transition-Labels ="%" [
;    let denom sum [[ts-num-events] of my-out-transitions] of end1
;    ifelse denom = 0 [
    ifelse not any? people [
      set label ""
    ]
    [
      set label (word (round (as-perc ts-num-events)) "%     ")
    ]
    stop
  ]
  if Transition-Labels ="Time-Param" [set label ifelse-value (is-number? ts-time-param) [(word "     " ts-time-param)] [""] stop]
  if Transition-Labels ="Rate" [set label (word "     " precision ts-rate 3) stop]
  if Transition-Labels ="Probability" [
    ifelse ts-rate = 0 [
      set label (word "     " "0%")
    ]
    [
      set label (word "     " round (100 * ts-rate / [ds-transition-rate] of end1) "%")
    ]
    stop
  ]
  set label ""
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to relabel-all-transitions
  foreach sort transitions [? ->
    ask ? [
      relabel-transition
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to relabel-all-disease-states
  ; NB: foreach sort: Avoid ask on multiple turtles, because it uses random numbers.
  foreach sort disease-states [? ->
    ask ? [
      relabel-disease-state
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to relabel-disease-state
  if disease-state-labels = "" [set label "" stop]
  if disease-state-labels = "Name" [set label (word ds-name "      ") stop]
  if disease-state-labels = "Name-And-Current" [set label (word ds-name ": " ds-cur-num "       ") stop]
  if disease-state-labels = "Name-And-Peak" [set label (word ds-name ": " ds-max-num "      ") stop]
  if disease-state-labels = "Current" [set label (word ds-cur-num "     ") stop]
  if disease-state-labels = "Current-And-Peak" [set label (word ds-cur-num " : " ds-max-num " on Day " (int time-as-days ds-max-at-time)) stop]
  set label ""
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Plots and stats

to calc-stats

  if locked-down? [set num-days-in-lockdown num-days-in-lockdown + 1] ; Careful! Must run consistently before lockdown switched on/off.

  if Calculate-Susceptibles-Net? [calc-susceptibles-network] ; Output. Not clear if this gives useful information. Comment out to save computer time?
  set deaths fput ([ds-cur-num] of dead) deaths

  ask transition-between Susceptible Exposed [
    if max-new-exposed < first ts-num-new-events [set max-new-exposed first ts-num-new-events]
  ]
  ask new-case-transition [
    if max-new-cases < first ts-num-new-events [set max-new-cases first ts-num-new-events]
  ]

  let cur-max max daily-new-deaths
  if max-new-deaths < cur-max [set max-new-deaths cur-max]


  ; Stats reported in the Lancet paper from the LSHTM model.

  ; Any data on "cases" is likely to miss many infected people.
  ; Some are missed because asymptomatic.
  ; Some are missed because they are not tested anywhere.
  ; Some are missed because they do not report at a hospital.
  ; Some are missed because the tests make false-negative errors (while making false positves too!).
  ; So to simulate the generation of "cases", we have some choices.

  ask new-case-transition [
    set num-new-cases first ts-num-new-events
    set total-cases ts-num-events
    set num-peak-week-cases max weekly-from-daily ts-num-new-events
    set num-weeks-to-peak-week-cases position-of-max weekly-from-daily ts-num-new-events
  ]

  set num-new-infected [first ts-num-new-events] of transition-between Susceptible Exposed
  set num-new-icu-beds [first ts-num-new-events] of transition-between Hospitalize h-icu
  if max-new-hospitalize < num-new-icu-beds + [first ts-num-new-events] of transition-between Hospitalize h-non-icu [
    set max-new-hospitalize num-new-icu-beds + [first ts-num-new-events] of transition-between Hospitalize h-non-icu
  ]

  set num-new-deaths first daily-new-deaths

  set total-deaths [ds-cur-num] of dead
  set num-peak-week-deaths max weekly-from-daily daily-new-deaths
  set num-peak-icu-beds [ds-max-num] of h-icu
  set num-peak-non-icu-beds [ds-max-num] of h-non-icu

  ; Number actually infected, not number showing symptoms, tested postive, or otherwise counted as cases.
  set total-infected [ts-num-events] of transition-between Susceptible Exposed
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report daily-new-deaths
  report reduce [[y z] -> (map [[a b] -> a + b] y z)] [ts-num-new-events] of [my-in-transitions] of dead
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report weekly-from-daily [given-list]
  ; Given a list of daily data (latest first)
  ; report a list of weekly data (earliest first).

  let num-days length given-list
  report n-values ceiling (num-days / 7) [w ->
    sum sublist (reverse given-list) (w * 7) min (list (length given-list) ((w + 1) * 7))
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report position-of-max [given-list]
  ; Given a list
  ; reports the first position with the largest value.
  let best-val max given-list
  report position best-val given-list
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-initial-plots
  set-current-plot "People By Disease State"
  setup-extra-plot-pens
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-extra-plot-pens
  clear-plot
  let pen-name "Susceptible"
;  set-current-plot-pen pen-name
;  set-plot-pen-color [c-color] of susceptible
  foreach compartments-sorted [cur-comp ->
    set pen-name [ds-name] of cur-comp
    create-temporary-plot-pen pen-name
    set-current-plot-pen pen-name
    set-plot-pen-color [ds-color] of cur-comp
    set-plot-pen-mode 0
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to do-compartments-plot
  let pen-name ""
  foreach compartments-sorted [cur-comp ->
    set pen-name [ds-name] of cur-comp
    set-current-plot-pen pen-name
    plotxy sim-days [as-perc ds-cur-num] of cur-comp
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to do-time-series-plots
  set-current-plot "People By Disease State"
  do-compartments-plot

  set-current-plot "Changes"
  let cur-num 0
  set-current-plot-pen "New Exposed"
  set cur-num [first ts-num-new-events] of transition-between Susceptible Exposed
  plotxy sim-days cur-num
;  if cur-num > 0 [plotxy sim-days log cur-num 10]
  set-current-plot-pen "New Clinical Cases"
  set cur-num [first ts-num-new-events] of transition-between I-Preclinical I-Clinical
  plotxy sim-days cur-num
;  if cur-num > 0 [plotxy sim-days log cur-num 10]
  set-current-plot-pen "New Hospitalize"
  set cur-num [first ts-num-new-events] of transition-between I-Clinical Hospitalize
  plotxy sim-days cur-num
;  if cur-num > 0 [plotxy sim-days log cur-num 10]
  set-current-plot-pen "New Dead"
  set cur-num [first ts-num-new-events] of transition-between H-ICU Dead + [first ts-num-new-events] of transition-between H-Non-ICU Dead
  plotxy sim-days cur-num
;  if cur-num > 0 [plotxy sim-days log cur-num 10]
  set-current-plot-pen "New Recovered"
  set cur-num [first ts-num-new-events] of transition-between H-ICU Recovered + [first ts-num-new-events] of transition-between H-Non-ICU Recovered + [first ts-num-new-events] of transition-between I-Subclinical Recovered ;+ [first ts-num-new-events] of transition-between I-Clinical Recovered
  plotxy sim-days cur-num
;  if cur-num > 0 [plotxy sim-days log cur-num 10]

  set-current-plot "Cumulative"
  set-current-plot-pen "Exposed"
  set cur-num [ts-num-events] of transition-between Susceptible Exposed
  plotxy sim-days cur-num
;  if cur-num > 0 [plotxy sim-days log cur-num 10]
  set-current-plot-pen "Hospitalize?"
  set cur-num [ts-num-events] of transition-between I-Clinical Hospitalize
  plotxy sim-days cur-num
;  if cur-num > 0 [plotxy sim-days log cur-num 10]
  set-current-plot-pen "Dead"
  set cur-num [ts-num-events] of transition-between H-ICU Dead + [ts-num-events] of transition-between H-Non-ICU Dead
  plotxy sim-days cur-num
;  if cur-num > 0 [plotxy sim-days log cur-num 10]
  set-current-plot-pen "Recovered"
  set cur-num [ts-num-events] of transition-between H-ICU Recovered + [ts-num-events] of transition-between H-Non-ICU Recovered + [ts-num-events] of transition-between I-Subclinical Recovered ;+ [ts-num-events] of transition-between I-Clinical Recovered
  plotxy sim-days cur-num
;  if cur-num > 0 [plotxy sim-days log cur-num 10]

  set-current-plot "Estimated R"
  set-current-plot-pen "R"
  plotxy sim-days estimated-R
  set-current-plot-pen "R = 1"
  plotxy sim-days 1

  if Calculate-Susceptibles-Net? [
    set-current-plot "Largest Component Size / Susceptibles"
    if num-susceptible > 0 [
      plotxy sim-days (largest-sus-net-component / num-susceptible)
    ]

    set-current-plot "Sus. Neighbours of Infectious"
    ifelse 0 < num-infectious [
      let list-of-num-susc-nghbrs [count susceptible-neighbours] of infectors
      set-current-plot-pen "Max"
      plotxy sim-days max list-of-num-susc-nghbrs
      set-current-plot-pen "Mean"
      plotxy sim-days mean list-of-num-susc-nghbrs
      set-current-plot-pen "Min"
      plotxy sim-days min list-of-num-susc-nghbrs
    ]
    [
      set-current-plot-pen "Max"
      plotxy sim-days 0
      set-current-plot-pen "Mean"
      plotxy sim-days 0
      set-current-plot-pen "Min"
      plotxy sim-days 0
    ]

    set-current-plot "Susceptibles Net Components"
    plotxy sim-days (num-sus-net-components)
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to do-histo-plot
  set-current-plot "People Histogram"
  clear-plot
  let histo-list runresult histo-x-axis
  if empty? histo-list [clear-plot stop]
  set-plot-pen-interval 1
  set-plot-x-range 0 (1 + int max histo-list)
  histogram histo-list

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report pop-rescalar
  report (count people) / 66400000
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report as-perc [given-count]
  ; Report numbers as percentages of the number of people.
  ; Useful for comparing outputs between sim runs with different population sizes.
  report 100 * given-count / count people
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report as-prop [given-count]
  ; Report numbers as proportions of the number of people.
  ; Useful for comparing outputs between sim runs with different population sizes.
  report given-count / count people
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report estimated-R
  ; Estimate made from last 5 days.
  ; Should replace this with correct calculation of Effective Reproduction Number (R).
  if social-interactions = "Activity-Locations" [
    if any? post-disease-people with [(-1 - 5 + int time-as-days sim-time) < (int time-as-days first first p-ds-history)] [
      report mean [p-num-infected] of post-disease-people with [(-1 - 5 + int time-as-days sim-time) < (int time-as-days first first p-ds-history)]
    ]
    report r0
  ]
  report r0
  ; R = Duration * Opportunity * Transmission-Chance * Susceptibles (DOTS)
  ; report weighted-time-infectious * mean_over_time_infectious(sum_age-groups(contacts from ag * susceptibles in ag)) * (susceptibility-given-r0 r0)
  ; report DOTS= R0 * S = R
  ; So for R=1, need herd immunity S = 1/R0.
  ; Analogy: Inverse of force of infection.
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Compute the Susceptibles-Susceptibles social network
;; Did epidemic end because we ran out of infectious people,
;; or because the Susceptibles-Network fragmented?
; Only for outputs.
; Not very useful? Therefore delete?

to calc-susceptibles-network
  ; Q. Does the epidemic halt because the network breaks up?
  let num-components 0
  let cur-component-sze 0
  let largest-component-size 0
  let stack []
  let cur-node nobody

  ask people with [susceptible?] [set p-component 0]
  foreach sort people with [susceptible?] [? ->

    ask ? [
      set stack fput self stack
      if p-component = 0 [
        if largest-component-size < cur-component-sze [
          set largest-component-size cur-component-sze
        ]
        set cur-component-sze 0
        set num-components num-components + 1
      ]
    ]
    while [not empty? stack] [
      set cur-node first stack
      set stack but-first stack
      ask cur-node [
        if p-component = 0 [
          set p-component num-components
          set cur-component-sze cur-component-sze + 1
          ask susceptible-neighbours with [p-component = 0] [
            set stack fput self stack
          ]
        ]
      ]
    ]
  ]
  if largest-component-size < cur-component-sze [
    set largest-component-size cur-component-sze
  ]

  set num-sus-net-components num-components
  set largest-sus-net-component largest-component-size
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report susceptible-neighbours
  ; Note the u in neighbours. To avoid confusing Netlogo with neighbors.
  ; To edit: Should exclude people currently in hospital.
  let return-set (turtle-set )
  if No-Isolation-In-Household? [
    ask out-habitant-neighbors [
      set return-set (turtle-set return-set (in-habitant-neighbors with [susceptible?]))
    ]
  ]
  if not p-off-school? [
    ask out-studying-neighbors [
      set return-set (turtle-set return-set (in-studying-neighbors with [susceptible? and not p-off-school?]))
    ]
  ]
  if not p-off-work? [
    ask out-staffing-neighbors [
      set return-set (turtle-set return-set (in-staffing-neighbors with [susceptible? and not p-off-work?]))
    ]
  ]
  if Go-To-Friend? [
    set return-set (turtle-set return-set (friendship-neighbors with [susceptible?]))
  ]
  report return-set with [self != myself]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
@#$#@#$#@
GRAPHICS-WINDOW
460
10
858
409
-1
-1
1.95
1
12
1
1
1
0
0
0
1
0
199
0
199
1
1
1
ticks
30.0

BUTTON
345
75
445
108
NIL
Setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
225
15
452
48
Perc-Patches-With-Household
Perc-Patches-With-Household
0
100
50.0
5
1
%
HORIZONTAL

MONITOR
185
735
335
780
NIL
Count People
17
1
11

MONITOR
10
1320
202
1365
Friendships Network Density (%)
100 * 2 * (count friendships) / (count people * ((count people) - 1))
2
1
11

MONITOR
225
1320
332
1365
Mean Friendships
2 * (count friendships) / (count people)
1
1
11

BUTTON
1230
985
1452
1018
Show Seed Infected Person's Links
ask seed-infected [ask my-links [set hidden? false]]\nask seed-infected [ask my-out-links [set hidden? false]]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

SLIDER
360
1005
577
1038
Base-Transmission-Chance
Base-Transmission-Chance
0
25
10.0
0.5
1
%
HORIZONTAL

BUTTON
260
270
322
303
Go
Go b-events-next-time
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

BUTTON
328
270
445
303
Go 1 Time Point
go b-events-next-time
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

PLOT
875
10
1220
220
People By Disease State
Days Since Start
People
0.0
1.0
0.0
1.0
true
true
"" ""
PENS
"Susceptible" 1.0 0 -10899396 true "" ""

BUTTON
350
310
445
343
Go 1 Week
go-1-week
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

SWITCH
360
735
565
768
Go-To-Work?
Go-To-Work?
0
1
-1000

SWITCH
360
770
565
803
Go-To-School?
Go-To-School?
0
1
-1000

SWITCH
360
805
565
838
Go-To-Friend?
Go-To-Friend?
0
1
-1000

SWITCH
360
840
565
873
No-Isolation-In-Household?
No-Isolation-In-Household?
0
1
-1000

BUTTON
1230
1020
1367
1053
Highlight Infectious
ask people with [infectious?] [set size size * 2]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

SWITCH
360
940
565
973
Parent-Looks-After-Child?
Parent-Looks-After-Child?
0
1
-1000

BUTTON
1230
1055
1397
1088
Highlight Not Susceptible
ask people with [not susceptible?] [set size size * 2]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

CHOOSER
260
350
445
395
Halt-When
Halt-When
"No Susceptible" "All Post-Infectious" "Infections Impossible" "6 Months" "9 Months" "2 Years"
5

TEXTBOX
10
10
225
45
Disease Decisions
24
0.0
1

TEXTBOX
10
45
203
64
(C) Christopher J Watts, 2021.
11
0.0
1

BUTTON
665
420
752
453
Change
change-view
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

SLIDER
360
875
597
908
Stay-At-Home-With-Symptoms?
Stay-At-Home-With-Symptoms?
0
100
75.0
5
1
%
HORIZONTAL

MONITOR
1230
30
1320
75
Susceptible
num-susceptible
17
1
11

MONITOR
1230
75
1320
120
Exposed
num-exposed
17
1
11

MONITOR
1230
435
1320
480
All Infectious
num-infectious
17
1
11

MONITOR
1230
480
1320
525
All Hospitalize
num-hospitalize
17
1
11

MONITOR
1230
345
1320
390
Recovered
num-recovered
17
1
11

MONITOR
1230
390
1320
435
Dead
num-dead
17
1
11

INPUTBOX
10
790
162
850
People-To-1-Hospital
100000.0
1
0
Number

INPUTBOX
10
850
162
910
People-To-1-Workplace
20.0
1
0
Number

INPUTBOX
10
910
162
970
People-To-1-Primary-School
1250.0
1
0
Number

INPUTBOX
10
970
162
1030
People-To-1-Secondary-School
5000.0
1
0
Number

MONITOR
185
795
335
840
NIL
Count Hospitals
17
1
11

MONITOR
185
855
335
900
NIL
Count Workplaces
17
1
11

MONITOR
185
915
335
960
Count Primary Schools
Count Schools with [l-max-age < 12]
17
1
11

MONITOR
185
975
335
1020
Count Secondary Schools
count schools with [l-max-age > 12]
17
1
11

MONITOR
185
685
335
730
NIL
Count Households
17
1
11

TEXTBOX
15
680
165
705
City Definition:
16
0.0
1

TEXTBOX
360
980
510
998
Disease Transitions:
16
0.0
1

TEXTBOX
365
675
515
700
People's Behaviour:
16
0.0
1

TEXTBOX
10
1050
160
1080
Friendship Network:
16
0.0
1

INPUTBOX
10
1095
162
1155
Friendship-Base-Param
0.0
1
0
Number

INPUTBOX
10
1155
162
1215
Friendship-Distance
0.25
1
0
Number

INPUTBOX
165
1095
335
1155
Friendship-Not-Same-Sex
0.5
1
0
Number

INPUTBOX
10
1215
162
1275
Friendship-Age-Gap
0.25
1
0
Number

INPUTBOX
165
1215
335
1275
Friendship-Not-Same-Workplace
0.5
1
0
Number

INPUTBOX
165
1155
335
1215
Friendship-Not-Same-School
0.75
1
0
Number

TEXTBOX
15
1075
300
1101
Prob(Friendship) = exp(-(Base + Sum(Param * Test)))
11
0.0
1

MONITOR
315
170
375
215
HH:MM
sim-time-as-hh-mm
2
1
11

MONITOR
380
170
445
215
NIL
Day
17
1
11

MONITOR
235
170
310
215
NIL
Sim-Time
17
1
11

PLOT
1235
735
1435
885
People Histogram
X-Axis
Count
0.0
1.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" ""

CHOOSER
1235
690
1637
735
Histo-X-Axis
Histo-X-Axis
"[count my-in-habitants] of households" "[count my-friendships] of people" "[count susceptible-neighbours] of people" "[time-as-days p-infection-time] of people with [not susceptible?]" "[p-num-infected] of people" "[p-sex] of post-disease-people" "[p-age] of post-disease-people" "[p-age-group] of people"
4

BUTTON
1445
775
1572
808
Update Histogram
do-histo-plot
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

BUTTON
345
420
497
453
Schedule 1 Infection
infect-one-susceptible
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

MONITOR
235
220
335
265
Days Since Start
time-as-days sim-time
2
1
11

MONITOR
340
220
445
265
Weeks Since Start
(time-as-days sim-time) / 7
1
1
11

BUTTON
1475
1020
1590
1053
Update Labels
relabel-all-disease-states\nrelabel-all-transitions
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

BUTTON
10
1280
72
1313
NIL
Setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
1475
950
1705
975
Disease State Transition View:
16
0.0
1

TEXTBOX
1230
955
1380
975
City View:
16
0.0
1

MONITOR
1550
140
1670
185
Workplaces
array:item loc-type-num-infections-here 2
17
1
11

MONITOR
1550
95
1670
140
Schools
array:item loc-type-num-infections-here 1
17
1
11

MONITOR
1550
50
1670
95
Households
array:item loc-type-num-infections-here 0
17
1
11

MONITOR
1550
185
1670
230
Hospitals
array:item loc-type-num-infections-here 3
17
1
11

TEXTBOX
1550
10
1675
50
Where Transmission Occur:
13
0.0
1

MONITOR
350
120
445
165
NIL
Count People
17
1
11

INPUTBOX
785
700
937
760
Seed-Setup
0.0
1
0
Number

INPUTBOX
785
760
937
820
Seed-Go
0.0
1
0
Number

TEXTBOX
785
670
1005
690
Random Number Seeds:
16
0.0
1

BUTTON
940
715
1002
748
Clear
set seed-setup 0
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
940
770
1002
803
Clear
set seed-go 0
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1005
715
1107
748
Use Previous
set seed-setup prev-seed-setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1005
770
1107
803
Use Previous
set seed-go prev-seed-go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
545
625
745
645
See below for more parameters.
13
0.0
1

MONITOR
1705
375
1842
420
Case Fatality Rate (%)
100 * total-deaths / total-cases
1
1
11

MONITOR
1775
460
1847
505
Current R
estimated-R
3
1
11

PLOT
1745
575
1945
725
Estimated R
Time (Days)
Estimate
0.0
1.0
0.0
1.0
true
false
"" ""
PENS
"R" 1.0 0 -16777216 true "" ""
"R = 1" 1.0 0 -7500403 true "" ""

PLOT
2235
195
2510
345
Largest Component Size / Susceptibles
Sim-Time (Days)
Standardized Metric
0.0
1.0
0.0
1.0
true
false
"" ""
PENS
"Largest Size" 1.0 0 -14730904 true "" ""

MONITOR
2460
140
2622
185
Size of Largest Component
largest-sus-net-component
1
1
11

MONITOR
2460
90
2552
135
# Components
num-sus-net-components
1
1
11

PLOT
2235
45
2440
195
Susceptibles Net Components
Sim-Time (Days)
Metric
0.0
1.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -2674135 true "" ""

TEXTBOX
2460
50
2665
96
Potential Transmission Network between Susceptibles:
13
0.0
1

PLOT
2235
350
2495
525
Sus. Neighbours of Infectious
Sim-Time (Days)
Sus. Neighbors
0.0
1.0
0.0
1.0
true
true
"" ""
PENS
"Max" 1.0 0 -14333415 true "" ""
"Mean" 1.0 0 -10899396 true "" ""
"Min" 1.0 0 -4399183 true "" ""

MONITOR
1905
425
2042
470
 Count Superspreaders
count superspreaders
17
1
11

MONITOR
1905
475
2137
520
% of Infections due to Superspreaders
100 * (sum [p-num-infected] of superspreaders) / total-infected
1
1
11

MONITOR
1905
525
2062
570
Worst Spreader's Infected
max-infected-by-one
17
1
11

BUTTON
1230
1090
1332
1123
Jiggle People
ask people [\nmove-to p-location\nset heading random 360\nfd 0.4\n]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
785
1190
892
1223
Print B-Events
print \"\"\nprint (word \"B-Events at sim-time=\" sim-time \", ticks=\" ticks)\nforeach b-events [? ->\nprint ?\n]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
260
310
342
343
Go 1 Day
go-1-day
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

CHOOSER
1595
980
1747
1025
Disease-State-Labels
Disease-State-Labels
"Name" "Name-And-Current" "Name-And-Peak" "Current" "Current-And-Peak" ""
0

MONITOR
560
420
660
465
NIL
Current-View
17
1
11

CHOOSER
1595
1030
1733
1075
Transition-Labels
Transition-Labels
"Count" "%" "Time-Param" "Rate" "Probability" ""
5

INPUTBOX
580
1000
732
1060
Contacts-Per-Hour
2.0
1
0
Number

SWITCH
785
880
932
913
Print-Scheduling?
Print-Scheduling?
1
1
-1000

SWITCH
785
915
987
948
Print-Processing-B-Events?
Print-Processing-B-Events?
1
1
-1000

TEXTBOX
785
850
935
870
Debugging Aids:
16
0.0
1

INPUTBOX
10
1415
163
1475
Base-Size
1.0
1
0
Number

BUTTON
785
1065
912
1098
Print p-ds-history
ask ifelse-value (empty? trace-turtles) [one-of people] [trace-turtles with [is-person? self]] [\nforeach reverse p-ds-history [\np -> show (word \n(precision (time-as-days first p) 1) \n\": \" \n([ds-name] of last p)\n)]]\n
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
925
1000
1142
1060
Trace-Turtles-Reporter
seed-infected
1
0
String (reporter)

SWITCH
785
1000
917
1033
Trace-Turtles?
Trace-Turtles?
1
1
-1000

TEXTBOX
925
1065
1165
1160
Enter code that reports a list, turtle-set, or individual turtle of breeds People and/or Locations.\nThe result will be stored in a sorted list trace-turtles.
11
0.0
1

TEXTBOX
360
915
510
933
Not currently in use by code:
11
0.0
1

MONITOR
945
865
1072
910
NIL
Sim-Stopping-Reason
17
1
11

CHOOSER
10
120
215
165
Population-Generator
Population-Generator
"Demographic Data"
0

INPUTBOX
10
730
162
790
People-To-1-Household
4.0
1
0
Number

MONITOR
900
1190
1002
1235
NIL
Length B-Events
17
1
11

CHOOSER
10
70
215
115
Disease-Model
Disease-Model
"LSHTM-Covid-19"
0

CHOOSER
10
170
215
215
Social-Interactions
Social-Interactions
"Contacts-Matrices"
0

INPUTBOX
10
490
155
550
R0
2.68651581109931
1
0
Number

BUTTON
785
1150
1017
1183
Print Age-Group's Force-Of-Infection
foreach n-values (array:length age-groups) [ag-id -> ag-id] [ag-id -> print (word ag-id \": \" force-of-infection ag-id age-groups-total-infectiousness age-groups-total-people)]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
1005
915
1137
948
Using-Profiler?
Using-Profiler?
1
1
-1000

CHOOSER
10
400
337
445
Seed-Infectors
Seed-Infectors
"n-Per-Day-For-t-Days-From-d 1 0 0" "n-Per-Day-For-t-Days-From-d 0.5 2 0" "n-Per-Day-For-t-Days-From-d 0.5 7 0" "n-Per-Day-For-t-Days-From-d 0.5 14 0" "n-Per-Day-For-t-Days-From-d 0.5 21 0" "n-Per-Day-For-t-Days-From-d 0.5 28 0" "n-Per-Day-For-t-Days-From-d 1 1 0" "n-Per-Day-For-t-Days-From-d 1 2 0" "n-Per-Day-For-t-Days-From-d 1 7 0" "n-Per-Day-For-t-Days-From-d 1 14 0" "n-Per-Day-For-t-Days-From-d 1 21 0" "n-Per-Day-For-t-Days-From-d 1 28 0" "n-Per-Day-For-t-Days-From-d 2 1 0" "n-Per-Day-For-t-Days-From-d 2 2 0" "n-Per-Day-For-t-Days-From-d 2 7 0" "n-Per-Day-For-t-Days-From-d 2 14 0" "n-Per-Day-For-t-Days-From-d 2 21 0" "n-Per-Day-For-t-Days-From-d 2 28 0" "n-Per-Day-For-t-Days-From-d 2 28 (random 7)" "n-Per-Day-For-t-Days-From-d 2 28 (random 21)" "n-Per-Day-For-t-Days-From-d 2 28 Seed-Start-Day"
20

BUTTON
1080
865
1142
898
NIL
Setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

CHOOSER
10
220
215
265
Intervention
Intervention
"Base" "School Closures" "Social Distancing" "Elderly Shielding" "Self-Isolation" "Combination"
0

INPUTBOX
270
490
370
550
Intervention-Day
27.0
1
0
Number

TEXTBOX
365
710
550
736
Use with Activity-Locations:
12
0.0
1

INPUTBOX
10
270
110
330
Start-Date
2020-1-29
1
0
String

TEXTBOX
10
335
120
390
YYYY-M-D\n(Used to convert school holiday dates.)
11
0.0
1

SLIDER
590
735
762
768
Weekend-Work
Weekend-Work
0
200
10.0
5
1
%
HORIZONTAL

SWITCH
590
770
767
803
Simulating-Weekends?
Simulating-Weekends?
1
1
-1000

TEXTBOX
15
555
190
590
R0 is used with contacts matrices to calculate Susceptibility.
11
0.0
1

TEXTBOX
167
1422
317
1463
Used for controlling default sizes of people, households, etc.
11
0.0
1

TEXTBOX
590
710
785
736
Use with Contacts-Matrices:
12
0.0
1

TEXTBOX
785
980
935
998
Trace-Turtles:
13
0.0
1

PLOT
875
225
1220
425
Changes
Days Since Start
New People
0.0
1.0
0.0
1.0
true
true
"" ""
PENS
"New Exposed" 1.0 0 -4079321 true "" ""
"New Clinical Cases" 1.0 0 -2674135 true "" ""
"New Hospitalize" 1.0 0 -13791810 true "" ""
"New Recovered" 1.0 0 -2064490 true "" ""
"New Dead" 1.0 0 -16777216 true "" ""

PLOT
875
430
1220
630
Cumulative
Days Since Start
People
0.0
1.0
0.0
1.0
true
true
"" ""
PENS
"Exposed" 1.0 0 -4079321 true "" ""
"Hospitalize?" 1.0 0 -13791810 true "" ""
"Recovered" 1.0 0 -2064490 true "" ""
"Dead" 1.0 0 -16777216 true "" ""

MONITOR
235
120
347
165
NIL
Count Households
17
1
11

TEXTBOX
1905
385
2190
426
\"Superspreaders\" have infected at least the threshold number of people.
11
0.0
1

TEXTBOX
2235
10
2490
40
Susceptibles Network Analysis:
16
0.0
1

TEXTBOX
1235
665
1475
701
User-Definable Histogram:
16
0.0
1

TEXTBOX
1905
355
2110
380
Superspreaders Analysis:
16
0.0
1

MONITOR
1230
120
1320
165
I-Preclinical
num-i-preclinical
17
1
11

MONITOR
1230
165
1320
210
I-Clinical
num-i-clinical
17
1
11

MONITOR
1230
210
1320
255
I-Subclinical
num-i-subclinical
17
1
11

MONITOR
1230
255
1320
300
H-Non-ICU
num-h-non-icu
17
1
11

MONITOR
1230
300
1320
345
H-ICU
num-h-icu
17
1
11

MONITOR
1320
30
1377
75
%
as-perc num-susceptible
1
1
11

MONITOR
1320
75
1377
120
%
as-perc num-exposed
1
1
11

MONITOR
1320
120
1377
165
%
as-perc num-i-preclinical
1
1
11

MONITOR
1320
165
1377
210
%
as-perc num-i-clinical
1
1
11

MONITOR
1320
210
1377
255
%
as-perc num-i-subclinical
1
1
11

MONITOR
1320
255
1377
300
%
as-perc num-h-non-icu
2
1
11

MONITOR
1320
300
1377
345
%
as-perc num-h-icu
1
1
11

MONITOR
1320
345
1377
390
%
as-perc num-recovered
1
1
11

MONITOR
1320
390
1377
435
%
as-perc num-dead
2
1
11

TEXTBOX
1230
10
1315
28
Current:
13
0.0
1

TEXTBOX
1385
60
1460
80
Peak:
13
0.0
1

MONITOR
1385
80
1465
125
Exposed
[ds-max-num] of Exposed
17
1
11

MONITOR
1385
125
1465
170
I-Preclinical
[ds-max-num] of i-preclinical
17
1
11

MONITOR
1385
170
1465
215
I-Clinical
[ds-max-num] of i-clinical
17
1
11

MONITOR
1385
215
1465
260
I-Subclinical
[ds-max-num] of i-subclinical
17
1
11

MONITOR
1385
290
1465
335
H-Non-ICU
Num-Peak-Non-ICU-Beds\n;[ds-max-num] of h-non-icu
17
1
11

MONITOR
1385
335
1465
380
H-ICU
Num-Peak-ICU-Beds\n;[ds-max-num] of h-icu
17
1
11

MONITOR
1465
80
1522
125
%
as-perc [ds-max-num] of Exposed
1
1
11

MONITOR
1465
125
1522
170
%
as-perc [ds-max-num] of i-preclinical
1
1
11

MONITOR
1465
170
1522
215
%
as-perc [ds-max-num] of i-clinical
1
1
11

MONITOR
1465
215
1522
260
%
as-perc [ds-max-num] of i-subclinical
1
1
11

MONITOR
1465
290
1522
335
%
as-perc Num-Peak-Non-ICU-Beds ;[ds-max-num] of h-non-icu
1
1
11

MONITOR
1465
335
1522
380
%
as-perc Num-Peak-ICU-Beds ;[ds-max-num] of h-icu
1
1
11

TEXTBOX
1555
280
1700
300
Peak Transition Rates:
13
0.0
1

MONITOR
1550
305
1640
350
S -> E
max-new-exposed
17
1
11

MONITOR
1550
350
1640
395
-> Hospitalize
max-new-hospitalize
17
1
11

MONITOR
1550
395
1640
440
-> Dead
max-new-deaths
17
1
11

MONITOR
1640
395
1697
440
%
as-perc max-new-deaths
2
1
11

MONITOR
1640
350
1697
395
%
as-perc max-new-hospitalize
1
1
11

MONITOR
1640
305
1697
350
%
as-perc max-new-exposed
1
1
11

MONITOR
1230
560
1320
605
Exposed
total-exposed
17
1
11

MONITOR
1230
605
1320
650
Hospitalize
total-hospitalize
17
1
11

TEXTBOX
1230
535
1380
553
Total To Date:
13
0.0
1

MONITOR
1320
560
1377
605
%
as-perc total-exposed
1
1
11

MONITOR
1320
605
1377
650
%
as-perc total-hospitalize
1
1
11

MONITOR
1385
10
1477
55
NIL
Locked-Down?
17
1
11

INPUTBOX
225
55
335
115
Population-Size
39697.0
1
0
Number

MONITOR
1685
10
1770
55
NIL
Total-Cases
17
1
11

MONITOR
1830
10
1915
55
NIL
Total-Deaths
17
1
11

MONITOR
1685
80
1860
125
NIL
Num-Peak-Week-Cases
17
1
11

MONITOR
1685
125
1860
170
NIL
Num-Peak-Week-Deaths
17
1
11

MONITOR
1685
170
1887
215
NIL
Num-Weeks-To-Peak-Week-Cases
17
1
11

TEXTBOX
1685
60
1835
78
Weekly Outputs:
13
0.0
1

TEXTBOX
1385
270
1505
288
Health Burden:
13
0.0
1

INPUTBOX
580
1130
732
1190
Lockdown-ICU-Level
1.0
1
0
Number

CHOOSER
355
1145
497
1190
Lockdown-Trigger
Lockdown-Trigger
"No-Lockdown" "Intervention-Day" "ICU-Level"
0

TEXTBOX
355
1130
570
1156
Lockdown-Trigger not yet implemented.
11
0.0
1

SWITCH
10
455
195
488
Generate-Random-R0?
Generate-Random-R0?
1
1
-1000

BUTTON
830
1245
992
1278
Print Random R0 values
\nprint n-values 11 [-> precision (random-normal 2.675739 0.5719293) 3]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
1390
570
1515
615
NIL
Num-New-ICU-Beds
17
1
11

MONITOR
1390
480
1515
525
NIL
Num-New-Infected
17
1
11

MONITOR
1390
525
1515
570
NIL
Num-New-Cases
17
1
11

MONITOR
1390
615
1515
660
NIL
Num-New-Deaths
17
1
11

TEXTBOX
1390
455
1540
473
Daily Counts:
13
0.0
1

MONITOR
1705
330
1867
375
Infection Fatality Rate (%)
100 * total-deaths / total-infected
1
1
11

MONITOR
1705
280
1832
325
New Cases (Ip -> Ic)
max-new-cases
17
1
11

MONITOR
1830
280
1887
325
%
as-perc max-new-cases
2
1
11

MONITOR
140
335
237
380
NIL
School-Holiday?
17
1
11

INPUTBOX
355
1370
910
1430
Input-Data-Folder
NIL
1
0
String

MONITOR
140
275
235
320
Date (Y-M-D)
date-string sim-time
17
1
11

BUTTON
1000
955
1142
988
Move to Households
foreach sort households [h -> ask h[setup-cur-household]]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
10
625
182
658
LSHTM-Run
LSHTM-Run
1
200
5.0
1
1
NIL
HORIZONTAL

BUTTON
185
625
267
658
Reset-R0
setup-LSHTM-R0-and-Intervention
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
10
585
187
618
Intervention-Shift
Intervention-Shift
-14
84
0.0
7
1
Days
HORIZONTAL

SWITCH
2535
205
2742
238
Calculate-Susceptibles-Net?
Calculate-Susceptibles-Net?
1
1
-1000

INPUTBOX
2045
410
2197
470
Superspreader-Threshold
8.0
1
0
Number

TEXTBOX
2070
525
2220
581
NB: If using contact matrices and force of infection, program does not count how many you infect.
11
0.0
1

MONITOR
1915
10
1972
55
%
as-perc Total-Deaths
2
1
11

MONITOR
1770
10
1827
55
%
as-perc total-cases
1
1
11

INPUTBOX
160
490
265
550
Seed-Start-Day
17.0
1
0
Number

INPUTBOX
205
555
340
615
Intervention-Duration
84.0
1
0
Number

MONITOR
1975
10
2037
55
NIL
Total-ICU
17
1
11

MONITOR
2035
10
2092
55
%
as-perc Total-ICU
2
1
11

SWITCH
830
1290
992
1323
Hide-Peoples-Links?
Hide-Peoples-Links?
0
1
-1000

BUTTON
355
595
497
628
Day to YYYY-MM-DD
let text-entry user-input \"Enter day number:\"\nif text-entry != false [\nuser-message (word \"Day number \" text-entry \" is \" (date-string (start-date-as-num + (1440 * runresult text-entry))) \" (YYYY-MM-DD).\")\n]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
355
555
497
588
YYYY-MM-DD to Day
let text-entry user-input \"Enter date in the format YYYY-MM-DD:\"\nif text-entry != false [\nuser-message (word text-entry \" is day number \" ((parsed-date text-entry) - start-date-as-num) \".\")\n]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
360
1325
447
1358
Set Folder
setup-input-folder\n
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
450
1325
547
1358
Clear Folder
set input-data-folder \"\"
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1030
1150
1197
1183
Print Symptomatic-Rates
print array:to-list symptomatic-rates
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

@#$#@#$#@
# DISEASE DECISIONS

Version 3.1. This program in NetLogo (C) Christopher J Watts, 2021.

## WHAT IS IT?

This program simulates a disease epidemic in a population. The representation of the disease is based on the model of Covid-19 transmission published in Davies et al. (2020) and developed at the LSHTM. The LSHTM model is an SEIR-type compartmental model, stochastic, and individual-based. It uses contact matrices, giving the contact rate between each pair of 16 age groups, to determine the number of people per day in age group j encountered by a person in age group i. These matrices were derived from the POLYMOD study (Mossong et al. 2008). Also input to the model is demographic data, stating how many people there are in each age group (also broken down by sex, though this model does not use that information).

The current program is a reproduction of the LSHTM model, using NetLogo 6 and some techniques belonging to an agent-based simulation approach. The main points of comparison and difference are:

* Both models sample from a probability distribution, for each person in each compartment, how long that person spends in that compartment before transitioning ("maturing") to another. (In this sense, the LSHTM model is "individual-based", unlike traditional differential equation models such as the SIR model of Kermack & McKendrick 1927).
* Both models sample how many people are flowing out of the Susceptible compartment in the current time step. (This is Binomially distrbuted.) The ABM, however, then has to sample also which people (agents) they are (which requires generating a random permutation, thus a lot more random numbers, requiring more computer time). Similar sampling occurs where the compartmental model contains branch points (e.g. symptomatic / asymptomatic routes, hospitalize / not, death / recovered).
* The ABM uses a discrete-event simulation engine. This should not have any consequences for model behaviour at present, but is used for future extensions.
* The ABM can represent agents as having other attributes, and much more heterogeneity than is used by the present model processes. Again, this is for future extensions. The LSHTM model only uses Age-Group and Disease-State to determine what happens to people.

## HOW IT WORKS

There are two views: "City" and "Disease States". Clicking the "Change View" button under the NetLogo World will switch between the two.

* In "City" view, every person and location is visible. For reproducing the LSHTM model, agents' coordinates and locations play no role in the model. (Contacts are based on the Contact Matrices, not on who you are with in your current location.) So we recommend running with the "view updates" tickbox unchecked to speed things up.
* In "Disease States" view, the compartmental structure is shown instead using a NetLogo turtle breed (diseae-states) for compartments, and directed links (transitions) for flow paths.

On running "Setup":

* Disease-states and transitions are created, and associated with processes to determine a person's transition routes and times, given their age group.
* Demographic data are loaded from the file "Demog.csv", households are created to locate people in, and people are created to match the demographics.
* Contact matrices are loaded from the files "cm_home.csv", "cm_work.csv", "cm_school.csv" and "cm_other.csv". The overall contact rates between each pair of age groups is the sum of these four matrices.
* Intervention start and end dates may be defined. During an intervention contact rates are a weighted sum of the four input matrices. (E.g. for "School Closures", the school matrix is multiplied by 0%.)
* School holidays may be scheduled. Dates are included in the code for the UK. During school holidays, school contact rates are 0.
* Seed infection events are scheduled. Starting on a given day d, there will be n people infected as seeds per day, for t days.

### Disease progression

The disease progresses through several stages. The LSHTM model contained a transmission model (SEI3R structure), and added a health burden model and a deaths model to run in parallel to it. Parallel models are underdesirable in an ABM, so we have attempted to merge them into the main structure.

* **Susceptible**: The normal state a person is in, susceptible to the disease being transmitted to them by anyone they encounter who has it.
* **Exposed**: A person has been infected, but is not yet infectious, nor displaying symptoms.
* **I-Preclinical**: A person can now infect others, and will become clinically ill, but does not yet display symptoms.
* **I-Clinical**: A person is infectious and nows displays symptoms such as coughing, a temperature, and aches, and might now alter their behaviour in response (e.g. staying home).
* **I-Subclinical**: A person is less infectious (50% of I-Preclinical and I-Clinical), and asymptomatic.
* **Hospitalize?**: A dummy state, in which it is decided whether or not a person is proceding directly to Recovered, or becoming so ill they need hospital treatment. 
* **H-Non-ICU**: After hospitalization, a person occupies a bed on a General Ward.
* **H-ICU**: After hospitalization, a person occupies a bed in an Intensive Care Unit.
* **Recovered**: A person is no longer ill or infectious. For now, they are immune to further infection.
* **Dead**: A person has died and performs no further activities.

People in the states I-Preclinical, I-Clinical, and I-Subclinical are capable of transmitting the disease to other people. (In the current version, nobody works at the hospital, and so there is no transmission there.)

Recovery can occur from the I-Subclinical, and H- states.

Death can occur only from the H- states.

We assume Recovered people do not lose their immunity. A simple extension would allow them to pass again to Susceptible after a given time period. For some other infectious diseases (the common cold and regular seasonal flu are examples), immunity can become irrelevant within months.

### The simulation run

The simulation starts at midnight with everyone in their homes (households). 

On running the Go procedure, the clock ("sim-time") is advanced to the time of the first event in the events list ("b-events"). All events scheduled for that time are processed.

In a seed event, a randomly chosen Susceptible person is made Exposed (infected by someone outside the represented population).

Infection events are scheduled for every 6-hour time step. In an infection event, the rate of people becoming infected in each age group is calculated ("the force of infection"). This is based on the contact rate with and proportion of the population infectious in each age group, multiplied by the chance of a contact with an infectious person resulting in infection ("Susceptibility"). Susceptibility is derived from the Basic Reproduction Number, R0.

On entering a new disease-state, a person's transition route out of it and their transition time are sampled. Transition events may be processed at times other than the 6-hourly time steps of infection events.

Each day at midnight, an update of output metrics is scheduled, so that, e.g., we know how many people died that day.

The simulation runs until a halting event. This could be, e.g., the clock reaching 2 years since the start.

## HOW TO USE IT

Click "Setup". Click one of the "Go" buttons.

## THINGS TO NOTICE

Providing there are sufficient numbers of seed infections, there should appear an epidemic. The number of new cases, and the numbers of people in each disease-state after Susceptible rise and fall. The total numbers of people Recovered or Dead rise in S-shaped curves.

## THINGS TO TRY

Try running with an Intervention. This alters the contact rates for a 12-week period. Also try different start dates and start shifts for Interventions. What happens to the curves?

## EXTENDING THE MODEL

### Extending the disease model

Suggestions include:

* Alternative disease models: Examine the procedure "setup-disease-LSHTM-Covid-19" to understand how to create compartments and transitions.
* Replacing the contact matrices: These are crude simplifications. They state mean rates, with no variations (no one can be a superspreader). Their use assumes contacts can be with anyone in the age group, with no preference for family, co-workers, schoolmates, friends, or neighbours. They do not state how long contacts are; 10 hours with a spouse at home is 1 contact, as is 30 seconds talking to the bus driver. Could an ABM be designed and calibrated that improves on this depiction of social contact?

## NETLOGO FEATURES

Note how the three-phase approach to discrete-event simulation (Tocher 1963) is implemented in NetLogo (procedures beginning "b-events-...").

## RELATED MODELS

The SIR model of Kermack & McKendrick (1927) was the original compartmental model. Many others have been proposed since (SEIR, SEIRD, SEIHRD, etc.)
The S-I model, based on the logistic function, is sometimes used as a model of the diffusion of innovations.

## CREDITS AND REFERENCES

The LSHTM model

The disease model is derived from one developed at the London School of Hygiene and Tropical Medicine (LSHTM), and described in the following paper, its supplementary materials, and the associated program posted on github (especially the files UK.R and corona.cpp).

Davies, Nicholas, Adam J Kucharski1, Rosalind M Eggo1, CMMID nCov working group & W John Edmunds (2020) "The effect of non-pharmaceutical interventions on COVID-19 cases, deaths and demand for hospital services in the UK: a modelling study". The Lancet. https://doi.org/10.1016/S2468-2667(20)30133-X

https://github.com/cmmid/covid-UK


### Other references

Kermack William Ogilvy and McKendrick A. G. (1927) A contribution to the mathematical theory of epidemicsProc. R. Soc. Lond. A115700–721
http://doi.org/10.1098/rspa.1927.0118

Tocher, K.D. (1963) "The Art of Simulation", English Universities Press.
Sterman, John (2000) "Business Dynamics: Systems Thinking and Modeling for a Complex World", McGraw-Hill Education.

@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

ball football
false
0
Polygon -7500403 false true 301 133 301 164 275 192 229 224 167 236 137 236 74 224 30 194 3 162 2 138 30 104 76 74 134 62 168 62 228 74 274 105
Polygon -7500403 true true 300 150 300 165 270 195 225 225 163 236 134 236 75 225 30 195 2 162 2 140 30 105 75 75 136 63 165 63 225 75 270 105 300 135
Line -16777216 false 300 155 5 155
Polygon -1 true false 28 193 28 107 51 91 51 209
Rectangle -1 true false 90 150 210 160
Rectangle -1 true false 198 141 205 170
Rectangle -1 true false 183 141 190 170
Rectangle -1 true false 168 141 175 170
Rectangle -1 true false 153 141 160 170
Rectangle -1 true false 138 141 145 170
Rectangle -1 true false 123 141 130 170
Rectangle -1 true false 108 141 115 170
Rectangle -1 true false 93 141 100 170
Polygon -1 true false 272 193 272 107 249 91 249 209

ball tennis
false
0
Circle -7500403 true true 30 30 240
Circle -7500403 false true 30 30 240
Polygon -16777216 true false 50 82 54 90 59 107 64 140 64 164 63 189 59 207 54 222 68 236 76 220 81 195 84 163 83 139 78 102 72 83 63 67
Polygon -16777216 true false 250 82 246 90 241 107 236 140 236 164 237 189 241 207 246 222 232 236 224 220 219 195 216 163 217 139 222 102 228 83 237 67
Polygon -1 true false 247 79 243 86 237 106 232 138 232 167 235 199 239 215 244 225 236 234 229 221 224 196 220 163 221 138 227 102 234 83 240 71
Polygon -1 true false 53 79 57 86 63 106 68 138 68 167 65 199 61 215 56 225 64 234 71 221 76 196 80 163 79 138 73 102 66 83 60 71

bottle
false
0
Circle -7500403 true true 90 240 60
Rectangle -1 true false 135 8 165 31
Line -7500403 true 123 30 175 30
Circle -7500403 true true 150 240 60
Rectangle -7500403 true true 90 105 210 270
Rectangle -7500403 true true 120 270 180 300
Circle -7500403 true true 90 45 120
Rectangle -7500403 true true 135 27 165 51

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bread
false
0
Polygon -16777216 true false 140 145 170 250 245 190 234 122 247 107 260 79 260 55 245 40 215 32 185 40 155 31 122 41 108 53 28 118 110 115 140 130
Polygon -7500403 true true 135 151 165 256 240 196 225 121 241 105 255 76 255 61 240 46 210 38 180 46 150 37 120 46 105 61 47 108 105 121 135 136
Polygon -1 true false 60 181 45 256 165 256 150 181 165 166 180 136 180 121 165 106 135 98 105 106 75 97 46 107 29 118 30 136 45 166 60 181
Polygon -16777216 false false 45 255 165 255 150 180 165 165 180 135 180 120 165 105 135 97 105 105 76 96 46 106 29 118 30 135 45 165 60 180
Line -16777216 false 165 255 239 195

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

building institution
false
0
Rectangle -7500403 true true 0 60 300 270
Rectangle -16777216 true false 130 196 168 256
Rectangle -16777216 false false 0 255 300 270
Polygon -7500403 true true 0 60 150 15 300 60
Polygon -16777216 false false 0 60 150 15 300 60
Circle -1 true false 135 26 30
Circle -16777216 false false 135 25 30
Rectangle -16777216 false false 0 60 300 75
Rectangle -16777216 false false 218 75 255 90
Rectangle -16777216 false false 218 240 255 255
Rectangle -16777216 false false 224 90 249 240
Rectangle -16777216 false false 45 75 82 90
Rectangle -16777216 false false 45 240 82 255
Rectangle -16777216 false false 51 90 76 240
Rectangle -16777216 false false 90 240 127 255
Rectangle -16777216 false false 90 75 127 90
Rectangle -16777216 false false 96 90 121 240
Rectangle -16777216 false false 179 90 204 240
Rectangle -16777216 false false 173 75 210 90
Rectangle -16777216 false false 173 240 210 255
Rectangle -16777216 false false 269 90 294 240
Rectangle -16777216 false false 263 75 300 90
Rectangle -16777216 false false 263 240 300 255
Rectangle -16777216 false false 0 240 37 255
Rectangle -16777216 false false 6 90 31 240
Rectangle -16777216 false false 0 75 37 90
Line -16777216 false 112 260 184 260
Line -16777216 false 105 265 196 265

building store
false
0
Rectangle -7500403 true true 30 45 45 240
Rectangle -16777216 false false 30 45 45 165
Rectangle -7500403 true true 15 165 285 255
Rectangle -16777216 true false 120 195 180 255
Line -7500403 true 150 195 150 255
Rectangle -16777216 true false 30 180 105 240
Rectangle -16777216 true false 195 180 270 240
Line -16777216 false 0 165 300 165
Polygon -7500403 true true 0 165 45 135 60 90 240 90 255 135 300 165
Rectangle -7500403 true true 0 0 75 45
Rectangle -16777216 false false 0 0 75 45

bus
false
0
Polygon -7500403 true true 15 206 15 150 15 120 30 105 270 105 285 120 285 135 285 206 270 210 30 210
Rectangle -16777216 true false 36 126 231 159
Line -7500403 false 60 135 60 165
Line -7500403 false 60 120 60 165
Line -7500403 false 90 120 90 165
Line -7500403 false 120 120 120 165
Line -7500403 false 150 120 150 165
Line -7500403 false 180 120 180 165
Line -7500403 false 210 120 210 165
Line -7500403 false 240 135 240 165
Rectangle -16777216 true false 15 174 285 182
Circle -16777216 true false 48 187 42
Rectangle -16777216 true false 240 127 276 205
Circle -16777216 true false 195 187 42
Line -7500403 false 257 120 257 207

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

factory
false
0
Rectangle -7500403 true true 76 194 285 270
Rectangle -7500403 true true 36 95 59 231
Rectangle -16777216 true false 90 210 270 240
Line -7500403 true 90 195 90 255
Line -7500403 true 120 195 120 255
Line -7500403 true 150 195 150 240
Line -7500403 true 180 195 180 255
Line -7500403 true 210 210 210 240
Line -7500403 true 240 210 240 240
Line -7500403 true 90 225 270 225
Circle -1 true false 37 73 32
Circle -1 true false 55 38 54
Circle -1 true false 96 21 42
Circle -1 true false 105 40 32
Circle -1 true false 129 19 42
Rectangle -7500403 true true 14 228 78 270

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

hospital
false
0
Rectangle -7500403 true true 15 30 285 255
Rectangle -16777216 true false 120 195 180 255
Line -7500403 true 150 195 150 255
Rectangle -16777216 true false 45 195 105 240
Rectangle -16777216 true false 210 195 270 240
Line -7500403 true 75 195 75 240
Line -7500403 true 240 195 240 240
Line -16777216 false 285 105 285 255
Line -16777216 false 15 180 285 180
Rectangle -16777216 true false 45 120 105 165
Rectangle -16777216 true false 210 120 270 165
Rectangle -1 true false 120 45 195 90
Line -7500403 true 75 120 75 165
Line -7500403 true 240 120 240 165
Line -16777216 false 15 105 285 105
Rectangle -2674135 true false 150 45 165 90
Rectangle -2674135 true false 135 60 180 75
Rectangle -16777216 true false 45 45 105 90
Rectangle -16777216 true false 210 45 270 90
Rectangle -16777216 true false 120 120 195 165
Line -7500403 true 75 120 75 165
Line -7500403 true 75 45 75 90
Line -7500403 true 240 45 240 90
Line -7500403 true 150 120 150 165

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

house ranch
false
0
Rectangle -7500403 true true 270 120 285 255
Rectangle -7500403 true true 15 180 270 255
Polygon -7500403 true true 0 180 300 180 240 135 60 135 0 180
Rectangle -16777216 true false 120 195 180 255
Line -7500403 true 150 195 150 255
Rectangle -16777216 true false 45 195 105 240
Rectangle -16777216 true false 195 195 255 240
Line -7500403 true 75 195 75 240
Line -7500403 true 225 195 225 240
Line -16777216 false 270 180 270 255
Line -16777216 false 0 180 300 180

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

train
false
0
Rectangle -7500403 true true 30 105 240 150
Polygon -7500403 true true 240 105 270 30 180 30 210 105
Polygon -7500403 true true 195 180 270 180 300 210 195 210
Circle -7500403 true true 0 165 90
Circle -7500403 true true 240 225 30
Circle -7500403 true true 90 165 90
Circle -7500403 true true 195 225 30
Rectangle -7500403 true true 0 30 105 150
Rectangle -16777216 true false 30 60 75 105
Polygon -7500403 true true 195 180 165 150 240 150 240 180
Rectangle -7500403 true true 135 75 165 105
Rectangle -7500403 true true 225 120 255 150
Rectangle -16777216 true false 30 203 150 218

train passenger car
false
0
Polygon -7500403 true true 15 206 15 150 15 135 30 120 270 120 285 135 285 150 285 206 270 210 30 210
Circle -16777216 true false 240 195 30
Circle -16777216 true false 210 195 30
Circle -16777216 true false 60 195 30
Circle -16777216 true false 30 195 30
Rectangle -16777216 true false 30 140 268 165
Line -7500403 true 60 135 60 165
Line -7500403 true 60 135 60 165
Line -7500403 true 90 135 90 165
Line -7500403 true 120 135 120 165
Line -7500403 true 150 135 150 165
Line -7500403 true 180 135 180 165
Line -7500403 true 210 135 210 165
Line -7500403 true 240 135 240 165
Rectangle -16777216 true false 5 195 19 207
Rectangle -16777216 true false 281 195 295 207
Rectangle -13345367 true false 15 165 285 173
Rectangle -2674135 true false 15 180 285 188

train passenger engine
false
0
Rectangle -7500403 true true 0 180 300 195
Polygon -7500403 true true 283 161 274 128 255 114 231 105 165 105 15 105 15 150 15 195 15 210 285 210
Circle -16777216 true false 17 195 30
Circle -16777216 true false 50 195 30
Circle -16777216 true false 220 195 30
Circle -16777216 true false 253 195 30
Rectangle -16777216 false false 0 195 300 180
Rectangle -1 true false 11 111 18 118
Rectangle -1 true false 270 129 277 136
Rectangle -16777216 true false 91 195 210 210
Rectangle -16777216 true false 1 180 10 195
Line -16777216 false 290 150 291 182
Rectangle -16777216 true false 165 90 195 90
Rectangle -16777216 true false 290 180 299 195
Polygon -13345367 true false 285 180 267 158 239 135 180 120 15 120 16 113 180 113 240 120 270 135 282 154
Polygon -2674135 true false 284 179 267 160 239 139 180 127 15 127 16 120 180 120 240 127 270 142 282 161
Rectangle -16777216 true false 210 115 254 135
Line -7500403 true 225 105 225 150
Line -7500403 true 240 105 240 150

tram
false
0
Polygon -7500403 true true 15 206 15 150 15 120 30 105 270 105 285 120 285 135 285 206 270 210 30 210
Rectangle -16777216 true false 165 126 231 165
Rectangle -16777216 true false 36 126 105 165
Line -7500403 true 60 135 60 165
Line -7500403 true 60 120 60 165
Line -7500403 true 90 120 90 165
Line -7500403 true 120 120 120 165
Line -7500403 true 180 120 180 165
Line -7500403 true 210 120 210 165
Line -7500403 true 240 135 240 165
Rectangle -16777216 true false 15 174 285 182
Circle -16777216 true false 48 187 42
Rectangle -16777216 true false 120 127 150 195
Circle -16777216 true false 195 187 42
Line -7500403 true 137 120 137 207
Line -7500403 true 60 105 45 90
Line -7500403 true 45 90 60 75
Line -7500403 true 60 75 75 90
Line -7500403 true 75 90 60 105
Line -7500403 true 240 105 255 90
Line -7500403 true 255 90 240 75
Line -7500403 true 240 75 225 90
Line -7500403 true 225 90 240 105
Rectangle -16777216 true false 255 120 285 165
Line -7500403 true 285 120 285 195
Line -7500403 true 0 225 300 225

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="experiment-LSHTM-Docking-50-R0s" repetitions="5" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>setup-LSHTM-R0-and-Intervention
setup</setup>
    <go>go-1-day</go>
    <metric>timer</metric>
    <metric>time-taken-setup</metric>
    <metric>time-taken-go</metric>
    <metric>prev-seed-setup</metric>
    <metric>prev-seed-go</metric>
    <metric>Input-Data-Folder</metric>
    <metric>R0</metric>
    <metric>Seed-Start-Day</metric>
    <metric>Intervention-Day</metric>
    <metric>array:to-list symptomatic-rates</metric>
    <metric>count patches</metric>
    <metric>count people</metric>
    <metric>count households</metric>
    <metric>sim-time</metric>
    <metric>last-midnight = sim-time</metric>
    <metric>ceiling time-as-days sim-time</metric>
    <metric>max-new-cases</metric>
    <metric>total-cases</metric>
    <metric>total-deaths</metric>
    <metric>total-icu</metric>
    <metric>num-peak-week-cases</metric>
    <metric>num-peak-week-deaths</metric>
    <metric>num-peak-icu-beds</metric>
    <metric>num-peak-non-icu-beds</metric>
    <metric>num-weeks-to-peak-week-cases</metric>
    <metric>num-days-in-lockdown</metric>
    <metric>total-infected</metric>
    <metric>[ds-cur-num] of susceptible</metric>
    <metric>[ds-cur-num] of exposed</metric>
    <metric>[ds-cur-num] of i-preclinical</metric>
    <metric>[ds-cur-num] of i-clinical</metric>
    <metric>[ds-cur-num] of i-subclinical</metric>
    <metric>[ds-cur-num] of hospitalize</metric>
    <metric>[ds-cur-num] of h-icu</metric>
    <metric>[ds-cur-num] of h-non-icu</metric>
    <metric>[ds-cur-num] of dead</metric>
    <metric>[ds-cur-num] of recovered</metric>
    <metric>[ds-max-num] of exposed</metric>
    <metric>[ds-max-num] of i-preclinical</metric>
    <metric>[ds-max-num] of i-clinical</metric>
    <metric>[ds-max-num] of i-subclinical</metric>
    <metric>[ds-max-at-time] of exposed</metric>
    <metric>[ds-max-at-time] of i-preclinical</metric>
    <metric>[ds-max-at-time] of i-clinical</metric>
    <metric>[ds-max-at-time] of i-subclinical</metric>
    <enumeratedValueSet variable="Population-Size">
      <value value="39697"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Model">
      <value value="&quot;LSHTM-Covid-19&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Population-Generator">
      <value value="&quot;Demographic Data&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-Interactions">
      <value value="&quot;Contacts-Matrices&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Seed-Infectors">
      <value value="&quot;n-Per-Day-For-t-Days-From-d 2 28 Seed-Start-Day&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Intervention">
      <value value="&quot;Base&quot;"/>
      <value value="&quot;School Closures&quot;"/>
      <value value="&quot;Social Distancing&quot;"/>
      <value value="&quot;Elderly Shielding&quot;"/>
      <value value="&quot;Self-Isolation&quot;"/>
      <value value="&quot;Combination&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Intervention-Shift">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Intervention-Duration">
      <value value="84"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Halt-When">
      <value value="&quot;2 Years&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="LSHTM-Run">
      <value value="5"/>
      <value value="6"/>
      <value value="8"/>
      <value value="9"/>
      <value value="11"/>
      <value value="12"/>
      <value value="14"/>
      <value value="17"/>
      <value value="19"/>
      <value value="25"/>
      <value value="34"/>
      <value value="38"/>
      <value value="40"/>
      <value value="41"/>
      <value value="42"/>
      <value value="51"/>
      <value value="52"/>
      <value value="56"/>
      <value value="66"/>
      <value value="70"/>
      <value value="75"/>
      <value value="82"/>
      <value value="88"/>
      <value value="91"/>
      <value value="96"/>
      <value value="97"/>
      <value value="105"/>
      <value value="110"/>
      <value value="113"/>
      <value value="114"/>
      <value value="116"/>
      <value value="119"/>
      <value value="120"/>
      <value value="126"/>
      <value value="127"/>
      <value value="129"/>
      <value value="138"/>
      <value value="143"/>
      <value value="145"/>
      <value value="146"/>
      <value value="153"/>
      <value value="154"/>
      <value value="159"/>
      <value value="164"/>
      <value value="176"/>
      <value value="185"/>
      <value value="188"/>
      <value value="190"/>
      <value value="192"/>
      <value value="193"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Generate-Random-R0?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Perc-Patches-With-Household">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Simulating-Weekends?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Weekend-Work">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Start-Date">
      <value value="&quot;2020-1-29&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="People-To-1-Hospital">
      <value value="100000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Print-Processing-B-Events?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Print-Scheduling?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Using-Profiler?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Seed-Setup">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Seed-Go">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-LSHTM-Docking-IVShift-50-R0s" repetitions="10" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>setup-LSHTM-R0-and-Intervention
setup</setup>
    <go>go-1-day</go>
    <metric>timer</metric>
    <metric>time-taken-setup</metric>
    <metric>time-taken-go</metric>
    <metric>prev-seed-setup</metric>
    <metric>prev-seed-go</metric>
    <metric>Input-Data-Folder</metric>
    <metric>R0</metric>
    <metric>Seed-Start-Day</metric>
    <metric>Intervention-Day</metric>
    <metric>array:to-list symptomatic-rates</metric>
    <metric>count patches</metric>
    <metric>count people</metric>
    <metric>count households</metric>
    <metric>sim-time</metric>
    <metric>last-midnight = sim-time</metric>
    <metric>ceiling time-as-days sim-time</metric>
    <metric>max-new-cases</metric>
    <metric>total-cases</metric>
    <metric>total-deaths</metric>
    <metric>total-icu</metric>
    <metric>num-peak-week-cases</metric>
    <metric>num-peak-week-deaths</metric>
    <metric>num-peak-icu-beds</metric>
    <metric>num-peak-non-icu-beds</metric>
    <metric>num-weeks-to-peak-week-cases</metric>
    <metric>num-days-in-lockdown</metric>
    <metric>total-infected</metric>
    <metric>[ds-cur-num] of susceptible</metric>
    <metric>[ds-cur-num] of exposed</metric>
    <metric>[ds-cur-num] of i-preclinical</metric>
    <metric>[ds-cur-num] of i-clinical</metric>
    <metric>[ds-cur-num] of i-subclinical</metric>
    <metric>[ds-cur-num] of hospitalize</metric>
    <metric>[ds-cur-num] of h-icu</metric>
    <metric>[ds-cur-num] of h-non-icu</metric>
    <metric>[ds-cur-num] of dead</metric>
    <metric>[ds-cur-num] of recovered</metric>
    <metric>[ds-max-num] of exposed</metric>
    <metric>[ds-max-num] of i-preclinical</metric>
    <metric>[ds-max-num] of i-clinical</metric>
    <metric>[ds-max-num] of i-subclinical</metric>
    <metric>[ds-max-at-time] of exposed</metric>
    <metric>[ds-max-at-time] of i-preclinical</metric>
    <metric>[ds-max-at-time] of i-clinical</metric>
    <metric>[ds-max-at-time] of i-subclinical</metric>
    <enumeratedValueSet variable="Population-Size">
      <value value="39697"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Model">
      <value value="&quot;LSHTM-Covid-19&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Population-Generator">
      <value value="&quot;Demographic Data&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-Interactions">
      <value value="&quot;Contacts-Matrices&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Seed-Infectors">
      <value value="&quot;n-Per-Day-For-t-Days-From-d 2 28 Seed-Start-Day&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Intervention">
      <value value="&quot;Combination&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Intervention-Shift">
      <value value="-7"/>
      <value value="0"/>
      <value value="7"/>
      <value value="14"/>
      <value value="21"/>
      <value value="28"/>
      <value value="35"/>
      <value value="42"/>
      <value value="56"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Intervention-Duration">
      <value value="84"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Halt-When">
      <value value="&quot;2 Years&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="LSHTM-Run">
      <value value="5"/>
      <value value="6"/>
      <value value="8"/>
      <value value="9"/>
      <value value="11"/>
      <value value="12"/>
      <value value="14"/>
      <value value="17"/>
      <value value="19"/>
      <value value="25"/>
      <value value="34"/>
      <value value="38"/>
      <value value="40"/>
      <value value="41"/>
      <value value="42"/>
      <value value="51"/>
      <value value="52"/>
      <value value="56"/>
      <value value="66"/>
      <value value="70"/>
      <value value="75"/>
      <value value="82"/>
      <value value="88"/>
      <value value="91"/>
      <value value="96"/>
      <value value="97"/>
      <value value="105"/>
      <value value="110"/>
      <value value="113"/>
      <value value="114"/>
      <value value="116"/>
      <value value="119"/>
      <value value="120"/>
      <value value="126"/>
      <value value="127"/>
      <value value="129"/>
      <value value="138"/>
      <value value="143"/>
      <value value="145"/>
      <value value="146"/>
      <value value="153"/>
      <value value="154"/>
      <value value="159"/>
      <value value="164"/>
      <value value="176"/>
      <value value="185"/>
      <value value="188"/>
      <value value="190"/>
      <value value="192"/>
      <value value="193"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Generate-Random-R0?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Perc-Patches-With-Household">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Simulating-Weekends?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Weekend-Work">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Start-Date">
      <value value="&quot;2020-1-29&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="People-To-1-Hospital">
      <value value="100000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Print-Processing-B-Events?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Print-Scheduling?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Using-Profiler?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Seed-Setup">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Seed-Go">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

transition
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 60 225
Line -7500403 true 150 150 240 225
@#$#@#$#@
1
@#$#@#$#@
