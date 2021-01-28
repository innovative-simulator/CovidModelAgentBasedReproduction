# A Covid Model's Agent-Based Reproduction

This is an agent-based simulation developed to reproduce the model of COVID-19 transmission and non-pharmaceutical interventions. The model being reproduced is that of Davies et al. (2020), and is available at

https://github.com/cmmid/covid-UK 

Our reproducing ABM is described in Chattoe-Brown et al. (2021, not peer reviewed), which discusses similarities and differences between the two models' mechanisms, and compares their behaviour.

All code on this git is (C) Christopher J. Watts (2021), and covered by GNU Public Licence 3.0.

## Subfolders in this git

### Rutland_Model

This contains the NetLogo model 

* __DiseaseDecisions.nlogo__

from which the simulation experiment "experiment-LSHTM-Docking-50-R0s" was run. 

Also in this folder are all data files required as input to the model (in .csv format):

* __Demog.csv__ : The population size for the county of Rutland, by age group.
* __cm_....csv__ : Contact matrices, based on the POLYMOD study, and rescaled for Rutland. There is one file for each of home, work, school, and other.
* __Runs_R0_SeedStart_Peakt.csv__ : So that our comparison of models focuses on mechanisms rather than parameters, we output from the LSHTM model values for R0, the start date of seed infections, and the date of peak cases in the Base scenario (from which intervention start dates are calculated). This file supplies them for each of 200 runs.
* __2-linelist_symp_fit_fIa0.5.csv__ : This is taken from the LSHTM model. It contains Markov Chain estimates for the joint posterior probability distribution for age-based symptomatic rates.

Two shell/batch scripts for running the NetLogo experiment in headless mode.

* __DiseaseDecisions_docking.bat__ : For Windows.
* __DiseaseDecisions_docking.sh__ : For Linux.

The scripts are included as guides only, and you may need to edit the pathnames in these for your own computer.

### Output_Analysis

There are three .qs data files:

* __X.qs__ : This is the data from the NetLogo BehaviourSpace experiment.
* __1-dynamics_Corrected_Rutland.qs__ : This is the data output from the LSHTM model, using our correction of the Seed Infections bug.
* __1-dynamics_Uncorrected_Rutland.qs__ : This is the data output from the LSHTM model, using the original, uncorrected version of UK.R.

The results from NetLogo and the LSHTM model are compared and the figures in the paper drawn in the script:

* __Cases_Analysis.R__

Additional files used in producing the analysis include:

* __Docking_X.R__ : We ran the NetLogo model on various computers. Each one generated a .csv file. This script combines the .csv files into a single data table and saves it as the file X.qs
* __UK_corrected.R__ : This is our corrected version of the LSHTM model's main script, "UK.R". The original contained a bug in its creation of seed infections for intervention scenarios. As well as copying seed infection schedules over from the Baseline scenario, our version also saves the 1-dynamics.qs file periodically, rather than just at the end of the experiment.
* __R0s_Select_50.xlsx__ : The analysis in Davies et al. is based on 200 runs, each with its own R0 parameter value. We selected 50 of these for our reproduction, evenly spread across the range. This .xlsx file contains our selection.


## References

Chattoe-Brown, E., Gilbert, N., Robertson, D. A., Watts, C. (2021, not peer reviewed) "Reproduction as a Means of Evaluating Policy Models: A Case Study of a COVID-19 Simulation".

Davies, Nicholas G., Kucharski, Adam J., Eggo, Rosalind M., Gimma, Amy and Edmunds, W. John on behalf of the Centre for the Mathematical Modelling of Infectious Diseases COVID-19 Working Group (2020) ‘Effects of Non-Pharmaceutical Interventions on COVID-19 Cases, Deaths, and Demand for Hospital Services in the UK: A Modelling Study’, Lancet Public Health, 5(7), July, pp. E375-E385. doi:10.1016/S2468-2667(20)30133-X

