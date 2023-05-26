########################################################################
########################################################################
## Load datasets from LSHTM model and DiseaseDecisions.nlogo experiment 
## Do plots for Paper R, comparing programs by R0
## (Also code for investigating how many repetitions (in NetLogo)
## should be run for each R0 value.
########################################################################
########################################################################


library(data.table)
library(qs)
library(tidyverse)
library(ggplot2)

# Set the working directory as required for your computer:
setwd("~\\CovidModelAgentBasedReproduction\\Output_Analysis")
#setwd("D:\\MyDocus\\Simulation\\NetLogo\\Diffusion\\DiseaseDecisions\\git\\CovidModelAgentBasedReproduction\\Output_Analysis")

########################################################################

# Data file from LSHTM model (Target Model):
# Decide which one you want to use.
#D <- qread("1-dynamics_Uncorrected_Rutland.qs")
D <- qread("1-dynamics_Corrected_Rutland.qs")

# Restrict analysis to one region
D <- D[run <= 200 & region=="England"]

# Add population size and rename field
DN <- D[compartment %in% list("S", "E", "Ip", "Is", "Ia", "R")]
D <- cbind(D, DN[, sum(value), by=list(scenario, run, t, region)])
D <- D[,list(scenario,run,t,compartment,region,value,Pop_Size=V1)]

# Add the R0 values for each run:
R0s <- data.table(read.csv("Runs_R0_SeedStart_Peakt.csv"), key="R0")
R0s <- R0s[,list(R0=R0),key=run]
setkey(D, run)
D <- R0s[D] # Left outer join
head(D)  # Check we still have data

########################################################################

# Read in the data from DiseaseDecisions.nlogo (Reproduction Model):
X <- qread("X.qs") # X comes from the Table csv files output by DiseaseDecisions.nlogo

# Check what you have loaded/created:
#head(D)
#unique(D[,compartment])
#unique(D[,R0])
#head(X)
#colnames(X)
#length(unique(X[,R0]))

# Use only the R0s that are in X:
Selected_R0s <- data.table(unique(X[,R0]))
setkey(D, R0)
dim(D)
D <- D[Selected_R0s]
dim(D)

# Mapping scenario field in D to that in X
# (Unnecessary now - the two programs have now been synched.)
scenarios <- data.table(scenario=unique(D[,scenario]), intervention=unique(D[,scenario]))		
scenarios

setkey(D, R0)
setkey(X, R0)

########################################################################

# Construct comparable data sets and plot them

# Peak Cases
plot_peak_cases <- function(scen="Base", iv_delay=0, ylim=c(0,5)) {
	DC <- D[scenario==scen & compartment=="cases" & region=="England", max(value),by=list(scenario, compartment, region, run, R0, Pop_Size)]
	DC <- DC[,list(run, R0=round(R0,5), Value=100 * V1 / Pop_Size)]
	XC <- X[Intervention==scenarios[scenario==scen,intervention] & Intervention.Shift==iv_delay,list(count.people, R0=round(R0,5), ABM_Value=100 * max.new.cases / count.people), keyby=max.new.cases]
	plot_comparison(XC, DC, ylim, ylab="Peak Cases as % of Population")
}

# Total Cases
plot_total_cases <- function(scen="Base", iv_delay=0, ylim=c(0,5)) {
	DC <- D[scenario==scen & compartment=="cases" & region=="England", sum(value),by=list(scenario, compartment, region, run, R0, Pop_Size)]
	DC <- DC[,list(run, R0=round(R0,5), Value=100 * V1 / Pop_Size)]
	XC <- X[Intervention==scenarios[scenario==scen,intervention] & Intervention.Shift==iv_delay,list(count.people, R0=round(R0,5), ABM_Value=100 * total.cases / count.people), keyby=max.new.cases]
	plot_comparison(XC, DC, ylim, ylab="Total Cases as % of Population")
}

# Min S
plot_final_S <- function(scen="Base", iv_delay=0, ylim=c(0,5)) {
	DC <- D[scenario==scen & compartment=="S" & region=="England", min(value),by=list(scenario, compartment, region, run, R0, Pop_Size)]
	DC <- DC[,list(run, R0=round(R0,5), Value=100 * V1 / Pop_Size)]
	XC <- X[Intervention==scenarios[scenario==scen,intervention] & Intervention.Shift==iv_delay,list(count.people, R0=round(R0,5), ABM_Value=100 * X.ds.cur.num..of.susceptible / count.people), keyby=max.new.cases]
	plot_comparison(XC, DC, ylim, ylab="Final Susceptibles as % of Population")
}

# Total Deaths
plot_total_deaths <- function(scen="Base", iv_delay=0, ylim=c(0,5)) {
	DC <- D[scenario==scen & compartment=="death_o" & region=="England", sum(value),by=list(scenario, compartment, region, run, R0, Pop_Size)]
	DC <- DC[,list(run, R0=round(R0,5), Value=100 * V1 / Pop_Size)]
	XC <- X[Intervention==scenarios[scenario==scen,intervention] & Intervention.Shift==iv_delay,list(count.people, R0=round(R0,5), ABM_Value=100 * total.deaths / count.people), keyby=max.new.cases]
	plot_comparison(XC, DC, ylim, ylab="Total Deaths as % of Population")
}

# (Total ICU can't be calculated from 1-dynamics)

# Peak ICU
plot_peak_icu <- function(scen="Base", iv_delay=0, ylim=c(0,5)) {
	DC <- D[scenario==scen & compartment=="icu_p" & region=="England", max(value),by=list(scenario, compartment, region, run, R0, Pop_Size)]
	DC <- DC[,list(run, R0=round(R0,5), Value=100 * V1 / Pop_Size)]
	XC <- X[Intervention==scenarios[scenario==scen,intervention] & Intervention.Shift==iv_delay,list(count.people, R0=round(R0,5), ABM_Value=100 * num.peak.icu.beds / count.people), keyby=max.new.cases]
	plot_comparison(XC, DC, ylim, ylab="Peak ICU Beds as % of Population")
}

# Peak non-ICU
plot_peak_nonicu <- function(scen="Base", iv_delay=0, ylim=c(0,5)) {
	DC <- D[scenario==scen & compartment=="nonicu_p" & region=="England", max(value),by=list(scenario, compartment, region, run, R0, Pop_Size)]
	DC <- DC[,list(run, R0=round(R0,5), Value=100 * V1 / Pop_Size)]
	XC <- X[Intervention==scenarios[scenario==scen,intervention] & Intervention.Shift==iv_delay,list(count.people, R0=round(R0,5), ABM_Value=100 * num.peak.non.icu.beds / count.people), keyby=max.new.cases]
	plot_comparison(XC, DC, ylim, ylab="Peak Non-ICU Beds as % of Population")
}

# Time to Peak Cases
plot_peak_t_cases <- function(scen="Base", iv_delay=0, ylim=c(0,20)) {
	DC <- D[compartment == "cases" & scenario==scen & region=="England", .(total_cases = sum(value)), by = .(t, run, R0)][,t[which.max(total_cases)],by=.(run,R0)];
	DC <- DC[,list(run, R0=round(R0,5), Value=V1/7)]
	XC <- X[Intervention==scenarios[scenario==scen,intervention] & Intervention.Shift==iv_delay,list(count.people, R0=round(R0,5), ABM_Value=num.weeks.to.peak.week.cases), keyby=max.new.cases]
	plot_comparison(XC, DC, ylim, ylab="Weeks to Peak in Cases")
}

# Plot comparable data sets
plot_comparison <- function(XC, DC, ylim=c(0,100), ylab="Output as % of Population") {
	Z <- XC[,list(Mean_Output=mean(ABM_Value), p0.025=quantile(ABM_Value, 0.025), p0.05=quantile(ABM_Value, 0.05), p0.10=quantile(ABM_Value, 0.10), p0.25=quantile(ABM_Value, 0.25), p0.50=quantile(ABM_Value, 0.50), p0.75=quantile(ABM_Value, 0.75), p0.90=quantile(ABM_Value, 0.90), p0.95=quantile(ABM_Value, 0.95), p0.975=quantile(ABM_Value, 0.975)), by=R0]
	xlim <- c(min(Z[,R0]), max(Z[,R0]))
	
	plot(x=xlim, y=ylim, type="n", xlab="R0", ylab=ylab, xlim=xlim, ylim=ylim)
	lines(x=Z[,R0], y=Z[,p0.90], type="p", pch=3, col="blue")
	lines(x=Z[,R0], y=Z[,p0.50], type="p", pch=3, col="green")
	lines(x=Z[,R0], y=Z[,p0.10], type="p", pch=3, col="red")
	lines(x=DC[,R0], y=DC[,Value], type="p", pch=4, col="black")
	legend("topright", legend=c("Reproduction 90%ile", "Reproduction 50%ile", "Reproduction 10%ile", "Target Model" ), pch=c(3,3,3,4), col=c("blue", "green", "red", "black"), inset = c(0.05, 0.05));
}

# Plot comparable data sets (NG's version)
plot_comparison <- function(XC, DC, ylim=c(0,100), ylab="Output as % of Population") {
	Z <- XC[,list(Mean_Output=mean(ABM_Value), Min_Output=min(ABM_Value), Max_Output=max(ABM_Value), p0.025=quantile(ABM_Value, 0.025), p0.05=quantile(ABM_Value, 0.05), p0.10=quantile(ABM_Value, 0.10), p0.25=quantile(ABM_Value, 0.25), p0.50=quantile(ABM_Value, 0.50), p0.75=quantile(ABM_Value, 0.75), p0.90=quantile(ABM_Value, 0.90), p0.95=quantile(ABM_Value, 0.95), p0.975=quantile(ABM_Value, 0.975)), by=R0]
	ggplot(Z) +
	theme_light(base_size = 14) +
	theme(plot.title = element_text(size=14), axis.title.x = element_text(size=14), axis.title.y = element_text(size=14)) +
	labs(x=expression(R[0]), y=ylab) +
	ylim(ylim) +
	geom_point(aes(x=DC[,R0], y=DC[,Value]), shape=16, colour="blue") +
	geom_point(aes(x=Z[,R0], y = Z[,p0.50]), shape=3) +
	geom_ribbon(aes(x=Z[,R0], ymin = Z[,p0.25], ymax = Z[,p0.75]), fill = "grey50", alpha=0.2) +
#	geom_ribbon(aes(x=Z[,R0], ymin = Z[,p0.10], ymax = Z[,p0.90]), fill = "grey55", alpha=0.2)
#	geom_ribbon(aes(x=Z[,R0], ymin = Z[,p0.05], ymax = Z[,p0.95]), fill = "grey60", alpha=0.2) +
	geom_ribbon(aes(x=Z[,R0], ymin = Z[,p0.025], ymax = Z[,p0.975]), fill = "grey60", alpha=0.2)
#	geom_ribbon(aes(x=Z[,R0], ymin = Z[,Min_Output], ymax = Z[,Max_Output]), fill = "grey70", alpha=0.2)
}

########################################################################

# Plots for Paper R:

plot_peak_cases(scen="Base", iv_delay=0, ylim=c(0,4))
plot_peak_cases(scen="School Closures", iv_delay=0, ylim=c(0,4))
plot_peak_cases(scen="Combination", iv_delay=0, ylim=c(0,4))

plot_total_cases(scen="Base", iv_delay=0, ylim=c(0,100))
plot_total_cases(scen="Combination", iv_delay=0, ylim=c(0,100))

plot_total_deaths(scen="Base", iv_delay=0, ylim=c(0,1.25))
plot_total_deaths(scen="Combination", iv_delay=0, ylim=c(0,1.25))

plot_peak_icu(scen="Base", iv_delay=0, ylim=c(0,2))
plot_peak_icu(scen="Combination", iv_delay=0, ylim=c(0,2))

plot_peak_t_cases(scen="Base", iv_delay=0, ylim=c(0,60))
plot_peak_t_cases(scen="Combination", iv_delay=0, ylim=c(0,60))

########################################################################
########################################################################

# Alternative plots, not in paper R

# Other scenarios:

plot_peak_cases(scen="Social Distancing", iv_delay=0, ylim=c(0,5))

plot_total_cases(scen="School Closures", iv_delay=0, ylim=c(0,100))
plot_total_cases(scen="Social Distancing", iv_delay=0, ylim=c(0,100))
#plot_total_cases(scen="Combination", iv_delay=14, ylim=c(0,100))

# Final Susceptibles (not used in paper)
plot_final_S(scen="Base", iv_delay=0, ylim=c(0,100))
plot_final_S(scen="School Closures", iv_delay=0, ylim=c(0,100))
plot_final_S(scen="Social Distancing", iv_delay=0, ylim=c(0,100))
plot_final_S(scen="Combination", iv_delay=0, ylim=c(0,100))

# Peak Non-ICU beds (only used it to diagnose problem in ICU)
plot_peak_nonicu(scen="Base", iv_delay=0, ylim=c(0,5))

########################################################################
########################################################################

# Repetitions chart

colnames(X)

pop <- 39697

plot_reps_peak_cases <- function(lo_R0=3.0, hi_R0=3.05, scen="Base") {
	x <- X[Intervention==scen & R0>lo_R0 & R0<hi_R0, .(R0=R0, run_number=X.run.number., ABM_value=max.new.cases)]
	DC <- D[scenario==scen & R0>lo_R0 & R0<hi_R0 & compartment=="cases" & region=="England", max(value),by=list(scenario, compartment, region, run, R0, Pop_Size)]
	DC <- DC[,list(run, R0=round(R0,5), Value=100 * V1 / pop)]
	Y <- repetitions_data(x)
	plot_repetitions(Y, DC)
}


# Check there's no correlation between value and data order.
#plot(x=x[,run_number],y=x[,ABM_value])

repetitions_data <- function(x) {
	n <- dim(x)[1]
	Y <- data.table(n=integer(), cum_mean=double(), cum_se=double(), cum_min=double(), cum_max=double(), cum_median=double(), p0.10=double(), p0.90=double(), p0.05=double(), p0.95=double(), p0.25=double(), p0.75=double(), p0.025=double(), p0.975=double())
	for (i in 2:n) {
		Y <- rbind(Y, x[1:i, .(
			n=i, 
			cum_mean=mean(100*ABM_value/pop), 
			cum_se=qt(p=0.975,df=i-1)*sd(100*ABM_value/pop)/sqrt(i), 
			cum_min=min(100*ABM_value/pop), 
			cum_max=max(100*ABM_value/pop),
			cum_median=median(100*ABM_value/pop), 
			p0.10=quantile(100*ABM_value/pop,0.10), 
			p0.90=quantile(100*ABM_value/pop,0.90), 
			p0.05=quantile(100*ABM_value/pop,0.05), 
			p0.95=quantile(100*ABM_value/pop,0.95), 
			p0.25=quantile(100*ABM_value/pop,0.25), 
			p0.75=quantile(100*ABM_value/pop,0.75), 
			p0.025=quantile(100*ABM_value/pop,0.025), 
			p0.975=quantile(100*ABM_value/pop,0.975)
		)])
	}
	return(Y);
}

# Median and min, max, and percentiles
plot_repetitions <- function(Y, DC) {
	ggplot(Y) +
	theme_light(base_size = 14) +
	theme(plot.title = element_text(size=14), axis.title.x = element_text(size=14), axis.title.y = element_text(size=14)) +
	labs(x=expression(Repetitions), y="% of Population", title="Cumulative median with min, max, \nand 2.5/25/75/97.5 percentiles versus Repetitions \nRed line shows value from single TM run") +
#	ylim(c(0,25)) +
	geom_line(aes(x=Y[,n], y=Y[,cum_median]), shape=16, colour="blue") +
	geom_line(aes(x=Y[,n], y=DC[, Value]), shape=16, colour="red") +
	geom_ribbon(aes(x=Y[,n], ymin = Y[,cum_min], ymax = Y[,cum_max]), fill = "grey60", alpha=0.2) + 
	geom_ribbon(aes(x=Y[,n], ymin = Y[,p0.025], ymax = Y[,p0.975]), fill = "grey50", alpha=0.2) + 
	geom_ribbon(aes(x=Y[,n], ymin = Y[,p0.25], ymax = Y[,p0.75]), fill = "grey30", alpha=0.2)
}

# Mean and 95% CI
plot_repetitions_mean <- function(Y, DC) {
	ggplot(Y) +
	theme_light(base_size = 14) +
	theme(plot.title = element_text(size=14), axis.title.x = element_text(size=14), axis.title.y = element_text(size=14)) +
	labs(x=expression(Repetitions), y="% of Population", title="Cumulative mean with 95% C.I. versus Repetitions\nRed line shows TM value") +
#	ylim(c(0,25)) +
	geom_line(aes(x=Y[,n], y=Y[,cum_mean]), shape=16, colour="blue") +
	geom_line(aes(x=Y[,n], y=DC[, Value]), shape=16, colour="red") +
	geom_ribbon(aes(x=Y[,n], ymin = Y[,cum_mean-cum_se], ymax = Y[,cum_mean+cum_se]), fill = "grey50", alpha=0.2)
}

lo_R0 <- 2.65
hi_R0 <- 2.7
#lo_R0 <- 3.0
#hi_R0 <- 3.05
#lo_R0 <- 3.2
#hi_R0 <- 3.25
scen <- "Base"
scen <- "Combination"

plot_reps_peak_cases(2.65, 2.70, "Base")
plot_reps_peak_cases(3.0, 3.05, "Base")
plot_reps_peak_cases(3.2, 3.25, "Base")

plot_reps_peak_cases(2.65, 2.70, "Combination")
plot_reps_peak_cases(3.0, 3.05, "Combination")
plot_reps_peak_cases(3.2, 3.25, "Combination")

########################################################################
########################################################################

# TM points outside the RM prediction interval

points_outside_interval <- function(XC, DC, lp = 0.025, up = 0.975) {
	Z <- XC[, list(lq=quantile(ABM_Value, lp), uq=quantile(ABM_Value, up)), by=R0]
	ZD <- merge(Z, DC)
	n <- sum(ZD[, ((Value<lq) | (Value>uq)) ])
	return(n)
}

# Time to peak cases: Base
DC <- D[scenario=="Base" & compartment=="cases" & region=="England", .(total_cases = sum(value)), by = .(t, run, R0)][,t[which.max(total_cases)],by=.(run,R0)];
DC <- DC[,list(run, R0=round(R0,5), Value=V1/7)]
XC <- X[Intervention==scenarios[scenario=="Base",intervention] & Intervention.Shift==0,list(count.people, R0=round(R0,5), ABM_Value=num.weeks.to.peak.week.cases), keyby=max.new.cases]

# Time to peak cases: Combination
DC <- D[scenario=="Combination" & compartment=="cases" & region=="England", .(total_cases = sum(value)), by = .(t, run, R0)][,t[which.max(total_cases)],by=.(run,R0)];
DC <- DC[,list(run, R0=round(R0,5), Value=V1/7)]
XC <- X[Intervention==scenarios[scenario=="Combination",intervention] & Intervention.Shift==0,list(count.people, R0=round(R0,5), ABM_Value=num.weeks.to.peak.week.cases), keyby=max.new.cases]

# Total cases: Base
DC <- D[scenario=="Base" & compartment=="cases" & region=="England", sum(value),by=list(scenario, compartment, region, run, R0, Pop_Size)]
DC <- DC[,list(run, R0=round(R0,5), Value=100 * V1 / Pop_Size)]
XC <- X[Intervention==scenarios[scenario=="Base",intervention] & Intervention.Shift==0,list(count.people, R0=round(R0,5), ABM_Value=100 * total.cases / count.people), keyby=max.new.cases]

# Peak cases: Base
DC <- D[scenario=="Base" & compartment=="cases" & region=="England", max(value),by=list(scenario, compartment, region, run, R0, Pop_Size)]
DC <- DC[,list(run, R0=round(R0,5), Value=100 * V1 / Pop_Size)]
XC <- X[Intervention==scenarios[scenario=="Base",intervention] & Intervention.Shift==0,list(count.people, R0=round(R0,5), ABM_Value=100 * max.new.cases / count.people), keyby=max.new.cases]

# Peak cases: Combination
DC <- D[scenario=="Combination" & compartment=="cases" & region=="England", max(value),by=list(scenario, compartment, region, run, R0, Pop_Size)]
DC <- DC[,list(run, R0=round(R0,5), Value=100 * V1 / Pop_Size)]
XC <- X[Intervention==scenarios[scenario=="Combination",intervention] & Intervention.Shift==0,list(count.people, R0=round(R0,5), ABM_Value=100 * max.new.cases / count.people), keyby=max.new.cases]

points_outside_interval(XC, DC, lp=0.025, up=0.975)
points_outside_interval(XC, DC, lp=0.05, up=0.95)
points_outside_interval(XC, DC, lp=0.25, up=0.75)

# Likelihood of getting this many points outside a 95% interval
n <- points_outside_interval(XC, DC, lp=0.025, up=0.975)
dbinom(n, size=50, prob=0.05) # Prob(Get n points outside interval)
1 - pbinom(n-1, size=50, prob=0.05) # Prob(Get >n-1 points outside interval)

# A reminder of what the Binomial Distribution looks like
plot(0:50, dbinom(0:50, size=50, prob=0.05),type='h', xlab="Successes", ylab="Probability", main="Binomial Distribution (n=50, p=0.05)", lwd=3)

# NB: Our prediction intervals (percentiles) are only estimates.
# The true percentiles from the population of RM runs may be different.
# We could "acquire" a sampling distribution for each R0 value's percentiles
# using bootstrap resampling.


########################################################################
########################################################################
