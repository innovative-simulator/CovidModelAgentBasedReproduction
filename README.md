# A Covid Model's Agent-Based Reproduction

This is an agent-based simulation developed to reproduce the model of COVID-19 transmission and non-pharmaceutical interventions. The model being reproduced is that of Davies et al. (2020), and is available at

https://github.com/cmmid/covid-UK 

Our reproducing ABM is described in Chattoe-Brown et al. (2021, not peer reviewed), which discusses similarities and differences between the two models' mechanisms, and compares their behaviour.

All code on this git is (C) Christopher J. Watts (2021), and covered by GNU Public Licence 3.0.

## Subfolders in this git

### The ABM

* ____Rutland_Model____ contains the NetLogo model DiseaseDecisions.nlogo from which the simulation experiment was run. Also in this folder are all data files required as input to the model (in .csv format), and scripts for running the NetLogo experiment in headless mode (.sh for Linux systems, .bat for Windows). You may need to edit the pathnames in these for your own computer.

### Data analysis

## References

Chattoe-Brown, E., Gilbert, N., Robertson, D. A., Watts, C. (2021, not peer reviewed) "Reproduction as a Means of Evaluating Policy Models: A Case Study of a COVID-19 Simulation".

Davies, Nicholas G., Kucharski, Adam J., Eggo, Rosalind M., Gimma, Amy and Edmunds, W. John on behalf of the Centre for the Mathematical Modelling of Infectious Diseases COVID-19 Working Group (2020) ‘Effects of Non-Pharmaceutical Interventions on COVID-19 Cases, Deaths, and Demand for Hospital Services in the UK: A Modelling Study’, Lancet Public Health, 5(7), July, pp. E375-E385. doi:10.1016/S2468-2667(20)30133-X

