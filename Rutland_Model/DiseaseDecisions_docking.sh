#!/bin/bash

echo "Running Netlogo experiment..."

/home/pi/Applications/Netlogo/NetLogo-6.1.1/netlogo-headless.sh \
	--model ./DiseaseDecisions.nlogo \
	--experiment experiment-LSHTM-Docking-50-R0s \
	--table ./L11-experiment-LSHTM-Docking-50-R0s-Rutland.csv \
	--threads 4

#cp ./CW-experiment-LSHTM-Docking-50-R0s-Rutland.csv ~/Shares/

echo "Done!"
