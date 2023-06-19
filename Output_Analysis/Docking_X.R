##############################################################################
# Loads .csv files containing output from NetLogo's BehaviorSpace experiment
# then combines them into one data table and outputs it to the file X.qs
##############################################################################

library(data.table)
library(qs)

setwd("Z:\\git\\CovidModelAgentBasedReproduction\\Output_Analysis")
setwd("D:\\Share\\git\\CovidModelAgentBasedReproduction\\Output_Analysis")

# Convert BehaviorSpace output, and add a field to identify the source later.
X <- cbind(as.data.table(read.csv("MyPC-experiment-LSHTM-Docking-50-R0s-Rutland.csv", header=TRUE, sep=",", skip=6)), source="L800_230601")

# Appending more files.
#X <- rbind(X,X2)
#rm(X2)

dim(X)

#colnames(X)
#head(X)

# .qs file will be much faster to load than .csv, and more compact.
qsave(X, "X.qs")

# Check what we have now.
unique(X[,R0])
unique(X[,"Intervention"])
test_R0 <- unique(X[,R0])[1]
test_R0
dim(X[R0==test_R0 & Intervention=="Base"])[1]
test_R0 <- unique(X[,R0])[50]
test_R0
dim(X[R0==test_R0 & Intervention=="Base"])[1]


