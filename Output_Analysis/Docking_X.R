##############################################################################
# Loads .csv files containing output from NetLogo's BehaviorSpace experiment
# then combines them into one data table and outputs it to the file X.qs
##############################################################################

library(data.table)
library(qs)

setwd("~\\DiseaseDecisions\\Docking\\All_Interventions")

# Use this if you ran from BehaviorSpace
X1 <- cbind(as.data.table(read.csv("DiseaseDecisions experiment-LSHTM-Docking-50-R0s-table-NG-20201021.csv", header=TRUE, sep=",", skip=6)), source="NG21")
X2 <- cbind(as.data.table(read.csv("ORS-experiment-LSHTM-Docking-50-R0s-Rutland.csv", header=TRUE, sep=",", skip=6)), source="ORS")
X3 <- cbind(as.data.table(read.csv("DiseaseDecisions experiment-LSHTM-Docking-50-R0s-table-NG-20201022.csv", header=TRUE, sep=",", skip=6)), source="NG22")
X4 <- cbind(as.data.table(read.csv("DiseaseDecisions experiment-LSHTM-Docking-50-R0s-table_NG1.csv", header=TRUE, sep=",", skip=6)), source="NG1")
X5 <- cbind(as.data.table(read.csv("DiseaseDecisions experiment-LSHTM-Docking-50-R0s-table_NG2.csv", header=TRUE, sep=",", skip=6)), source="NG2")
X6 <- cbind(as.data.table(read.csv("DiseaseDecisions experiment-LSHTM-Docking-50-R0s-table_NG3.csv", header=TRUE, sep=",", skip=6)), source="NG3")
X7 <- cbind(as.data.table(read.csv("L11a-experiment-LSHTM-Docking-50-R0s-Rutland.csv", header=TRUE, sep=",", skip=6)), source="L11a")
X8 <- cbind(as.data.table(read.csv("L900_DiseaseDecisions experiment-LSHTM-Docking-50-R0s-table.csv", header=TRUE, sep=",", skip=6)), source="L900")
X9 <- cbind(as.data.table(read.csv("ORSa-experiment-LSHTM-Docking-50-R0s-Rutland.csv", header=TRUE, sep=",", skip=6)), source="ORSa")
X10 <- cbind(as.data.table(read.csv("L8unf-experiment-LSHTM-Docking-50-R0s-Rutland.csv", header=TRUE, sep=",", skip=6)), source="L8")
X11 <- cbind(as.data.table(read.csv("L4unf-experiment-LSHTM-Docking-50-R0s-Rutland.csv", header=TRUE, sep=",", skip=6)), source="L4")
X12 <- cbind(as.data.table(read.csv("L8unf2-experiment-LSHTM-Docking-50-R0s-Rutland.csv", header=TRUE, sep=",", skip=6)), source="L8")
X13 <- cbind(as.data.table(read.csv("L6-experiment-LSHTM-Docking-50-R0s-Rutland.csv", header=TRUE, sep=",", skip=6)), source="L6")
X14 <- cbind(as.data.table(read.csv("L11b-experiment-LSHTM-Docking-50-R0s-Rutland.csv", header=TRUE, sep=",", skip=6)), source="L11b")
X15 <- cbind(as.data.table(read.csv("DiseaseDecisions experiment-LSHTM-Docking-50-R0s-table_NG4.csv", header=TRUE, sep=",", skip=6)), source="NG4")
X16 <- cbind(as.data.table(read.csv("L4-experiment-LSHTM-Docking-50-R0s-Rutland.csv", header=TRUE, sep=",", skip=6)), source="L4")
X17 <- cbind(as.data.table(read.csv("DiseaseDecisions experiment-LSHTM-Docking-50-R0s-table_NG5.csv", header=TRUE, sep=",", skip=6)), source="NG5")
X18 <- cbind(as.data.table(read.csv("L900b-experiment-LSHTM-Docking-50-R0s-Rutland.csv", header=TRUE, sep=",", skip=6)), source="L900b")



X <- X1
X <- rbind(X,X2)
X <- rbind(X,X3)
X <- rbind(X,X4)
X <- rbind(X,X5)
X <- rbind(X,X6)
X <- rbind(X,X7)
X <- rbind(X,X8)
X <- rbind(X,X9)
X <- rbind(X,X10)
X <- rbind(X,X11)
X <- rbind(X,X12)
X <- rbind(X,X13)
X <- rbind(X,X14)
X <- rbind(X,X15)
X <- rbind(X,X16)
X <- rbind(X,X17)
X <- rbind(X,X18)
rm(X1)
rm(X2)
rm(X3)
rm(X4)
rm(X5)
rm(X6)
rm(X7)
rm(X8)
rm(X9)
rm(X10)
rm(X11)
rm(X12)
rm(X13)
rm(X14)
rm(X15)
rm(X16)
rm(X17)
rm(X18)


dim(X)

#colnames(X)
#head(X)

qsave(X, "X.qs")

unique(X[,R0])
unique(X[,"Intervention"])
test_R0 <- unique(X[,R0])[1]
test_R0
dim(X[R0==test_R0 & Intervention=="Base"])[1]
test_R0 <- unique(X[,R0])[50]
test_R0
dim(X[R0==test_R0 & Intervention=="Base"])[1]


