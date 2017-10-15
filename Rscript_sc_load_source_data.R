###########################################################
# Data analysis for enzymes activities at Haibei station
#
# Authors: Xin Jing
# First created: May 1, 2016

# Contact Xin Jing (jingxin0123@gmail.com) 
###########################################################
# remove objects in global environment
rm(list = ls()) 

###########################################################
# load library
library(dplyr) # Calls: select
library(reshape2) # Calls: melt
library(ggplot2) # Calls: ggplot, ggsave

###########################################################
# read table
enz.dat <- read.csv("./data/data_processing/enz_dat_short_comm.csv")
str(enz.dat)
summary(enz.dat)

###########################################################
# remove outlier
enz.dat$PER[enz.dat$PER < 35] <- NA
enz.dat$MBC[enz.dat$MBC < 200] <- NA

###########################################################
# calculate specific enzyme activities
enz.dat$BGspe <- enz.dat$BG/enz.dat$MBC
enz.dat$CBspe <- enz.dat$CB/enz.dat$MBC
enz.dat$NAGspe <- enz.dat$NAG/enz.dat$MBC
enz.dat$PHOSspe <- enz.dat$PHOS/enz.dat$MBC
enz.dat$LAPspe <- enz.dat$LAP/enz.dat$MBC
enz.dat$POXspe <- enz.dat$POX/enz.dat$MBC
enz.dat$PERspe <- enz.dat$PER/enz.dat$MBC

###########################################################
# calculate ecoenzymatic ratios of C, N and P
enz.dat$LogC <- log(enz.dat$BG)
enz.dat$LogN <- log(enz.dat$NAG + enz.dat$LAP)
enz.dat$LogP <- log(enz.dat$PHOS)

enz.dat$CNratio <- enz.dat$LogC/enz.dat$LogN
enz.dat$CPratio <- enz.dat$LogC/enz.dat$LogP
enz.dat$NPratio <- enz.dat$LogN/enz.dat$LogP

###########################################################
# End of the script
###########################################################