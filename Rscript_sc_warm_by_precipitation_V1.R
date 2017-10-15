###########################################################
# Data analysis for enzymes activities at Haibei station
#
###########################################################

# Authors: Xin Jing
# First created: May 1, 2016

# Please contact Xin Jing (jingxin0123@gmail.com) before utilizing any part of the script
###########################################################
# load source data
source("./R/Rscript_sc_load_source_data.R")

# select subset data
enz.dat <- subset(enz.dat, exp == 1)
enz.dat <- droplevels(enz.dat)
enz.dat <- enz.dat %>%
  select(block, treat, BG, CB, NAG, LAP, PHOS, POX, PER,
         BGspe, CBspe, NAGspe, LAPspe, PHOSspe, POXspe, PERspe,
         CNratio, CPratio, NPratio)
str(enz.dat)

# reshape data
df <- melt(enz.dat, id.vars = c("block", "treat"))
str(df)
df$block <- factor(df$block, ordered = TRUE,
                   levels = c("A", "D", "E", "F"))
df$treat <- factor(df$treat, ordered = TRUE,
                   levels = c("control", "NWDP", "NWIP", "SWNP", "SWDP", "SWIP", "WWNP"))
levels(df$variable) <- c("BG", "CB", "NAG", "LAP", "PHOS", "POX", "PER",
                         "BGmbc", "CBmbc", "NAGmbc", "LAPmbc", "APmbc", "POXmbc", "PERmbc",
                         "C:N ratio", "C:P ratio", "N:P ratio")

# t test
nenz <- levels(df$variable)
ntreat <- c("NWDP", "NWIP", "SWNP", "SWDP", "SWIP")
# ntreat <- c("NWDP", "NWIP", "SWNP", "SWDP", "SWIP", "WWNP")

for (i in nenz) {
  out = subset(df, variable == i)
  g1 <- subset(out, treat == 'control')
  for (j in ntreat) {
    g2 <- subset(out, treat == j)
    dif <- (g2$value - g1$value)/g1$value*100
    # cat(i, j, dif, "\n")
    ttest <- t.test(dif, na.action = na.omit)
    t.val <- ttest$statistic
    p.val <- ttest$p.value
    meandif <- mean((g2$value - g1$value)/g1$value*100, na.rm = TRUE)
    yminus <- ttest$conf.int[1]
    yplus <- ttest$conf.int[2]
    cat(i, j, t.val, p.val, meandif, yminus, yplus, "\n") # print p.value, mean and se of the relative change
  }
}
# store the results in ./data/data_processing/enz_dat_short_comm_WP_V1.csv

###########################################################
# End of the script
###########################################################