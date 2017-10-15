###########################################################
# Data analysis for enzymes activities at Haibei station
#
#
# Author: Xin Jing
# First created: May 1, 2016
#
# Contact Xin Jing (jingxin0123@gmail.com) 
###########################################################
rm(list = ls())

# load library
library(ggplot2) # Calls: ggplot
library(grid) # Calls: grid.newpage, grid.draw

# load data
enz.dat <- read.csv("./data/data_processing/enz_dat_short_comm_V3.csv")

# drop levels
enz.dat <- subset(enz.dat, treat != "WWNP")
enz.dat <- droplevels(enz.dat)

# reorder the data
levels(enz.dat$exp)
levels(enz.dat$enz.type)
enz.dat$enzyme <- factor(enz.dat$enzyme, ordered = TRUE,
                         levels = c("BG", "CB", "NAG", "LAP", "PHOS", "POX", "PER",
                                    "BGmbc", "CBmbc", "NAGmbc", "LAPmbc", "APmbc", "POXmbc", "PERmbc",
                                    "C:N", "C:P", "N:P"))
levels(enz.dat$enzyme) <- c("BG", "CB", "NAG", "LAP", "AP", "POX", "PER",
                            "BGmbc", "CBmbc", "NAGmbc", "LAPmbc", "APmbc", "POXmbc", "PERmbc",
                            "C:N ratio", "C:P ratio", "N:P ratio")
enz.dat$treat.type <- factor(enz.dat$treat.type, ordered = TRUE,
                             levels = c("WP", "NutA"))
levels(enz.dat$treat.type) <- c("Experiment 1", "Experiment 2")
enz.dat$treat <- factor(enz.dat$treat, ordered = TRUE,
                        levels = c("SWIP", "SWDP", "SWNP", "NWIP","NWDP",
                                   "N100P50", "P50", "N100", "N50", "N25"))
levels(enz.dat$treat) <- c("T+ P+", "T+ P-", "T+ P0", "T0 P+", "T0 P-", 
                           "N100 P50", "P50", "N100", "N50", "N25")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# enzyme activites
theme_set(theme_minimal(base_size = 16))
p1 <- ggplot(data = subset(enz.dat, enz.type == "enz" 
                           & treat.type == "Experiment 1")) + 
  geom_point(aes(x = treat, y = mean), shape = 1, size = 1.8) +
  geom_errorbar(aes(x = treat, ymin = yminus, ymax = yplus, 
                    width = 0.08), size = 0.25) +
  coord_cartesian(ylim = c(-200, 200)) +
  scale_y_continuous(breaks = c(-150, 0, 150)) +
  facet_wrap(facets = "enzyme", nrow = 1) +
  coord_flip() + 
  geom_hline(yintercept = 0, lty = 2, size = 0.2) +
  # ylim(-165, 165) +
  xlab("Experiment 1") +
  ylab("") +
  theme(axis.text.x = element_text(size = 9))

p2 <- ggplot(data = subset(enz.dat, enz.type == "enz" 
                           & treat.type == "Experiment 2")) +
  geom_point(aes(x = treat, y = mean), shape = 1, size = 1.8) +
  geom_errorbar(aes(x = treat, ymin = yminus, ymax = yplus, 
                    width = 0.08), size = 0.25) +
  coord_cartesian(ylim = c(-200, 200)) +
  scale_y_continuous(breaks = c(-50, 0, 50)) +
  facet_wrap(facets = "enzyme", nrow = 1) +
  coord_flip() + geom_hline(yintercept = 0, lty = 2, size = 0.2) +
  # ylim(-165, 165) +
  xlab("Experiment 2") +
  ylab("\nRelative change (% mean ± 95%CI)") +
  theme(axis.text.x = element_text(size = 9))

grid.newpage()
pdf("./outs/enz_dat_short_comm_enz.pdf", width = 8, height = 7)
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# specific activities
theme_set(theme_minimal(base_size = 16))
p3 <- ggplot(data = subset(enz.dat, enz.type == "spe" 
                           & treat.type == "Experiment 1")) +
  geom_point(aes(x = treat, y = mean), shape = 1, size = 1.8) +
  geom_errorbar(aes(x = treat, ymin = yminus, ymax = yplus, 
                    width = 0.08), size = 0.25) +
  coord_cartesian(ylim = c(-200, 200)) +
  scale_y_continuous(breaks = c(-200, 0, 200)) +
  facet_wrap(facets = "enzyme", nrow = 1) +
  coord_flip() + geom_hline(yintercept = 0, lty = 2, size = 0.2) +
  # ylim(-150, 230) +
  xlab("Experiment 1") +
  ylab("") +
  theme(axis.text.x = element_text(size = 8))

p4 <- ggplot(data = subset(enz.dat, enz.type == "spe" 
                           & treat.type == "Experiment 2")) + 
  geom_point(aes(x = treat, y = mean), shape = 1, size = 1.8) +
  geom_errorbar(aes(x = treat, ymin = yminus, ymax = yplus, 
                    width = 0.08), size = 0.25) +
  coord_cartesian(ylim = c(-200, 200)) +
  scale_y_continuous(breaks = c(-150, 0, 150)) +
  facet_wrap(facets = "enzyme", nrow = 1) +
  coord_flip() + geom_hline(yintercept = 0, lty = 2, size = 0.2) +
  # ylim(-150, 230) +
  xlab("Experiment 2") +
  ylab("\nRelative change (% mean ± 95%CI)") +
  theme(axis.text.x = element_text(size = 8))

grid.newpage()
pdf("./outs/enz_dat_short_comm_spe.pdf", width = 8, height = 7)
grid.draw(rbind(ggplotGrob(p3), ggplotGrob(p4), size = "last"))
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# stoichiometry
theme_set(theme_minimal(base_size = 16))
p5 <- ggplot(data = subset(enz.dat, enz.type == "sto" 
                           & treat.type == "Experiment 1")) + 
  geom_point(aes(x = treat, y = mean), shape = 1, size = 1.8) +
  geom_errorbar(aes(x = treat, ymin = yminus, ymax = yplus, 
                    width = 0.08), size = 0.25) +
  facet_wrap(facets = "enzyme", nrow = 1) +
  coord_flip() + geom_hline(yintercept = 0, lty = 2, size = 0.2) +
  ylim(-20, 18) +
  xlab("Experiment 1") +
  ylab("") +
  theme(axis.text.x = element_text(size = 9))

p6 <- ggplot(data = subset(enz.dat, enz.type == "sto" 
                           & treat.type == "Experiment 2")) + 
  geom_point(aes(x = treat, y = mean), shape = 1, size = 1.8) +
  geom_errorbar(aes(x = treat, ymin = yminus, ymax = yplus, 
                    width = 0.08), size = 0.25) +
  facet_wrap(facets = "enzyme", nrow = 1) +
  coord_flip() + geom_hline(yintercept = 0, lty = 2, size = 0.2) +
  ylim(-20, 18) +
  xlab("Experiment 2") +
  ylab("\nRelative change (% mean ± 95%CI)") +
  theme(axis.text.x = element_text(size = 9))

grid.newpage()
pdf("./outs/enz_dat_short_comm_sto.pdf", width = 8, height = 7)
grid.draw(rbind(ggplotGrob(p5), ggplotGrob(p6), size = "last"))
dev.off()

###########################################################
# End of the script
###########################################################