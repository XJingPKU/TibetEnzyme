###########################################################
# Meta analysis: environmental changes on soil enzyme activity
#
###########################################################

# Authors: Xin Jing
# First created: Oct. 29, 2017

# Please contact Xin Jing (jingxin0123@gmail.com) before utilizing any part of the script
###########################################################
options(useFancyQuotes = FALSE)
rm(list = ls())

# load library
library(plyr)
library(dplyr)
library(ggplot2)

# load source data
enz.meta <- read.csv("./data/data_processing/meta_data_enz_apline_grasslands.csv")

# data summary
str(enz.meta)
summary(enz.meta)
head(enz.meta)
names(enz.meta)

# clean data
enz.meta <- enz.meta[enz.meta$enz.func != "others", ]

# calculate lnRR and samlping variance for each lnRR
enz.meta$yi <- with(enz.meta, log(treat/CK))
enz.meta$vi <- with(enz.meta, (1/treat.N)*(treat.sd2/treat)^2 + (1/CK.N)*(CK.sd2/CK)^2)

# calculate raw lnRR based on several factors

ddply(enz.meta, c("Expt", "Expt.treat"), function(x) {
  df <- data.frame(value = mean(x[, "yi"]),
  n.obs = length(x[, "yi"]),
  lci = mean(x[, "yi"]) - qt(0.975, length(x[, "yi"]))*sd(x[, "yi"])/sqrt(length(x[, "yi"])),
  uci = mean(x[, "yi"]) + qt(0.975, length(x[, "yi"]))*sd(x[, "yi"])/sqrt(length(x[, "yi"]))
  )
})

df1 <- ddply(enz.meta, c("Expt", "Expt.treat", "enz.func"), function(x) {
  df <- data.frame(value = mean(x[, "yi"]),
                   n.obs = length(x[, "yi"]),
                   lci = mean(x[, "yi"]) - qt(0.975, length(x[, "yi"]))*sd(x[, "yi"])/sqrt(length(x[, "yi"])),
                   uci = mean(x[, "yi"]) + qt(0.975, length(x[, "yi"]))*sd(x[, "yi"])/sqrt(length(x[, "yi"]))
  )
})
df1 %>% 
  mutate(Expt = factor(Expt, levels = c("Experiment 1", "Experiment 2")),
         Expt.treat = factor(Expt.treat, levels = c("warm", "N", "P", "NP"),
                             labels = c("Warming", "N", "P", "NP")),
         enz.func = factor(enz.func, levels = c("P cycling", "N cycling", "C cycling"))) %>% 
  ggplot(aes(x = enz.func, y = value, ymin = lci, ymax = uci)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey70") +
  geom_point(size = 4.5, shape = 1) +
  geom_errorbar(width = 0.08, size = 0.25) +
  geom_text(aes(x = enz.func, y = value, 
                label = n.obs, hjust = 0.45, vjust = 2.3)) +
  facet_grid( ~ Expt.treat) +
  labs(x = "", y = "\nlnRR (mean ± 95% CI)") +
  theme_minimal(base_size = 16) +
  theme(strip.background = element_blank()) +
  coord_flip()
ggsave("./outs/meta_expt_enz_func.pdf", width = 6.8, height = 4.5)

###########################################################
# End of the script
###########################################################