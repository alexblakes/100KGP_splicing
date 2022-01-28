library(ggplot2)
library(ggseqlogo)
library(dplyr)

setwd("~/re_gecip/machine_learning/AlexBlakes/near_splice/write_up/paper/code/branchpoints/plotting/")

df <- read.csv("../stats/pwm.tsv", sep="\t")

dfa <- select(df[df$region == "Acceptor", ], -1:-2)
dfd <- select(df[df$region == "Donor", ], -1:-2)
dfb <- select(df[df$region == "Branchpoint", ], -1:-2)

a <- as.matrix(t(dfa))
d <- as.matrix(t(dfd))
b <- as.matrix(t(dfb))

pwms <- list(b,a,d)

ggplot() + 
  geom_logo(pwms) +
  theme_logo() +
  facet_wrap(~seq_group, ncol=1)
  
