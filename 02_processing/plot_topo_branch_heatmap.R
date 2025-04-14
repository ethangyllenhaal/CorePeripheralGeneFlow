##########
# By: Ethan Gyllenhaal
# Updated 10 April 2025
#
# R script used for making heatmaps of various simulated values
# First plots topology by parameter value, then plots branch length metrics
########

# load libraries
library(tidyverse)
library(gridExtra)

# set up directories
setwd("/path/to/dir/")
raw_input <- read.csv("coresim_output_branch.csv")

# looad in input 
input <- raw_input %>%
  # count is the sum of all topologies
  mutate(count = Anest + Bnest + ABnest + Cmono + Other + Non.monophyletic + Poor.resolution) %>%
  # nres is the sum of resolved topologies (note that Cmono is not always well resolved at the base)
  mutate(nres = Anest + Bnest + ABnest + Cmono + Other) %>%
  # not used in final figure, but the number of correct topological patterns
  mutate(true = case_when(history=="ABC" ~ Cmono/count,
                             history=="CAB" ~ ABnest/count,
                             history=="ACB" ~ Bnest/count)) %>%
  # make x and y names, don't think I used this variable for this project
  mutate(x = paste(disp, time, sep="_")) %>%
  mutate(y = paste(ne, int, sep="_")) %>%
  # output the summary category
  mutate(summary = case_when((Non.monophyletic + Poor.resolution)/count > 0.4 ~ "poor_resolution", # + Poor.resolution for bootstrap version
                             Anest/nres > 0.6 ~ "A_nested",
                             Bnest/nres > 0.6 ~ "B_nested",
                             ABnest/nres > 0.6 ~ "A_and_B_nested",
                             Cmono/nres > 0.6 ~ "C_monophyly", 
                             .default = "mixed")) %>%
  # output the proportion of trees that match the summary category
  mutate(summary_val = case_when((Non.monophyletic + Poor.resolution)/count > 0.4 ~ "", # + Poor.resolution for bootstrap version
                             Anest/nres > 0.6 ~ as.character(100*Anest/count),
                             Bnest/nres > 0.6 ~ as.character(100*Bnest/count),
                             ABnest/nres > 0.6 ~ as.character(100*ABnest/count),
                             Cmono/nres > 0.6 ~ as.character(100*Cmono/count), 
                             .default = "")) %>%
  # make branch length ratios
  # note that these are skipping dividing the branch length values by the count because it cancels out
  # each a/b/core value is the sum from the summary program, though, not mean
  rowwise() %>%
  mutate(core_ratio = core/total) %>%
  mutate(a_ratio = a/total) %>%
  mutate(b_ratio = b/total) %>%
  mutate(core_per = core/min(c(a,b))) %>%
  mutate(core_per_1 = if_else(core_per>1, 0, 1)) %>% # effectively a boolean to highlight values greater than 1 
  arrange(disp, ne, time, int) # arrange for plotting

# color assignments, note that A_nested is never recovered
color_assignments <- c("A_nested" = "lightblue3", "B_nested" = "green3", "A_and_B_nested" = "blue", "C_monophyly" = "orange", "poor_resolution" = "black")

# Heatmap assignment plots
# Used for Figure 2 and S1 (S1 is from an old version of phylogeny classification)
## NOTE ###
## The X axis for dispersal is actually 10x the number of migrants per generation
## For conveniently parsing files I took  out the 0.
# each one of these is the same, so only annotating the first
abc <- ggplot(filter(input, history=="ABC"), aes(x=factor(disp), y=factor(ne))) +
  geom_tile(aes(fill=summary), color = "white") + # tile plot of summary
  geom_text(aes(label=summary_val), size=3) + # percent of trees matching summary
  facet_wrap(vars(factor(int), factor(time)), nrow=3) + # wrapping
  scale_fill_manual(values = color_assignments) + # color by assignments
  # theme takes out facet wrap info
  # if you want to check this, run only code above this 
  # (I recommend you do this before interpretting)
  theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.background = element_blank(), panel.spacing = unit(0.2, "lines")) +
  ggtitle("Core last")

acb <- ggplot(filter(input, history=="ACB"), aes(x=factor(disp), y=factor(ne))) +
  geom_tile(aes(fill=summary), color = "white") + 
  facet_wrap(vars(factor(int), factor(time)), nrow=3) +
  geom_text(aes(label=summary_val), size=3) +
  scale_fill_manual(values = color_assignments) +
  theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.background = element_blank(), panel.spacing = unit(0.2, "lines")) +
  ggtitle("A first, then core")

cab <- ggplot(filter(input, history=="CAB"), aes(x=factor(disp), y=factor(ne))) +
  geom_tile(aes(fill=summary), color = "white") + 
  facet_wrap(vars(factor(int), factor(time)), nrow=3) +
  geom_text(aes(label=summary_val), size=3) +
  scale_fill_manual(values = color_assignments) +
  theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.background = element_blank(), panel.spacing = unit(0.2, "lines")) +
  ggtitle("Core first")

# plot all 3 together
grid.arrange(abc + theme(legend.position="none"), acb+ theme(legend.position="none"), cab+ theme(legend.position="none"), nrow=1)


# Branch length plot, Figure S2
# plots the ratio of mean core branch length (of the 4 core pops) to the minimum peripheral pop branch length
# all based on one sample per population
# As with above, only describing the first one
# no alpha for abc because never >1
acb_s2 <- ggplot(filter(input, history=="ACB"), aes(x=factor(disp), y=factor(ne))) + 
  # tile plot of ratio, add alpha aesthetic to highlight values greater than 1
  geom_tile(aes(fill=core_per, alpha=core_per_1), color = "white") + 
  # set alpha scale
  scale_alpha(range = c(0.5,1)) +
  # wrap as above
  facet_wrap(vars(factor(int), factor(time)), nrow=3) +
  # ensure same viridis scale for all 3
  scale_fill_viridis_c(limits = c(0.35,1.5)) +
  # as with above, theme call takes out facet wrap info
  theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.background = element_blank(), panel.spacing = unit(0.2, "lines")) +
  ggtitle("A first, then core")
abc_s2 <- ggplot(filter(input, history=="ABC"), aes(x=factor(disp), y=factor(ne))) + 
  geom_tile(aes(fill=core_per), color = "white") + 
  facet_wrap(vars(factor(int), factor(time)), nrow=3) +
  scale_fill_viridis_c(limits = c(0.35,1.5)) +
  theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.background = element_blank(), panel.spacing = unit(0.2, "lines")) +
  ggtitle("Core last")
cab_s2 <- ggplot(filter(input, history=="CAB"), aes(x=factor(disp), y=factor(ne))) + 
  geom_tile(aes(fill=core_per, alpha=core_per_1), color = "white") + 
  scale_alpha(range = c(0.5,1)) +
  facet_wrap(vars(factor(int), factor(time)), nrow=3) +
  scale_fill_viridis_c(limits = c(0.35,1.5)) +
  theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.background = element_blank(), panel.spacing = unit(0.2, "lines")) +
  ggtitle("Core first")

grid.arrange(abc_s2, acb_s2, cab_s2, nrow=3) # grab legend from here
grid.arrange(abc_s2 + theme(legend.position="none"), acb_s2 + theme(legend.position="none"), cab_s2 + theme(legend.position="none"), nrow=1)

# Bonus violin plots, included because this is where I made input for Figure 3A before sending it to Luke

# make input for these plots/for Luke's script
violin_input <- input %>%
  # select branch length columns
  select(disp, ne, time, int, a, b, core, total, name, history, core_per, count) %>%
  # get the mean of each pattern, DO NOT USE THE RAW VALUE, it is the sum of each replicate matching the pattern
  mutate(mean_a = a/count) %>%
  mutate(mean_b = b/count) %>%
  mutate(mean_core = core/count)

# output/input CSV to make sure it's usable
write.csv(violin_input, file="violin_full.csv", quote=F)
violin_input <- read.csv("violin_full.csv", header=T)

# quick violin plots, not in paper but in short:
# history (i.e., underlying tree) is color and divide across X axis
# each plot is the mean branch length for a, b, and core
# colors is not really used, but lets us split up gene flow levels (a lil ggplot cheat)
# NOTE: each point is the mean of a SET of parameter combinations
# it is NOT an individual simulation
# this aims to show the general pattern across simulations, likely a better match for variable empricial data
violin_a <- ggplot(data = violin_input, aes(x=history, fill=history, y=mean_a, colors=factor(disp*0.1))) +
  geom_violin(position = position_dodge(0.9)) +
  geom_beeswarm(dodge.width=0.9, corral="wrap", corral.width = 0.2, size=0.7, color="gray30", alpha=1) +
  stat_summary(fun = "mean", position = position_dodge(0.9), shape=16) + ylim(0.065,0.2)+
  scale_fill_manual(values =c("orange", "green3", "blue")) + ggtitle("A")
violin_a

violin_b <- ggplot(data = violin_input, aes(x=history, fill=history, y=mean_b, colors=factor(disp*0.1))) +
  geom_violin(position = position_dodge(0.9)) +
  geom_beeswarm(dodge.width=0.9, corral="wrap", corral.width = 0.2, size=0.7, color="gray30", alpha=1) +
  stat_summary(fun = "mean", position = position_dodge(0.9), shape=16) + ylim(0.065,0.2)+
  scale_fill_manual(values =c("orange", "green3", "blue")) + ggtitle("B")
violin_b

violin_core <- ggplot(data = violin_input, aes(x=history, fill=history, y=mean_core, colors=factor(disp*0.1))) +
  geom_violin(position = position_dodge(0.9)) +
  geom_beeswarm(dodge.width=0.9, corral="wrap", corral.width = 0.2, size=0.7, color="gray30", alpha=1) +
  stat_summary(fun = "mean", position = position_dodge(0.9), shape=16) + ylim(0.065,0.2) +
  scale_fill_manual(values =c("orange", "green3", "blue")) + ggtitle("C")
violin_core

grid.arrange(violin_a+ theme(legend.position="none"), violin_b+ theme(legend.position="none"), violin_core+ theme(legend.position="none"), nrow=1)
