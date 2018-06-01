## POAP Experiment Plots
source("plottingUtils.R")

dat = read_csv("../Results/portExp__8675309000.csv_full_200.csv")
dat <- clean_data(dat)
dat.sum <- summarize_dat(dat)


##### 
## The Tau0 plot
#####
g <- dat.sum %>% filter(isBayes, Method != "SAA") %>%
  ggplot(aes(n, avgTau0, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=downTau0, ymax=upTau0), position=pd)

#Create a graph of Tau0
g <- make_pretty(g) + scale_y_continuous() + 
  ylab("Tau0") + ylim(0, 10) +
  theme(legend.position=c(.8, .7))

ggsave("../Results/Tau0Plot_temp.pdf", 
         g, width=3.25, height=3.25, units="in")

#########
## Regularization Plot
########
g <- dat.sum %>% filter(!isBayes, ! Method %in% c("SAA", "FullInfo")) %>%
  ggplot(aes(n, avgTau0, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=downTau0, ymax=upTau0), position=pd)

#Create a graph of Tau0
g <- make_pretty(g) + scale_y_continuous() + 
  ylab("Gamma") + 
  theme(legend.position=c(.8, .6))

ggsave("../Results/GammaPlot_temp.pdf", 
       g, width=3.25, height=3.25, units="in")



g <- make_pretty(g) + scale_y_continuous()
  theme(legend.direction = "horizontal",
        legend.position=c(.5, .1),
        legend.justification = "center") +
  guides(shape=guide_legend(nrow=2, byrow=TRUE)) + 
  ylab("(%) of Full-Info") + 
  scale_y_continuous(labels=scales::percent, 
                     limits = c(.70, .97),
                     breaks=seq(.75, .95, .1)) + 
  scale_x_log10(labels=scales::comma) 

g
