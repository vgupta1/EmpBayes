### Three part CLT Analysis and Plots
source("plottingUtils.R")
library(stringr)

dat = read_csv("../Results/POAPCLT_plot___8675309.csv_full_100.csv")
dist_type = "exponential"

dat <- clean_data(dat)
dat.sum <- summarize_clt_dat(dat)

pd = position_dodge(.5)
#Just the Bayes ones
g <- filter(dat.sum, isBayes) %>%
  ggplot(aes(N, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd)

g <- make_pretty(g) + 
  scale_x_continuous() + 
  scale_y_continuous(labels=scales::percent, limits=c(.775, .9)) +
  theme(legend.justification = "center", 
        legend.position=c(.5, .1), 
        legend.direction="horizontal") + 
  ylab("(%) of Full-Info") + xlab("S")
g

spath = str_c("../../../../OR Submission_1/Figures/POAPCLT_Bayes_",
              dist_type, ".pdf")
ggsave(spath, 
       g, width=3.25, height=3.25, units="in")


#Just the REg ones
g <- filter(dat.sum, isReg) %>%
  ggplot(aes(N, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd)

g <- make_pretty(g) + 
  theme(legend.justification = "center", 
        legend.position= c(.5, .1), 
        legend.direction="horizontal") + 
  scale_x_continuous() +
  scale_y_continuous(labels=scales::percent, limits=c(.75, .95)) +
  ylab("(%) of Full-Info") + xlab("S")
g

spath = str_c("../../../../OR Submission_1/Figures/POAPCLT_Reg_",
              dist_type, ".pdf")
ggsave(spath, 
       g, width=3.25, height=3.25, units="in")


