### Large Sample Analysis
source("plottingUtils.R")
library(stringr)

dist_type = "exponential"
dist_type = "t"
dist_type = "pareto"
dist_type = "uniform"
dist_type = "bernoulli"

dat_path = str_c("../Results/POAPLargeSample___", 
                 dist_type, "_8675309.csv_full_200.csv")
dat = read_csv(dat_path)

dat <- clean_data(dat)
dat.sum <- summarize_clt_dat(dat)

pd = position_dodge(.5)
#Just the Bayes ones
g <- filter(dat.sum, isBayes, 
                      !Method %in% c("HO_EB", "LOO_EB")) %>%
  ggplot(aes(N, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd)

g <- make_pretty(g, TRUE) + 
  scale_x_continuous() + 
  scale_y_continuous(labels=scales::percent) +
  theme(legend.justification = "center", 
        legend.position=c(.5, .1), 
        legend.direction="horizontal") + 
  ylab("(%) of Full-Info") + xlab("S")
g

###
#
ggsave(str_c("../../../../MS Revision/Figures/LS_", dist_type, "_plot.pdf"), g, width=3.25, height=3.25, units="in")



###Just the Reg ones
g <- filter(dat.sum, isReg) %>%
  ggplot(aes(N, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd)

g <- make_pretty(g, FALSE) + 
  theme(legend.justification = "center", 
        legend.position= c(.5, .1), 
        legend.direction="horizontal") + 
  scale_x_continuous() +
  scale_y_continuous(labels=scales::percent) +
  ylab("(%) of Full-Info") + xlab("S")
g

ggsave(str_c("../../../../MS Revision/Figures/LSReg_", dist_type, "_plot.pdf"), g, width=3.25, height=3.25, units="in")




#### 
#something seems funny with regopt
#try plotting the gammas...
filter(dat.sum, Method %in% c("OracleReg", "SteinReg", "K5_Reg")) %>%
  ggplot(aes(N, avgTau0, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=avgTau0 - stdTau0, ymax=avgTau0 + stdTau0), position=pd)





spath = str_c("../../../EmpBayes/POAPCLT_REG_Pres_", dist_type, ".pdf")

spath = str_c("../../../../OR Submission_1/Figures/POAPCLT_Reg_",
              dist_type, ".pdf")

#Fiddle with limits to get legend to fit
# #exponenntial
# g <- g + scale_y_continuous(labels=scales::percent,
#                             limits= c(.78, .95))

#uniform
g <- g + scale_y_continuous(labels=scales::percent,
                            limits= c(.775, .93))
#pareto
g <- g + scale_y_continuous(labels=scales::percent,
                            limits= c(.8, .95))
g

ggsave(spath, 
       g, width=3.25, height=3.25, units="in")

ggsave(spath,
       g, width = 5.35, height=9.04, units="in")

