### Three part CLT Analysis and Plots
source("plottingUtils.R")
library(stringr)

dist_type = "exponential"
dist_type = "t"
dist_type = "pareto"
dist_type = "uniform"
dist_type = "bernoulli"
 
dat_path = str_c("../Results/POAPCLT_plot___", 
              dist_type, "_8675309.csv_full_200.csv")
dat = read_csv(dat_path)

dat <- clean_data(dat)
dat.sum <- summarize_clt_dat(dat)

pd = position_dodge(.5)
#Just the Bayes ones
g <- filter(dat.sum, isBayes) %>%
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
g <- g + theme(legend.position = c(.5, .9))
ggsave("student_t_emp_bayes.pdf", g, width=4, height=4, units="in")

ggsave("uniform_emp_bayes.pdf", g, width=4, height=4, units="in")


###

spath = str_c("../../../../OR Submission_1/Figures/POAPCLT_Bayes_",
              dist_type, ".pdf")

##fiddle with the limits to amke the legend fit
#exponential
g <- g + scale_y_continuous(labels=scales::percent,
                            limits= c(.78, .9))

## uniform
g <- g + scale_y_continuous(labels=scales::percent, 
                            limits = c(.77, .9))

## pareto
g <- g + scale_y_continuous(labels=scales::percent, 
                            limits=c( .8, .93))

ggsave(spath, 
       g, width=3.25, height=3.25, units="in")


#Just the REg ones
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

