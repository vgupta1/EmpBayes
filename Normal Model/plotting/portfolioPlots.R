## POAP Experiment Plots
source("plottingUtils.R")

#old run with bug fix
dat = read_csv("../Results/portExpOrig2_8675309.csv_full_200.csv")

#New run for paper with larger Gammas
dat = read_csv("../Results/portExp2Big100__8675309.csv_full_200.csv")

###
dat <- clean_data(dat)
dat.sum <- summarize_dat(dat)
dat.sum2 <- summarize_dat(dat, FALSE)

## plot for main text
g <- dat.sum %>% filter(isBayes, 
                        !Method %in%c("EB_MM", "HO_EB", "LOO_EB") ) %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd)

g <- make_pretty(g, TRUE) + 
    theme(legend.direction = "horizontal",
        legend.position=c(.5, .1),
        legend.justification = "center") +
  guides(shape=guide_legend(nrow=2, byrow=TRUE)) + 
  ylab("(%) of Full-Info") + 
  scale_y_continuous(labels=scales::percent, 
                     limits = c(.70, .97),
                    breaks=seq(.75, .95, .1)) + 
  scale_x_log10(labels=scales::comma) +

ggsave("../../../../MS Revision/Figures/portPerfBayesb.pdf",
       g, width=3.25, height=3.25, units="in")

## The absolute performance version
g <- dat.sum2 %>% filter(isBayes, 
                         !Method %in%c("EB_MM", "HO_EB", "LOO_EB") ) %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd)

g <- make_pretty(g, TRUE) + 
  theme(legend.direction = "horizontal",
        legend.position=c(.5, .1),
        legend.justification = "center") +
  guides(shape=guide_legend(nrow=2, byrow=TRUE)) + 
  ylab("(%) of Full-Info") + 
  scale_y_continuous() + 
  scale_x_log10(labels=scales::comma) 

g <- g + scale_y_continuous() + 
  ylab("Absolute Performance")
ggsave("../../../../MS Revision/Figures/portPerfBayes_Absb.pdf",
       g, width=3.25, height=3.25, units="in")

################
#Std error version for Reviewer
g <- dat.sum %>% filter(isBayes, 
                        !Method %in%c("EB_MM", "HO_EB", "LOO_EB") ) %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=avg - std/sqrt(200), ymax=avg+ std/sqrt(200)), position=pd)

g <- make_pretty(g, TRUE) + 
  theme(legend.direction = "horizontal",
        legend.position=c(.5, .1),
        legend.justification = "center") +
  guides(shape=guide_legend(nrow=2, byrow=TRUE)) + 
  ylab("(%) of Full-Info") + 
  scale_y_continuous(labels=scales::percent, 
                     limits = c(.70, .97),
                     breaks=seq(.75, .95, .1)) + 
  scale_x_log10(labels=scales::comma)
  
  ggsave("../../../../MS Revision/Figures/portPerfBayes_stdErr.pdf",
         g, width=3.25, height=3.25, units="in")






##############

#####
###Create a Bigger version for appendix with all the methods
g <- dat.sum %>% filter(isBayes) %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd)

g <- make_pretty(g, TRUE) + 
  theme(
    legend.direction = "horizontal", 
    legend.position=c(.5, .1), 
    legend.justification = "center") +
  guides(shape=guide_legend(nrow=2, byrow=TRUE)) + 
  ylab("(%) of Full-Info") 

g
ggsave("../../../../MS Revision/Figures/portPerfBayesBig.pdf",
       g, width=6.5, height=3.25, units="in")

#################
#Regularization versions
########
g <- filter(dat.sum, isReg, 
            !Method %in% c("FWRO_Eps_.05","HO_Reg", "LOO_Reg")) %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd)

g <- make_pretty(g, FALSE) + 
  theme(legend.direction = "horizontal", 
        legend.position=c(.5, .1), 
        legend.justification =  "center") +
  guides(shape=guide_legend(nrow=2, byrow=TRUE)) + 
  ylab("(%) of Full-Info") +
  scale_y_continuous(labels=scales::percent, limits =c(.7, .97),
                    breaks=seq(.75, .95, .1))
g

ggsave("../../../../MS Revision/Figures/portPerfRegb.pdf",
       g, width=3.25, height=3.25, units="in")

####
# The absolute perf version
#################
#Regularization versions
########
g <- filter(dat.sum2, isReg, 
            !Method %in% c("FWRO_Eps_.05", "LOO_Reg", "HO_Reg")) %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd)

g <- make_pretty(g, FALSE) + 
  theme(legend.direction = "horizontal", 
        legend.position=c(.5, .1), 
        legend.justification =  "center") +
  guides(shape=guide_legend(nrow=2, byrow=TRUE)) + 
  ylab("(%) of Full-Info") +
  scale_y_continuous()
g

g + scale_y_continuous() +
  ylab("Scaled Performance")

ggsave("../../../../MS Revision/Figures/portPerfReg_Abs.pdf", 
       g, width=3.25, height=3.25, units="in")



### Make a bigger version for appendix
g <- dat.sum %>% filter(isReg) %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd)

g <- make_pretty(g, FALSE) +  
  theme(legend.direction = "horizontal", 
        legend.position=c(.5, .1), 
        legend.justification =  "center")+
  guides(shape=guide_legend(nrow=2, byrow=TRUE)) + 
  ylab("(%) of Full-Info") +
  scale_y_continuous(labels=scales::percent, limits =c(.7, .97), 
                     breaks=seq(.75, .95, .1))
g

ggsave("../../../../MS Revision/Figures/portPerfRegBig.pdf", 
       g, width=6.5, height=3.25, units="in")


### sill version with std. err
g <- filter(dat.sum, isReg, 
            !Method %in% c("FWRO_Eps_.05","HO_Reg", "LOO_Reg")) %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=avg - std/sqrt(200), ymax=avg + std/sqrt(200)), position=pd)

g <- make_pretty(g, FALSE) + 
  theme(legend.direction = "horizontal", 
        legend.position=c(.5, .1), 
        legend.justification =  "center") +
  guides(shape=guide_legend(nrow=2, byrow=TRUE)) + 
  ylab("(%) of Full-Info") +
  scale_y_continuous(labels=scales::percent, limits =c(.7, .97),
                     breaks=seq(.75, .95, .1))
g

ggsave("../../../../MS Revision/Figures/portPerfReg_stderr.pdf",
       g, width=3.25, height=3.25, units="in")




######
# Similar plots for variance?
##########
gstdBayes <- filter(dat.sum, isBayes, Method != "EB_MM") %>%
  ggplot(aes(n, std, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) 

gstdBayes <- make_pretty(gstdBayes, TRUE) + 
  theme(legend.position=c(.23, .2), 
        legend.justification = "center", 
        legend.direction="vertical")+
  ylab("Std. Dev of (%) of Full-Info") + 
  scale_y_log10(label=scales::percent) + 
  guides(shape=guide_legend(ncol=2, byrow=FALSE))

ggsave("../../../../MS Revision/Figures/portVarBayes.pdf", 
       gstdBayes, width=3.25, height=3.25, units="in")

gstdReg <- filter(dat.sum, isReg, Method!= "RO_Eps_.05") %>%
  ggplot(aes(n, std, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd)

gstdReg <- make_pretty(gstdReg, FALSE) + 
  theme(legend.position=c(.23, .2), 
        legend.justification = "center", 
        legend.direction="vertical")+
        scale_y_log10(labels=scales::percent) + 
        ylab("Std. Dev of (%) of Full-Info")+ 
  guides(shape=guide_legend(ncol=2, byrow=FALSE))

ggsave("../../../../MS Revision/Figures/portVarReg.pdf", 
       gstdReg, width=3.25, height=3.25, units="in")


### Investigate the computational time dimension
dat.sum %>% 
  arrange(desc(n), desc(avgTime), Method) %>% 
  select(n, Method, avgTime, upTime) %>%
  View(.)
  
#Present average time relative to SAA
dat.sum<- dat.sum %>% ungroup()
t <- dat.sum %>% filter(Method == "SAA") %>%
  select(n, avgTime)
dat.sum <- inner_join(dat.sum, t, by=c("n")) %>%
  rename(SAATime = avgTime.y, 
         avgTime = avgTime.x)  %>%
  mutate(TimeRatio = avgTime/SAATime)

dat.sum %>% filter(Method %in% "SteinReg") %>% arrange(n) %>% select(avgTime, TimeRatio)

#Try to generate a plot with n along the top, and methods along the side
dat.sum %>% filter(Method %in% c("BoxStein", "EB_MLE", "SURE_MSE", "SteinReg", "LOO", "RO_Eps_.01", "SAA")) %>%
  select(n, Label, TimeRatio) %>%
  spread(Label, TimeRatio) %>%
  write_csv("../Results/port_avg_time_to_solve.csv")

dat.sum %>% filter(Method %in% c("SAA")) %>%
  select(n, Label, avgTime) %>%
  write_csv("../Results/port_avg_time_saa.csv")



#############
##Plots in Tau0
############
#Convergence in n for the interesting ones
g <- dat.sum %>% filter(isBayes, 
                        !Label %in% c("SAA", "EB LOO", "EB SURE")) %>%
  ggplot(aes(n, avgTau0, color=Label, group=Label, shape=Label, linetype=Label)) + 
  geom_point(position=pd) + geom_line(position=pd) +
  #geom_errorbar(aes(ymin=downTau0, ymax=upTau0), position=pd)+ 
  theme_bw() + scale_x_log10()

g <- make_pretty(g, TRUE) + 
  ylab("Tau") + scale_y_continuous()+ 
  theme(legend.position=c(.8, .7))

ggsave("../../../../MS Revision/Figures/convInTau.pdf", 
       g, width=3.25, height=3.25, units="in")

g <- dat %>% filter(log2(n) == 17, 
               isBayes, 
                !Method %in%c("SAA", "EB_MM") ) %>%
  ggplot(aes(x=Label, y=tau0, color=Method, group=Label)) + 
  geom_boxplot() + 
  theme_minimal(base_size=10) + 
  theme(legend.position="none", 
      text=element_text(family="Times New Roman")) + 
  xlab("") + ylab("Tau")

ggsave("../../../../MS Revision/Figures/boxPlotTau.pdf", 
       g, width=3.25, height=3.25, units="in")

dat.tau <- dat %>% filter(log2(n) == 17, 
                          Method %in% c("K5_EB", "HO_EB", "LOO_EB", "OR", "BoxStein")) %>%
            select(Run, Method, tau0) %>%
              spread(Method, tau0)
dat.tau %>% mutate_at(vars(-OR, -Run), 
                      funs(. > dat.tau$OR)) %>%
    summarise_at(vars(-OR, -Run), mean)

### REgularization versions
#Convergence in n for the interesting ones
g<- dat.sum %>% filter( !isBayes, 
          Method %in% c("K5_Reg", "HO_Reg", "LOO_Reg", "OracleReg", "SteinReg") ) %>% 
  ggplot(aes(n, avgTau0, color=Label, group=Label)) + 
  geom_point(position=pd) + geom_line(position=pd) +
#  geom_errorbar(aes(ymin=downTau0, ymax=upTau0), position=pd)+ 
  theme_bw() + scale_x_log10()
g <- make_pretty(g, FALSE) + 
  guides(color=guide_legend(nrow=2, byrow=TRUE))  + 
  ylab("Gamma") +   scale_y_continuous()+
  theme(legend.position=c(.5, .9), 
        legend.direction = "horizontal")
g

# ggsave("../../../../MS Revision//Figures/oldConvOfGamma.pdf", 
#        g, width=3.25, height=3.25, units="in")
ggsave("../../../../MS Revision//Figures/ConvInGamma.pdf", 
       g, width=3.25, height=3.25, units="in")

g <- dat %>% filter(log2(n) == 17, 
               Method %in% c("K5_Reg", "HO_Reg", "LOO_Reg", "OracleReg", "SteinReg")) %>%
  ggplot(aes(x=Label, y=tau0, color=Method, group=Label)) + 
  geom_boxplot() + 
  theme_minimal(base_size=10) + 
  theme(legend.position="none", 
        text=element_text(family="Times New Roman")) + 
  xlab("") + ylab("Gamma")

ggsave("../../../../MS Revision//Figures/boxPlotGamma.pdf", 
       g, width=3.25, height=3.25, units="in")

