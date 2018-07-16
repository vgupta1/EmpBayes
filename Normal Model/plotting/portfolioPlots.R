## POAP Experiment Plots
source("plottingUtils.R")

###
#Old data from original publication.  Needs to be regenerated.  
dat = read_csv("../Results/portExp__8675309000.csv_full_200.csv")
dat = read_csv("../Results/portExp__8675309000.csv_full_100.csv")


dat <- clean_data(dat)
dat.sum <- summarize_dat(dat)
dat.sum2 <- summarize_dat(dat, FALSE)

## plot for main text
g <- dat.sum %>% filter(isBayes, Method != "EB_MM") %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd)

g <- make_pretty(g) + 
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
ggsave("../../../../OR Submission_1/Figures/portPerfBayes.pdf", 
       g, width=3.25, height=3.25, units="in")


## The absolute performance version
g <- dat.sum2 %>% filter(isBayes, Method != "EB_MM") %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd)

g <- make_pretty(g) + 
  theme(legend.direction = "horizontal",
        legend.position=c(.5, .1),
        legend.justification = "center") +
  guides(shape=guide_legend(nrow=2, byrow=TRUE)) + 
  ylab("(%) of Full-Info") + 
  scale_y_continuous(labels=scales::percent, 
                     limits = c(.70, .97),
                     breaks=seq(.75, .95, .1)) + 
  scale_x_log10(labels=scales::comma) 

g <- g + scale_y_continuous() + 
  ylab("Scaled Performance")
ggsave("../../../../OR Submission_1/Figures/portPerfBayes_Abs.pdf", 
       g, width=3.25, height=3.25, units="in")

##############



#####
###Create a Bigger version for appendix with all the methods
g <- dat.sum %>% filter(isBayes) %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd)

g <- make_pretty(g) + 
  theme(
    legend.direction = "horizontal", 
    legend.position=c(.5, .1), 
    legend.justification = "center") +
  guides(shape=guide_legend(nrow=2, byrow=TRUE)) + 
  ylab("(%) of Full-Info") 

g
ggsave("../../../../OR Submission_1/Figures/portPerfBayesBig.pdf", 
       g, width=6.5, height=3.25, units="in")

#################
#Regularization versions
########
g <- filter(dat.sum, isReg, Method != "RO_Eps_.05") %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd)

g <- make_pretty(g) + 
  theme(legend.direction = "horizontal", 
        legend.position=c(.5, .1), 
        legend.justification =  "center") +
  guides(shape=guide_legend(nrow=2, byrow=TRUE)) + 
  ylab("(%) of Full-Info") +
  scale_y_continuous(labels=scales::percent, limits =c(.7, .97),
                    breaks=seq(.75, .95, .1))
g

ggsave("../../../../OR Submission_1/Figures/portPerfReg.pdf", 
       g, width=3.25, height=3.25, units="in")

####
# The absolute perf version
#################
#Regularization versions
########
g <- filter(dat.sum2, isReg, Method != "RO_Eps_.05") %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd)

g <- make_pretty(g) + 
  theme(legend.direction = "horizontal", 
        legend.position=c(.5, .1), 
        legend.justification =  "center") +
  guides(shape=guide_legend(nrow=2, byrow=TRUE)) + 
  ylab("(%) of Full-Info") +
  scale_y_continuous(labels=scales::percent, limits =c(.7, .97),
                     breaks=seq(.75, .95, .1))
g

g + scale_y_continuous() +
  ylab("Scaled Performance")

ggsave("../../../../OR Submission_1/Figures/portPerfReg_Abs.pdf", 
       g, width=3.25, height=3.25, units="in")



### Make a bigger version for appendix
g <- dat.sum %>% filter(isReg) %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd)

g <- make_pretty(g) +  
  theme(legend.direction = "horizontal", 
        legend.position=c(.5, .1), 
        legend.justification =  "center")+
  guides(shape=guide_legend(nrow=2, byrow=TRUE)) + 
  ylab("(%) of Full-Info") +
  scale_y_continuous(labels=scales::percent, limits =c(.7, .97), 
                     breaks=seq(.75, .95, .1))
g

ggsave("../../../../OR Submission_1/Figures/portPerfRegBig.pdf", 
       g, width=6.5, height=3.25, units="in")

######
# Similar plots for variance?
##########
gstdBayes <- filter(dat.sum, isBayes, Method != "EB_MM") %>%
  ggplot(aes(n, std, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) 

gstdBayes <- make_pretty(gstdBayes) + 
  theme(legend.position=c(.23, .3), 
        legend.justification = "center", 
        legend.direction="vertical")+
  ylab("Std. Dev of (%) of Full-Info") + 
  scale_y_log10(label=scales::percent)

ggsave("../../../../OR Submission_1/Figures/portVarBayes.pdf", 
       gstdBayes, width=3.25, height=3.25, units="in")


gstdReg <- filter(dat.sum, isReg, Method!= "RO_Eps_.05") %>%
  ggplot(aes(n, std, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd)

gstdReg <- make_pretty(gstdReg) + 
  theme(legend.position=c(.23, .3), 
        legend.justification = "center", 
        legend.direction="vertical")+
        scale_y_log10(labels=scales::percent) + 
        ylab("Std. Dev of (%) of Full-Info")

ggsave("../../../../OR Submission_1/Figures/portVarReg.pdf", 
       gstdReg, width=3.25, height=3.25, units="in")


### Investigate the computational time dimension
dat.sum %>% filter(Method %in% c("FullInfo", "BoxStein", "SURE_MSE", "SAA")
                   )%>%
  arrange(desc(n), Method) %>% select(n, Method, avgTime, upTime)
  

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




