## Messing around with the port Simulation Experiment
library(tidyverse)
library(ggplot2)
library(forcats)
library(showtext)
font_add("Times New Roman", "Times New Roman.ttf")
showtext_auto()
font = "Times New Roman"

#Uses updated experimetnal values
dat = read_csv("../Results/portExp__8675309000.csv_full_200.csv")


dat$Method = as.factor(dat$Method)
##Reorder and relabel factors
dat <- mutate(dat, Method = fct_relevel(Method, 
                                        "FullInfo", "OR", "BoxStein", "DiracStein", 
                                        "OR_MSE", "EB_MLE", "EB_MM", "SURE_MSE", 
                                        "OracleReg", "SteinReg", "SteinRegBnded",
                                        "LOO", "RO_Eps_.01", "RO_Eps_.05", "RO_Eps_.1", 
                                        "SAA"))

#Create a new column of "pretty" labels
dat <- mutate(dat, Label = fct_recode(Method, 
                                      `EB OPT` = "BoxStein", 
                                      `EB Exact` = "DiracStein", 
                                      `EB MLE` = "EB_MLE", 
                                      `EB MM` = "EB_MM", 
                                      `MSE OR`= "OR_MSE", 
                                      `EB OR` = "OR", 
                                      `EB SURE` = "SURE_MSE", 
                                      `Stein Reg` = "SteinReg",
                                      `Stein Bnd` = "SteinRegBnded",
                                      `Reg OR` = "OracleReg", 
                                      `RO_.01` = "RO_Eps_.01", 
                                      `RO_.05` = "RO_Eps_.05", 
                                      `RO_.10` = "RO_Eps_.1")
)


#Re-express everything as a fraction of Oracle value
t <- dat %>% filter(Method == "FullInfo") %>%
  select(Run, n, thetaVal)
dat <- inner_join(dat, t, by=c("Run", "n")) %>%
  rename(thetaVal = thetaVal.x, 
         FullInfo = thetaVal.y)  %>%
  mutate(Ratio = thetaVal/FullInfo)

dat.sum <- dat %>% group_by(n, Method, Label) %>%
  summarise(avg = mean(Ratio), 
            avgTau0 = mean(tau0), 
            std  = sd(Ratio), 
            stdTau0 = sd(tau0), 
            up = quantile(Ratio, .9), 
            down = quantile(Ratio, .1), 
            avgTime = mean(time), 
            upTime = quantile(time, .9), 
            downTime = quantile(time, .1))

##What are typical values of Tau OR and TauHat?
dat.sum %>% filter(Method %in% c("OR", "BoxStein")) %>%
  arrange(n, Method) %>% select(n, Method, avgTau0)
  

## plot some shit
pd = position_dodge(.3)
g <- dat.sum %>% filter(Method %in% c("SAA","BoxStein", "SURE_MSE", "EB_MLE", "OR")
                  ) %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd) + 
  theme_minimal(base_size=11) +
  theme(legend.title=element_blank(), 
        text=element_text(family=font), 
        legend.direction = "horizontal", 
        legend.position=c(.5, .1), 
        legend.justification = "center")+
  guides(shape=guide_legend(nrow=2, byrow=TRUE)) + 
  scale_x_log10(labels=scales::comma) + 
  ylab("(%) of Full-Info") +
  scale_y_continuous(labels=scales::percent, limits = c(.70, .97),
                     breaks=seq(.75, .95, .1))
  
g
ggsave("../../../../OR Submission_1/Figures/portPerfBayes.pdf", 
       g, width=3.25, height=3.25, units="in")


###Create a Bigger version for appendix with all the methods
g <- dat.sum %>% filter(Method %in% c("SAA","BoxStein", "SURE_MSE", "EB_MLE", "OR", "EB_MM")
) %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd) + 
  theme_minimal(base_size=11) +
  theme(legend.title=element_blank(), 
        text=element_text(family=font), 
        legend.direction = "horizontal", 
        legend.position=c(.5, .1), 
        legend.justification = "center")+
  guides(shape=guide_legend(nrow=2, byrow=TRUE)) + 
  scale_x_log10(labels=scales::comma) + 
  ylab("(%) of Full-Info") +
  scale_y_continuous(labels=scales::percent, limits = c(.70, .97),
                     breaks=seq(.75, .95, .1))

g
ggsave("../../../../OR Submission_1/Figures/portPerfBayesBig.pdf", 
       g, width=6.5, height=3.25, units="in")


#Regularization version

pd = position_dodge(.3)
g1 <- dat.sum %>% filter(Method %in% c("SAA", "SteinReg", "SteinRegBnded", "LOO", 
                                      "OracleReg", "RO_Eps_.01")) %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd) + 
  theme_minimal(base_size=11) +
  theme(legend.title=element_blank(), 
        text=element_text(family=font), 
        legend.direction = "horizontal", 
        legend.position=c(.5, .1), 
        legend.justification =  "center")+
  guides(shape=guide_legend(nrow=2, byrow=TRUE)) + 
  scale_x_log10(labels=scales::comma) + 
  ylab("(%) of Full-Info") +
  scale_y_continuous(labels=scales::percent, limits =c(.7, .97), 
                    breaks=seq(.75, .95, .1))
g1

ggsave("../../../../OR Submission_1/Figures/portPerfReg.pdf", 
       g1, width=3.25, height=3.25, units="in")

### Make a bigger version for appendix
pd = position_dodge(.3)
g1 <- dat.sum %>% filter(Method %in% c("SAA", "SteinReg", "SteinRegBnded", "LOO", 
                                       "OracleReg", "RO_Eps_.01", "RO_Eps_.05", "RO_Eps_.1")) %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd) + 
  theme_minimal(base_size=11) +
  theme(legend.title=element_blank(), 
        text=element_text(family=font), 
        legend.direction = "horizontal", 
        legend.position=c(.5, .1), 
        legend.justification =  "center")+
  guides(shape=guide_legend(nrow=2, byrow=TRUE)) + 
  scale_x_log10(labels=scales::comma) + 
  ylab("(%) of Full-Info") +
  scale_y_continuous(labels=scales::percent, limits =c(.7, .97), 
                     breaks=seq(.75, .95, .1))
g1

ggsave("../../../../OR Submission_1/Figures/portPerfRegBig.pdf", 
       g1, width=6.5, height=3.25, units="in")

######
# Similar plots for variance?

gstdBayes <- dat.sum %>% filter(Method %in% c("SAA","BoxStein", "SURE_MSE", "EB_MLE", "EB_MM")
) %>%
  ggplot(aes(n, std, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  theme_minimal(base_size=11) +
  theme(legend.title=element_blank(), 
        text=element_text(family=font), 
        legend.position=c(.23, .3), 
        legend.justification = "center", 
        legend.direction="vertical")+
  scale_x_log10(labels=scales::comma) + scale_y_log10(labels=scales::percent) + 
  ylab("(%) of Full-Info")

ggsave("../../../../OR Submission_1/Figures/portVarBayes.pdf", 
       gstdBayes, width=3.25, height=3.25, units="in")


gstdReg <- dat.sum %>% 
  filter(Method %in% c("SAA", "SteinReg", "LOO", 
                       "RO_Eps_.01", "RO_Eps_.1")) %>%
  ggplot(aes(n, std, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  theme_minimal(base_size=11) +
  theme(legend.title=element_blank(), 
        text=element_text(family=font), 
        legend.position=c(.23, .3), 
        legend.justification = "center", 
        legend.direction="vertical")+
  scale_x_log10(labels=scales::comma) + scale_y_log10(labels=scales::percent) + 
  ylab("(%) of Full-Info")

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




