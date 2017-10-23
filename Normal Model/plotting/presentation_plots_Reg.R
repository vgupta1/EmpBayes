### Special version for the presentation
### renames steinBnd and trims n to make things slower.

library(tidyverse)
library(ggplot2)
library(forcats)
library(extrafont)
font = "Times New Roman"

dat = read_csv("../Results/3Part_plot__3part_0.0_0.1_20.001_1.0_1.0_20.001_0.3_4.0_20.001_8675309000.csv_full_200.csv")

dat$Method = as.factor(dat$Method)

#Reorder the levels in place to make plots consistent.
#This may need to be supplemented with specific color specs
dat <- mutate(dat, Method = fct_relevel(Method, 
                                        "FullInfo", "OR", "BoxStein", "DiracStein", 
                                        "OR_MSE", "EB_MLE", "EB_MM", "SURE_MSE", 
                                        "OracleReg", "SteinReg", "SteinRegBnded",
                                        "LOO", "RO_Eps_.01", "RO_Eps_.05", "RO_Eps_.1", 
                                        "SAA"))

#Create a new column of "pretty" labels
dat <- mutate(dat, Label = fct_recode(Method, 
                                      `EB OPT` = "BoxStein", 
                                      `EB OPT Dirac` = "DiracStein", 
                                      `EB MLE` = "EB_MLE", 
                                      `EB MM` = "EB_MM", 
                                      `MSE OR`= "OR_MSE", 
                                      `EB OR` = "OR", 
                                      `EB SURE` = "SURE_MSE", 
                                      `Stein Reg Old` = "SteinReg",
                                      `Stein Reg` = "SteinRegBnded",
                                      `Reg OR` = "OracleReg", 
                                      `RO_1%` = "RO_Eps_.01", 
                                      `RO_5%` = "RO_Eps_.05", 
                                      `RO_10%` = "RO_Eps_.1")
)


#Re-express everything as a fraction of Oracle value
t <- dat %>% filter(Method == "FullInfo") %>%
  select(Run, n, thetaVal)
dat <- inner_join(dat, t, by=c("Run", "n")) %>%
  rename(thetaVal = thetaVal.x, 
         FullInfo = thetaVal.y)  %>%
  mutate(Ratio = thetaVal/FullInfo)

#Summarize the data across runs
dat.sum <- dat %>% group_by(n, Method, Label) %>%
  summarise(avg = mean(Ratio), 
            avgTau0 = mean(tau0), 
            std  = sd(Ratio), 
            stdTau0 = sd(tau0), 
            up = quantile(Ratio, .9), 
            down = quantile(Ratio, .1), 
            upTau0 = quantile(tau0, .9), 
            downTau0 = quantile(tau0, .1))



pd = position_dodge(.3)
g <- dat.sum %>% filter(Method %in% c("SteinRegBnded", "OracleReg", 
                                      "LOO", "SAA", "RO_Eps_.01"), 
                        n < 1200 
) %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd) + 
  theme_minimal(base_size=16) +
  theme(legend.position = c(.55, .9), 
        legend.title=element_blank(), 
        text=element_text(family=font), 
        legend.direction = "horizontal", 
        legend.text=element_text(size=16),
        legend.justification = "center") + 
#  guides(shape=guide_legend(nrow=2, byrow=TRUE)) +
  scale_x_log10(labels=scales::comma) + 
  scale_y_continuous(labels=scales::percent, limits=c(0, 1.05)) +
  ylab("(%) of Full-Info")
g

ggsave("../../../../Presentations/3PartReg_Results.pdf", 
       g, width=9, height=5.3, units="in")


############
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
                                      `Stein Old` = "SteinReg",
                                      `Stein Reg` = "SteinRegBnded",
                                      `Reg OR` = "OracleReg", 
                                      `RO_1%` = "RO_Eps_.01", 
                                      `RO_5%` = "RO_Eps_.05", 
                                      `RO_10%` = "RO_Eps_.1")
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

pd = position_dodge(.3)
g1 <- dat.sum %>% filter(Method %in% c("SAA", "SteinRegBnded", "LOO", 
                                       "OracleReg", "RO_Eps_.01"), 
                         n < 1200) %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd) + 
  theme_minimal(base_size=16) +
  theme(legend.title=element_blank(),
        legend.text=element_text(size=16),
        text=element_text(family=font), 
        legend.direction = "horizontal", 
        legend.position=c(.5, .1), 
        legend.justification =  "center")+
  #guides(shape=guide_legend(nrow=2, byrow=TRUE)) + 
  scale_x_log10(labels=scales::comma) + 
  ylab("(%) of Full-Info") +
  scale_y_continuous(labels=scales::percent, limits =c(.7, .97), 
                     breaks=seq(.75, .95, .1))
g1

ggsave("../../../../Presentations/PortReg_Results.pdf", 
       g1, width=9, height=5.3, units="in")







