### Special version for the presentation
### renames steinBnd and trims n to make things cleaner.
library(tidyverse)
library(ggplot2)
library(forcats)
library(extrafont)
font = "Times New Roman"

dat = read_csv("../Results/3Part_plot__3part_0.0_0.1_20.001_1.0_1.0_20.001_0.3_4.0_20.001_8675309000.csv_full_200.csv")
dat = read_csv("../Results/3Part_plot_presentation___3part_0.0_0.1_20.001_1.0_1.0_20.001_0.3_4.0_20.001_8675309000.csv_full_80.csv")
dat$Method = as.factor(dat$Method)

#Reorder the levels in place to make plots consistent.
#This may need to be supplemented with specific color specs
dat <- mutate(dat, Method = fct_relevel(Method, 
                                        "FullInfo", "OR", "BoxStein", "DiracStein", 
                                        "OR_MSE", "EB_MLE", "EB_MM", "SURE_MSE", 
                                        "OracleReg", "OracleReg_1", "OracleReg_5", 
                                        "SteinReg", "SteinReg_1", "SteinReg_5",
                                        "LOO", "LOO_1", "LOO_5", 
                                        "RO_Eps_.01", "RO_Eps_.05", "RO_Eps_.1", 
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
                                      `Stein Reg .01` = "SteinReg",
                                      `Stein Reg 1` = "SteinReg_1", 
                                      `Stein Reg` = "SteinReg_5",
                                      `Oracle` = "OracleReg",
                                      `Oracle 1` = "OracleReg_1",
                                      `Oracle 5` = "OracleReg_5",
                                      `RO_1%` = "RO_Eps_.01", 
                                      `RO_5%` = "RO_Eps_.05", 
                                      `RO_10%` = "RO_Eps_.1", 
                                      `LOO .01` = "LOO", 
                                      `LOO`="LOO_1")
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



pd = position_dodge(.1)
g <- dat.sum %>% filter(Method %in% c("OR", 
                                      "EB_MLE", "EB_MM", 
                                      "SAA", "BoxStein", 
                                      "SURE_MSE")
) %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd) + 
  theme_minimal(base_size=12) +
  theme(legend.position = c(.7, .3), 
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

ggsave("../../../../Presentations/3PartEB_Results.png", 
       g, width=9.04, height=5.35, units="in")


############
dat = read_csv("../Results/portExp__8675309000.csv_full_200.csv")
dat = read_csv("../Results/portExp_presentation__8675309000.csv_full_80.csv")
dat$Method = as.factor(dat$Method)
##Reorder and relabel factors
dat <- mutate(dat, Method = fct_relevel(Method, 
                                        "FullInfo", "OR", "BoxStein", "DiracStein", 
                                        "OR_MSE", "EB_MLE", "EB_MM", "SURE_MSE", 
                                        "OracleReg", "OracleReg_1", "OracleReg_5", 
                                        "SteinReg", "SteinReg_1", "SteinReg_5",
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
                                      `Stein Reg` = "SteinReg",
                                      `Stein Reg 1` = "SteinReg_1", 
                                      `Stein Reg 5` = "SteinReg_5",
                                      `Oracle` = "OracleReg",
                                      `Oracle 1` = "OracleReg_1",
                                      `Oracle 5` = "OracleReg_5",
                                      `RO_1%` = "RO_Eps_.01", 
                                      `RO_5%` = "RO_Eps_.05", 
                                      `RO_10%` = "RO_Eps_.1", 
                                      `LOO .01` = "LOO", 
                                      `LOO 1`="LOO_1", 
                                      `LOO 5` = "LOO_5")
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

pd = position_dodge(.1)
g1 <- dat.sum %>% filter(Method %in% c("SAA","LOO_5", 
                                       "SteinReg_5", 
                                       "OracleReg", "RO_Eps_.01")
) %>%
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
  scale_y_continuous(labels=scales::percent)
g1

g1 <- dat.sum %>% filter(Method %in% c("OR", 
                                       "EB_MLE", "EB_MM", 
                                       "SAA", "BoxStein", 
                                       "SURE_MSE")
) %>%
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
  scale_y_continuous(labels=scales::percent)
g1




ggsave("../../../../Presentations/EBPortReg_Results.png", 
       g1, width=9.04, height=5.35, units="in")


###############
dat = read_csv("../Results/POAPCLT_plot_presentation__8675309.csv_full_80.csv")

dat <- mutate(dat, Method = fct_relevel(Method, 
                                        "FullInfo", "OR", "BoxStein", "DiracStein", 
                                        "OR_MSE", "EB_MLE", "EB_MM", "SURE_MSE", 
                                        "OracleReg", "OracleReg_1", "OracleReg_5", 
                                        "SteinReg", "SteinReg_1", "SteinReg_5",
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
                                      `Stein Reg .01` = "SteinReg",
                                      `Stein Reg 1` = "SteinReg_1", 
                                      `Stein Reg 5` = "SteinReg_5",
                                      `Oracle` = "OracleReg",
                                      `Oracle 1` = "OracleReg_1",
                                      `Oracle 5` = "OracleReg_5",
                                      `RO_1%` = "RO_Eps_.01", 
                                      `RO_5%` = "RO_Eps_.05", 
                                      `RO_10%` = "RO_Eps_.1")
)


#Re-express everything as a fraction of Oracle value
t <- dat %>% filter(Method == "FullInfo") %>%
  select(Run, N, thetaVal)
dat <- inner_join(dat, t, by=c("Run", "N")) %>%
  rename(thetaVal = thetaVal.x, 
         FullInfo = thetaVal.y)  %>%
  mutate(Ratio = thetaVal/FullInfo)

dat.sum <- dat %>% group_by(N, Method, Label) %>%
  summarise(avg = mean(Ratio), 
            avgTau0 = mean(tau0), 
            std  = sd(Ratio), 
            stdTau0 = sd(tau0), 
            up = quantile(Ratio, .9), 
            down = quantile(Ratio, .1), 
            avgTime = mean(time), 
            upTime = quantile(time, .9), 
            downTime = quantile(time, .1))



g <- dat.sum %>% filter(Method %in% c("OracleReg", "LOO_1", 
                                      "SteinReg_5", 
                                      "RO_Eps_.01")) %>%
  ggplot(aes(N, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd) + 
  theme_minimal() + 
  theme(legend.title=element_blank(), 
        text=element_text(family=font), 
        legend.justification = "center", 
        legend.position= c(.5, .1), 
        legend.direction="horizontal") + 
  scale_y_continuous(labels=scales::percent) +
  ylab("(%) of Full-Info") + xlab("S")
g

ggsave("../../../../Presentations/CLTPort_Results.png", 
       g, width=9, height=5.3, units="in")




