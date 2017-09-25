### Three part Analysis and Plots for the Reg
library(tidyverse)
library(ggplot2)
library(forcats)
library(extrafont)
#font = "Times New Roman"
font = "CM Roman"

dat = read_csv("../Results/3Part_plot__3part_0.0_0.1_20.001_1.0_1.0_20.001_0.3_4.0_20.001_8675309000.csv_full_200.csv")
spath = "../../../../OR Submission_1/Figures/3PartReg20.pdf"

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
                                      `EB Opt` = "BoxStein", 
                                      `EB Opt Dirac` = "DiracStein", 
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

#Summarize the data across runs
dat.sum <- dat %>% group_by(n, Method) %>%
  summarise(avg = mean(Ratio), 
            avgTau0 = mean(tau0), 
            std  = sd(Ratio), 
            stdTau0 = sd(tau0), 
            up = quantile(Ratio, .9), 
            down = quantile(Ratio, .1))

dat.sum <- mutate(dat.sum, 
                  Label = fct_recode(Method, 
                                     `EB Opt` = "BoxStein", 
                                     `EB Opt Dirac` = "DiracStein", 
                                     `EB MLE` = "EB_MLE", 
                                     `EB MM` = "EB_MM", 
                                     `MSE OR`= "OR_MSE", 
                                     `EB OR` = "OR", 
                                     `EB SURE` = "SURE_MSE", 
                                     `Stein Reg` = "SteinReg",
                                     `Stein Reg` = "SteinReg",
                                     `Stein Bnd` = "SteinRegBnded",
                                     `Reg OR` = "OracleReg", 
                                     `RO_.01` = "RO_Eps_.01", 
                                     `RO_.05` = "RO_Eps_.05", 
                                     `RO_.10` = "RO_Eps_.1")
)

pd = position_dodge(.3)
g <- dat.sum %>% filter(Method %in% c("SteinReg", "SteinRegBnded", "OracleReg", 
                                      "LOO", "RO_Eps_.01", "RO_Eps_.05")
                        ) %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd) + 
  theme_minimal(base_size=11) +
  theme(legend.position = c(.5, .1), 
        legend.title=element_blank(), 
        text=element_text(family=font), 
        legend.direction = "horizontal", 
        legend.justification = "center") + 
  scale_x_log10(labels=scales::comma) + 
  scale_y_continuous(labels=scales::percent, limits=c(-.15, 1)) +
  ylab("(%) of Full-Info")
g

ggsave(spath, 
       g, width=4, height=4, units="in")


g <- dat.sum %>% filter(Method %in% c("SteinReg", "SteinRegBnded", "OracleReg", 
                                      "LOO")
) %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd) + 
  theme_minimal(base_size=11) +
  theme(legend.position = c(.5, .1), 
        legend.title=element_blank(), 
        text=element_text(family=font), 
        legend.direction = "horizontal", 
        legend.justification = "center") + 
  scale_x_log10(labels=scales::comma) + 
  scale_y_continuous(labels=scales::percent, limits=c(-.15, 1)) +
  ylab("(%) of Full-Info")






#############
#What is a typical value of Gamma?
g1 <- dat %>% filter(Method %in% c("SteinReg", "SteinRegBnded", "OracleReg", "LOO"), 
                     n == 8192) %>%
  ggplot(aes(tau0)) + 
  geom_density(aes(group=Label, fill=Label), 
               alpha=.5, linetype="blank") + 
  xlab("Fitted tau") + ylab("") + 
  theme_minimal(base_size=11) + 
  theme(text=element_text(family=font),
        legend.title=element_blank(), 
        legend.position=c(.8, .8))

