### Special version for the presentation
library(tidyverse)
library(ggplot2)
library(forcats)
library(extrafont)
font = "Times New Roman"

dat = read_csv("../Results/2PartDensity_plot___8675309000.csv_full_100.csv")
dat$Method = as.factor(dat$Method)

#Reorder the levels in place to make plots consistent.
#This is customized for the presentation.... Needs to be made more consistent.
dat <- mutate(dat, Method = fct_relevel(Method, 
                                        "FullInfo", "OR", "BoxStein",  
                                        "OracleReg", "OracleReg_5", 
                                        "SteinReg", "SteinReg_5",
                                        "LOO", "LOO_5", 
                                        "FWRO_Eps_.1", 
                                        "SAA"))

#Create a new column of "pretty" labels
dat <- mutate(dat, Label = fct_recode(Method, 
                                      `EB OPT` = "BoxStein", 
                                      `EB OR` = "OR", 
                                      `Stein Reg` = "SteinReg", 
                                      `Oracle` = "OracleReg",
                                      `Oracle 5` = "OracleReg_5",
                                      `RO_1%` = "FWRO_Eps_.1", 
                                      `LOO` = "LOO")
)

#Re-express everything as a fraction of FullInfo value
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
g <- dat.sum %>% filter(Method %in% c("OracleReg", 
                                      "SteinReg",
                                      "LOO", 
                                      "SAA", "FWRO_Eps_.1")
                  ) %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd) + 
  theme_minimal(base_size=12) +
  theme(legend.position = c(.55, .9), 
        legend.title=element_blank(), 
        text=element_text(family=font), 
        legend.direction = "horizontal", 
        legend.text=element_text(size=16),
        legend.justification = "center") + 
  scale_x_log10(labels=scales::comma) + 
  scale_y_continuous(labels=scales::percent, limits=c(-1, 1)) +
  ylab("(%) of Full-Info")
g

ggsave("../../../../Presentations/BetterLoo.png",
       g, width=9.04, height=5.35, units="in")





