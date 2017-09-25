library(tidyverse)
library(ggplot2)
library(forcats)
library(extrafont)
font = "Times New Roman"
font = "CM Roman"



###SAA Plot stuff
#Swithced to using this one to include the regularization on the plot
dat = read_csv("../Results/SAA_Plot__OddEven_1.0_2.1_8675309000.csv_full_204.csv")
# dat = read_csv("../Results/SAA_Plot__OddEven_1.0_2.0_8675309000.csv_full_200.csv")
#spath = "../../../../OR Submission_1/Figures/SAA1.pdf"

#dat = read_csv("../Results/SAA_Plot__OddEven_1.0_5.0_8675309000.csv_full_200.csv")
# spath = "../../../../OR Submission_1/Figures/SAA2.pdf"

#TEMP not meant to be saved
#Being used to get a sense of hte oracle regularization parameter
#dat = read_csv("../Results/SAA_Plot__OddEven_1.0_2.1_8675309000.csv_full_204.csv")

dat$Method = as.factor(dat$Method)
dat<-  mutate(dat, Method = fct_recode(Method, `EB Opt`="BoxStein", 
                                                 `Stein Reg` = "SteinReg"), 
                    Method = fct_relevel(Method, "SAA", after = Inf))

dat.sum <- dat %>%
  group_by(n, Method) %>%
  summarise(avg = mean(thetaVal) / .1, 
            up = quantile(thetaVal, .9) / .1, 
            down=quantile(thetaVal, .1) / .1)


#########
## Temp only
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

pd = position_dodge(.2)
dat.sum %>% filter(Method %in% c("SAA", "EB Opt", "OR_MSE", "OracleReg", "SteinReg")) %>%
  ggplot(aes(n, avg, group=Method, color=Method)) + 
  geom_point(aes(shape=Method), position=pd) + 
  geom_line(aes(linetype=Method), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd) + 
  theme_minimal(base_size=11) +
  theme(legend.position = c(.8, .2), 
        legend.title=element_blank()) + 
  scale_x_log10(labels=scales::comma) + 
  scale_y_continuous(labels=scales::percent, limits=c(.35, 1)) +
  ylab("(%) of Full-Info")

names(dat)
dat %>% filter(Method=="OracleReg", n > 20000) %>%
  ggplot(aes(tau0, fill=as.factor(n), group=as.factor(n))) + 
  geom_density(alpha = .5, linetype="blank")

filter(dat, Method=="OracleReg") %>% group_by(n) %>%
  summarise(max = max(tau0), min = min(tau0), 
            low_bound = mean(tau0 < .1 + 1e-10))

dat %>% filter(Method =="OracleReg") %>%
  ggplot(aes(n, tau0, group=n)) + 
  geom_boxplot() + scale_x_log10()

# End Temp
##############



pd = position_dodge(.2)
g <- dat.sum %>% filter(Method %in% c("SAA", "EB Opt", "Stein Reg")) %>%
  ggplot(aes(n, avg, group=Method, color=Method)) + 
  geom_point(aes(shape=Method), position=pd) + 
  geom_line(aes(linetype=Method), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd) + 
  theme_minimal(base_size=12) +
  theme(legend.position = c(.8, .2), 
        legend.title=element_blank(), 
        text=element_text(family=font)) + 
  scale_x_log10(labels=scales::comma) + 
  scale_y_continuous(labels=scales::percent, limits=c(.35, 1)) +
  ylab("(%) of Full-Info")

ggsave(spath, 
       g, width=3.25, height=3.25, units="in")
