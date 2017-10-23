### Three part CLT Analysis and Plots
library(tidyverse)
library(ggplot2)
library(forcats)

library(extrafont)
font = "Times New Roman"
font = "CM Roman"

#dat = read_csv("../Results/3PartCLT_plot__3partCLT_131072_8675309.csv_full_8.csv")
dat = read_csv("../Results/POAPCLT_plot__123456.csv_full_8.csv")
#spath = "../../../../OR Submission_1/Figures/3PartCLT.pdf"

dat$Method = as.factor(dat$Method)

#Reorder the levels in place to make plots consistent.
#This may need to be supplemented with specific color specs
dat <- mutate(dat, Method = fct_relevel(Method, 
                                        "FullInfo", "OR", "BoxStein", "DiracStein", 
                                        "OR_MSE", "EB_MLE", "EB_MM", "SURE_MSE", "SAA"))

#Create a new column of "pretty" labels
dat <- mutate(dat, Label = fct_recode(Method, 
                                      `EB OPT` = "BoxStein", 
                                      `EB OPT Dirac` = "DiracStein", 
                                      `EB MLE` = "EB_MLE", 
                                      `EB MM` = "EB_MM", 
                                      `MSE OR`= "OR_MSE", 
                                      `EB OR` = "OR", 
                                      `EB SURE` = "SURE_MSE")
)

#Re-express everything as a fraction of full-info
t <- dat %>% filter(Method == "FullInfo") %>%
  select(Run, N, thetaVal)

dat <- inner_join(dat, t, by=c("Run", "N")) %>%
  rename(thetaVal = thetaVal.x, 
         FullInfo = thetaVal.y)  %>%
  mutate(Ratio = thetaVal/FullInfo)


#Summarize the data across runs
dat.sum <- dat %>% group_by(N, Method, Label) %>%
  summarise(avg = mean(Ratio), 
            avgTau0 = mean(tau0), 
            std  = sd(Ratio), 
            stdTau0 = sd(tau0), 
            up = quantile(Ratio, .9), 
            down = quantile(Ratio, .1))

pd = position_dodge(.5)
#Just the Bayes ones
g <- dat.sum %>% filter(!Method %in% c("FullInfo", "DiracStein", "OR_MSE", 
                                       "LOO", "OracleReg", "RO_Eps_.1", "RO_Eps_.05", "RO_Eps_.01", "SteinReg", "SteinRegBnded")) %>%
  ggplot(aes(N, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd) + 
  theme_minimal() + 
  theme(legend.title=element_blank(), 
        text=element_text(family=font), 
        legend.justification = "center", 
        legend.position=c(.5, .1), 
        legend.direction="horizontal") + 
  scale_y_continuous(labels=scales::percent, limits=c(.775, .9)) +
  ylab("(%) of Full-Info") + xlab("S")
g

ggsave("../../../../OR Submission_1/Figures/POAPCLT_Bayes.pdf", 
       g, width=3.25, height=3.25, units="in")

#Just the REg ones
g <- dat.sum %>% filter(!Method %in% c("FullInfo", "DiracStein", "OR_MSE", 
                                       "OR", "BoxStein", "OR_MSE", "EB_MLE", "EB_MM", "SURE_MSE", "RO_Eps_.1", "RO_Eps_.01")) %>%
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
  scale_y_continuous(labels=scales::percent, limits=c(.75, .95)) +
  ylab("(%) of Full-Info") + xlab("S")
g
ggsave("../../../../OR Submission_1/Figures/POAPCLT_Reg.pdf", 
       g, width=3.25, height=3.25, units="in")


