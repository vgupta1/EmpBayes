###  analyzing the empirical bayes Knapsack stuff
library(dplyr)
library(ggplot2)
library(reshape2)

setwd("~/Dropbox/Empirical Bayes/Experiments/EmpBayes/Normal Model/")
dat = read.csv("results_10_8675309.csv")
dat = read.csv("test2_results_0.5_0.2_8675309.csv")
dat = read.csv("test2_results_0.5_3_8675309.csv")
dat = read.csv("test1_resultsb_1_8675309.csv")
dat = read.csv("exp_results_Exp_1_1_8675309.csv")

##redone with new harness




##top level summary
dat %>% group_by(n, Method) %>%
  summarise(avg=mean(thetaVal))

filter(dat, n>65000, Method != "FullInfo") %>%
  group_by(Method) %>%
  summarise(avg=mean(thetaVal), 
            low=quantile(thetaVal,.05), 
            high=quantile(thetaVal, .95)) %>%
  arrange(desc(avg)) %>%
  ggplot(aes(x=Method, y=avg, color=Method), data=.) + 
  geom_errorbar(aes(ymin=low, ymax=high)) + 
  theme(legend.position="none")

#progression in n
pos = position_dodge(.2)
dat %>% filter(n>1e3) %>%
  filter(! Method %in% c("FullInfo")) %>%
#  filter(Method %in% c(""NaiveZ", "MLE","MM", "OracleZ", "Rescaled")) %>%
  group_by(n, Method) %>%
  summarise(avg=mean(thetaVal), 
            low=quantile(thetaVal,.05), 
            high=quantile(thetaVal, .95)) %>%
  ggplot(aes(x=n, y=avg, color=Method, group=Method), 
         data=., position=pos) + 
  geom_point(position=pos) +
  geom_errorbar(aes(ymin=low, ymax=high), position=pos) + 
  scale_x_log10() + 
  theme(legend.title=element_blank())
  
dat %>% filter(Run == 6) %>%
  filter(Method %in% c("MLE", "MM", "OracleZ", "Rescaled")) %>%
  ggplot(aes(x=n, y=thetaVal, color=Method, group=Method), data=.) + 
  geom_point() + geom_line() + 
  scale_x_log10()

##Convergence of taus?
dat %>% 
  filter(Method %in% c("MLE", "MM", "OracleZ", "ZZ", "Rescaled")) %>%
ggplot(aes(x=n, y=tau0, color=Method, group=Method), data=.) + 
  geom_smooth(aes(fill=Method)) +
  scale_x_log10()  + 
  theme_bw()


########
### Understanding resaling diffs
#######
dat = read.csv("test3_results_tau_scale_8675309.csv")
dat = read.csv("test4_results_tau_scale_0.5_0.2_8675309.csv")
dat = read.csv("test4_results_tau_scale_0.5_3_8675309.csv")



dat$FScale = factor(dat$Scale)
dat %>% group_by(n, FScale) %>%
  summarise(avg= mean(thetaVal)) %>%
  ggplot(aes(x=n, y=avg, color=FScale), data=.) + 
  geom_point() + geom_line() +
  theme_bw() +
  scale_x_log10()

##Does the wrong tau eventually converge?
pos = position_dodge(.2)
dat %>% 
  group_by(n, FScale) %>%
  summarise(avg=mean(thetaVal), 
            low=quantile(thetaVal,.05), 
            high=quantile(thetaVal, .95)) %>%
  ggplot(aes(x=n, y=avg, color=FScale), 
         data=., position=pos) + 
  geom_point() + geom_line() +
  scale_x_log10() + 
  theme(legend.title=element_blank())


  


######
##Everyone should be dominated path by path by OracleZ
dat.sum <- dcast(dat, Run + n ~ Method, value.var="thetaVal")

dat.sum %>% filter(n>65000) %>%
  ggplot(aes(x=OracleZ, y=OracleX), data=.) + 
  geom_point() + geom_abline()


#Further tests
#What does our new method unscaled look like?
#What does scaling but using xs look like?
#Where does TauX and Rescaled TauX Converge to?
#Does Rao-Blackkwelization help in estimating taux?
#How big is the interval?  anyway to check?





