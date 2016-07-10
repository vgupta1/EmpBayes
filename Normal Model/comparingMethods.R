###  analyzing the empirical bayes Knapsack stuff
library(dplyr)
library(ggplot2)
library(reshape2)

setwd("~/Dropbox/Empirical Bayes/Experiments/EmpBayes/Normal Model/")

##redone with new harness
dat = read.csv("Results/testResults_Gamma_1_0.5_8675309.csv")
dat = read.csv("Results/testResults_Gamma_1_2_8675309.csv")
dat = read.csv("Results/testResults_Uniform_-1_1_8675309.csv")
dat = read.csv("Results/testResults_Uniform_-3_3_8675309.csv")
dat = read.csv("Results/testResults_Beta_0.5_0.5_8675309.csv")

dat = read.csv("Results/test_gaussian_means_-3_2_8675309.csv")
dat = read.csv("Results/test_gaussian_means_0_2_8675309.csv")
dat = read.csv("temp_0_2_8675309.csv")


#limit yourself to interesting metohds
method_bad = c("MLE", "FullInfo", "NaiveZ", "NaiveX", "OracleX", "Rescaled", "Rescaled_2.0", "Rescaled_4.0")

##top level summary
dat %>% group_by(Method) %>% filter(n>4000, ! Method %in% method_bad) %>%
  summarise(avg=mean(thetaVal)) %>% arrange(desc(avg))

filter(dat, n>65000, !Method %in% method_bad) %>%
  group_by(Method) %>%
  summarise(avg=mean(thetaVal), 
            low=quantile(thetaVal,.05), 
            high=quantile(thetaVal, .95)) %>%
  arrange(desc(avg)) %>%
  ggplot(aes(x=Method, y=avg, color=Method), data=.) + 
  geom_errorbar(aes(ymin=low, ymax=high)) + 
  geom_point() +
  theme(legend.position="none")

#progression in n
pos = position_dodge(.2)
dat %>%filter(n >2000) %>%
  filter(! Method %in% c("FullInfo", "MLE", "ZZ", "NaiveX", "EmpBayesX", 
                         "OracleX", "NaiveZ", "MM") ) %>%
  group_by(n, Method) %>%
  summarise(avg=mean(thetaVal), 
            low=quantile(thetaVal,.05), 
            high=quantile(thetaVal, .95)) %>%
  ggplot(aes(x=n, y=avg, color=Method, group=Method), 
         data=., position=pos) + 
  geom_point(aes(shape=Method), position=pos) +
  geom_errorbar(aes(ymin=low, ymax=high), position=pos) + 
  scale_x_log10() + 
  theme(legend.title=element_blank())
  
##Convergence of taus?
dat %>% 
  filter(! Method %in% c("FullInfo", "MLE", "ZZ", "NaiveX", "EmpBayesX", 
                         "OracleX", "NaiveZ", "MM") ) %>%
  filter(! Method %in% c("Rescaled_0.125", "Rescaled_0.25")) %>%
  filter(n < 1e3) %>%
  ggplot(aes(x=n, y=tau0, color=Method, group=Method), data=.) + 
  geom_smooth(se=FALSE) +
  scale_x_log10()  + 
  theme_bw()

dat %>% 
  filter(! Method %in% c("FullInfo", "MLE", "ZZ", "NaiveX", "EmpBayesX", 
                         "OracleX", "NaiveZ", "MM") ) %>%
  group_by(Method, n) %>%
  summarize(mean(tau0)) %>% dcast(n ~ Method) %>% select(n, OracleZ, Rescaled)  


###################
# Sanity Checks
###################
# Looking specifically at the gaussian case
dat.small <- filter(dat, Method %in% c("OracleZ", "Bayes"))
dat.small %>%
  dcast(Run + n ~ Method, value.var = "thetaVal") %>%
  group_by(n) %>% mutate(Diff=OracleZ - Bayes, err = Diff/Bayes) %>%
  summarise( mindiff= min(Diff), 
             maxdiff=max(Diff), 
             indx = which(Diff == min(Diff)), 
             maxerr = max(err))

dat.small %>% filter(n == 512) %>%
  dcast(Run ~ Method, value.var = "thetaVal") %>% 
  ggplot(aes(x=Bayes, y=OracleZ), data=.) + geom_point()
  
##pick one place where things go wrong...
dat.small %>% filter(n == 512) %>%











pos = position_dodge(.2)
dat %>%filter(n >2000) %>%
  filter(! Method %in% method_bad) %>%
  filter( Method != "EmpBayesX") %>%
  group_by(n, Method) %>%
  summarise(avg=mean(thetaVal), 
            low=quantile(thetaVal,.05), 
            high=quantile(thetaVal, .95)) %>%
  ggplot(aes(x=n, y=avg, color=Method, group=Method), 
         data=., position=pos) + 
  geom_point(position=pos) + geom_line()
  scale_x_log10() + 
  theme(legend.title=element_blank())




dat %>% filter(Run == 6) %>%
  filter(Method %in% c("MLE", "MM", "OracleZ", "Rescaled")) %>%
  ggplot(aes(x=n, y=thetaVal, color=Method, group=Method), data=.) + 
  geom_point() + geom_line() + 
  scale_x_log10()



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





