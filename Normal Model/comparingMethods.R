###  analyzing the empirical bayes Knapsack stuff
library(dplyr)
library(ggplot2)
library(reshape2)

setwd("~/Dropbox/Empirical Bayes/Experiments/EmpBayes/Normal Model/")
dat = read.csv("results_1_8675309.csv")

##top level summary
dat %>% group_by(n, Method) %>%
  summarise(avg=mean(thetaVal))

filter(dat, n==65536) %>%
  group_by(Method) %>%
  summarise(avg=mean(thetaVal), 
            low=quantile(thetaVal,.05), 
            high=quantile(thetaVal, .95)) %>%
  arrange(desc(avg)) %>%
  ggplot(aes(x=Method, y=avg, color=Method), data=.) + 
  geom_errorbar(aes(ymin=low, ymax=high)) + 
  theme(legend.position="none")

pos = position_dodge(.2)
dat %>% 
  filter(Method %in% c("Bayes", "MLE","MM", "OracleZ", "Rescaled")) %>%
  group_by(n, Method) %>%
  summarise(avg=mean(thetaVal), 
            low=quantile(thetaVal,.05), 
            high=quantile(thetaVal, .95)) %>%
  ggplot(aes(x=n, y=avg, color=Method, group=Method), 
         data=., position=pos) + 
  geom_point(position=pos) +
  geom_errorbar(aes(ymin=low, ymax=high), position=pos) + 
  scale_x_log10() + 
  theme_bw() + theme(legend.title=element_blank())
  
  
##Sanity Checks
##Tau MLE, and TauMM be converging to 1.  
dat %>% 
  filter(Method %in% c("MLE", "MM", "EmpBayesX", "OracleX", "OracleZ", "Rescaled")) %>%
ggplot(aes(x=n, y=tau0, color=Method, group=Method), data=.) + 
  geom_smooth() +
  scale_x_log10() + ylim(0, 5) + 
  theme_bw()


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







##some spot checks
ggplot(aes(x=n, y=thetaVal, color=Method, group=Method), 
       data=dat) + 
  geom_point() + geom_line()

dat %>% filter(n==100000, Run==1) %>%
  select(Run, Method, thetaVal) %>%
  arrange(desc(thetaVal))

#compare across for a singel n
filter(dat, n==25600) %>%
ggplot(aes(x=Method, y=thetaVal, color=Method, group=Method), data=.) + 
  geom_boxplot() + theme(legend.position="none")

dat %>% group_by(n, Method) %>%
  summarise(avg = mean(thetaVal), med=median(thetaVal), 
            low = quantile(thetaVal, .05), up = quantile(thetaVal, .95)) %>%
  arrange(n, desc(avg))

#what's happenign with taus?
dat %>% filter(Method=="EmpBayesX") %>%
  ggplot(aes(x=n, y=tau0), data=.) + 
  geom_point() + 
  scale_x_log10() + 
  geom_smooth()+ ylim(0, 3)




##NaiveX's Yval compare to NaiveX's thetaval?  Would expect unbiased estimate.
t_dat = filter(dat, Method=="NaiveX") %>%
  select(n, Run, YVal, thetaVal) %>%
  mutate(diff = YVal - thetaVal)

t_dat %>% group_by(n) %>%
  summarise(avg = mean(diff), std = sd(diff)/10, 
            low = avg-std, high=avg+std) %>%
  select(n, avg, low, high)

  ggplot(aes(x=YVal, y=thetaVal), data=.) + 
  geom_point() + geom_abline()
  
  
  
  melt(id.vars = c("n", "Run"), measure.vars=c("thetaVal", "YVal")) %>%
  cast(Run + n ~ variable) %>%
  head()



##Does NaiveZ perform better than Naive X?

##Q's for full data-set.
# How do EmpBayesX and EmpBayesAvg compare?




###Sanity Checks
##Some checks: OracleX >= EmpBayesX path by path for thetaval.  Should be nearly identical
filter(dat, Method %in% c("EmpBayesX", "OracleX")) %>%
  select(Run, n, Method, thetaVal) %>%
  dcast(formula=Run + n ~ Method) %>%
  ggplot(aes(x=EmpBayesX, y=OracleX, color=factor(n)), data=.) + 
  geom_point() + 
  geom_abline()

##EmpBayesX >= NaiveX path by path for yval... should hold approx for thetaval
filter(dat, Method %in% c("EmpBayesX", "NaiveX")) %>%
  select(Run, n, Method, thetaVal) %>%
  dcast(formula=Run + n ~ Method) %>%
  ggplot(aes(x=EmpBayesX, y=NaiveX, color=factor(n)), data=.) + 
  geom_point() + 
  geom_abline()


