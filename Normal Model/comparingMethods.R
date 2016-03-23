###  analyzing the empirical bayes Knapsack stuff
library(dplyr)
library(ggplot2)
library(reshape2)

setwd("~/Dropbox/Empirical Bayes/Experiments/Normal Model/")
dat = read.csv("test1_results_1_8675309.csv")

##reorder so that methods in a sensible order
dat$Method = factor(dat$Method, levels=c("NaiveX", "EmpBayesAvg", "EmpBayesX", "NaiveZ", "OracleX", "OracleZ", "FullInfo"))

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


