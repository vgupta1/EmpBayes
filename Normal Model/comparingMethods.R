###  analyzing the empirical bayes Knapsack stuff
library(dplyr)
library(ggplot2)
library(reshape2)

setwd("~/Dropbox/Empirical Bayes/Experiments/Normal Model/")
dat = read.csv("test1_results_1_8675309.csv")

#spot check that things look sensible

ggplot(aes(x=n, y=thetaVal, color=Method, group=Method), data=dat) + 
  geom_point() + scale_x_log10() + geom_smooth()

dat %>% group_by(Method) %>%
  summarise(avg = mean(thetaVal), med=median(thetaVal), 
            low = quantile(thetaVal, .05), up = quantile(thetaVal, .95)) %>%
  arrange(desc(avg))

##Some checks: OracleX >= EmpBayesX path by path for thetaval.  Should be nearly identical
filter(dat, Method %in% c("EmpBayesX", "OracleX")) %>%
  select(Run, n, Method, thetaVal) %>%
  dcast(formula=Run + n ~ Method) %>%
  ggplot(aes(x=EmpBayesX, y=OracleX, color=factor(n)), data=.) + 
  geom_point() + 
  geom_abline()

##EmpBayesX >= NaiveX path by path for yval... should hold approx for thetaval
filter(dat, Method %in% c("EmpBayesX", "NaiveX")) %>%
  select(Run, n, Method, YVal) %>%
  dcast(formula=Run + n ~ Method) %>%
  ggplot(aes(x=EmpBayesX, y=NaiveX, color=factor(n)), data=.) + 
  geom_point() + 
  geom_abline()

##Still something seems buggy up...
dat.sum <- dcast(dat,formula = Run + n ~ Method, value.var="YVal")
filter(dat.sum, EmpBayesX < NaiveX ) %>%
  head()



##NaiveX's Yval compare to NaiveX's thetaval?  Would expect unbiased estimate.


##Does NaiveZ perform better than Naive X?

##Q's for full data-set.
# How do EmpBayesX and EmpBayesAvg compare?

