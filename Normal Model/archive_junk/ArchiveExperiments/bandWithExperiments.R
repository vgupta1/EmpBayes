library(dplyr)
library(ggplot2)
library(reshape2)

setwd("~/Dropbox/Empirical Bayes/Experiments/EmpBayes/Normal Model/")

dat = read.csv("bandwidth2.csv")

filter(dat, Method == "Regularization") %>%
  group_by(bandwidth) %>% 
  summarise(avg=mean(thetaVal), 
            upval = mean(thetaVal) + sd(thetaVal), downval = mean(thetaVal) - sd(thetaVal)) %>% 
  ggplot(aes(x=bandwidth, y=avg)) + 
  geom_line() + geom_point() + geom_errorbar(aes(ymin=downval, ymax=upval))


filter(dat, Method == "Gauss") %>%
  group_by(bandwidth) %>% 
  summarise(avg=mean(thetaVal), 
            upval = mean(thetaVal) + sd(thetaVal), downval = mean(thetaVal) - sd(thetaVal)) %>% 
  ggplot(aes(x=bandwidth, y=avg)) + 
  geom_line() + geom_point() + geom_errorbar(aes(ymin=downval, ymax=upval)) + 
  xlim(-.5, 0)


filter(dat, Method == "ImpulseStein") %>%
  group_by(bandwidth) %>% 
  ggplot(aes(x=bandwidth, y=tau0)) + 
  geom_smooth()

