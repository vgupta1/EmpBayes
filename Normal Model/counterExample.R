## CounterExample ot the TauXY/2Policy

library(ggplot2)
library(plyr)
library(dplyr)

setwd("~/Dropbox/Empirical Bayes/Experiments/EmpBayes/Normal Model/")

dat1 = read.csv("UniformLikelihood_10.0_1.0_15.0_14.1_5164174290.csv")
dat2 = read.csv("UniformLikelihood_10.0_1.0_15.0_14.1_5167462266.csv")
dat3 = read.csv("UniformLikelihood_10.0_1.0_15.0_14.1_8675309.csv")

dat2 <- mutate(dat2, Run = Run + 10)
dat3 <- mutate(dat3, Run = Run + 20)

dat <- rbind(dat1, dat2, dat3)
pos = position_dodge(.2)
dat %>% filter(n>2000, Method %in% c("OracleZ", "Rescaled")) %>%
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

