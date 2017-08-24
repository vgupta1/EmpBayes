###  analyzing the empirical bayes Knapsack stuff
library(ggplot2)
library(reshape2)
library(dplyr)
library(stringr)

setwd("~/Dropbox/Empirical Bayes/Experiments/EmpBayes/Normal Model/")

dat = read.csv("gaussian_2.0_0.0_2.0__parallel_results.csv")
dat = read.csv("threePart_0.01_0.01_1.5_0.1_8675309")

dat = read.csv("_CLTExp_4_2.0_4_parallel_results.csv")
dat = read.csv("_CLTExp_16_2.0_16_parallel_results.csv")


##Check the stability of the CV Estimators 
dat %>% filter(Method %in% c("CV_Shrink", "CV_Shrink2")) %>% 
  group_by(n, Method) %>% summarise(mean(tau0 < 1e-6))


#########
#Rescale the results so that we are looking at percentages
dat.scale <- dcast(dat, Run + n ~ Method, value.var="thetaVal")
dat.scale <- melt(dat.scale, id.vars = c("Run", "n", "OracleZ"), 
                              variable.name="Method", 
                              value.name = "thetaVal") %>%
    mutate( relPerf = thetaVal / OracleZ) 
head(dat.scale)

###
#progression in n for
###
pos = position_dodge(.2)
methods = c("SAA", "MLE", "CV_Shrink", "Rescaled", "MM", "ExactStein", "OracleReg", "SURE", "GaussO", "Box")
methods = c("BoxSAK", "GaussOSAK", "SincSAK", "ExactStein")


dat.scale %>%
  filter(Method %in% methods) %>%
  dplyr::group_by(n, Method) %>%
  dplyr::summarise(avg=mean(relPerf), 
            low=quantile(relPerf,.05), 
            high=quantile(relPerf, .95)) %>%
  ggplot(aes(x=n, y=avg, color=Method, group=Method), 
         data=., position=pos) + 
  geom_point(position=pos) +
  geom_errorbar(aes(ymin=low, ymax=high), position=pos) + 
    geom_line(position=pos) +
  scale_x_log10() + 
  theme_minimal(base_size=14) +
  theme(legend.title=element_blank())

###############
####
# Notes on the oddEven Case with tau_odd = .01 : oddEven_0.01_0.01_8675309.csv
####
# MLE does essentially perfectly...  How can we fool it better?  
#    Idea:  Need to construct example where "right" amount of shrinkage is distinct from variability in thetas
# Both CV and CV_shrink2 do very well.  (CV_shrink does better).
#   We expected them to undershrink.  Check the values. 

# Understanding the choice of kernel
# For this case it seems the kernel methods all fail for large n!


## Notes on the CLT Experiment
## 1.  Need a better examle that distinguishes between SAA and others
## 2.  Why does oracle regularization tank so hard?



####
#Studying Tau Convergence
####
methods = c("Box", "Gauss", "ExactStein", "OracleZ", "BoxS", "GaussS", "SincS")
dat.tau <- dcast(dat, Run + n ~ Method, value.var="tau0")
dat.tau <- melt(dat.tau, id.vars = c("Run", "n"), 
                  variable.name="Method", 
                  value.name = "tau0")  
head(dat.tau)

dat.tau %>%
  filter(Method %in% methods) %>%
  dplyr::group_by(n, Method) %>%
  dplyr::summarise(avg=mean(tau0), 
                   low=quantile(tau0,.05), 
                   high=quantile(tau0, .95)) %>%
  ggplot(aes(x=n, y=avg, color=Method, group=Method), 
         data=., position=pos) + 
  geom_point(position=pos) +
#  geom_errorbar(aes(ymin=low, ymax=high), position=pos) + 
  geom_line(position=pos) +
  scale_x_log10() + 
  theme_minimal(base_size=14) +
  theme(legend.title=element_blank())



#Study how many observations you need for a uniform before a CLT kicks in







#####
##
filter(dat, Method %in% c("OracleZ", "ExactStein", "ImpulseStein", "Rescaled")) %>%
  group_by(n, Method) %>%
  summarize(avg = mean(thetaVal))

filter(dat, n == 131072) %>%
  dcast(Run ~ Method, value.var="thetaVal") %>%
  ggplot(aes(x=ImpulseStein, y=ExactStein), data=.) + geom_point() + 
  geom_abline()


