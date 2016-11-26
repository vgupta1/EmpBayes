###  analyzing the empirical bayes Knapsack stuff
library(ggplot2)
library(reshape2)
library(dplyr)


setwd("~/Dropbox/Empirical Bayes/Experiments/EmpBayes/Normal Model/")

##Ridge regression is unstable for some reason
##for now throw away the bad ones.
dat <- filter(dat, Method!="Ridge" | thetaVal > 1e-6)

#########
#Rescale the results so that we are looking at percentages
dat.scale <- dcast(dat, Run + n ~ Method, value.var="thetaVal")
dat.scale <- melt(dat.scale, id.vars = c("Run", "n", "OracleZ")) %>%
    mutate( thetaVal = value / OracleZ) 
names(dat.scale)[4] <- "Method"
head(dat.scale)

###
#progression in n
###
pos = position_dodge(.2)
dat.scale %>%
#  filter(!Method %in% c("FullInfo", "OracleX", "NaiveX", "EmpBayesX", "Regularization_.01", "Regularization_.1", "MLE", "MM")) %>%
  filter(Method %in% c("ExactStein", "Regularization_1","Ridge", "MLE", "NaiveZ")) %>%
#  filter(Method %in% c("Ridge", "RidgeHalf", "MM")) %>%
#  filter(Method %in% c("Box", "Gauss", "ExactStein", "Sinc", "Regularization_1", "Ridge", "MLE", "NaiveZ")) %>%
  dplyr::group_by(n, Method) %>%
  dplyr::summarise(avg=mean(thetaVal), 
            low=quantile(thetaVal,.05), 
            high=quantile(thetaVal, .95)) %>%
  ggplot(aes(x=n, y=avg, color=Method, group=Method), 
         data=., position=pos) + 
#  geom_point(aes(shape=Method), position=pos) +
  geom_point(position=pos) +
  geom_errorbar(aes(ymin=low, ymax=high), position=pos) + 
   geom_line(position=pos) +
  scale_x_log10() + 
  theme_minimal(base_size=14) +
  theme(legend.title=element_blank())





## generate the pdf plots
filter(dat, Method == "Gauss") %>%
  ggplot(aes(x=thetaVal, group=factor(n), fill=factor(n))) + 
  geom_density(linetype="blank", alpha = .5) + 
  theme_minimal()

ggplot()



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


##Convergence of taus?
dat %>% 
  filter(Method %in% c("Bayes", "ExactStein", "ImpulseStein_3", "ImpulseStein_45", "ImpulseStein_6", "ImpulseStein_9")) %>%
  ggplot(aes(x=n, y=tau0, color=Method, group=Method), data=.) + 
  geom_smooth() +
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
             maxerr = max(err), 
             avgerr = mean(err))

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





