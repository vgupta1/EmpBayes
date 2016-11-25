##When does CLT kick in?
library(ggplot2)
library(reshape2)
library(dplyr)

setwd("~/Dropbox/Empirical Bayes/Experiments/EmpBayes/Normal Model/")

dat1 = read.csv("CLT_4_parallel_results.csv")
dat1$N = 4
dat2 = read.csv("CLT_8_parallel_results.csv")
dat2$N = 8
dat3 = read.csv("CLT_16_parallel_results.csv")
dat3$N = 16
dat4 = read.csv("CLT_32_parallel_results.csv")
dat4$N = 32
dat5 = read.csv("CLT_64_parallel_results.csv")
dat5$N = 64

dat = do.call("rbind", list(dat1, dat2, dat3, dat4, dat5))

#look at the rescaled thing
dat.scale <- dcast(dat, Run + n + N~ Method, value.var="thetaVal")
dat.scale <- melt(dat.scale, id.vars = c("Run", "n", "N", "OracleZ")) %>%
  mutate( thetaVal = value / OracleZ) 
names(dat.scale)[5] <- "Method"
head(dat.scale)

#just create a boxplot for different values of exact stein
filter(dat.scale, Method == "ExactStein", n==256) %>%
  ggplot(aes(x=factor(N), y=thetaVal)) + 
  geom_boxplot()


#better idea... do a line for each value of $N$
dat.scale %>% filter(Method=="ExactStein") %>%
  dplyr::group_by(n, Method, N) %>%
  dplyr::summarize(avg = mean(thetaVal), 
            upval = quantile(thetaVal, .95), 
            downval = quantile(thetaVal, .05)) %>%
  ggplot(aes(x=n, y=avg, color=factor(N))) +
  geom_point() +geom_line() + scale_x_log10()


dat.scale %>% filter(Method=="ExactStein") %>%
  dplyr::group_by(n, Method, N) %>%
  dplyr::summarize(avg = mean(thetaVal), 
                   upval = quantile(thetaVal, .95), 
                   downval = quantile(thetaVal, .05)) %>%
  ggplot(aes(x=n, y=avg, color=factor(N))) +
  geom_point() +geom_line() + scale_x_log10()


