####
# Generate simulated data according to the model of Pani et. al
# Used by julia script to run experimetns
###
library(tidyverse)
library(ggplot2)
library(copula)
library(purrr)

##Parameters taken from personal correspondence with M. Sahin
gumbCop <- gumbelCopula(2.)
dist <- mvdc(gumbCop, margins=c("cauchy", "lnorm"), 
                    paramMargins=list(list(location=7.958527, scale = 12.208889), 
                                      list(meanlog=2.205216, sdlog=1.430539))
             )

#Generate excess data, and save down first 2^17 elements
set.seed(8675309)
n = 5000000
dat <- rMvdc(n, dist)
dat <- as_tibble(dat)
names(dat)<- c("Beta0", "Beta1")

#Filtering parameters obtained from personel correspondence
dat<- dat %>% mutate(exp_term = exp(-Beta0/Beta1)) %>%
  filter( -700 <= Beta0, Beta0 <= 100, 
          0.5 <= Beta1, Beta1 <= 800, 
          exp_term <= 20)

#Create the phis and the c's
#Since only one bid level, use middle "viable" bid.
dat <- dat %>% mutate(cmax = Beta1/.1, 
             c = cmax/.5,
              psi = exp(-Beta0/Beta1),
              phi = Beta0 + Beta1 * log(c + psi))

dat <- mutate(dat, thetas = phi) %>%
        select(c, thetas)

# Rescale for convenience so numbers closer to unity
# Constants chosen "by hand" but does not affect final computation
#dat <- dat[1:2^21, ]
dat <- dat %>% mutate(thetas = thetas/200, c = c/200)

dat %>% ggplot(aes(thetas/c)) + geom_density()
summary(dat$thetas/dat$c)

### 
# Generate vs 
dat <- dat %>% mutate(ratio = thetas/c)
dat <- dat %>% mutate(ratio_rank = rank(ratio)/nrow(dat))

#Goals:
 #Low ratio items should have enough noise to dominate top 10% of noisy ratio distribution
 #Med ratio items should have high precision (remember, these are to fool the other methods)
 #High ratio items shoudl enough noise to to seem worse than medium ratio items

#simplest.  Try a triangular shape...  
#Fix the endpoints to make the top/bottom 10% do what you want.  
#Save down a couple of placement of peaks
dat <- dat %>% mutate(type = ifelse(ratio_rank < .33, "Low", 
                                    ifelse(ratio_rank < .66, "Med", "High")))

dat %>% group_by(type) %>% 
  summarise(avgTheta = mean(thetas), 
            avgc = mean(c), 
            avgratio = mean(ratio)
            )

# #do a linear interpolation?
# v_fun <- approxfun(c(0, .5, 1), 
#             c(.1, 1, .1))
# t = seq(0, 1, by=.05)
# qplot(x=t, y=map_dbl(t, v_fun)) + geom_point()

#Try a piecewise constant approximation

dat2 <- dat %>% mutate(vs = ifelse(type=="Low", .1,
                            ifelse(type=="Med", 10, 
                                                8)),  #8
                       muhat = thetas + rnorm(nrow(dat2))/ sqrt(vs), 
                       ratio_noisy = muhat/c, 
                       rank_noisy = rank(ratio_noisy)/nrow(dat2), 
                       type_noisy = ifelse(rank_noisy < .33, "Low", 
                                           ifelse(rank_noisy < .66, "Med", "High")))


#Now compute tthe confusion
#table(dat2$type, dat2$type_noisy)

dat2 %>% filter(type_noisy == "High") %>% group_by(type) %>%
  summarise(num = n()) %>%
  mutate(prop = num/sum(num)) %>%
  select(-num)


##Current Best
#c(.01,10, 1))
#  35   38  27

### write down a samplet to try out
dat2 %>% select(thetas, vs, c) %>%
  write_csv("../Results/param_portExp_mtn1.csv")

ggplot(dat2, aes(ratio_rank, vs)) + geom_point()


#sanity check, compute the values for SAA, random and and full-info
mean(dat$ratio)
quantile(dat$ratio, .8)

#What is the true distribution among SAA
dat %>% filter(ratio_noisy > quantile(dat$ratio_noisy, .8)) %>%
  group_by(type) %>% summarise(n()/26214)

dat %>% filter(ratio_noisy > quantile(dat$ratio_noisy, .8)) %>%
  summarise(mean(ratio))


dat %>%
  filter(ratio_noisy > 0) %>%
  ggplot(aes(ratio, ratio_noisy, color=type)) + 
    geom_point(alpha = .5, size = 1) + 
  xlim(0, 4)

  
# dat %>% select(thetas, vs, c) %>%
#   write_csv("../Results/param_portExp_Linear_.5.csv")

dat %>% select(thetas, vs, c) %>%
  write_csv("../Results/param_big_portExp_Linear_.5.csv")

### Read it in, scale it down and save it back
dat <- read_csv("../Results/param_portExp_Linear_.5.csv")

dat %>% mutate(thetas = thetas / 10, c = c/10) %>%
  write_csv("../Results/param_portExp_scaled_.5.csv")




