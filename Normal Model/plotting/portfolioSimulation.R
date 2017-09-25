####
# Generate simulated data according to the model of Pani et. al
# Used by julia script to run experimetns
###
library(tidyverse)
library(ggplot2)
library(copula)

##Parameters taken from personal correspondence with M. Sahin
gumbCop <- gumbelCopula(2.)
dist <- mvdc(gumbCop, margins=c("cauchy", "lnorm"), 
                    paramMargins=list(list(location=7.958527, scale = 12.208889), 
                                      list(meanlog=2.205216, sdlog=1.430539))
             )

#Generate excess data, and save down first 2^17 elements
set.seed(8675309)
n = 1000000
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
dat <- dat[1:2^17, ]
dat <- dat %>% mutate(thetas = thetas/200, c = c/200)

dat %>% ggplot(aes(thetas/c)) + geom_density()
summary(dat$thetas/dat$c)

### 
# Generate vs according to a variety of models
# all have the features
# low value/weight ratio -> very low precision
# medium value/weight ratio -> high precision
# high value/weight ratio -> medium precision
dat <- dat %>% mutate(ratio = thetas/c)
dat %>% ggplot(aes(ratio)) + geom_density()
dat <- dat %>% mutate(ratio_rank = rank(ratio)/nrow(dat))

#Goals:
 #Low ratio items should have enough noise to dominate top 10% of noisy ratio distribution
 #Med ratio items should have high precision (remember, these are to fool the other methods)
 #High ratio items shoudl enough noise to to seem worse than medium ratio items

#simplest.  Try a triangular shape...  
#Fix the endpoints to make the top/bottom 10% do what you want.  
#Save down a couple of placement of peaks
dat <- dat %>% mutate(type = ifelse(ratio_rank < .1, "Low", 
                                    ifelse(ratio_rank < .9, "Med", "High")))

dat %>% group_by(type) %>% 
  summarise(avgTheta = mean(thetas), 
            avgc = mean(c), 
            avgratio = mean(ratio)
            )
#Need low sigma to be about .2 = .6 * .46 - .073
#Double it (.4) to give room for high items.  So typical low item ratio  = 1.0

#high sigma can't make them beat low sigma often. (i.e. 2 stdev)
#(1.3 + 2 * sigma)/ 2.95 < 1... yields Sigma .825
#Try a .4 (for symmetry) for now

#do a linear interpolation?
v_fun <- approxfun(c(0, .2, 1), 
            c(1/.4^2, 1./01^2, 1/.4^2))

dat <- dat %>% mutate(vs = v_fun(ratio_rank))
ggplot(dat, aes(x=ratio, y=vs)) + geom_point()

## compute a realization of muhat
dat <- dat %>% mutate(muhat = thetas + rnorm(nrow(dat)/ sqrt(vs)))
dat <- dat %>% mutate(ratio_noisy = muhat/c) 

#sanity check, compute the values for SAA, random and and full-info
mean(dat$ratio)
quantile(dat$ratio, .95)

#What is the true distribution among SAA
dat %>% filter(ratio_noisy > quantile(dat$ratio_noisy, .95)) %>%
  group_by(type) %>% summarise(n()/6554)
  
# dat %>% select(thetas, vs, c) %>%
#   write_csv("../Results/param_portExp_Linear_.5.csv")

dat %>% select(thetas, vs, c) %>%
  write_csv("../Results/param_portExp_Linear_.2.csv")


