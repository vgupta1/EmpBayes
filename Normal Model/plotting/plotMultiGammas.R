###
# Plots for the Small-Data Linear Opt Paper
## #Convergence in Gamma
library(ggplot2)
library(stringr)
library(tidyverse)
library(extrafont) #for latex fonts
#font = "Times New Roman"
font = "CM Roman"

spath = "../../../../OR Submission_1/Figures/CompareConvByGamma.pdf"
dat = read_csv("../Results/OddEven_multi_gamma_8675309.csv")



dat <- dat %>% spread(Method, val) %>%
  mutate(Diff = abs(Oracle - Estimate))

dat.sum <- dat %>% group_by(n, Gamma) %>%
  summarise(avg = mean(Diff), 
            up = quantile(Diff, .9),
            down = quantile(Diff, .1))

#Code the labels by hand out of laziness
dat.sum <- mutate(dat.sum, GammaFactor = as.factor(Gamma))
# mylabs = c(expression(paste(Gamma, " = .01")), 
#            expression(paste(Gamma, " = .1")), 
#            expression(paste(Gamma, " = 1")), 
#            expression(paste(Gamma, " = 10"))) 
mylabs = levels(dat.sum$GammaFactor)


pd = position_dodge(.2)
g <- dat.sum %>% filter(n >= 400) %>%
  ggplot(aes(n, avg, group=Gamma, color=GammaFactor)) + 
  geom_point(aes(shape=GammaFactor), position=pd) + 
  geom_line(aes(linetype=GammaFactor), position=pd) +
  geom_errorbar(aes(ymin=down, ymax=up), position=pd) + 
  theme_minimal(base_size=12) +
  theme(legend.position = c(.7, .7), 
        legend.title=element_blank(), 
        text=element_text(family=font)) + 
  scale_x_log10(labels=scales::comma) + ylab("Abs Difference") +
  scale_color_discrete(labels = mylabs) +
  scale_shape_discrete(labels = mylabs) +
  scale_linetype_discrete(labels=mylabs)

g

ggsave(spath, 
       g, width=3.25, height=3.25, units="in")




