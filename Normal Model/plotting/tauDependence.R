## The TauDependence Plot
library(tidyverse)
library(ggplot2)
library(showtext)

font_add("Times New Roman", "Times New Roman.ttf")
showtext_auto()

# g <- qplot( x= seq(1:10), y = rnorm(10)) + 
#   geom_point() +
#   theme(text=element_text(family="Times New Roman")) + 
#   xlab(expression(tau))
# 
# g
# ggsave("temp.pdf", g)

###  Single plot with dependence in tau
dat = read_csv("../Results/tauDependence.csv")
g<- dat %>% 
  ggplot(aes(tau, thetaVal)) + 
  geom_point() + geom_line() + 
  theme_minimal(base_size=11) +
  theme(text=element_text(family="Times New Roman", size=11)) + 
  xlab(expression(tau)) + ylab("Objective") +
  ylim(0, .05)

ggsave("../../../../OR Submission_1/Figures/3PartTauDep.pdf", 
       g, width=3.25, height=3.25, units="in")
