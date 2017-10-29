## The TauDependence Plot
library(tidyverse)
library(ggplot2)
library(showtext)

font_add("Times New Roman", "Times New Roman.ttf")
showtext_auto()

###  Single plot with dependence in tau
dat = read_csv("../Results/tauDepependence_131072_8675309.csv")
dat2 = read_csv("../Results/tauDepependence_512_8675309.csv")
dat = rbind(dat, dat2)
dat$n = as.factor(dat$n)

g<- dat %>% 
  ggplot(aes(tauVal, RelFullInfo, linetype=n, color=n)) + 
  geom_line() + 
  theme_minimal(base_size=10) +
  theme(text=element_text(family="Times New Roman"), 
        legend.position=c(.8, .2)) + 
  xlab(expression(tau)) + ylab("(%) of Full-Info") + 
  scale_y_continuous(labels=scales::percent)
g

ggsave("../../../../OR Submission_1/Figures/3PartTauDep.pdf", 
       g, width=3, height=3, units="in")
