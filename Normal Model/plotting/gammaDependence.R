#####  Single Plot with dependence in Gamma
library(tidyverse)
library(ggplot2)
library(forcats)
font = "Times New Roman"
dat = read_csv("../Results/GammaDependence.csv")
g<- dat %>% 
  ggplot(aes(Gamma, thetaVal)) + 
  geom_point() + geom_line() + 
  theme_minimal(base_size=10) +
  theme(text=element_text(family=font, size=10)) + 
  xlab(expression(Gamma)) + ylab("(%) of Full-Info")

ggsave("../../../../OR Submission_1/Figures/3PartGammaDep.pdf", 
       g, width=3, height=3, units="in")


dat %>% filter(Method == "OracleReg") %>%
  ggplot(aes(tau0)) + 
  geom_density()

dat %>% filter(Method == "OracleReg") %>%
  group_by(n) %>%
  summarise(max = max(tau0), 
            min = min(tau0), 
            avg = mean(tau0), 
            bnd_up = mean(tau0 >= 19.999), 
            bnd_low = mean(tau0 <= .1), 
  )

