library(tidyverse)
library(ggplot2)
library(stringr)

library(extrafont) #for latex fonts
#font = "Times New Roman"
font = "CM Roman"


dat = read_csv("../Results/Unif_plot.csv")
dat$Method = as.factor(dat$Method)

#mark the true one as solid... everyone else is dashed
dat <- dat %>% mutate(isApprox = !Method == "Target")

g1<- dat %>% filter(!str_detect(Method, "NonU.*")) %>%
  ggplot(aes(x_val, val, group=Method, color=Method)) +
  geom_line(aes(linetype=isApprox)) + 
  xlab("") + ylab("") + 
  theme(legend.position="none") + 
  theme_minimal(base_size=12) + 
  theme(text=element_text(family=font), 
        legend.position="none") + 
  ylim(-10, 30)

g2<- dat %>% filter(str_detect(Method, "NonU.*") | Method=="Target") %>%
  ggplot(aes(x_val, val, group=Method, color=Method)) +
  geom_line(aes(linetype=isApprox)) + 
  xlab("") + ylab("") + 
  theme(legend.position="none") + 
  theme_minimal(base_size=12) + 
  theme(text=element_text(family=font), 
        legend.position="none") + 
  ylim(-10, 30)

ggsave("../../Normal Model/../../../OR Submission_1/Figures/UniformApprox1.pdf", 
       g1, width=3.25, height=3.25, units="in")

ggsave("../../Normal Model/../../../OR Submission_1/Figures/UniformApprox2.pdf", 
       g2, width=3.25, height=3.25, units="in")
