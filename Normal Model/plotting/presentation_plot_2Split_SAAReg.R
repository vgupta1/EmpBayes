### A presentation Plot
#Show how the SAA example with 2 types for regularization works

library(tidyverse)
library(ggplot2)
library(forcats)
library(showtext)
# font_add("Times New Roman", "Times New Roman.ttf")
# showtext_auto()
# font = "Times New Roman"

library(gridExtra)

dat = read_csv("../Results/2PartDensity_Reg.csv")
dat <- dat %>% mutate(Type=as.factor(thetas), 
                      Type = fct_recode(Type, 
                                        Low="-1.75", High="0.15")
)

##First create plot for baseline
g <- dat %>%
  ggplot(aes(muhat, fill = Type)) +
  geom_density(alpha = .5, linetype = "blank") +
  geom_vline(xintercept = quantile(dat$muhat, .9),
             linetype = "dotted") +
  xlab("") + ylab("") +
  theme_minimal(base_size = 10) +
  theme(
    legend.title = element_blank(),
    legend.position = c(.2, .8),
    panel.grid = element_blank()
  )
g

ggsave("../../../../Presentations/2PartDensity_Pres.pdf", 
       g, height = 3.7, width=5, units="in")

