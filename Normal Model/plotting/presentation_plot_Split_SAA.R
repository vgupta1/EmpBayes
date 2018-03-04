### A presentation Plot
#Show how the SAA example with 3 types of items decomposes into mixture pieces.  

library(tidyverse)
library(ggplot2)
library(forcats)
library(showtext)
font_add("Times New Roman", "Times New Roman.ttf")
showtext_auto()
font = "Times New Roman"


library(gridExtra)
dat = read_csv("../Results/3PartDensities_131072_20.00001.csv")
dat <- dat %>% mutate(Type=as.factor(mu), 
                      Type = fct_recode(Type, 
                                        Low="0", Med="0.3", High="1")
)

##First create plot for baseline
g <- dat %>%
  ggplot(aes(muhat, fill = Type)) +
  geom_density(alpha = .5, linetype = "blank") +
  geom_vline(xintercept = quantile(dat$muhat, .99),
             linetype = "dotted") +
  xlab("") + ylab("") +
  theme_minimal(base_size = 10, base_family = font) +
  theme(
    legend.title = element_blank(),
    legend.position = c(.2, .8),
    panel.grid = element_blank()
  )
g

g2 <- dat %>%
  ggplot(aes(muhat)) +
  geom_density(alpha=.5, fill="grey", linetype="blank") + 
  geom_vline(xintercept = quantile(dat$muhat, .99), 
             linetype="dotted") +
  xlab("") + ylab("") + 
  theme_minimal(base_size=11, base_family=font)

ggsave("../../../../Presentations/3PartDensity_Pres.pdf", 
       g, height = 3.7, width=5, units="in")

ggsave("../../../../Presentations/3PartDensity_grey_Pres.pdf", 
       g2, height = 3.7, width=5, units="in")

