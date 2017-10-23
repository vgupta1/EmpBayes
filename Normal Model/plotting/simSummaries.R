library(tidyverse)
library(ggplot2)
library(gridExtra)
library(showtext)
font_add("Times New Roman", "Times New Roman.ttf")
showtext_auto()
font = "Times New Roman"

dat = read_csv("../Results/param_portExp_mtn1.csv")

##Trim it down to ake it possible to graph things
dat <- dat[1:2^(17), ]

g<- ggplot(dat, aes(thetas)) + 
  geom_histogram(aes(y=..density..), binwidth=.5) +
  geom_rug() + ylab("") +
  theme_minimal(base_family = "Times New Roman") + 
  xlab(expression(mu))

## Make an inlaid table for the grob
theta_table <- tribble(
  ~Stat, ~Value, 
  "Min", .001, 
  "Max", 31.7, 
  "Mean",.787, 
  "Median", .370
)


g<- g + annotation_custom(tableGrob(theta_table, rows=NULL, cols=NULL,
                                    theme=ttheme_minimal(base_size=8)), 
                          xmin=15, xmax=30, ymin=.4, ymax=.75)

ggsave("../../../../OR Submission_1/Figures/portMuPlot.png", g,
       width=3.25, height=3.25, units="in", dpi=600)

##############
##Repeat the process for c
g<- ggplot(dat, aes(c)) + 
  geom_histogram(aes(y=..density..), binwidth=.5) +
  geom_rug() + ylab("") +
  theme_minimal(base_size=11, base_family = "Times New Roman") + 
  xlab(expression(c))

## Make an inlaid table for the grob
cost_table <- tribble(
  ~Stat, ~Value, 
  "Min",  .050, 
  "Max",  65.9, 
  "Mean", 2.15, 
  "Median", 1.16
)


g<- g + annotation_custom(tableGrob(cost_table, rows=NULL, cols=NULL,
                                    theme=ttheme_minimal(base_size=8)), 
                          xmin=40, xmax=65, ymin=.3, ymax=.5)

ggsave("../../../../OR Submission_1/Figures/portCPlot.png", g,
       width=3.25, height=3.25, units="in", dpi=600)


####################
### Ratios
dat <- mutate(dat, ratio = thetas/c)

g<- ggplot(dat, aes(ratio)) + 
  geom_histogram(aes(y=..density..), binwidth=.05) +
  geom_rug() + ylab("") +
  theme_minimal(base_size=11, base_family = "Times New Roman") + 
  xlab(expression(mu/c)) + 
  xlim(0, 4)

## Make an inlaid table for the grob
ratio_table <- tribble(
  ~Stat, ~Value, 
  "Min",  .022, 
  "Max",  6.92, 
  "Mean", .335, 
  "Median", .335
)

g<- g + annotation_custom(tableGrob(ratio_table, rows=NULL, cols=NULL,
                                    theme=ttheme_minimal(base_size=8)), 
                          xmin=2, xmax=3.5, ymin=3, ymax=5)

ggsave("../../../../OR Submission_1/Figures/portRatioPlot.png", g,
       width=3.25, height=3.25, units="in", dpi=600)


#Nu versus ratio
g <- ggplot(dat, aes(ratio, vs)) + 
  geom_point(size=.5) + 
  theme_minimal(base_size=11, base_family = "Times New Roman") + 
  xlim(0, 4) + 
  ylab(expression(nu)) + xlab(expression(mu/c))

ggsave("../../../../OR Submission_1/Figures/portNuPlot.png", g,
       width=3.25, height=3.25, units="in", dpi=600)


