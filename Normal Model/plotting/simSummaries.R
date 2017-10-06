library(tidyverse)
library(ggplot2)
library(gridExtra)
library(showtext)
font_add("Times New Roman", "Times New Roman.ttf")
showtext_auto()
font = "Times New Roman"

dat = read_csv("../Results/param_portExp_Linear_.5.csv")

g<- ggplot(dat, aes(thetas)) + 
  geom_histogram(aes(y=..density..), binwidth=.5) +
  geom_rug() + ylab("") +
  theme_minimal(base_size=11, base_family = "Times New Roman") + 
  xlab(expression(mu))

## Make an inlaid table for the grob
theta_table <- tribble(
  ~Stat, ~Value, 
  "Min", .001, 
  "Max", 37.6, 
  "Mean",.792, 
  "Median", .368
)


g<- g + annotation_custom(tableGrob(theta_table, rows=NULL, cols=NULL,
                                    theme=ttheme_minimal(base_size=8)), 
                          xmin=20, xmax=35, ymin=.4, ymax=.75)

ggsave("../../../../OR Submission_1/Figures/portMuPlot.pdf", g,
       width=3.25, height=3.25, units="in")

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
  "Max",  77.5, 
  "Mean", 2.16, 
  "Median", 1.154
)


g<- g + annotation_custom(tableGrob(cost_table, rows=NULL, cols=NULL,
                                    theme=ttheme_minimal(base_size=8)), 
                          xmin=40, xmax=65, ymin=.3, ymax=.5)

ggsave("../../../../OR Submission_1/Figures/portCPlot.pdf", g,
       width=3.25, height=3.25, units="in")


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
  "Min",  .025, 
  "Max",  8.21, 
  "Mean", .335, 
  "Median", .335
)

g<- g + annotation_custom(tableGrob(ratio_table, rows=NULL, cols=NULL,
                                    theme=ttheme_minimal(base_size=8)), 
                          xmin=2, xmax=3.5, ymin=3, ymax=5)

ggsave("../../../../OR Submission_1/Figures/portRatioPlot.pdf", g,
       width=3.25, height=3.25, units="in")


#Nu versus ratio
g <- ggplot(dat, aes(ratio, vs)) + 
  geom_point(size=.5) + 
  theme_minimal(base_size=11, base_family = "Times New Roman") + 
  xlim(0, 4) + 
  ylab(expression(nu)) + xlab(expression(mu/c))

ggsave("../../../../OR Submission_1/Figures/portNuPlot.pdf", g,
       width=3.25, height=3.25, units="in")


