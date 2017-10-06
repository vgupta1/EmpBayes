#Plots to illustrate converegence ot normal distribution
#For our simulation method

library(tidyverse)
library(ggplot2)

library(showtext)
font_add("Times New Roman", "Times New Roman.ttf")
showtext_auto()
font = "Times New Roman"



## Try to do the normal generation phase
sim<- function(N){
  n = 100000
  x <- matrix(rexp(n * N) - 1, n, N) %>% rowSums()
  x <- x / sqrt(N)
  return(x)
}

library(purrr)
N_grid = c(1, 5, 10, 15, 20)
t <- N_grid %>% 
  map(sim) %>%
  data.frame() 
names(t) <- N_grid
dat <- t %>% gather(N, x, -N)

library(forcats)
dat <- mutate(dat, Label = as.factor(N))
levels(dat$Label )
dat$Label <- fct_relevel(dat$Label, "5", after=1)

g <- ggplot(dat, aes(x=x, group=N, 
                color=as.factor(Label), 
                linetype=as.factor(Label)
                )
      ) + 
  geom_density() + 
  theme_minimal(base_size=11, base_family="Times New Roman") + 
  theme(legend.position = c(.7, .7), 
        legend.title=element_blank()) + 
  ylab("") + xlab("") + xlim(-4, 6)

ggsave("../../../../OR Submission_1/Figures/NormalConv.pdf", 
       g, width=3.25, height=3.25, units="in")


