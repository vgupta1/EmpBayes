## Plot for Thm: LOO Validation is Not Best in Class
library(tidyverse)
library(ggplot2)
library(showtext)

font_add("Times New Roman", "Times New Roman.ttf")
showtext_auto()

dat = read_csv("../Results/Thm41.csv", col_names=FALSE)
head(dat)
names(dat) <- c("Gamma", "Oracle", "LOO")

#Create an auxiliary tibble to place the labels
dat_labels = tribble(
            ~Gamma,   ~Perf,   ~Text, 
            5.75,        .008,  "Oracle", 
            5.75,        -.003, "LOO"
)


g <- dat %>% gather(Method, Perf, -Gamma) %>%
  ggplot(aes(Gamma, Perf)) + 
  geom_line(aes(color=Method)) + 
  theme_minimal(base_size=8) + 
  geom_text(data=dat_labels, aes(label=Text), 
            family="Times New Roman", size=3) +
  theme(text=element_text(family="Times New Roman"), 
        legend.position="none", 
        legend.title=element_blank()) + 
  scale_y_continuous(labels=scales::percent) + 
  ylab(" (%) of Full-Info") + 
  xlab(expression(Gamma))

ggsave("../../../../OR Submission_1/Figures/LooThmPic.pdf", 
       g, width=2.5, height=2.5, units="in")







