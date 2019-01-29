## Plot for Thm: LOO Validation is Not Best in Class
library(tidyverse)
library(ggplot2)
library(showtext)

font_add("Times New Roman", "Times New Roman.ttf")
showtext_auto()

#Updated for the new loo theorem
#dat = read_csv("../Results/Thm41.csv", col_names=FALSE)
dat = read_csv("../Results/UpdatedKFoldThmPlot.csv")

head(dat)
#make pretty labels and reorder the factor
dat<- dat %>% gather(Method, Perf, -Gamma) %>%
        mutate(Method = fct_relevel(Method, c("K2", "K5", "K10", "OR")))



# #Create an auxiliary tibble to place the labels
# dat_labels = tribble(
#             ~Gamma,   ~Perf,   ~Text, 
#             5.75,        .008,  "Oracle", 
#             5.75,        -.003, "LOO"
# )

g <- dat %>%
  ggplot(aes(Gamma, Perf)) +
  geom_line(aes(color=Method, linetype=Method)) +
  theme_minimal(base_size=8)
# g <-g + 
#   geom_text(data=dat_labels, aes(label=Text), 
#             family="Times New Roman", size=3) 
g<- g +
  theme(text=element_text(family="Times New Roman"), 
        legend.position=c(.8, .3), 
        legend.title=element_blank()) + 
  scale_y_continuous(labels=scales::percent) + 
  ylab("Performance") + 
  xlab(expression(Gamma)) + 
  ylim(-.006, .0012) + xlim(0, 25)
g

# ggsave("../../../../OR Submission_1/Figures/LooThmPic.pdf", 
#        g, width=2.5, height=2.5, units="in")
ggsave("../../../../MS Revision/Figures/KFoldThmPic.pdf", 
       g, width=2.5, height=2.5, units="in")








