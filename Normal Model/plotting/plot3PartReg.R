### Three part Analysis and Plots for the Reg
library(tidyverse)
library(ggplot2)
library(forcats)
library(extrafont)
font = "Times New Roman"
font = "CM Roman"

dat = read_csv("../Results/3Part_plot__3part_0.0_0.1_20.001_1.0_1.0_20.001_0.3_4.0_20.001_8675309000.csv_full_200.csv")

#Temp smalll run to try and fix LOO
dat = read_csv("../Results/3Part_plot__3part_0.0_0.1_20.001_1.0_1.0_20.001_0.3_4.0_20.001_8675309000.csv_full_20.csv")

spath = "../../../../OR Submission_1/Figures/3PartReg20.pdf"

3
ggsave(spath, 
       g, width=3, height=3, units="in")


##Now a big version for the appendix
g <- dat.sum %>% filter(Method %in% c("SteinReg", "SteinRegBnded", "OracleReg", 
                                      "LOO", "SAA", "RO_Eps_.05", "RO_Eps_.01", "RO_Eps_.1")
) %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd) + 
  theme_minimal(base_size=11) +
  theme(legend.position = c(.8, .3), 
        legend.title=element_blank(), 
        text=element_text(family=font), 
        legend.direction = "horizontal", 
        #legend.text=element_text(size=8),
        legend.justification = "center") + 
  guides(shape=guide_legend(nrow=3, byrow=TRUE)) +
  scale_x_log10(labels=scales::comma) + 
  scale_y_continuous(labels=scales::percent, limits=c(0, 1)) +
  ylab("(%) of Full-Info")
g

ggsave("../../../../OR Submission_1/Figures/3PartReg20Big.pdf", 
       g, width=6.5, height=3.25, units="in")




#############
#What is a typical value of Gamma?
g1 <- dat %>% filter(Method %in% c("SteinReg", "SteinRegBnded", "OracleReg", "LOO"), 
                     n == 8192) %>%
  ggplot(aes(tau0)) + 
  geom_density(aes(group=Label, fill=Label), 
               alpha=.5, linetype="blank") + 
  xlab("Fitted tau") + ylab("") + 
  theme_minimal(base_size=11) + 
  theme(text=element_text(family=font),
        legend.title=element_blank(), 
        legend.position=c(.8, .8))

