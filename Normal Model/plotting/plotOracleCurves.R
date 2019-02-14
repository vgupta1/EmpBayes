#####
#  Small plotting functions for the oracle curves in POAP Experiment
#####
library(tidyverse)
library(ggplot2)
library(forcats)
library(showtext)
font_add("Times New Roman", "Times New Roman.ttf")
showtext_auto()

dat = read_csv("../Results/RegOracleCurve.csv")

#use the full-info value to rescale everyone
zstar = dat$FullInfo[1]
dat <- dat %>% mutate_at(vars(-Gamma, -FullInfo), 
                  funs(./zstar))
#hold onto the oracle values for safekeeping
oracle_vals <- dat$Oracle

#melt it, change the labels
dat <- dat %>% gather(Method, Perf, -Gamma) %>%
  mutate(Method = fct_recode(Method,
      `Full-Info` = "FullInfo",
      `Reg OR` = "Oracle",
      `Reg OPT` = "Stein",
      `Reg LOO` = "LOO",
      `Reg HO` = "HO",
      `Reg K5` = "K5"
    )
  ) 

#build an auxiliary dataset for the optimal points
dat.maxGammas <- dat %>%
  group_by(Method) %>%
  summarise(index = which.max(Perf), 
            Gamma = Gamma[index], 
            EstValue = Perf[index], 
            Value = oracle_vals[index]) %>%
  filter(Method != "Full-Info")

g <- dat %>% filter(Method == "Reg OR") %>%
  ggplot(aes(Gamma, Perf)) + 
  geom_line() + 
  geom_point(aes(Gamma, Value, color=Method, shape=Method), 
             data=dat.maxGammas) + 
  theme_minimal(base_size=10) + 
  theme(legend.title=element_blank(), 
        legend.position = c(.8, .4), 
        text=element_text(family="Times New Roman")) +
  scale_y_continuous(labels=scales::percent)+ 
  ylab("% of Full-Info")
  

#By hand remap the colors/shapes to match the originals
library(scales)
Reg_colors<- hue_pal()(9)
names(Reg_colors) <- c("SAA", "Reg OR", "Reg OPT", "RO 1%", "RO 5%", "Reg LOO", "Reg HO", "Reg K5", "XX")
Reg_shapes = c(16, 17, 15,  3,  7,  8, 5, 6, 0)
names(Reg_shapes) <- c("SAA", "Reg OR", "Reg OPT", "RO 1%", "RO 5%", "Reg LOO", "Reg HO", "Reg K5", "XX")

g <- g + scale_color_manual(values=Reg_colors) + scale_shape_manual(values=Reg_shapes)
  
ggsave("../../../../MS Revision/Figures/OracleCurveReg.pdf", 
       g, width=3.25, height=3.25, units="in")

### now create a similar plot but with the many curves
g2 <- dat %>% filter(Method != "Full-Info") %>%
  ggplot(aes(Gamma, Perf, color=Method, linetype=Method)) + 
    geom_line() + 
  theme_minimal(base_size=10) + 
  theme(legend.title=element_blank(), 
        legend.position = c(.8, .4), 
        text=element_text(family="Times New Roman")) +
  scale_y_continuous(labels=scales::percent)+ 
  ylab("% of Full-Info")
g


ggsave("../../../../MS Revision/Figures/CompareCurvesReg.pdf", 
       g2, width=3.25, height=3.25, units="in")

#################
####Now repeat for EB Method
####
dat = read_csv("../Results/EBOracleCurve.csv")

#use the full-info value to rescale everyone
zstar = dat$FullInfo[1]
dat <- dat %>% mutate_at(vars(-Tau, -FullInfo), 
                         funs(./zstar))

#hold onto the oracle values for safekeeping
oracle_vals <- dat$Oracle

#melt it, change the labels
dat <- dat %>% gather(Method, Perf, -Tau) %>%
  mutate(Method = fct_recode(Method,
                             `Full-Info` = "FullInfo",
                             `EB OR` = "Oracle", 
                             `EB OPT` = "Stein", 
                             `EB LOO` = "LOO", 
                             `EB K5` = "K5", 
                             `EB HO` = "HO" 
                            )
          ) 


#build an auxiliary dataset for the optimal points
dat.maxTaus <- dat %>%
  group_by(Method) %>%
  summarise(index = which.max(Perf), 
            Tau = Tau[index], 
            EstValue = Perf[index], 
            Value = oracle_vals[index]) %>%
  filter(Method != "Full-Info")

g<- dat %>% filter(Method == "EB OR") %>%
  ggplot(aes(Tau, Perf)) + 
  geom_line() + 
  geom_point(aes(Tau, Value, color=Method, shape=Method), 
             data=dat.maxTaus) + 
  theme_minimal(base_size=10) + 
  theme(legend.title=element_blank(), 
        legend.position = c(.8, .4), 
        text=element_text(family="Times New Roman")) +
  scale_y_continuous(labels=scales::percent)+ 
  ylab("% of Full-Info")

EB_shapes = c(16, 17, 15,  3,  7,  8, 5, 6, 0)
names(EB_shapes) <- c("SAA", "EB OR", "EB OPT", "EB MLE", "EB MM", "EB SURE", "EB LOO", "EB K5", "EB HO")
EB_colors <- hue_pal()(9)
names(EB_colors) <- c("SAA", "EB OR", "EB OPT", "EB MLE", "EB MM", "EB SURE", "EB LOO", "EB K5", "EB HO")
EB_linetypes = linetype_pal()(9)
names(EB_linetypes) <- c("SAA", "EB OR", "EB OPT", "EB MLE", "EB MM", "EB SURE", "EB LOO", "EB K5", "EB HO")


g <- g + scale_color_manual(values=EB_colors) + scale_shape_manual(values=EB_shapes)

ggsave("../../../../MS Revision/Figures/oracleCurveEB.pdf", 
       g, width=3.25, height=3.25, units="in")


### now create a similar plot but with the many curves
g2 <- dat %>% filter(Method != "Full-Info") %>%
  ggplot(aes(Tau, Perf, color=Method, linetype=Method)) + 
  geom_line() + 
  theme_minimal(base_size=10) + 
  theme(legend.title=element_blank(), 
        legend.position = c(.8, .4), 
        text=element_text(family="Times New Roman")) +
  scale_y_continuous(labels=scales::percent)+ 
  ylab("% of Full-Info") + 
  + scale_color_manual(values=EB_colors) + 
  scale_shape_manual(values=EB_shapes)+ 
  scale_linetype_manual(values=EB_linetypes)

ggsave("../../../../MS Revision/Figures/CompareCurvesEb.pdf", 
       g2, width=3.25, height=3.25, units="in")


    