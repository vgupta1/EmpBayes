### Three part CLT Analysis and Plots
library(tidyverse)
library(ggplot2)
library(forcats)

library(extrafont)
font = "Times New Roman"
font = "CM Roman"

dat = read_csv("../Results/3PartCLT_plot__3partCLT_65536_8675309.csv_full_200.csv")
#spath = "../../../../OR Submission_1/Figures/3PartCLT.pdf"

dat$Method = as.factor(dat$Method)

#Reorder the levels in place to make plots consistent.
#This may need to be supplemented with specific color specs
dat <- mutate(dat, Method = fct_relevel(Method, 
                                        "FullInfo", "OR", "BoxStein", "DiracStein", 
                                        "OR_MSE", "EB_MLE", "EB_MM", "SURE_MSE", "SAA"))

#Create a new column of "pretty" labels
dat <- mutate(dat, Label = fct_recode(Method, 
                                      `EB OPT` = "BoxStein", 
                                      `EB OPT Dirac` = "DiracStein", 
                                      `EB MLE` = "EB_MLE", 
                                      `EB MM` = "EB_MM", 
                                      `MSE OR`= "OR_MSE", 
                                      `EB OR` = "OR", 
                                      `EB SURE` = "SURE_MSE")
)

#Re-express everything as a fraction of full-info
t <- dat %>% filter(Method == "FullInfo") %>%
  select(Run, N, thetaVal)

dat <- inner_join(dat, t, by=c("Run", "N")) %>%
  rename(thetaVal = thetaVal.x, 
         FullInfo = thetaVal.y)  %>%
  mutate(Ratio = thetaVal/FullInfo)


#Summarize the data across runs
dat.sum <- dat %>% group_by(N, Method, Label) %>%
  summarise(avg = mean(Ratio), 
            avgTau0 = mean(tau0), 
            std  = sd(Ratio), 
            stdTau0 = sd(tau0), 
            up = quantile(Ratio, .9), 
            down = quantile(Ratio, .1))

pd = position_dodge(.3)
g <- dat.sum %>% filter(Method %in% c("BoxStein", "OR")) %>%
  ggplot(aes(N, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd) 

+ 
  theme_minimal(base_size=11) +
  theme(legend.position = c(.5, .1), 
        legend.title=element_blank(), 
        text=element_text(family=font), 
        legend.direction = "horizontal", 
        legend.justification = "center") + 
  scale_x_log10(labels=scales::comma) + 
  scale_y_continuous(labels=scales::percent, limits=c(-.15, 1)) +
  ylab("(%) of Full-Info")
g

ggsave(spath, 
       g, width=3.25, height=3.25, units="in")


#############
#is it because people are overshrinking?
g1 <- dat %>% filter(Method %in% c("BoxStein", "EB_MM", "EB_MLE", "SURE_MSE", "OR"), 
                     n == 131072) %>%
  ggplot(aes(tau0)) + 
  geom_density(aes(group=Label, fill=Label), 
               alpha=.5, linetype="blank") + 
  xlab("Fitted tau") + ylab("") + 
  theme_minimal(base_size=11) + 
  theme(text=element_text(family=font),
        legend.title=element_blank(), 
        legend.position=c(.8, .8))

##Create an auxiliary table to do the labeling
dat.labels <- tribble(
  ~Method, ~x, ~y, 
  "EB OR", .9, 15, 
  "EB OPT", 1.1, 4.8, 
  "EB SURE", 1.8, 11, 
  "EB MM", 3, 3, 
  "EB MLE", 5, 6
)

g1 <- g1 + theme(legend.position="none") + 
  geom_text(data=dat.labels, 
            aes(x, y, label=Method), 
            family=font, 
            size=3, hjust="c", vjust="m")

ggsave("../../../../OR Submission_1/Figures/3PartTaus.pdf", 
       g1, width=3.25, height=3.25, units="in")


######
##Plot the relative densities
#######
dat = read_csv("../Results/3PartDensities_131072_20.00001.csv")
dat <- dat %>% mutate(Type=as.factor(mu), 
                      Type = fct_recode(Type, 
                                        Low="0", Med="0.3", High="1")
)

##First create plot for baseline
g <- dat %>%
  ggplot(aes(muhat, fill=Type)) +
  geom_density(linetype="blank", alpha=.5) + 
  geom_vline(xintercept = quantile(dat$muhat, .95), 
             linetype="dotted") +
  xlab("") + ylab("") + 
  theme_minimal(base_size=11) +
  theme(legend.title=element_blank(), 
        text = element_text(family=font), 
        legend.position=c(.2,.8))

props_table <-
  tribble(
    ~Type, ~Percent,
    "Low",  "94%",     
    "Med",  "0%",
    "High", "6%"
  )

g<- g + annotation_custom(tableGrob(props_table, rows=NULL, 
                                    theme=ttheme_minimal(base_size=8)), 
                          xmin=5, xmax=12, ymin=.5, ymax=.75)

ggsave("../../../../OR Submission_1/Figures/3PartDensitySAA.pdf", 
       g, width=3.25, height=3.25, units="in")

####
#Next plot tauOR
g <- dat %>%
  ggplot(aes(rOR, fill=Type)) +
  geom_density(linetype="blank", alpha=.5) + 
  geom_vline(xintercept = quantile(dat$rOR, .95), 
             linetype="dotted") +
  xlab("") + ylab("") + 
  theme_minimal(base_size=11) +
  theme(legend.title=element_blank(), 
        text = element_text(family=font), 
        legend.position=c(.2,.8))

props_table <-
  tribble(
    ~Type, ~Percent,
    "Low",  "1.7%",     
    "Med",  "4.4%",
    "High", "94.0%"
  )

g<- g + annotation_custom(tableGrob(props_table, rows=NULL, 
                                    theme=ttheme_minimal(base_size=8)), 
                          xmin=1.7, xmax=3, ymin=.62, ymax=.83)

ggsave("../../../../OR Submission_1/Figures/3PartDensityOr.pdf", 
       g, width=3.25, height=3.25, units="in")


###  Single plot with dependence in tau
dat = read_csv("../Results/tauDependence.csv")
g<- dat %>% 
  ggplot(aes(tau, thetaVal)) + 
  geom_point() + geom_line() + 
  theme_minimal(base_size=11) +
  theme(text=element_text(family=font size=11)) + 
  xlab("Tau") + ylab("Obj.")

ggsave("../../../../OR Submission_1/Figures/3PartTauDep.pdf", 
       g, width=3.25, height=3.25, units="in")


#####  Single Plot with dependence in Gamma
dat = read_csv("../Results/GammaDependence.csv")
g<- dat %>% 
  ggplot(aes(Gamma, thetaVal)) + 
  geom_point() + geom_line() + 
  theme_minimal(base_size=11) +
  theme(text=element_text(family=font, size=11)) + 
  xlab("Gamma") + ylab("Obj.")

ggsave("../../../../OR Submission_1/Figures/3PartGammaDep.pdf", 
       g, width=3.25, height=3.25, units="in")




