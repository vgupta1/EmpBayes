### Three part Analysis and Plots for the EB methods

source("plottingUtils.R")

#VG Data updated after fixing bug in SteinEB Feb 2019
#dat = read_csv("../Results/3Part_plot___3part_0.0_0.1_20.001_1.0_1.0_20.001_0.3_4.0_20.001_8675309000.csv_full_200.csv")
dat = read_csv("../Results/3Partb__8675309.csv_full_200.csv")

dat <- clean_data(dat)
dat.sum <- summarize_dat(dat)

#A small version for the paper
g <- filter(dat.sum, isBayes, 
            !Method %in% c("EB_MM", "K5_EB", "LOO_EB", "HO_EB") )%>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd) 
g <- make_pretty(g) + 
  theme(legend.position = c(.65, .3), 
        legend.text=element_text(size=8),
        legend.direction = "horizontal", 
        legend.justification = "center") + 
  ylab("(%) of Full-Info") +
  guides(shape=guide_legend(nrow=3, byrow=TRUE))
g

ggsave("../../../../MS Revision/Figures/3Part20b.pdf", 
       g, width=3, height=3, units="in")

#A Bigger version for the appendix
g <- filter(dat.sum, isBayes, 
            !Method %in% c("LOO_EB", "K5_EB", "HO_EB")) %>%
  ggplot(aes(n, avg, group=Label, color=Label)) + 
  geom_point(aes(shape=Label), position=pd) + 
  geom_line(aes(linetype=Label), position=pd) + 
  geom_errorbar(aes(ymin=down, ymax=up), position=pd)
g <- make_pretty(g) + 
  theme(legend.position = c(.7, .4), 
        legend.direction = "horizontal", 
        legend.justification = "center") + 
  ylab("(%) of Full-Info")
g

ggsave("../../../../MS Revision/Figures/3Part20Bigb.pdf",
       g, width=6, height=3.25, units="in")


#############
## densities of the fitted taus
g1 <- filter(dat, isBayes, n == 131072, 
             !Method %in% c("LOO_EB", "K5_EB", "HO_EB", "SAA")) %>%
  ggplot(aes(tau0)) + 
  geom_density(aes(group=Label, fill=Label), 
               alpha=.5, linetype="blank")

g1 <- make_pretty(g1) + 
  scale_x_continuous()+ scale_y_continuous() + 
  theme(legend.position="none") + 
  xlab("Fitted tau") + ylab("")

##Create an auxiliary table to do the labeling
dat.labels <- tribble(
  ~Method, ~x, ~y, 
  "EB OR", .9, 15, 
  "EB OPT", 1.1, 4.8, 
  "EB SURE", 1.8, 11, 
  "EB MM", 3, 3, 
  "EB MLE", 4.8, 6
  )

g1 <- g1 + geom_text(data=dat.labels, 
              aes(x, y, label=Method), 
              family="Times New Roman", 
              size=2.5, hjust="c", vjust="m")
g1 <- g1 + xlab(expression(tau))


ggsave("../../../../MS Revision/Figures/3PartTausb.pdf", 
       g1, width=2.275, height=2.275, units="in")


######
##Plot the relative densities
#######
library(gridExtra)
dat = read_csv("../Results/3PartDensities_131072_20.00001.csv")
dat <- dat %>% mutate(Type=as.factor(mu), 
                      Type = fct_recode(Type, 
                                 Low="0", Med="0.3", High="1")
                      )

##First create plot for baseline
g <- dat %>%
  ggplot(aes(muhat, fill=Type, linetype=Type, color=Type)) +
  geom_density(alpha=.5) + 
  geom_vline(xintercept = quantile(dat$muhat, .95), 
             linetype="dotted") +
  xlab("") + ylab("") + 
  theme_minimal(base_size=10) +
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
                      xmin=6.5, xmax=11.5, ymin=.5, ymax=.75)

ggsave("../../../../OR Submission_1/Figures/3PartDensitySAA.pdf", 
       g, width=2.75, height=2.75, units="in")

####
#Next plot tauOR
g <- dat %>%
  ggplot(aes(rOR, fill=Type, linetype=Type, color=Type)) +
  geom_density(alpha=.5) + 
  geom_vline(xintercept = quantile(dat$rOR, .95), 
             linetype="dotted") +
  xlab("") + ylab("") + 
  theme_minimal(base_size=10) +
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
                          xmin=2, xmax=3, ymin=.62, ymax=.83)

ggsave("../../../../OR Submission_1/Figures/3PartDensityOr.pdf", 
       g, width=2.75, height=2.75, units="in")



###Finally do TauMLE
g <- dat %>%
  ggplot(aes(rMLE, fill=Type, linetype=Type, color=Type)) +
  geom_density(alpha=.5) + 
  geom_vline(xintercept = quantile(dat$rMLE, .95), 
             linetype="dotted") +
  xlab("") + ylab("") + 
  theme_minimal(base_size=10) +
  theme(legend.title=element_blank(), 
        text = element_text(family=font), 
        legend.position=c(.2,.8))

props_table <-
  tribble(
    ~Type, ~Percent,
    "Low",  "0.0%",     
    "Med",  "59.4%",
    "High", "40.6%"
  )

g<- g + annotation_custom(tableGrob(props_table, rows=NULL, 
                                    theme=ttheme_minimal(base_size=8)), 
                          xmin=.7, xmax=1, ymin=4, ymax=6)

ggsave("../../../../OR Submission_1/Figures/3PartDensityMLE.pdf", 
       g, width=2.75, height=2.75, units="in")


## Back track out the tau values for hte paper
dat %>% mutate(tauOR = vs/rOR  * (muhat-rOR), 
               tauMLE =vs/rMLE * (muhat-rMLE), 
               tauMM = vs/rMM  * (muhat-rMM)) %>% 
  select(muhat, tauOR, tauMLE, tauMM) %>% head()


####





