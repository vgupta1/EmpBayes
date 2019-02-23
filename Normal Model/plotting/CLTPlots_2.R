### New Robustness to Non-Normality plots
source("plottingUtils.R")
library(stringr)
library(purrr)
library(forcats)
library(modelr)
library(scam)

###
# First load up the various TV plots and take a look
###
dat_unif = read_csv("../Results/TV_plot_uniform2.csv_uniform_8675309.csv")
dat_unif$dist = "uniform"
dat_expon = read_csv("../Results/TV_plot_exponential2.csv_exponential_8675309.csv")
dat_expon$dist = "exponential"
# dat_bern = read_csv("../Results/TV_plot_bernoulli.csv_bernoulli_8675309.csv")
# dat_bern$dist = "bernoulli"
dat_pareto = read_csv("../Results/TV_plot_pareto2.csv_pareto_8675309.csv")
dat_pareto$dist = "pareto"
dat_t = read_csv("../Results/TV_plot_StudentT2.csv_t_8675309.csv")
dat_t$dist = "t"

dat.tv = rbind(dat_unif, dat_expon, dat_pareto, dat_t) #dat_bern
rm(dat_unif, dat_expon, dat_bern, dat_pareto, dat_t)

### Only subgaussian things
dat_unif = read_csv("../Results/TV_plot_uniform2.csv_uniform_8675309.csv")
dat_unif$dist = "uniform"
dat_bern = read_csv("../Results/TV_bernoulli_bernoulli_8675309.csv")
dat_bern$dist = "bernoulli"
dat_beta = read_csv("../Results/TV_Beta_beta_8675309.csv")
dat_beta$dist= "beta"
dat_beta2 = read_csv("../Results/TV_Beta2_beta2_8675309.csv")
dat_beta2$dist= "beta2"
dat_beta5 = read_csv("../Results/TV_Beta5_beta5_8675309.csv")
dat_beta5$dist= "beta5"
dat_arcsin = read_csv("../Results/TV_arcsin_arcsine_8675309.csv")
dat_arcsin$dist = "arcsine"

dat.tv = rbind(dat_unif, dat_bern, dat_beta, dat_beta2, dat_beta5, dat_arcsin)


dat.tv <- mutate(dat.tv, dist = as.factor(dist))
dat.tv <- dat.tv %>% mutate(Labels = fct_recode(dist, Exponential="exponential", Pareto="pareto", Uniform="uniform", 'Student-t'="t"))

g <- dat.tv %>%
  ggplot(aes(S, TV, group=Labels, color=Labels)) + 
  geom_point(aes(shape=Labels))  + geom_line(aes(linetype=Labels))+
  theme_minimal() + 
  theme(legend.title=element_blank(), 
        legend.position=c(.8, .8)) 

ggsave("../../../../MS Revision/Figures/TVConvergence.pdf",
       g, width=3.25, height=3.25, units="in")

### 
# Play a little bit to see if you can't smooth things out
###
fit <- function(df) { loess(TV ~S, data=df)}
fit2 <- function(df){
  scam(TV ~ s(S, k =8, bs = "mpd"), data = df)
}

dat.tv %>% filter(dist=="beta") %>%
  fit2(.)

dat.tv <- 
  dat.tv %>% group_by(dist) %>% 
  nest() %>%
  mutate(smoothModel = map(data, fit), 
         preds = map2(data, smoothModel, add_predictions), 
         monotoneModel = map(data, fit2), 
         predMon = map2(data, monotoneModel, add_predictions)) %>%
  select(-smoothModel, -monotoneModel) %>%
  unnest() %>%
  select(-S1, -TV1, -Labels1, -S2, -TV2, -Labels2)
  
dat.tv %>% filter(dist == "beta5") %>%
  ggplot(aes(S, pred, group=Labels, color=Labels)) + 
  geom_line(aes(linetype=Labels))+
  geom_point(aes(y=TV), size=1, color="green") +
  geom_line(aes(y=pred1), linetype="dotted", color="black") +
      theme_minimal() + 
  theme(legend.title=element_blank(), 
        legend.position=c(.8, .8)) 

#beta has some non-monotonicity problems
#beta2
#Beta 

###
#For each distribution load in data-set.
#Strip to just relative performance and combine.
load_dist_dat <- function(dist_type){
  dat_path = str_c("../Results/POAPCLT_plot___", 
                   dist_type, "_8675309.csv_full_200.csv")
  dat = read_csv(dat_path)
  
  dat <- clean_data(dat)
  dat.sum <- summarize_clt_dat(dat)
  names(dat.sum)[1] <- "S"
  
  dat.sum <- dat.sum %>% ungroup() %>%
    filter(Method %in%
             c("OR", "BoxStein", "OracleReg", "SteinReg")) %>%
    select(S, Method, Label, avg, up, down, std)
  dat.sum$dist = dist_type
  return(dat.sum)
}
#dists <- c("exponential", "t", "pareto", "uniform")  #drop bernoulli for now

dists <- c("beta2", "beta5", "beta", "uniform")
dat.dists <- map(dists, load_dist_dat)
dat.dists <- do.call(rbind, dat.dists)
dat.dists <- mutate(dat.dists, dist = as.factor(dist))

#merge with prevous to get relevant TV
dat.dists<- inner_join(dat.dists, dat.tv, by=c("S", "dist")) %>%
            mutate(dist = as.factor(dist))  

# form the diff of Oracle and Stein 
dat.dists %>% ungroup() %>% 
  filter(Method %in% c("OR", "BoxStein")) %>%
  select(S, TV, Method, dist, avg, pred, pred1) %>%
  spread(Method, avg) %>%
  mutate(diff = OR-BoxStein)  %>%
  ggplot(aes(pred1, diff, group=dist, color=dist)) + 
  geom_point() + geom_line() +
  theme_minimal() + 
  theme(legend.title=element_blank(), legend.position = "top") + 
  scale_y_continuous(labels=scales::percent) + 
  scale_x_reverse()

dat.dists %>% ungroup() %>% 
  filter(Method %in% c("OracleReg", "SteinReg")) %>%
  select(S, TV, Method, dist, avg, pred, pred1) %>%
  spread(Method, avg) %>%
  mutate(diff = OracleReg-SteinReg)  %>%
  ggplot(aes(pred1, diff, group=dist, color=dist)) + 
  geom_point() + geom_line() + 
  theme_minimal() + 
  theme(legend.title=element_blank(), legend.position = "top") + 
  scale_y_continuous(labels=scales::percent) + 
  scale_x_reverse()



