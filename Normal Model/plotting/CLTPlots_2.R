### New Robustness to Non-Normality plots
source("plottingUtils.R")
library(stringr)
library(purrr)

###
# First load up the various TV plots and take a look
###
dat_unif = read_csv("../Results/TV_plot_uniform_small.csv_uniform_8675309.csv")
dat_unif$dist = "uniform"
dat_expon = read_csv("../Results/TV_plot_exponential_small.csv_exponential_8675309.csv")
dat_expon$dist = "exponential"
dat_bern = read_csv("../Results/TV_plot_bernoulli.csv_bernoulli_8675309.csv")
dat_bern$dist = "bernoulli"
dat_pareto = read_csv("../Results/TV_plot_pareto_small.csv_pareto_8675309.csv")
dat_pareto$dist = "pareto"
dat_t = read_csv("../Results/TV_plot_StudentT_small.csv_t_8675309.csv")
dat_t$dist = "t"

dat.tv = rbind(dat_unif, dat_expon, dat_bern, dat_pareto, dat_t)
rm(dat_unif, dat_expon, dat_bern, dat_pareto, dat_t)
dat.tv <- mutate(dat.tv, dist = as.factor(dist))

dat.tv %>% filter(dist != "bernoulli") %>%
  ggplot(aes(S, TV, group=dist, color=dist)) + 
  geom_point()  + geom_smooth()+
  theme_minimal()

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
dists <- c("exponential", "t", "pareto", "uniform")  #drop bernoulli for now
dat.dists <- map(dists, load_dist_dat)
dat.dists <- do.call(rbind, dat.dists)
dat.dists <- mutate(dat.dists, dist = as.factor(dist))

#merge with prevous to get relevant TV
dat.dists<- inner_join(dat.dists, dat.tv, by=c("S", "dist")) %>%
            mutate(dist = as.factor(dist))  

# form the ratio of Oracle and Stein 
dat.dists %>% ungroup() %>% 
  filter(Method %in% c("OR", "BoxStein")) %>%
  select(S, TV, Method, dist, avg) %>%
  spread(Method, avg) %>%
  mutate(diff = OR-BoxStein)  %>%
  ggplot(aes(TV, diff, group=dist, color=dist)) + 
  geom_point() + geom_line() + 
  theme_minimal() + 
  theme(legend.title=element_blank(), legend.position = "top") + 
  scale_y_continuous(labels=scales::percent) + 
  scale_x_reverse()
