###  Standard code for setting up the plots in a consistent way.  
library(tidyverse)
library(ggplot2)
library(forcats)
library(showtext)
font_add("Times New Roman", "Times New Roman.ttf")
showtext_auto()

pd = position_dodge(.3)

#VG fix this up after final run. 
clean_data <- function(dat){
  #Drop the stuff you don't want to deal with for simplicty
  dat <- dat %>% filter(!Method %in% c("DiracStein", "LOO_5", "LOO_10",
                                       "OR_MSE", "OracleReg_5", "OracleReg_10",
                                       "RO_Eps_.1", "SteinReg_5", "SteinReg_10")
                        ) %>%
                  mutate(Method = factor(Method), 
                         Method = fct_relevel(Method, 
                                              "FullInfo", "OR", 
                                              "BoxStein", "EB_MLE", "SURE_MSE", 
                                              "OracleReg", 
                                              "SteinReg", "LOO", "RO_Eps_.01", 
                                              "SAA", 
                                              "EB_MM", "RO_Eps_.05")
                           )

  ## Now add the clean labels
  dat <- dat %>% mutate(Label = fct_recode(Method, 
                              `Full-Info` = "FullInfo",
                              `EB OR` = "OR", 
                              `EB OPT` = "BoxStein", 
                              `EB MLE` = "EB_MLE", 
                              `EB MM` = "EB_MM", 
                              `EB SURE` = "SURE_MSE", 
                              `Reg OR` = "OracleReg", 
                              `Reg OPT` = "SteinReg",
                              `RO 1%` = "RO_Eps_.01", 
                              `RO 5%` = "RO_Eps_.05")
                 )
  #add a filter for bayes vs reg
  dat <- mutate(dat, isBayes = Method %in% c("OR", "BoxStein", "EB_MLE", "EB_MM", "SURE_MSE", "SAA"), 
                     isReg = Method %in% c("OracleReg", "SteinReg", "LOO", "RO_Eps_.01", "RO_Eps_.05", "SAA")
                )
    return(dat)  
}

##Does the most common things to plots... additional tuning
#often necessary
make_pretty <- function(plt){
  plt <- plt + theme_minimal(base_size=10) + 
    theme(legend.title=element_blank(), 
          text=element_text(family="Times New Roman")) + 
    scale_y_continuous(labels=scales::percent) + 
    scale_x_log10(labels=scales::comma) 
    return(plt)
}

##most general summarizing for data... 
##
summarize_dat <- function(dat, ratio=TRUE){
  #Re-express everything as a fraction of Oracle value
    t <- dat %>% filter(Method == "FullInfo") %>%
    select(Run, n, thetaVal)
  dat <- inner_join(dat, t, by=c("Run", "n")) %>%
    rename(thetaVal = thetaVal.x, 
           FullInfo = thetaVal.y)  
  if (ratio) {
    dat <- mutate(dat, Ratio=thetaVal/FullInfo)
  }else
  {
    dat <- mutate(dat, Ratio=thetaVal)
  }

  dat.sum <- dat %>% group_by(n, Method, Label, isBayes, isReg) %>%
    summarise(avg = mean(Ratio), 
              avgTau0 = mean(tau0), 
              std  = sd(Ratio), 
              stdTau0 = sd(tau0), 
              up = quantile(Ratio, .9), 
              down = quantile(Ratio, .1), 
              upTau0 = quantile(tau0, .9), 
              downTau0 = quantile(tau0, .1), 
              avgTime = mean(time), 
              upTime = quantile(time, .9), 
              downTime = quantile(time, .1))
  
}

  
summarize_clt_dat <- function(dat){
  #Re-express everything as a fraction of full-info
  t <- dat %>% filter(Method == "FullInfo") %>%
    select(Run, N, thetaVal)
  
  dat <- inner_join(dat, t, by=c("Run", "N")) %>%
    rename(thetaVal = thetaVal.x, 
           FullInfo = thetaVal.y)  %>%
    mutate(Ratio = thetaVal/FullInfo)
  
  
  #Summarize the data across runs
  dat.sum <- dat %>% group_by(N, Method, Label, isBayes, isReg) %>%
    summarise(avg = mean(Ratio), 
              avgTau0 = mean(tau0), 
              std  = sd(Ratio), 
              stdTau0 = sd(tau0), 
              up = quantile(Ratio, .9), 
              down = quantile(Ratio, .1))
  return(dat.sum)
}


