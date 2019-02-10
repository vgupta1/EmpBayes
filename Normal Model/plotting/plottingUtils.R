###  Standard code for setting up the plots in a consistent way.  
library(tidyverse)
library(ggplot2)
library(forcats)
library(scales)
library(showtext)
font_add("Times New Roman", "Times New Roman.ttf")
showtext_auto()

pd = position_dodge(.3)

clean_data <- function(dat){
  #Drop the stuff you don't want to deal with for simplicty
  dat <- dat %>% filter(!Method %in% c("DiracStein",
                                       "OR_MSE", "OracleReg_5", 
                                       "FWRO_Eps_.1", "SteinReg_5")
                        ) %>%
                  mutate(Method = factor(Method), 
                         Method = fct_relevel(Method, 
                                              "FullInfo", "OR", "SAA", 
                                              "BoxStein", "EB_MLE", "SURE_MSE", "K5_EB",
                                              "EB_MM", "HO_EB", "LOO_EB", 
                                              "OracleReg", 
                                              "SteinReg", "FWRO_Eps_.01", "K5_Reg",  
                                              "FWRO_Eps_.05", "HO_Reg", "LOO_Reg")
                           )

  ## Now add the clean labels
  dat <- dat %>% mutate(Label = fct_recode(Method, 
                              `Full-Info` = "FullInfo",
                              `EB OR` = "OR", 
                              `EB OPT` = "BoxStein", 
                              `EB MLE` = "EB_MLE", 
                              `EB MM` = "EB_MM", 
                              `EB SURE` = "SURE_MSE", 
                              `EB LOO` = "LOO_EB", 
                              `EB K5` = "K5_EB", 
                              `EB HO` = "HO_EB", 
                              `Reg OR` = "OracleReg", 
                              `Reg OPT` = "SteinReg",
                              `RO 1%` = "FWRO_Eps_.01", 
                              `RO 5%` = "FWRO_Eps_.05", 
                              `Reg LOO` = "LOO_Reg", 
                              `Reg HO` = "HO_Reg", 
                              `Reg K5`= "K5_Reg")
                 )
  #add a filter for bayes vs reg
  dat <- mutate(dat, isBayes = Method %in% c("OR", "BoxStein", "EB_MLE", "EB_MM", "SURE_MSE", "SAA", "HO_EB", "LOO_EB", "K5_EB"), 
                     isReg = Method %in% c("OracleReg", "SteinReg", "FWRO_Eps_.01", "FWRO_Eps_.05", "SAA", "HO_Reg", "LOO_Reg", "K5_Reg")
                )
    return(dat)  
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
              down = quantile(Ratio, .1), 
              avgTime = mean(time), 
              stdTime = sd(time), 
              upTime = quantile(time, .9), 
              downTime = quantile(time, .1))
  return(dat.sum)
}

##Does the most common things to plots... additional tuning
#often necessary
make_pretty <- function(plt, isEB){
  plt <- plt + theme_minimal(base_size=10) + 
    theme(legend.title=element_blank(), 
          text=element_text(family="Times New Roman")) + 
    scale_y_continuous(labels=scales::percent) + 
    scale_x_log10(labels=scales::comma) 

  if(isEB){
    plt <- plt + scale_shape_manual(values=EB_shapes) + 
      scale_color_manual(values=EB_colors) + 
      scale_linetype_manual(values=EB_linetypes)
  }else{
    plt <- plt + scale_shape_manual(values=Reg_shapes) + 
      scale_color_manual(values=Reg_colors) + 
      scale_linetype_manual(values=Reg_linetypes)
  }
      return(plt)
}

#####
#Most plots separate Eb and Reg.
#manually assign the colors, shapes and linetypes for consistency
#####
#Shapes
EB_shapes = c(16, 17, 15,  3,  7,  8, 5, 6, 0)
names(EB_shapes) <- c("SAA", "EB OR", "EB OPT", "EB MLE", "EB MM", "EB SURE", "EB LOO", "EB K5", "EB HO")
Reg_shapes = c(16, 17, 15,  3,  7,  8, 5, 6, 0)
names(Reg_shapes) <- c("SAA", "Reg OR", "Reg OPT", "RO 1%", "RO 5%", "Reg LOO", "Reg HO", "Reg K5", "XX")

#This is how you will set these
#scale_shape_manual(values=shapes) + scale_color_manual(values=colors)

#Colors
EB_colors <- hue_pal()(9)
Reg_colors<- hue_pal()(9)
names(EB_colors) <- c("SAA", "EB OR", "EB OPT", "EB MLE", "EB MM", "EB SURE", "EB LOO", "EB K5", "EB HO")
names(Reg_colors) <- c("SAA", "Reg OR", "Reg OPT", "RO 1%", "RO 5%", "Reg LOO", "Reg HO", "Reg K5", "XX")

#linetypes
EB_linetypes = linetype_pal()(9)
Reg_linetypes = linetype_pal()(9)
names(EB_linetypes) <- c("SAA", "EB OR", "EB OPT", "EB MLE", "EB MM", "EB SURE", "EB LOO", "EB K5", "EB HO")
names(Reg_linetypes) <- c("SAA", "Reg OR", "Reg OPT", "RO 1%", "RO 5%", "Reg LOO", "Reg HO", "Reg K5", "XX")
