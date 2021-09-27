##################################################
## Supplementary Simulation: Mixed sens/spec    ##
##################################################
## ****; Sept 12, 2020
## X subject to sensitivity, Y to specificity
## and vice versa

library(ggplot2)
library(tidyverse)
library(patchwork)
library(wesanderson)
library(ggpattern)
library(MASS)
setwd("~/Dropbox/ASD & EMRs/Code/Plots")


####################
# Define Functions #

## expit function
expit <- function(x){
  exp(x)/(1+exp(x))
}

## logit
logit <- function(x){
  log(x/(1-x))
}

## function to generate data under DAG Type I
genI <- function(n=2000,      ## sample size
                 Xprev=0.5,   ## prevalence of X
                 Yprev=0.5,   ## baseline (unexposed) prevalence of Y
                 b1=0,        ## logOR for X~Y; default is null
                 meanN00=2,   ## mean no. of visits when X=0, Y=0 
                 meanN01=2,   ## mean no. of visits when X=0, Y=1
                 meanN10=2,   ## mean no. of visits when X=1, Y=0
                 meanN11=2,   ## mean no. of visits when X=1, Y=1
                 sensX=0.8,   ## sensitivity of X
                 specX=0.99,  ## specificity of X (1.0 --> no overdiagnosing due to ip)
                 sensY=0.8,   ## sensitivity of X
                 specY=0.99,  ## specificity of X (1.0 --> no overdiagnosing due to ip)
                 negbinTheta=1
){
  ## set baseline
  b0 <- logit(Yprev)
  
  ## generate true X,Y
  X <- rbinom(n,1,Xprev)
  Y <- rbinom(n,1,expit(b0+b1*X))
  
  ## generate number of visits
  # Nvis <- 1+rpois(n, -1 +    ## -1 since we add a minimum 1 to the no. visits
  #                   meanN00*(1-X)*(1-Y)+
  #                   meanN01*(1-X)*Y+
  #                   meanN10*X*(1-Y)+
  #                   meanN11*X*Y) 
  Nvis <- 1+rnegbin(n, -1 +    ## -1 since we add a minimum 1 to the no. visits
                      meanN00*(1-X)*(1-Y)+
                      meanN01*(1-X)*Y+
                      meanN10*X*(1-Y)+
                      meanN11*X*Y,
                    negbinTheta)
  
  ## generate observed Xobs,Yobs
  Xobs <- as.numeric(rbinom(n,Nvis,sensX*X + (1-specX)*(1-X))>0)
  Yobs <- as.numeric(rbinom(n,Nvis,sensY*Y + (1-specY)*(1-Y))>0)
  
  
  ## fit models
  beta <- coef(glm(Yobs~Xobs,family=binomial))
  betaN <- coef(glm(Yobs~Xobs+Nvis,family=binomial))
  
  ## save
  res <- list(b0=c(beta[1],betaN[1]),b1=c(beta[2],betaN[2]))
  names(res[[1]]) <- names(res[[2]]) <- c("naive","adjusted")
  return(res)
}





## function to run simulation I
simI <- function(R=4000,      ## no. of iterations
                 n=2000,      ## sample size
                 Xprev=0.5,   ## prevalence of X
                 Yprev=0.5,   ## baseline (unexposed) prevalence of Y
                 b1=0,        ## logOR for X~Y; default is null
                 meanN00=2,   ## mean no. of visits when X=0, Y=0 
                 meanN01=2,   ## mean no. of visits when X=0, Y=1
                 meanN10=2,   ## mean no. of visits when X=1, Y=0
                 meanN11=2,   ## mean no. of visits when X=1, Y=1
                 sensX=0.8,   ## sensitivity of X
                 specX=0.99,  ## specificity of X (1.0 --> no overdiagnosing due to ip)
                 sensY=0.8,   ## sensitivity of X
                 specY=0.99,  ## specificity of X (1.0 --> no overdiagnosing due to ip)
                 negbinTheta=1
){
  
  ## initialize
  df0 <- df1 <- c()
  
  ## loop over R iterations
  for(rr in 1:R){
    
    ## generate data and fit models
    est <- genI(n=n,Xprev=Xprev,Yprev=Yprev,b1=b1,meanN00=meanN00,meanN01=meanN01,meanN10=meanN10,meanN11=meanN11,sensX=sensX,specX=specX,sensY=sensY,specY=specY,negbinTheta=negbinTheta)
    
    ## collect estimates
    df0 <- rbind(df0,est$b0)
    df1 <- rbind(df1,est$b1)
  }
  
  ## return results
  return(list(b0=df0,b1=df1))
}




## function to plot distributions of results 
make_box <- function(res_list,setting_list,title="",true_val=0,slope=TRUE,ylim=c(-0.5,0.5)){
  df <- c()
  for(rr in 1:length(res_list)){
    if(slope==TRUE){
      res <- gather(tibble(data.frame(res_list[[rr]]$b1)),"Est","b1") 
    }else{
      res <- gather(tibble(data.frame(res_list[[rr]]$b0)),"Est","b1") 
    }
    names(res)[2] <- "b"
    res$setting <- setting_list[rr]
    res$rr <- rr
    df <- rbind(df,res)
  }
  df$Est[df$Est=="naive"] <- " Naive"
  df$Est[df$Est=="adjusted"] <- "Adjusted"
  
  ## color the dif
  df$col <- "1"
  df$col[df$setting %in% c("(3)","(4)","(13)","(24)","(12)")] <- "2" ## these are non dif
  ## stripes for collider bias
  df$pat <- "1"
  df$pat[df$setting %in% c("(12)","(123)","(124)","(1234)")] <- "2" ## these are blank
  
  both <- ggplot(data=df,aes(x=reorder(setting,rr),y=b,fill=Est))+
    geom_boxplot(outlier.size = 0.5)+
    scale_x_discrete(labels=set_names_parse) +
    scale_fill_manual(values=c("gray75","gray98"))+
    geom_abline(intercept=true_val,slope=0,linetype="dashed")+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    ylim(ylim[1],ylim[2])+
    labs(x="", y="log OR")+
    ggtitle(title)
  
  naive <- ggplot(data=df[df$Est==" Naive",],aes(x=reorder(setting,rr),y=b))+
    # geom_boxplot(outlier.size = 0.5,fill="gray98")+
    geom_boxplot_pattern(aes(fill=col,pattern_spacing=pat),pattern="stripe",pattern_fill="white",outlier.size = 0.5,alpha=0.5)+
    scale_fill_manual(values=c(wes_red,"white"),guide=F)+
    scale_pattern_spacing_manual(values=c(1,0.015),guide=F)+
    scale_x_discrete(labels=set_names_parse) +
    # scale_fill_manual(values=c("gray75","gray98"))+
    geom_abline(intercept=true_val,slope=0,linetype="dashed")+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    ylim(ylim[1],ylim[2])+
    labs(x="", y="log OR")+
    ggtitle(paste0(title," - Naive") )
  
  adjusted <- ggplot(data=df[df$Est=="Adjusted",],aes(x=reorder(setting,rr),y=b))+
    # geom_boxplot(outlier.size = 0.5,fill="gray75")+
    geom_boxplot_pattern(aes(fill=col,pattern_spacing=pat),pattern="stripe",pattern_fill="white",outlier.size = 0.5,alpha=0.5)+
    scale_fill_manual(values=c(wes_red,"white"),guide=F)+
    scale_pattern_spacing_manual(values=c(1,0.015),guide=F)+
    scale_x_discrete(labels=set_names_parse) +
    # scale_fill_manual(values=c("gray75","gray98"))+
    geom_abline(intercept=true_val,slope=0,linetype="dashed")+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    ylim(ylim[1],ylim[2])+
    labs(x="", y="log OR")+
    ggtitle(paste0(title," - Adjusted") )
  
  vert <- ggplot(data=df,aes(x=reorder(setting,-rr),y=b))+
    # geom_boxplot(outlier.size = 0.5,fill="gray98")+
    ## ADDED COLOR AND  PATTERN
    geom_boxplot_pattern(aes(fill=col,pattern_spacing=pat),pattern="stripe",pattern_fill="white",outlier.size = 0.5,alpha=0.5)+
    scale_fill_manual(values=c(wes_red,"white"),guide=F)+
    scale_pattern_spacing_manual(values=c(1,0.015),guide=F)+
    scale_x_discrete(labels=set_names_parse[sort(1:length(res_list),decreasing=T)]) +
    # scale_fill_manual(values=c("gray75","gray98"))+
    geom_abline(intercept=true_val,slope=0,linetype="dashed")+
    theme_bw() +
    ylim(ylim[1],ylim[2])+
    labs(x="", y="log OR")+
    coord_flip()+
    facet_wrap(~ Est)+
    theme(strip.background = element_rect(color="white", fill="white", size=0, linetype="solid"),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    ggtitle(paste0(title) )
  
  
  
  return(list(both=both,naive=naive,adjusted=adjusted,vert=vert))
}


####################
# Run Simulation I #
nn <- 2000
RR <- 1000 #2000
xprop <- 0.5
yprop <- 0.25
set_names <- c("(3)","(4)","(13)","(24)","(12)","(34)","(23)","(14)","(134)","(234)","(123)","(124)","(1234)")


### (a) NULL; CROEN SIZE: MIXED
### 
resMIXEDa_1200 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=1   ,specX=1  ,sensY=1   ,specY=1  )
resMIXEDa_0030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=1.6,meanN11=1.6,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resMIXEDa_1030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resMIXEDa_0230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=1.6,meanN11=2.3,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resMIXEDa_1230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resMIXEDa_0004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=1.6,meanN11=1.6,sensX=1   ,specX=1  ,sensY=1.0 ,specY=0.99)
resMIXEDa_1004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=1   ,specX=1  ,sensY=1.0 ,specY=0.99)
resMIXEDa_0204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=1.6,meanN11=2.3,sensX=1   ,specX=1  ,sensY=1.0 ,specY=0.99)
resMIXEDa_1204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=1   ,specX=1  ,sensY=1.0 ,specY=0.99)
resMIXEDa_0034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=1.6,meanN11=1.6,sensX=0.75,specX=1.0,sensY=1.0 ,specY=0.99)
resMIXEDa_1034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=0.75,specX=1.0,sensY=1.0 ,specY=0.99)
resMIXEDa_0234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=1.6,meanN11=2.3,sensX=0.75,specX=1.0,sensY=1.0 ,specY=0.99)
resMIXEDa_1234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=0.75,specX=1.0,sensY=1.0 ,specY=0.99)

resMIXEDa_box <- make_box(list(resMIXEDa_0030,resMIXEDa_0004,resMIXEDa_1030,resMIXEDa_0204,resMIXEDa_1200,resMIXEDa_0034,resMIXEDa_0230,resMIXEDa_1004,resMIXEDa_1034,resMIXEDa_0234,resMIXEDa_1230,resMIXEDa_1204,resMIXEDa_1234),
                      set_names,"X-Sens, Y-Spec; Low Frequency",true_val=0,ylim=c(-0.5,0.5))
ggsave("resMIXEDa_box_both.pdf",    plot=resMIXEDa_box$both)
ggsave("resMIXEDa_box_naive.pdf",   plot=resMIXEDa_box$naive)
ggsave("resMIXEDa_box_adjusted.pdf",plot=resMIXEDa_box$adjusted)



### (c) NULL; ISRAELI SIZE; MIXED
### 
resMIXEDc_1200 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=1   ,specY=1  )
resMIXEDc_0030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resMIXEDc_1030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resMIXEDc_0230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resMIXEDc_1230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resMIXEDc_0004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1   ,specX=1  ,sensY=1.0 ,specY=0.99)
resMIXEDc_1004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1   ,specX=1  ,sensY=1.0 ,specY=0.99)
resMIXEDc_0204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1   ,specX=1  ,sensY=1.0 ,specY=0.99)
resMIXEDc_1204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=1.0 ,specY=0.99)
resMIXEDc_0034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=0.75,specX=1.0,sensY=1.0 ,specY=0.99)
resMIXEDc_1034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=0.75,specX=1.0,sensY=1.0 ,specY=0.99)
resMIXEDc_0234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=0.75,specX=1.0,sensY=1.0 ,specY=0.99)
resMIXEDc_1234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=0.75,specX=1.0,sensY=1.0 ,specY=0.99)

resMIXEDc_box <- make_box(list(resMIXEDc_0030,resMIXEDc_0004,resMIXEDc_1030,resMIXEDc_0204,resMIXEDc_1200,resMIXEDc_0034,resMIXEDc_0230,resMIXEDc_1004,resMIXEDc_1034,resMIXEDc_0234,resMIXEDc_1230,resMIXEDc_1204,resMIXEDc_1234),
                      set_names,"X-Sens, Y-Spec; High Frequency",true_val=0,ylim=c(-0.75,0.75))
ggsave("resMIXEDc_box_both.pdf",    plot=resMIXEDc_box$both)
ggsave("resMIXEDc_box_naive.pdf",   plot=resMIXEDc_box$naive)
ggsave("resMIXEDc_box_adjusted.pdf",plot=resMIXEDc_box$adjusted)


#########################



### (a) Non-NULL; CROEN SIZE; MIXED
### 
resNNMIXEDa_1200 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=1   ,specX=1  ,sensY=1   ,specY=1  )
resNNMIXEDa_0030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=1.6,meanN11=1.6,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resNNMIXEDa_1030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resNNMIXEDa_0230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=1.6,meanN11=2.3,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resNNMIXEDa_1230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resNNMIXEDa_0004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=1.6,meanN11=1.6,sensX=1   ,specX=1  ,sensY=1.0 ,specY=0.99)
resNNMIXEDa_1004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=1   ,specX=1  ,sensY=1.0 ,specY=0.99)
resNNMIXEDa_0204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=1.6,meanN11=2.3,sensX=1   ,specX=1  ,sensY=1.0 ,specY=0.99)
resNNMIXEDa_1204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=1   ,specX=1  ,sensY=1.0 ,specY=0.99)
resNNMIXEDa_0034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=1.6,meanN11=1.6,sensX=0.75,specX=1.0,sensY=1.0 ,specY=0.99)
resNNMIXEDa_1034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=0.75,specX=1.0,sensY=1.0 ,specY=0.99)
resNNMIXEDa_0234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=1.6,meanN11=2.3,sensX=0.75,specX=1.0,sensY=1.0 ,specY=0.99)
resNNMIXEDa_1234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=0.75,specX=1.0,sensY=1.0 ,specY=0.99)

resNNMIXEDa_box <- make_box(list(resNNMIXEDa_0030,resNNMIXEDa_0004,resNNMIXEDa_1030,resNNMIXEDa_0204,resNNMIXEDa_1200,resNNMIXEDa_0034,resNNMIXEDa_0230,resNNMIXEDa_1004,resNNMIXEDa_1034,resNNMIXEDa_0234,resNNMIXEDa_1230,resNNMIXEDa_1204,resNNMIXEDa_1234),
                        set_names,"X-Sens, Y-Spec; Low Frequency",true_val=0.5,ylim=c(-0.0,1.0))
ggsave("resNNMIXEDa_box_both.pdf",    plot=resNNMIXEDa_box$both)
ggsave("resNNMIXEDa_box_naive.pdf",   plot=resNNMIXEDa_box$naive)
ggsave("resNNMIXEDa_box_adjusted.pdf",plot=resNNMIXEDa_box$adjusted)



### (c) Non-NULL; ISRAELI SIZE; MIXED
### 
resNNMIXEDc_1200 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=1   ,specY=1  )
resNNMIXEDc_0030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resNNMIXEDc_1030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resNNMIXEDc_0230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resNNMIXEDc_1230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resNNMIXEDc_0004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1   ,specX=1  ,sensY=1.0 ,specY=0.99)
resNNMIXEDc_1004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1   ,specX=1  ,sensY=1.0 ,specY=0.99)
resNNMIXEDc_0204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1   ,specX=1  ,sensY=1.0 ,specY=0.99)
resNNMIXEDc_1204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=1.0 ,specY=0.99)
resNNMIXEDc_0034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=0.75,specX=1.0,sensY=1.0 ,specY=0.99)
resNNMIXEDc_1034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=0.75,specX=1.0,sensY=1.0 ,specY=0.99)
resNNMIXEDc_0234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=0.75,specX=1.0,sensY=1.0 ,specY=0.99)
resNNMIXEDc_1234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=0.75,specX=1.0,sensY=1.0 ,specY=0.99)

resNNMIXEDc_box <- make_box(list(resNNMIXEDc_0030,resNNMIXEDc_0004,resNNMIXEDc_1030,resNNMIXEDc_0204,resNNMIXEDc_1200,resNNMIXEDc_0034,resNNMIXEDc_0230,resNNMIXEDc_1004,resNNMIXEDc_1034,resNNMIXEDc_0234,resNNMIXEDc_1230,resNNMIXEDc_1204,resNNMIXEDc_1234),
                        set_names,"X-Sens, Y-Spec; High Frequency",true_val=0.5,ylim=c(-0.0,1.0))
ggsave("resNNMIXEDc_box_both.pdf",    plot=resNNMIXEDc_box$both)
ggsave("resNNMIXEDc_box_naive.pdf",   plot=resNNMIXEDc_box$naive)
ggsave("resNNMIXEDc_box_adjusted.pdf",plot=resNNMIXEDc_box$adjusted)



##################
## Make figures
fig2 <- (resMIXEDa_box$vert)/ 
  (resMIXEDc_box$vert) + plot_annotation(tag_levels = 'A')
ggsave("fig2_mixed.pdf",plot=fig2,width=7,height=8)

fig3 <- (resNNMIXEDa_box$vert)/ 
  (resNNMIXEDc_box$vert) + plot_annotation(tag_levels = 'A')
ggsave("fig3_mixed.pdf",plot=fig3,width=7,height=8)





##########################
## X-spec, Y-sens


### (a) NULL; CROEN SIZE: MIXED
### 
resMIXED2a_1200 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=1   ,specX=1  ,sensY=1   ,specY=1  )
resMIXED2a_0030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=1.6,meanN11=1.6,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resMIXED2a_1030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resMIXED2a_0230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=1.6,meanN11=2.3,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resMIXED2a_1230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resMIXED2a_0004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=1.6,meanN11=1.6,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resMIXED2a_1004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resMIXED2a_0204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=1.6,meanN11=2.3,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resMIXED2a_1204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resMIXED2a_0034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=1.6,meanN11=1.6,sensX=1.0,specX=0.99,sensY=0.75,specY=1.0)
resMIXED2a_1034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=1.0,specX=0.99,sensY=0.75,specY=1.0)
resMIXED2a_0234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=1.6,meanN11=2.3,sensX=1.0,specX=0.99,sensY=0.75,specY=1.0)
resMIXED2a_1234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=1.0,specX=0.99,sensY=0.75,specY=1.0)

resMIXED2a_box <- make_box(list(resMIXED2a_0030,resMIXED2a_0004,resMIXED2a_1030,resMIXED2a_0204,resMIXED2a_1200,resMIXED2a_0034,resMIXED2a_0230,resMIXED2a_1004,resMIXED2a_1034,resMIXED2a_0234,resMIXED2a_1230,resMIXED2a_1204,resMIXED2a_1234),
                          set_names,"X-Spec, Y-Sens; Low Frequency",true_val=0,ylim=c(-0.5,0.5))
ggsave("resMIXED2a_box_both.pdf",    plot=resMIXED2a_box$both)
ggsave("resMIXED2a_box_naive.pdf",   plot=resMIXED2a_box$naive)
ggsave("resMIXED2a_box_adjusted.pdf",plot=resMIXED2a_box$adjusted)



### (c) NULL; ISRAELI SIZE; MIXED2
### 
resMIXED2c_1200 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=1   ,specY=1  )
resMIXED2c_0030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resMIXED2c_1030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resMIXED2c_0230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resMIXED2c_1230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resMIXED2c_0004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resMIXED2c_1004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resMIXED2c_0204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resMIXED2c_1204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resMIXED2c_0034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1.0,specX=0.99,sensY=0.75,specY=1.0)
resMIXED2c_1034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1.0,specX=0.99,sensY=0.75,specY=1.0)
resMIXED2c_0234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1.0,specX=0.99,sensY=0.75,specY=1.0)
resMIXED2c_1234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1.0,specX=0.99,sensY=0.75,specY=1.0)

resMIXED2c_box <- make_box(list(resMIXED2c_0030,resMIXED2c_0004,resMIXED2c_1030,resMIXED2c_0204,resMIXED2c_1200,resMIXED2c_0034,resMIXED2c_0230,resMIXED2c_1004,resMIXED2c_1034,resMIXED2c_0234,resMIXED2c_1230,resMIXED2c_1204,resMIXED2c_1234),
                          set_names,"X-Spec, Y-Sens; High Frequency",true_val=0,ylim=c(-0.75,0.75))
ggsave("resMIXED2c_box_both.pdf",    plot=resMIXED2c_box$both)
ggsave("resMIXED2c_box_naive.pdf",   plot=resMIXED2c_box$naive)
ggsave("resMIXED2c_box_adjusted.pdf",plot=resMIXED2c_box$adjusted)


#########################



### (a) Non-NULL; CROEN SIZE; MIXED2
### 
resNNMIXED2a_1200 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=1   ,specX=1  ,sensY=1   ,specY=1  )
resNNMIXED2a_0030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=1.6,meanN11=1.6,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resNNMIXED2a_1030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resNNMIXED2a_0230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=1.6,meanN11=2.3,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resNNMIXED2a_1230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resNNMIXED2a_0004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=1.6,meanN11=1.6,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resNNMIXED2a_1004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resNNMIXED2a_0204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=1.6,meanN11=2.3,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resNNMIXED2a_1204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resNNMIXED2a_0034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=1.6,meanN11=1.6,sensX=1.0,specX=0.99,sensY=0.75,specY=1.0)
resNNMIXED2a_1034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=1.0,specX=0.99,sensY=0.75,specY=1.0)
resNNMIXED2a_0234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=1.6,meanN11=2.3,sensX=1.0,specX=0.99,sensY=0.75,specY=1.0)
resNNMIXED2a_1234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=1.0,specX=0.99,sensY=0.75,specY=1.0)

resNNMIXED2a_box <- make_box(list(resNNMIXED2a_0030,resNNMIXED2a_0004,resNNMIXED2a_1030,resNNMIXED2a_0204,resNNMIXED2a_1200,resNNMIXED2a_0034,resNNMIXED2a_0230,resNNMIXED2a_1004,resNNMIXED2a_1034,resNNMIXED2a_0234,resNNMIXED2a_1230,resNNMIXED2a_1204,resNNMIXED2a_1234),
                            set_names,"X-Spec, Y-Sens; Low Frequency",true_val=0.5,ylim=c(-0.0,1.0))
ggsave("resNNMIXED2a_box_both.pdf",    plot=resNNMIXED2a_box$both)
ggsave("resNNMIXED2a_box_naive.pdf",   plot=resNNMIXED2a_box$naive)
ggsave("resNNMIXED2a_box_adjusted.pdf",plot=resNNMIXED2a_box$adjusted)



### (c) Non-NULL; ISRAELI SIZE; MIXED2
### 
resNNMIXED2c_1200 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=1   ,specY=1  )
resNNMIXED2c_0030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resNNMIXED2c_1030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resNNMIXED2c_0230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resNNMIXED2c_1230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resNNMIXED2c_0004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resNNMIXED2c_1004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resNNMIXED2c_0204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resNNMIXED2c_1204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resNNMIXED2c_0034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1.0,specX=0.99,sensY=0.75,specY=1.0)
resNNMIXED2c_1034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1.0,specX=0.99,sensY=0.75,specY=1.0)
resNNMIXED2c_0234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1.0,specX=0.99,sensY=0.75,specY=1.0)
resNNMIXED2c_1234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1.0,specX=0.99,sensY=0.75,specY=1.0)

resNNMIXED2c_box <- make_box(list(resNNMIXED2c_0030,resNNMIXED2c_0004,resNNMIXED2c_1030,resNNMIXED2c_0204,resNNMIXED2c_1200,resNNMIXED2c_0034,resNNMIXED2c_0230,resNNMIXED2c_1004,resNNMIXED2c_1034,resNNMIXED2c_0234,resNNMIXED2c_1230,resNNMIXED2c_1204,resNNMIXED2c_1234),
                            set_names,"X-Spec, Y-Sens; High Frequency",true_val=0.5,ylim=c(-0.0,1.0))
ggsave("resNNMIXED2c_box_both.pdf",    plot=resNNMIXED2c_box$both)
ggsave("resNNMIXED2c_box_naive.pdf",   plot=resNNMIXED2c_box$naive)
ggsave("resNNMIXED2c_box_adjusted.pdf",plot=resNNMIXED2c_box$adjusted)



##################
## Make figures
fig2_2 <- (resMIXED2a_box$vert)/ 
  (resMIXED2c_box$vert) + plot_annotation(tag_levels = 'A')
ggsave("fig2_mixed2.pdf",plot=fig2_2,width=7,height=8)

fig3_2 <- (resNNMIXED2a_box$vert)/ 
  (resNNMIXED2c_box$vert) + plot_annotation(tag_levels = 'A')
ggsave("fig3_mixed2.pdf",plot=fig3_2,width=7,height=8)

