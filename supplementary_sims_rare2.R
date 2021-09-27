#################################################################
## Supplementary Simulation: Rare Outcome, adjusted for N^2    ##
#################################################################
## ****; Sept 14, 2020
## Rare Outcome-- wasnt fully adjusted for by N
## here we adjusted for N^2 as well

library(ggplot2)
library(tidyverse)
library(patchwork)
library(MASS)
library(truncnorm)
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
genI2 <- function(n=2000,      ## sample size
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
                 specY=0.99,   ## specificity of X (1.0 --> no overdiagnosing due to ip)
                 negbinTheta=1 ## theta parameter for negative binomial visit count
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
                    negbinTheta) ## 1
  
  ## generate observed Xobs,Yobs
  Xobs <- as.numeric(rbinom(n,Nvis,sensX*X + (1-specX)*(1-X))>0)
  Yobs <- as.numeric(rbinom(n,Nvis,sensY*Y + (1-specY)*(1-Y))>0)
  
  Nvis2 <- Nvis^2-mean(Nvis^2)
  Nvis3 <- Nvis^3-mean(Nvis^3)
  logNvis <- log(Nvis)
  ## fit models
  beta <- coef(glm(Yobs~Xobs,family=binomial))
  betaN <- coef(glm(Yobs~Xobs+Nvis+Nvis2,family=binomial))
  
  ## save
  res <- list(b0=c(beta[1],betaN[1]),b1=c(beta[2],betaN[2]))
  names(res[[1]]) <- names(res[[2]]) <- c("naive","adjusted")
  return(res)
}


## function to run simulation I
simI2 <- function(R=2000,      ## no. of iterations
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
    est <- genI2(n=n,Xprev=Xprev,Yprev=Yprev,b1=b1,meanN00=meanN00,meanN01=meanN01,meanN10=meanN10,meanN11=meanN11,sensX=sensX,specX=specX,sensY=sensY,specY=specY,negbinTheta=negbinTheta)
    
    ## collect estimates
    df0 <- rbind(df0,est$b0)
    df1 <- rbind(df1,est$b1)
  }
  
  ## return results
  return(list(b0=df0,b1=df1))
}


## rare outcome
set.seed(1000) ## set seed

resRare2Id_1200 <- simI2(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=1   ,specY=1  )
resRare2Id_0030 <- simI2(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resRare2Id_1030 <- simI2(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resRare2Id_0230 <- simI2(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resRare2Id_1230 <- simI2(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resRare2Id_0004 <- simI2(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resRare2Id_1004 <- simI2(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resRare2Id_0204 <- simI2(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resRare2Id_1204 <- simI2(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resRare2Id_0034 <- simI2(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resRare2Id_1034 <- simI2(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resRare2Id_0234 <- simI2(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resRare2Id_1234 <- simI2(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)

resRare2Id_box <- make_box(list(resRare2Id_0030,resRare2Id_0004,resRare2Id_1030,resRare2Id_0204,resRare2Id_1200,resRare2Id_0034,resRare2Id_0230,resRare2Id_1004,resRare2Id_1034,resRare2Id_0234,resRare2Id_1230,resRare2Id_1204,resRare2Id_1234),
                          set_names,"Spec=0.99; Hi-Freq; Rare outcome; N^2",true_val=0,ylim=c(-1,1))
ggsave("resRare2Id_box_both.pdf",    plot=resRare2Id_box$both)
ggsave("resRare2Id_box_naive.pdf",   plot=resRare2Id_box$naive)
ggsave("resRare2Id_box_adjusted.pdf",plot=resRare2Id_box$adjusted)


# figRare2 <- resRare2Id_box$naive+resRare2Id_box$adjusted + plot_annotation(tag_levels = 'A')
ggsave("figRare2.pdf",plot=resRare2Id_box$vert,width=7,height=4)


## rare exposure
set.seed(1000) ## set seed

resRareX2Id_1200 <- simI2(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=1   ,specY=1  )
resRareX2Id_0030 <- simI2(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resRareX2Id_1030 <- simI2(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resRareX2Id_0230 <- simI2(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resRareX2Id_1230 <- simI2(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resRareX2Id_0004 <- simI2(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resRareX2Id_1004 <- simI2(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resRareX2Id_0204 <- simI2(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resRareX2Id_1204 <- simI2(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resRareX2Id_0034 <- simI2(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resRareX2Id_1034 <- simI2(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resRareX2Id_0234 <- simI2(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resRareX2Id_1234 <- simI2(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)

resRareX2Id_box <- make_box(list(resRareX2Id_0030,resRareX2Id_0004,resRareX2Id_1030,resRareX2Id_0204,resRareX2Id_1200,resRareX2Id_0034,resRareX2Id_0230,resRareX2Id_1004,resRareX2Id_1034,resRareX2Id_0234,resRareX2Id_1230,resRareX2Id_1204,resRareX2Id_1234),
                           set_names,"Spec=0.99; Hi-Freq; Rare exposure; N^2",true_val=0,ylim=c(-1,1))
ggsave("resRareX2Id_box_both.pdf",    plot=resRareX2Id_box$both)
ggsave("resRareX2Id_box_naive.pdf",   plot=resRareX2Id_box$naive)
ggsave("resRareX2Id_box_adjusted.pdf",plot=resRareX2Id_box$adjusted)


# figRareX2 <- resRareX2Id_box$naive+resRareX2Id_box$adjusted + plot_annotation(tag_levels = 'A')
ggsave("figRareX2.pdf",plot=resRareX2Id_box$vert,width=7,height=4)

