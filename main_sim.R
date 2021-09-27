##################################################
## Informative Presence: Preliminary Simulation ##
##################################################
## ****; Aug 13, 2020

library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggpattern)
library(wesanderson)
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

# colors
wes_red <- wes_palette(n=5, name="Darjeeling1")[1]
wes_green <- wes_palette(n=5, name="Darjeeling1")[2]
wes_gold <- wes_palette(n=5, name="Darjeeling1")[3]
wes_blue <- wes_palette(n=5, name="Darjeeling1")[5]

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
  
  
  ## fit models
  beta <- coef(glm(Yobs~Xobs,family=binomial))
  betaN <- coef(glm(Yobs~Xobs+Nvis,family=binomial))

  ## save
  res <- list(b0=c(beta[1],betaN[1]),b1=c(beta[2],betaN[2]))
  names(res[[1]]) <- names(res[[2]]) <- c("naive","adjusted")
  return(res)
}




## function to run simulation I
simI <- function(R=2000,      ## no. of iterations
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
    scale_x_discrete(labels=set_names_parse2[sort(1:length(res_list),decreasing=T)]) +
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


## Plot distribution of visits in scenario (1,4)
plotvisits <- function(n=2000,      ## sample size
                       Xprev=0.1,   ## prevalence of X
                       Yprev=0.3,   ## baseline (unexposed) prevalence of Y
                       b1=0,        ## logOR for X~Y; default is null
                       meanN0=2,   ## mean no. of visits when X=0 
                       meanN1=4,   ## mean no. of visits when X=1
                       negbinTheta=1
){
  ## set baseline
  b0 <- logit(Yprev)
  
  ## generate true X
  X <- rbinom(n,1,Xprev)
  
  ## generate number of visits
  # Nvis <- 1+rpois(n, -1 +    ## -1 since we add a minimum 1 to the no. visits
  #                   meanN0*(1-X)+
  #                   meanN1*X) 
  Nvis <- 1+rnegbin(n, -1 +    ## -1 since we add a minimum 1 to the no. visits
                      meanN0*(1-X)+
                      meanN1*X,
                    negbinTheta)
  
  df <- data.frame(X=X,Nvis=Nvis)
  pp <- ggplot(data=df,aes(x=as.factor(X),y=Nvis))+
    geom_boxplot(outlier.size = 0.5)+
    # scale_fill_manual(values=c("gray75","gray98"))+
    # geom_abline(intercept=true_val,slope=0,linetype="dashed")+
    theme_classic() +
    ylim(0,100)+
    labs(x="X", y="Number of Visits")+
    ggtitle("")
  
  
  
  ## save
  res <- list(plot=pp,corr=cor(X,Nvis))
  return(res)
}

####################
# Run Simulation I #
nn <- 2000
RR <- 2000 #2000
xprop <- 0.5
yprop <- 0.25
set_names <- c("(3)","(4)","(13)",
               "(24)","(12)","(34)",
               "(23)","(14)","(134)",
               "(234)","(123)","(124)","(1234)")
set_names_parse <- c(expression(paste("N\u2192",X^{"o"})),expression(paste("N\u2192",Y^{"o"})),expression(paste("X\u2192N\u2192",X^{"o"})),
                     expression(paste("Y\u2192N\u2192",Y^{"o"})),expression(paste("X,Y\u2192N")),expression(paste("N\u2192",X^{"o"},",",Y^{"o"})),
                     expression(paste("Y\u2192N\u2192",X^{"o"})),expression(paste("X\u2192N\u2192",Y^{"o"})),expression(paste("X\u2192N\u2192",X^{"o"},",",Y^{"o"})),
                     expression(paste("Y\u2192N\u2192",X^{"o"},",",Y^{"o"})),expression(paste("X,Y\u2192N\u2192",X^{"o"})),expression(paste("X,Y\u2192N\u2192",Y^{"o"})),
                     expression(paste("X,Y\u2192N\u2192",X^{"o"},",",Y^{"o"})) )
set_names_parse <- c(expression(paste(N%->%X^{o})),expression(paste(N%->%Y^{o})),expression(paste(X%->%N%->%X^{o})),
                     expression(paste(Y%->%N%->%Y^{o})),expression(paste(X,",",Y%->%N)),expression(paste(N%->%X^{o},",",Y^{o})),
                     expression(paste(Y%->%N%->%X^{o})),expression(paste(X%->%N%->%Y^{o})),expression(paste(X%->%N%->%X^{o},",",Y^{o})),
                     expression(paste(Y%->%N%->%X^{o},",",Y^{o})),expression(paste(X,",",Y%->%N%->%X^{o})),expression(paste(X,",",Y%->%N%->%Y^{o})),
                     expression(paste(X,",",Y%->%N%->%X^{o},",",Y^{o})) )
set_names_parse2 <- c(expression(paste(N%->%X^{o},"              (3)")),
                      expression(paste(N%->%Y^{o},"              (4)")),
                      expression(paste(X%->%N%->%X^{o},"           (1,3)")),
                      expression(paste(Y%->%N%->%Y^{o},"           (2,4)")),
                      expression(paste(X,",",Y%->%N,"                    (1,2)")),
                      expression(paste(N%->%X^{o},",",Y^{o},"      (3,4)")),
                      expression(paste(Y%->%N%->%X^{o},"           (2,3)")),    
                      expression(paste(X%->%N%->%Y^{o},"           (1,4)")),      
                      expression(paste(X%->%N%->%X^{o},",",Y^{o},"   (1,3,4)")),
                      expression(paste(Y%->%N%->%X^{o},",",Y^{o},"   (2,3,4)")),
                      expression(paste(X,",",Y%->%N%->%X^{o},"        (1,2,3)")), 
                      expression(paste(X,",",Y%->%N%->%Y^{o},"        (1,2,4)")),
                      expression(paste(X,",",Y%->%N%->%X^{o},",",Y^{o},"(1,2,3,4)")) )


## set seed for sim I
set.seed(1000)


### (a) NULL; CROEN SIZE; SENSITIVITY
### 
resIa_1200 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=1   ,specX=1  ,sensY=1   ,specY=1  )
resIa_0030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=1.6,meanN11=1.6,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resIa_1030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resIa_0230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=1.6,meanN11=2.3,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resIa_1230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resIa_0004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=1.6,meanN11=1.6,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resIa_1004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resIa_0204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=1.6,meanN11=2.3,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resIa_1204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resIa_0034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=1.6,meanN11=1.6,sensX=0.75,specX=1.0,sensY=0.75,specY=1.0)
resIa_1034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=0.75,specX=1.0,sensY=0.75,specY=1.0)
resIa_0234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=1.6,meanN11=2.3,sensX=0.75,specX=1.0,sensY=0.75,specY=1.0)
resIa_1234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=0.75,specX=1.0,sensY=0.75,specY=1.0)

resIa_box <- make_box(list(resIa_0030,resIa_0004,resIa_1030,resIa_0204,resIa_1200,resIa_0034,resIa_0230,resIa_1004,resIa_1034,resIa_0234,resIa_1230,resIa_1204,resIa_1234),
                      set_names,"Sens=0.75; Low Frequency",true_val=0,ylim=c(-0.5,0.5))
ggsave("resIa_box_both.pdf",    plot=resIa_box$both)
ggsave("resIa_box_naive.pdf",   plot=resIa_box$naive)
ggsave("resIa_box_adjusted.pdf",plot=resIa_box$adjusted)


### (b) NULL; CROEN SIZE; SPECIFICITY
### 
resIb_1200 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=1   ,specX=1  ,sensY=1   ,specY=1  )
resIb_0030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=1.6,meanN11=1.6,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resIb_1030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resIb_0230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=1.6,meanN11=2.3,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resIb_1230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resIb_0004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=1.6,meanN11=1.6,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resIb_1004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resIb_0204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=1.6,meanN11=2.3,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resIb_1204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resIb_0034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=1.6,meanN11=1.6,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resIb_1034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resIb_0234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=1.6,meanN11=2.3,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resIb_1234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)

resIb_box <- make_box(list(resIb_0030,resIb_0004,resIb_1030,resIb_0204,resIb_1200,resIb_0034,resIb_0230,resIb_1004,resIb_1034,resIb_0234,resIb_1230,resIb_1204,resIb_1234),
                     set_names,"Spec=0.99; Low Frequency",true_val=0,ylim=c(-0.5,0.5))
ggsave("resIb_box_both.pdf",    plot=resIb_box$both)
ggsave("resIb_box_naive.pdf",   plot=resIb_box$naive)
ggsave("resIb_box_adjusted.pdf",plot=resIb_box$adjusted)

### (c) NULL; ISRAELI SIZE; SENSITIVITY
### 
resIc_1200 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=1   ,specY=1  )
resIc_0030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resIc_1030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resIc_0230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resIc_1230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resIc_0004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resIc_1004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resIc_0204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resIc_1204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resIc_0034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=0.75,specX=1.0,sensY=0.75,specY=1.0)
resIc_1034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=0.75,specX=1.0,sensY=0.75,specY=1.0)
resIc_0234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=0.75,specX=1.0,sensY=0.75,specY=1.0)
resIc_1234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=0.75,specX=1.0,sensY=0.75,specY=1.0)

resIc_box <- make_box(list(resIc_0030,resIc_0004,resIc_1030,resIc_0204,resIc_1200,resIc_0034,resIc_0230,resIc_1004,resIc_1034,resIc_0234,resIc_1230,resIc_1204,resIc_1234),
                      set_names,"Sens=0.75; High Frequency",true_val=0,ylim=c(-0.75,0.75))
ggsave("resIc_box_both.pdf",    plot=resIc_box$both)
ggsave("resIc_box_naive.pdf",   plot=resIc_box$naive)
ggsave("resIc_box_adjusted.pdf",plot=resIc_box$adjusted)

### (d) NULL; ISRAELI SIZE; SPECIFICITY
### 
resId_1200 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=1   ,specY=1  )
resId_0030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resId_1030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resId_0230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resId_1230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resId_0004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resId_1004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resId_0204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resId_1204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resId_0034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resId_1034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resId_0234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resId_1234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)

resId_box <- make_box(list(resId_0030,resId_0004,resId_1030,resId_0204,resId_1200,resId_0034,resId_0230,resId_1004,resId_1034,resId_0234,resId_1230,resId_1204,resId_1234),
                      set_names,"Spec=0.99; High Frequency",true_val=0,ylim=c(-0.75,0.75))
ggsave("resId_box_both.pdf",    plot=resId_box$both)
ggsave("resId_box_naive.pdf",   plot=resId_box$naive)
ggsave("resId_box_adjusted.pdf",plot=resId_box$adjusted)


#########################



### (a) Non-NULL; CROEN SIZE; SENSITIVITY
### 
resNNIa_1200 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=1   ,specX=1  ,sensY=1   ,specY=1  )
resNNIa_0030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=1.6,meanN11=1.6,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resNNIa_1030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resNNIa_0230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=1.6,meanN11=2.3,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resNNIa_1230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resNNIa_0004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=1.6,meanN11=1.6,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resNNIa_1004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resNNIa_0204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=1.6,meanN11=2.3,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resNNIa_1204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resNNIa_0034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=1.6,meanN11=1.6,sensX=0.75,specX=1.0,sensY=0.75,specY=1.0)
resNNIa_1034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=0.75,specX=1.0,sensY=0.75,specY=1.0)
resNNIa_0234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=1.6,meanN11=2.3,sensX=0.75,specX=1.0,sensY=0.75,specY=1.0)
resNNIa_1234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=0.75,specX=1.0,sensY=0.75,specY=1.0)

resNNIa_box <- make_box(list(resNNIa_0030,resNNIa_0004,resNNIa_1030,resNNIa_0204,resNNIa_1200,resNNIa_0034,resNNIa_0230,resNNIa_1004,resNNIa_1034,resNNIa_0234,resNNIa_1230,resNNIa_1204,resNNIa_1234),
                      set_names,"Sens=0.75; Low Frequency",true_val=0.5,ylim=c(-0.0,1.0))
ggsave("resNNIa_box_both.pdf",    plot=resNNIa_box$both)
ggsave("resNNIa_box_naive.pdf",   plot=resNNIa_box$naive)
ggsave("resNNIa_box_adjusted.pdf",plot=resNNIa_box$adjusted)

### (b) Non-NULL; CROEN SIZE; SPECIFICITY
### 
resNNIb_1200 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=1   ,specX=1  ,sensY=1   ,specY=1  )
resNNIb_0030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=1.6,meanN11=1.6,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resNNIb_1030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resNNIb_0230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=1.6,meanN11=2.3,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resNNIb_1230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resNNIb_0004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=1.6,meanN11=1.6,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resNNIb_1004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resNNIb_0204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=1.6,meanN11=2.3,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resNNIb_1204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resNNIb_0034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=1.6,meanN11=1.6,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resNNIb_1034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resNNIb_0234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=1.6,meanN11=2.3,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resNNIb_1234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=2.3,meanN10=2.3,meanN11=3.6,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)

resNNIb_box <- make_box(list(resNNIb_0030,resNNIb_0004,resNNIb_1030,resNNIb_0204,resNNIb_1200,resNNIb_0034,resNNIb_0230,resNNIb_1004,resNNIb_1034,resNNIb_0234,resNNIb_1230,resNNIb_1204,resNNIb_1234),
                      set_names,"Spec=0.99; Low Frequency",true_val=0.5,ylim=c(-0.0,1.0))
ggsave("resNNIb_box_both.pdf",    plot=resNNIb_box$both)
ggsave("resNNIb_box_naive.pdf",   plot=resNNIb_box$naive)
ggsave("resNNIb_box_adjusted.pdf",plot=resNNIb_box$adjusted)

### (c) Non-NULL; ISRAELI SIZE; SENSITIVITY
### 
resNNIc_1200 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=1   ,specY=1  )
resNNIc_0030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resNNIc_1030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resNNIc_0230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resNNIc_1230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=0.75,specX=1.0,sensY=1   ,specY=1  )
resNNIc_0004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resNNIc_1004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resNNIc_0204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resNNIc_1204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=0.75,specY=1.0)
resNNIc_0034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=0.75,specX=1.0,sensY=0.75,specY=1.0)
resNNIc_1034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=0.75,specX=1.0,sensY=0.75,specY=1.0)
resNNIc_0234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=0.75,specX=1.0,sensY=0.75,specY=1.0)
resNNIc_1234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=0.75,specX=1.0,sensY=0.75,specY=1.0)

resNNIc_box <- make_box(list(resNNIc_0030,resNNIc_0004,resNNIc_1030,resNNIc_0204,resNNIc_1200,resNNIc_0034,resNNIc_0230,resNNIc_1004,resNNIc_1034,resNNIc_0234,resNNIc_1230,resNNIc_1204,resNNIc_1234),
                      set_names,"Sens=0.75; High Frequency",true_val=0.5,ylim=c(-0.0,1.0))
ggsave("resNNIc_box_both.pdf",    plot=resNNIc_box$both)
ggsave("resNNIc_box_naive.pdf",   plot=resNNIc_box$naive)
ggsave("resNNIc_box_adjusted.pdf",plot=resNNIc_box$adjusted)

### (d) Non-NULL; ISRAELI SIZE; SPECIFICITY
### 
resNNId_1200 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=1   ,specY=1  )
resNNId_0030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resNNId_1030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resNNId_0230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resNNId_1230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resNNId_0004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resNNId_1004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resNNId_0204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resNNId_1204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resNNId_0034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resNNId_1034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resNNId_0234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resNNId_1234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)

resNNId_box <- make_box(list(resNNId_0030,resNNId_0004,resNNId_1030,resNNId_0204,resNNId_1200,resNNId_0034,resNNId_0230,resNNId_1004,resNNId_1034,resNNId_0234,resNNId_1230,resNNId_1204,resNNId_1234),
                      set_names,"Spec=0.99; High Frequency",true_val=0.5,ylim=c(-0.0,1.0))
ggsave("resNNId_box_both.pdf",    plot=resNNId_box$both)
ggsave("resNNId_box_naive.pdf",   plot=resNNId_box$naive)
ggsave("resNNId_box_adjusted.pdf",plot=resNNId_box$adjusted)




##################
## Make figures
# fig2 <- (resIa_box$naive + resIa_box$adjusted)/ 
#   (resIb_box$naive + resIb_box$adjusted)/ 
#   (resIc_box$naive + resIc_box$adjusted)/ 
#   (resId_box$naive + resId_box$adjusted) + plot_annotation(tag_levels = 'A')
# ggsave("fig2.pdf",plot=fig2,width=9,height=12)
# 
# fig3 <- (resNNIa_box$naive + resNNIa_box$adjusted)/ 
#   (resNNIb_box$naive + resNNIb_box$adjusted)/ 
#   (resNNIc_box$naive + resNNIc_box$adjusted)/ 
#   (resNNId_box$naive + resNNId_box$adjusted) + plot_annotation(tag_levels = 'A')
# ggsave("fig3.pdf",plot=fig3,width=9,height=12)


fig2 <- (resIa_box$vert|resIc_box$vert)/
  (resIb_box$vert|resId_box$vert) + 
  plot_annotation(tag_levels = 'A')
ggsave("fig2.pdf",plot=fig2,width=12,height=9)


fig3 <- (resNNIa_box$vert|resNNIc_box$vert)/
  (resNNIb_box$vert|resNNId_box$vert) + 
  plot_annotation(tag_levels = 'A')
ggsave("fig3.pdf",plot=fig3,width=12,height=9)


## add other labels (lines)
fig2_lines <- (resIa_box$vert|resIc_box$vert)/
  (resIb_box$vert|resId_box$vert) + 
  plot_annotation(tag_levels = 'A')
ggsave("fig2_lines.pdf",plot=fig2_lines,width=12,height=9)


fig3_lines <- (resNNIa_box$vert|resNNIc_box$vert)/
  (resNNIb_box$vert|resNNId_box$vert) + 
  plot_annotation(tag_levels = 'A')
ggsave("fig3_lines.pdf",plot=fig3_lines,width=12,height=9)











#############################
# multiple-comorbidity sim ##

## function to generate data under (2,3) with multiple X's
gen14multi <- function(n=2000,      ## sample size
                       nX=5,        ## number of conditions
                       Xprev=0.1,   ## prevalence of X
                       Yprev=0.1,   ## baseline (unexposed) prevalence of Y
                       b1=0,        ## logOR for X~Y; default is null
                       meanN0=2,   ## mean no. of visits when X=0 
                       meanN1=4,   ## mean no. of visits when X=1
                       sensX=0.8,   ## sensitivity of X
                       specX=0.99,  ## specificity of X (1.0 --> no overdiagnosing due to ip)
                       sensY=0.8,   ## sensitivity of X
                       specY=0.99   ## specificity of X (1.0 --> no overdiagnosing due to ip)
){
  ## set baseline
  b0 <- logit(Yprev)
  
  ## generate true X,Y
  X <- rbinom(n,1,Xprev)
  Y <- c()
  for(i in 1:nX){
    Y <- cbind(Y,rbinom(n,1,expit(b0+b1*X)))
  }
  Y <- as.matrix(Y)
  
  ## generate number of visits
  # Nvis <- 1+rpois(n, -1 +    ## -1 since we add a minimum 1 to the no. visits
  #                   meanN0*(1-X)+
  #                   meanN1*X) 
  Nvis <- 1+rnegbin(n, -1 +    ## -1 since we add a minimum 1 to the no. visits
                      meanN0*(1-X)+
                      meanN1*X,
                    1)
  
  ## generate observed Xobs,Yobs
  Xobs <- as.numeric(rbinom(n,Nvis,sensX*X + (1-specX)*(1-X))>0)
  Yobs <- c()
  for(i in 1:nX){
    Yobs <- cbind(Yobs,rbinom(n,1,as.numeric(rbinom(n,Nvis,sensY*Y[,i] + (1-specY)*(1-Y[,i]))>0)))
  }
  Yobs <- as.matrix(Yobs)
  
  ## fit models
  pval <- pvalN <- c()
  for(i in 1:nX){
    pval <- c(pval,coef(summary(glm(Yobs[,i]~Xobs,family=binomial)))[2,4])
    pvalN <- c(pvalN,coef(summary(glm(Yobs[,i]~Xobs+Nvis,family=binomial)))[2,4])
  }
  
  
  
  ## save
  res <- list(pval=pval,pvalN=pvalN)
  return(res)
}



## function to run simulation I
sim14multi <- function(R=2000,      ## no. of iterations
                       n=2000,      ## sample size
                       nX=5,        ## number of conditions/covariates
                       Xprev=0.1,   ## prevalence of X
                       Yprev=0.1,   ## baseline (unexposed) prevalence of Y
                       b1=0,        ## logOR for X~Y; default is null
                       meanN0=1.6,   ## mean no. of visits when X=0 
                       meanN1=2.3,   ## mean no. of visits when X=1
                       sensX=0.7,   ## sensitivity of X
                       specX=1.0,  ## specificity of X (1.0 --> no overdiagnosing due to ip)
                       sensY=1,   ## sensitivity of X
                       specY=1   ## specificity of X (1.0 --> no overdiagnosing due to ip)
){
  
  ## initialize
  pvals <- pvalsN <- c()
  
  ## loop over R iterations
  for(rr in 1:R){
    
    ## generate data and fit models
    est <- gen14multi(n=n,nX=nX,Xprev=Xprev,Yprev=Yprev,b1=b1,meanN0=meanN0,meanN1=meanN1,sensX=sensX,specX=specX,sensY=sensY,specY=specY)
    
    ## collect estimates
    pvals <- rbind(pvals,est$pval)
    pvalsN <- rbind(pvalsN,est$pvalN)
  }
  
  
  ## individual type I error
  t1e <- apply(pvals<0.05,2,mean)
  t1eN <- apply(pvalsN<0.05,2,mean)
  
  ## familywise error rate
  fwer <- mean(apply(pvals,1,function(x) sum(x<0.05)>0))
  fwer_bonf <- mean(apply(pvals,1,function(x) sum(x<0.05/nX)>0))
  # adjusted for N
  fwerN <- mean(apply(pvalsN,1,function(x) sum(x<0.05)>0))
  fwer_bonfN <- mean(apply(pvalsN,1,function(x) sum(x<0.05/nX)>0))
  
  ## number of false positives
  num_fp <- apply(pvals<0.05,1,sum)
  num_fp_bonf <- apply(pvals<0.05/nX,1,sum)
  num_fpN <- apply(pvalsN<0.05,1,sum)
  num_fp_bonfN <- apply(pvalsN<0.05/nX,1,sum)
  
  ## return results
  return(list(fwer=fwer,
              fwer_bonf=fwer_bonf,
              fwerN=fwerN,
              fwer_bonfN=fwer_bonfN,
              t1e=t1e,
              t1eN=t1eN,
              num_fp=num_fp,
              num_fp_bonf=num_fp_bonf,
              num_fpN=num_fpN,
              num_fp_bonfN=num_fp_bonfN))
}



## set seed for multiple-comorb sim
set.seed(1000)

### simulate Type I error
sim14_spec <- c(0.990,0.992,0.994,0.996,0.998)#c(0.985,0.9875,0.990,0.9925,0.995)
sim14 <- list()
for(ss in 1:length(sim14_spec)){
  sim14[[ss]] <- sim14multi(R=2000,n=2000,nX=10,Xprev=0.1,Yprev=0.3,meanN0=10,meanN1=20,sensX=1,specX=1,sensY=1.0,specY=sim14_spec[ss])

}


## data for figure
df14 <- data.frame(fwer=      100*sapply(sim14,function(x) x$fwer),
                   fwer_bonf= 100*sapply(sim14,function(x) x$fwer_bonf),
                   fwerN=     100*sapply(sim14,function(x) x$fwerN),
                   fwer_bonfN=100*sapply(sim14,function(x) x$fwer_bonfN),
                   t1e=       100*sapply(sim14,function(x) mean(x$t1e)), 
                   t1eN=      100*sapply(sim14,function(x) mean(x$t1eN)), 
                   spec=      100*sim14_spec )
df14 <- df14 %>% gather(type,val,fwer:t1eN)
df14$col <- "Adjusted"
df14$col[df14$type %in%c("fwer","fwer_bonf","t1e")] <- "Unadjusted"
df14$lty <- "FWER"
df14$lty[df14$type %in%c("fwer_bonf","fwer_bonfN")] <- "FWER-Bonf"
df14$lty[df14$type %in%c("t1e","t1eN")] <- "T1E"
df14$ltybonf <- "Naive"
df14$ltybonf[df14$lty=="FWER-Bonf"] <- "Bonf"

fig14_T1E <- ggplot(data=df14[df14$lty=="T1E",],aes(x=spec,y=val,group=type,col=col))+
  geom_line(size=1)+
  scale_color_manual(values=c("gray75","black"))+ ## ,guide=FALSE
  # scale_linetype_manual(values=c(2,1,3))+ ## ,guide=FALSE
  labs(color = "", linetype="")+
  theme_classic() +
  ylim(0,100)+
  labs(x="Specificity (%)", y="%")+
  ggtitle("Type-I Error")
ggsave("fig14_T1E.pdf",plot=fig14_T1E,width=5,height=3)


fig14_FWER <- ggplot(data=df14[df14$lty=="FWER",],aes(x=spec,y=val,group=type,col=col))+
  geom_line(size=1)+
  scale_color_manual(values=c("gray75","black"))+ ## ,guide=FALSE
  # scale_linetype_manual(values=c(2,1,3))+ ## ,guide=FALSE
  labs(color = "", linetype="")+
  theme_classic() +
  ylim(0,100)+
  labs(x="Specificity (%)", y="%")+
  ggtitle("")
ggsave("fig14_FWER.pdf",plot=fig14_FWER,width=5,height=3)


fig14_FWER_bonf <- ggplot(data=df14[df14$lty=="FWER-Bonf",],aes(x=spec,y=val,group=type,col=col))+
  geom_line(size=1)+
  scale_color_manual(values=c("gray75","black"))+ ## ,guide=FALSE
  # scale_linetype_manual(values=c(2,1,3))+ ## ,guide=FALSE
  labs(color = "", linetype="")+
  theme_classic() +
  ylim(0,100)+
  labs(x="Specificity (%)", y="%")+
  ggtitle("")
ggsave("fig14_FWER_bonf.pdf",plot=fig14_FWER_bonf,width=5,height=3)


## both
fig14_FWER_both <- ggplot(data=df14[df14$lty=="FWER"|df14$lty=="FWER-Bonf",],aes(x=spec,y=val,group=type,linetype=ltybonf,col=col))+
  geom_line(size=1)+
  scale_color_manual(values=c("gray75","black"))+ ## ,guide=FALSE
  scale_linetype_manual(values=c(1,3))+ ## ,guide=FALSE
  labs(color = "", linetype="")+
  theme_classic() +
  labs(x="Specificity (%)", y="%")+
  ggtitle("Family-Wise Error Rate")
ggsave("fig14_FWER_both.pdf",plot=fig14_FWER_both,width=5,height=3)

## ALL
fig14 <- ggplot(data=df14,aes(x=spec,y=val,group=type,linetype=lty,col=col))+
  geom_line(size=1)+
  scale_color_manual(values=c("gray75","black"))+ ## ,guide=FALSE
  scale_linetype_manual(values=c(2,1,3))+ ## ,guide=FALSE
  labs(color = "", linetype="")+
  theme_classic() +
  labs(x="Specificity (%)", y="%")+
  ggtitle("")
ggsave("fig14.pdf",plot=fig14,width=5,height=3)


## data for second figure
df14_num <- data.frame(num_fp=sapply(sim14,function(x) (mean(x$num_fp))),
                       num_fp_bonf=sapply(sim14,function(x) (mean(x$num_fp_bonf))),
                       num_fpN=sapply(sim14,function(x) (mean(x$num_fpN))),
                       num_fp_bonfN=sapply(sim14,function(x) (mean(x$num_fp_bonfN))),
                       spec=100*sim14_spec )
df14_num <- df14_num %>% gather(type,val,num_fp:num_fp_bonfN)
df14_num$col <- "Adjusted"
df14_num$col[df14_num$type %in%c("num_fp","num_fp_bonf")] <- "Unadjusted"
df14_num$ltybonf <- "Naive"
df14_num$ltybonf[df14_num$type %in%c("num_fp_bonf","num_fp_bonfN")] <- "Bonf"

fig14_num_naive <- ggplot(data=df14_num[df14_num$ltybonf=="Naive",],aes(x=spec,y=val,group=type,col=col))+ #,linetype=lty
  geom_line(size=1)+
  scale_color_manual(values=c("gray75","black"))+ ## ,guide=FALSE
  # scale_linetype_manual(values=c(1,3))+ ## ,guide=FALSE
  labs(color = "")+ #, linetype=""
  theme_classic() +
  ylim(0,5)+
  labs(x="Specificity (%)", y="Avg. No. of False Positives")+
  ggtitle("")
ggsave("fig14_num_naive.pdf",plot=fig14_num_naive,width=5,height=3)

fig14_num_bonf <- ggplot(data=df14_num[df14_num$ltybonf=="Bonf",],aes(x=spec,y=val,group=type,col=col))+ #,linetype=lty
  geom_line(size=1)+
  scale_color_manual(values=c("gray75","black"))+ ## ,guide=FALSE
  # scale_linetype_manual(values=c(1,3))+ ## ,guide=FALSE
  labs(color = "")+ #, linetype=""
  theme_classic() +
  ylim(0,5)+
  labs(x="Specificity (%)", y="Avg. No. of False Positives")+
  ggtitle("")
ggsave("fig14_num_bonf.pdf",plot=fig14_num_bonf,width=5,height=3)


## both
fig14_num <- ggplot(data=df14_num,aes(x=spec,y=val,group=type,linetype=ltybonf,col=col))+ #,linetype=lty
  geom_line(size=1)+
  scale_color_manual(values=c("gray75","black"))+ ## ,guide=FALSE
  scale_linetype_manual(values=c(1,3))+ ## ,guide=FALSE
  labs(color = "",linetype="")+ #, linetype=""
  theme_classic() +
  labs(x="Specificity (%)", y="Avg. No. of False Positives")+
  ggtitle("Average Number of False Positives")
ggsave("fig14_num.pdf",plot=fig14_num,width=5,height=3)

## all 5 plots
fig14_T1E+fig14_FWER+fig14_FWER_bonf+fig14_num_naive+fig14_num_bonf


# ## all in 2
# fig_both <- fig14+fig14_num+ plot_annotation(tag_levels = 'A')
# ggsave("fig14_both.pdf",plot=fig_both,width=10,height=4)

## all in 3
fig_multiple <- fig14_T1E+fig14_FWER_both+fig14_num+ plot_annotation(tag_levels = 'A')
ggsave("fig_multiple.pdf",plot=fig_multiple,width=13,height=4)



## visits
pp_visits <- plotvisits(n=40000,Xprev=0.1,Yprev=0.3,meanN0=10,meanN1=20)
ggsave("boxplot_visits.pdf",plot=pp_visits$plot,width=4,height=4)













######################################
## Varying N~X correlation
## set seed for supplementary sims
set.seed(1000)

### (d) NULL; ISRAELI SIZE; SPECIFICITY
### 
simItheta <- function(theta){
  
  resThetaId_1200 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=1   ,specY=1  ,negbinTheta = theta)
  resThetaId_0030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1.0,specX=0.99,sensY=1   ,specY=1  ,negbinTheta = theta)
  resThetaId_1030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1.0,specX=0.99,sensY=1   ,specY=1  ,negbinTheta = theta)
  resThetaId_0230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1.0,specX=0.99,sensY=1   ,specY=1  ,negbinTheta = theta)
  resThetaId_1230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1.0,specX=0.99,sensY=1   ,specY=1  ,negbinTheta = theta)
  resThetaId_0004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99,negbinTheta = theta)
  resThetaId_1004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99,negbinTheta = theta)
  resThetaId_0204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99,negbinTheta = theta)
  resThetaId_1204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99,negbinTheta = theta)
  resThetaId_0034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99,negbinTheta = theta)
  resThetaId_1034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99,negbinTheta = theta)
  resThetaId_0234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99,negbinTheta = theta)
  resThetaId_1234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99,negbinTheta = theta)
  
  return(list(resThetaId_0030,resThetaId_0004,resThetaId_1030,resThetaId_0204,resThetaId_1200,resThetaId_0034,resThetaId_0230,resThetaId_1004,resThetaId_1034,resThetaId_0234,resThetaId_1230,resThetaId_1204,resThetaId_1234))
}

resTheta07Id<- simItheta(0.7)
resTheta10Id<- simItheta(1.0)
resTheta15Id<- simItheta(1.5)
resTheta20Id<- simItheta(2.0)
resTheta25Id<- simItheta(2.5)
resTheta30Id<- simItheta(3.0)

resTheta07Id_box <- make_box(resTheta07Id,set_names,"Spec=0.99; Hi-Freq; Theta=0.7",true_val=0,ylim=c(-0.75,0.75))
resTheta10Id_box <- make_box(resTheta10Id,set_names,"Spec=0.99; Hi-Freq; Theta=1.0",true_val=0,ylim=c(-0.75,0.75))
resTheta15Id_box <- make_box(resTheta15Id,set_names,"Spec=0.99; Hi-Freq; Theta=1.5",true_val=0,ylim=c(-0.75,0.75))
resTheta20Id_box <- make_box(resTheta20Id,set_names,"Spec=0.99; Hi-Freq; Theta=2.0",true_val=0,ylim=c(-0.75,0.75))
resTheta25Id_box <- make_box(resTheta25Id,set_names,"Spec=0.99; Hi-Freq; Theta=2.5",true_val=0,ylim=c(-0.75,0.75))
resTheta30Id_box <- make_box(resTheta30Id,set_names,"Spec=0.99; Hi-Freq; Theta=3.0",true_val=0,ylim=c(-0.75,0.75))

resTheta1 <- (resTheta07Id_box$vert)/
  (resTheta10Id_box$vert)/ 
  (resTheta15Id_box$vert)#+
  # plot_annotation(tag_levels = 'A')
ggsave("resTheta1.pdf",plot=resTheta1,width=7,height=12)
resTheta2 <- (resTheta20Id_box$vert)/
  (resTheta25Id_box$vert)/
  (resTheta30Id_box$vert)#+
# plot_annotation(tag_levels = 'A')
ggsave("resTheta2.pdf",plot=resTheta2,width=7,height=12)

## old version:
# pp <- (resTheta07Id_box$naive+resTheta07Id_box$adjusted)/
#       (resTheta10Id_box$naive+resTheta10Id_box$adjusted)/
#       (resTheta15Id_box$naive+resTheta15Id_box$adjusted)/
#       (resTheta20Id_box$naive+resTheta20Id_box$adjusted)/
#       (resTheta25Id_box$naive+resTheta25Id_box$adjusted)/
#       (resTheta30Id_box$naive+resTheta30Id_box$adjusted)
# 
# 
# ggsave("resTheta07Id_box.pdf",    plot=(resTheta07Id_box$naive+resTheta07Id_box$adjusted),width=10,height=4)
# ggsave("resTheta10Id_box.pdf",    plot=(resTheta10Id_box$naive+resTheta10Id_box$adjusted),width=10,height=4)
# ggsave("resTheta15Id_box.pdf",    plot=(resTheta15Id_box$naive+resTheta15Id_box$adjusted),width=10,height=4)
# ggsave("resTheta20Id_box.pdf",    plot=(resTheta20Id_box$naive+resTheta20Id_box$adjusted),width=10,height=4)
# ggsave("resTheta25Id_box.pdf",    plot=(resTheta25Id_box$naive+resTheta25Id_box$adjusted),width=10,height=4)
# ggsave("resTheta30Id_box.pdf",    plot=(resTheta30Id_box$naive+resTheta30Id_box$adjusted),width=10,height=4)

### distributions of visits under different theta
## corr=c(0.27 0.31 0.36 0.41 0.44 0.47)
pp_Theta07visits <- plotvisits(n=40000,Xprev=xprop,Yprev=yprop,meanN0=10,meanN1=20,negbinTheta = 0.7)
pp_Theta10visits <- plotvisits(n=40000,Xprev=xprop,Yprev=yprop,meanN0=10,meanN1=20,negbinTheta = 1.0)
pp_Theta15visits <- plotvisits(n=40000,Xprev=xprop,Yprev=yprop,meanN0=10,meanN1=20,negbinTheta = 1.5)
pp_Theta20visits <- plotvisits(n=40000,Xprev=xprop,Yprev=yprop,meanN0=10,meanN1=20,negbinTheta = 2.0)
pp_Theta25visits <- plotvisits(n=40000,Xprev=xprop,Yprev=yprop,meanN0=10,meanN1=20,negbinTheta = 2.5)
pp_Theta30visits <- plotvisits(n=40000,Xprev=xprop,Yprev=yprop,meanN0=10,meanN1=20,negbinTheta = 3.0)
ggsave("boxplot_Theta07visits.pdf",plot=pp_Theta07visits$plot,width=4,height=4)
ggsave("boxplot_Theta10visits.pdf",plot=pp_Theta10visits$plot,width=4,height=4)
ggsave("boxplot_Theta15visits.pdf",plot=pp_Theta15visits$plot,width=4,height=4)
ggsave("boxplot_Theta20visits.pdf",plot=pp_Theta20visits$plot,width=4,height=4)
ggsave("boxplot_Theta25visits.pdf",plot=pp_Theta25visits$plot,width=4,height=4)
ggsave("boxplot_Theta30visits.pdf",plot=pp_Theta30visits$plot,width=4,height=4)

plot_Thetavisits <- (pp_Theta07visits$plot+pp_Theta10visits$plot+pp_Theta15visits$plot+pp_Theta20visits$plot+pp_Theta25visits$plot+pp_Theta30visits$plot)+ plot_layout(ncol = 6)+plot_annotation(tag_levels = 'A')
ggsave("plot_Thetavisits.pdf",plot=plot_Thetavisits,width=12,height=3)

df_theta <- data.frame(theta=c(0.7, 1.0, 1.5, 2.0, 2.5, 3.0),
                       corr=c(pp_Theta07visits$corr,pp_Theta10visits$corr,pp_Theta15visits$corr,pp_Theta20visits$corr,pp_Theta25visits$corr,pp_Theta30visits$corr) )
figTheta <- ggplot(data=df_theta,aes(x=theta,y=corr))+ #,linetype=lty
  geom_line(size=1)+
  theme_classic() +
  labs(x="Theta", y="Corr(X,N)")+
  ggtitle("")
ggsave("figTheta.pdf",plot=figTheta,width=4,height=4)


##############################
## POISSON distribution for N
## function to generate data under Poisson distribution for N (corr=0.80)
genIpois <- function(n=2000,      ## sample size
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
                 specY=0.99   ## specificity of X (1.0 --> no overdiagnosing due to ip)
){
  ## set baseline
  b0 <- logit(Yprev)
  
  ## generate true X,Y
  X <- rbinom(n,1,Xprev)
  Y <- rbinom(n,1,expit(b0+b1*X))
  
  ## generate number of visits
  Nvis <- 1+rpois(n, -1 +    ## -1 since we add a minimum 1 to the no. visits
                    meanN00*(1-X)*(1-Y)+
                    meanN01*(1-X)*Y+
                    meanN10*X*(1-Y)+
                    meanN11*X*Y)

  
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





## function to run simulation under Poisson distribution for N (corr=0.80)
simIpois <- function(R=2000,      ## no. of iterations
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
                 specY=0.99  ## specificity of X (1.0 --> no overdiagnosing due to ip)
){
  
  ## initialize
  df0 <- df1 <- c()
  
  ## loop over R iterations
  for(rr in 1:R){
    
    ## generate data and fit models
    est <- genIpois(n=n,Xprev=Xprev,Yprev=Yprev,b1=b1,meanN00=meanN00,meanN01=meanN01,meanN10=meanN10,meanN11=meanN11,sensX=sensX,specX=specX,sensY=sensY,specY=specY)
    
    ## collect estimates
    df0 <- rbind(df0,est$b0)
    df1 <- rbind(df1,est$b1)
  }
  
  ## return results
  return(list(b0=df0,b1=df1))
}


resPoisId_1200 <- simIpois(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=1   ,specY=1  )
resPoisId_0030 <- simIpois(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resPoisId_1030 <- simIpois(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resPoisId_0230 <- simIpois(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resPoisId_1230 <- simIpois(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resPoisId_0004 <- simIpois(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resPoisId_1004 <- simIpois(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resPoisId_0204 <- simIpois(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resPoisId_1204 <- simIpois(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resPoisId_0034 <- simIpois(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resPoisId_1034 <- simIpois(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resPoisId_0234 <- simIpois(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resPoisId_1234 <- simIpois(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)

resPoisId_box <- make_box(list(resPoisId_0030,resPoisId_0004,resPoisId_1030,resPoisId_0204,resPoisId_1200,resPoisId_0034,resPoisId_0230,resPoisId_1004,resPoisId_1034,resPoisId_0234,resPoisId_1230,resPoisId_1204,resPoisId_1234),
                      set_names,"Spec=0.99; Hi-Freq; Poisson",true_val=0,ylim=c(-4.5,0.75))

ggsave("resPoisId_box_both.pdf",    plot=resPoisId_box$both)
ggsave("resPoisId_box_naive.pdf",   plot=resPoisId_box$naive)
ggsave("resPoisId_box_adjusted.pdf",plot=resPoisId_box$adjusted)

# figPois <- resPoisId_box$naive+resPoisId_box$adjusted+ plot_annotation(tag_levels = 'A')
ggsave("resPoisId_box.pdf",plot=resPoisId_box$vert,width=7,height=4)



## function to plot poisson visits in (1,4)
plotPoisvisits <- function(n=4000,      ## sample size
                       Xprev=0.1,   ## prevalence of X
                       Yprev=0.3,   ## baseline (unexposed) prevalence of Y
                       b1=0,        ## logOR for X~Y; default is null
                       meanN0=2,   ## mean no. of visits when X=0 
                       meanN1=4   ## mean no. of visits when X=1
){
  ## set baseline
  b0 <- logit(Yprev)
  
  ## generate true X
  X <- rbinom(n,1,Xprev)
  
  ## generate number of visits
  Nvis <- 1+rpois(n, -1 +    ## -1 since we add a minimum 1 to the no. visits
                    meanN0*(1-X)+
                    meanN1*X)

  
  df <- data.frame(X=X,Nvis=Nvis)
  pp <- ggplot(data=df,aes(x=as.factor(X),y=Nvis))+
    geom_boxplot(outlier.size = 0.5)+
    # scale_fill_manual(values=c("gray75","gray98"))+
    # geom_abline(intercept=true_val,slope=0,linetype="dashed")+
    theme_classic() +
    ylim(0,100)+
    labs(x="X", y="Number of Visits")+
    ggtitle("")
  
  
  
  ## save
  res <- list(plot=pp,corr=cor(X,Nvis))
  return(res)
}

pp_Poisvisits <- plotPoisvisits(n=40000,Xprev=0.5,Yprev=0.25,meanN0=10,meanN1=20) ## corr=0.69
ggsave("boxplot_Poisvisits.pdf",plot=pp_Poisvisits$plot,width=4,height=4)
















######################
## Rare outcome sim ##

### (d) NULL; ISRAELI SIZE; SPECIFICITY
### 
resRareId_1200 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=1   ,specY=1  )
resRareId_0030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resRareId_1030 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resRareId_0230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resRareId_1230 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resRareId_0004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resRareId_1004 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resRareId_0204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resRareId_1204 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resRareId_0034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resRareId_1034 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resRareId_0234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resRareId_1234 <- simI(R=RR,n=nn,Xprev=xprop,Yprev=0.05,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)

resRareId_box <- make_box(list(resRareId_0030,resRareId_0004,resRareId_1030,resRareId_0204,resRareId_1200,resRareId_0034,resRareId_0230,resRareId_1004,resRareId_1034,resRareId_0234,resRareId_1230,resRareId_1204,resRareId_1234),
                      set_names,"Spec=0.99; Hi-Freq; Rare Outcome",true_val=0,ylim=c(-1,1))
ggsave("resRareId_box_both.pdf",    plot=resRareId_box$both)
ggsave("resRareId_box_naive.pdf",   plot=resRareId_box$naive)
ggsave("resRareId_box_adjusted.pdf",plot=resRareId_box$adjusted)


# figRare <- resRareId_box$naive+resRareId_box$adjusted + plot_annotation(tag_levels = 'A')
ggsave("figRare.pdf",plot=resRareId_box$vert,width=7,height=4)

######################
## Rare exposure sim ##

### (d) NULL; ISRAELI SIZE; SPECIFICITY
### 
resRareXId_1200 <- simI(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=1   ,specY=1  )
resRareXId_0030 <- simI(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resRareXId_1030 <- simI(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resRareXId_0230 <- simI(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resRareXId_1230 <- simI(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1.0,specX=0.99,sensY=1   ,specY=1  )
resRareXId_0004 <- simI(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resRareXId_1004 <- simI(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resRareXId_0204 <- simI(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resRareXId_1204 <- simI(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1   ,specX=1  ,sensY=1.0,specY=0.99)
resRareXId_0034 <- simI(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=10,meanN11=10,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resRareXId_1034 <- simI(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resRareXId_0234 <- simI(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=10,meanN11=25,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)
resRareXId_1234 <- simI(R=RR,n=nn,Xprev=0.05,Yprev=yprop,b1=0,meanN00=10,meanN01=25,meanN10=25,meanN11=30,sensX=1.0,specX=0.99,sensY=1.0,specY=0.99)

resRareXId_box <- make_box(list(resRareXId_0030,resRareXId_0004,resRareXId_1030,resRareXId_0204,resRareXId_1200,resRareXId_0034,resRareXId_0230,resRareXId_1004,resRareXId_1034,resRareXId_0234,resRareXId_1230,resRareXId_1204,resRareXId_1234),
                          set_names,"Spec=0.99; Hi-Freq; Rare Exposure",true_val=0,ylim=c(-1,1))
ggsave("resRareXId_box_both.pdf",    plot=resRareXId_box$both)
ggsave("resRareXId_box_naive.pdf",   plot=resRareXId_box$naive)
ggsave("resRareXId_box_adjusted.pdf",plot=resRareXId_box$adjusted)


# figRareX <- resRareXId_box$naive+resRareXId_box$adjusted + plot_annotation(tag_levels = 'A')
ggsave("figRareX.pdf",plot=resRareXId_box$vert,width=7,height=4)
