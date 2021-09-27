##################################################
## Supplementary Simulation: Mixed sens/spec    ##
##################################################
## ****; April 2021
## Running more mixed simulations based on wider grid of sens/spec (for revision)
## Based on X-->N-->Yo (Lines 1,4; same as multiple comorb sim)


## need to first run functions from main_sim.R

library(reshape2)

## grid of sens/spec
sens_grid <- seq(0.2,0.8,length=4)
spec_grid <- seq(0.98,0.998,length=5)

## initialize list of results for sens
resIa_grid <- vector(mode = "list", length = length(sens_grid))
resIc_grid <- vector(mode = "list", length = length(sens_grid))
resNNIa_grid <- vector(mode = "list", length = length(sens_grid))
resNNIc_grid <- vector(mode = "list", length = length(sens_grid))
for(ll_sens in 1:(length(sens_grid))){ ## loop over sens
  
  ## initialize list of results for spec
  resIa_grid[[ll_sens]] <- vector(mode = "list", length = length(spec_grid))
  resIc_grid[[ll_sens]] <- vector(mode = "list", length = length(spec_grid))
  resNNIa_grid[[ll_sens]] <- vector(mode = "list", length = length(spec_grid))
  resNNIc_grid[[ll_sens]] <- vector(mode = "list", length = length(spec_grid))
  
  for(ll_spec in 1:(length(spec_grid))){ ## loop over spec
    print(c(ll_sens,ll_spec)) 
    
    ### (a) NULL; CROEN SIZE; 
    resIa_grid[[ll_sens]][[ll_spec]] <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=1   ,specX=1  ,sensY=sens_grid[ll_sens],specY=spec_grid[ll_spec])
    ### (c) NULL; ISRAELI SIZE; 
    resIc_grid[[ll_sens]][[ll_spec]] <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1   ,specX=1  ,sensY=sens_grid[ll_sens],specY=spec_grid[ll_spec])

    ### (a) Non-NULL; CROEN SIZE; 
    resNNIa_grid[[ll_sens]][[ll_spec]] <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=1.6,meanN01=1.6,meanN10=2.3,meanN11=2.3,sensX=1   ,specX=1  ,sensY=sens_grid[ll_sens],specY=spec_grid[ll_spec])
    ### (c) Non-NULL; ISRAELI SIZE; 
    resNNIc_grid[[ll_sens]][[ll_spec]] <- simI(R=RR,n=nn,Xprev=xprop,Yprev=yprop,b1=0.5,meanN00=10,meanN01=10,meanN10=25,meanN11=25,sensX=1   ,specX=1  ,sensY=sens_grid[ll_sens],specY=spec_grid[ll_spec])
    
  }
}

taba_grid <- matrix(NA,nrow=length(sens_grid),ncol=length(spec_grid))
tabc_grid <- matrix(NA,nrow=length(sens_grid),ncol=length(spec_grid))
tabNNa_grid <- matrix(NA,nrow=length(sens_grid),ncol=length(spec_grid))
tabNNc_grid <- matrix(NA,nrow=length(sens_grid),ncol=length(spec_grid))
for(ll_sens in 1:(length(sens_grid))){ ## loop over sens
  for(ll_spec in 1:(length(spec_grid))){ ## loop over spec

    taba_grid[ll_sens,ll_spec] <- median(resIa_grid[[ll_sens]][[ll_spec]]$b1[,1])
    tabc_grid[ll_sens,ll_spec] <- median(resIc_grid[[ll_sens]][[ll_spec]]$b1[,1])
    tabNNa_grid[ll_sens,ll_spec] <- median(resNNIa_grid[[ll_sens]][[ll_spec]]$b1[,1])-0.5
    tabNNc_grid[ll_sens,ll_spec] <- median(resNNIc_grid[[ll_sens]][[ll_spec]]$b1[,1])-0.5
    
  }
}
rownames(taba_grid) <- rownames(tabc_grid) <- rownames(tabNNa_grid) <- rownames(tabNNc_grid) <- sens_grid
colnames(taba_grid) <- colnames(tabc_grid) <- colnames(tabNNa_grid) <- colnames(tabNNc_grid) <- spec_grid

heatmapa <- ggplot(data = melt(taba_grid), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  geom_text(aes(Var1, Var2, label = round(value,2)), color = "white") +
  guides(fill=FALSE)+
  scale_x_continuous(breaks=sens_grid)+
  scale_y_continuous(breaks=spec_grid)+
  labs(x = "Sensitivity", y = "Specificity")+
  theme(panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank())
ggsave("heatmapa.pdf",    plot=heatmapa)

heatmapc <- ggplot(data = melt(tabc_grid), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  geom_text(aes(Var1, Var2, label = round(value,2)), color = "white") +
  guides(fill=FALSE)+
  scale_x_continuous(breaks=sens_grid)+
  scale_y_continuous(breaks=spec_grid)+
  labs(x = "Sensitivity", y = "Specificity")+
  theme(panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank())
ggsave("heatmapc.pdf",    plot=heatmapc)

heatmapNNa <- ggplot(data = melt(tabNNa_grid), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  geom_text(aes(Var1, Var2, label = round(value,2)), color = "white") +
  guides(fill=FALSE)+
  scale_x_continuous(breaks=sens_grid)+
  scale_y_continuous(breaks=spec_grid)+
  labs(x = "Sensitivity", y = "Specificity")+
  theme(panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank())
ggsave("heatmapNNa.pdf",    plot=heatmapNNa)

heatmapNNc <- ggplot(data = melt(tabNNc_grid), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  geom_text(aes(Var1, Var2, label = round(value,2)), color = "white") +
  guides(fill=FALSE)+
  scale_x_continuous(breaks=sens_grid)+
  scale_y_continuous(breaks=spec_grid)+
  labs(x = "Sensitivity", y = "Specificity")+
  theme(panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank())
ggsave("heatmapNNc.pdf",    plot=heatmapNNc)
