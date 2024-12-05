

##########################################################################################
## Obtain niche overlap and fitness difference using IPM
##########################################################################################

## Get vital rate predictions for each treatment and posterior draw
treatments <- expand.grid(nHetero=allN,nCon=allN,country=c('N','B'),temp=allTemp,clone=NA)
l <- array(NA,dim=c(nrow(treatments),ndraws))
dN <- array(NA,dim=c(n,n,2,length(allTemp),ndraws))

vr <- buildVitalrates(allsizes,
                      treatments=treatments,
                      ndraws=ndraws)

## Build IPM for each posterior draw
for (i in 1:nrow(treatments)) {

    incF <- vr$df$nHetero == treatments$nH[i] &
        vr$df$nCon == treatments$nC[i] &
        vr$df$country == treatments$country[i] &
        vr$df$temp == treatments$temp[i]
    
    for (j in 1:ndraws) {
        
        ipm <- createIPM(allsizes,
                         P=vr$P[j,incF],
                         growth=vr$growth[j,incF],
                         R1=vr$R1[j,incF],
                         R2=vr$R2[j,incF],
                         R3=vr$R3[j,incF],
                         R4=vr$R4[j,incF],
                         offsize=vr$offsize[j,incF],
                         nc=nc)
        
        l[i,j] <- getLambda(ipm$ipm)

        ## Save in dN array to ease plotting
        if (treatments$country[i] == 'B') {
            dN[which(allN == treatments$nC[i]),
               which(allN == treatments$nH[i]),1,
               which(allTemp == treatments$temp[i]),j] <- l[i,j] ## Belgian
            
        } else {
            dN[which(allN == treatments$nH[i]),
               which(allN == treatments$nC[i]),2,
               which(allTemp == treatments$temp[i]),j] <- l[i,j] ## Norwegian
        }

        cat('\r',i,'/',nrow(treatments),'-',j,'/',ndraws)

    }  
}

saveRDS(dN,file='Results/coexistence.rds')


## Mean predictions
dNmean <- array(NA,dim=c(n,n,2,length(allTemp)))
l <- rep(NA,nrow(treatments))

for (i in 1:nrow(treatments)) {

    incF <- vr$df$nHetero == treatments$nH[i] &
        vr$df$nCon == treatments$nC[i] &
        vr$df$country == treatments$country[i] &
        vr$df$temp == treatments$temp[i]
    
    
    ipm <- createIPM(allsizes,
                     P=colMeans(vr$P[,incF]),
                     growth=colMeans(vr$growth[,incF]),
                     R1=colMeans(vr$R1[,incF]),
                     R2=colMeans(vr$R2[,incF]),
                     R3=colMeans(vr$R3[,incF]),
                     R4=colMeans(vr$R4[,incF]),
                     offsize=colMeans(vr$offsize[,incF]),
                     nc=nc)
    
    l[i] <- getLambda(ipm$ipm)

    ## Save in dN array to ease plotting
    if (treatments$country[i] == 'B') {
        dNmean[which(allN == treatments$nC[i]),
               which(allN == treatments$nH[i]),1,
               which(allTemp == treatments$temp[i])] <- l[i] ## Belgian
        
    } else {
        dNmean[which(allN == treatments$nH[i]),
               which(allN == treatments$nC[i]),2,
               which(allTemp == treatments$temp[i])] <- l[i] ## Norwegian
    }

    cat('\r',i,'/',nrow(treatments))

}  

saveRDS(dNmean,file='Results/coexistenceMean.rds')
