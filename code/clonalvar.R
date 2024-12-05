
#########################################################################################
## Zero-growth isoclines for each clonal combination
#########################################################################################

## all  B-N combinations and temperatures
treatments <- expand.grid(nHetero=allN,nCon=allN,country=c('N','B'))
treatments <- do.call("rbind", replicate(6, treatments, simplify = FALSE)) ## for each clone
treatments$clone <- rep(c(1,7,2,8,3,9,4,10,5,11,6,12),each=36)

treatments <- do.call("rbind", replicate(length(allTemp),
                                         treatments, simplify = FALSE)) ## for each temperature
treatments$temp <- rep(allTemp,each=432)

vr <- buildVitalrates(allsizes,
                      treatments=treatments,
                      ndraws=ndraws)

## Mean predictions
dNmean <- array(NA,dim=c(n,n,2,length(allTemp),12)) ## 12 clones
l <- rep(NA,nrow(treatments))

## average over posterior draws
P <- colMeans(vr$P)
growth <- colMeans(vr$growth)
R1 <- colMeans(vr$R1)
R2 <- colMeans(vr$R2)
R3 <- colMeans(vr$R3)
R4 <- colMeans(vr$R4)
offsize <- colMeans(vr$offsize)


## Build IPM and get lambda for each treatment
for (i in 1:nrow(treatments)) {

    incF <- vr$df$nHetero == treatments$nH[i] &
        vr$df$nCon == treatments$nC[i] &
        vr$df$country == treatments$country[i] &
        vr$df$temp == treatments$temp[i] &
        vr$df$clone == treatments$clone[i]
    
    
    ipm <- createIPM(allsizes,
                     P=P[incF],
                     growth=growth[incF],
                     R1=R1[incF],
                     R2=R2[incF],
                     R3=R3[incF],
                     R4=R4[incF],
                     offsize=offsize[incF],
                     nc=nc)
    
    l[i] <- getLambda(ipm$ipm)

    ## Save in dN array to ease plotting
    if (treatments$country[i] == 'B') {
        dNmean[which(allN == treatments$nC[i]),
               which(allN == treatments$nH[i]),1,
               which(allTemp == treatments$temp[i]),
               treatments$clone[i]] <- l[i] ## Belgian
        
    } else {
        dNmean[which(allN == treatments$nH[i]),
               which(allN == treatments$nC[i]),2,
               which(allTemp == treatments$temp[i]),
               treatments$clone[i]] <- l[i] ## Norwegian
    }

    cat('\r',i,'/',nrow(treatments))

}  


saveRDS(dNmean,file=paste0('Results/clonalvariation.rds'))
