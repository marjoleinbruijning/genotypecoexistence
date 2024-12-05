

##########################################################################################
## Vital rate decomposition
##########################################################################################

changeVR <- list(NA,'surv','growth','eggs','repr','clutch','sexoff','offsize')

ind <- rep(list(rep('incF',7)),8)
for (z in 2:8) {ind[[z]][z-1] <- 'incF2'}

dN <- array(NA,dim=c(n,n,2,length(allTemp),length(changeVR),ndraws))
treatments <- expand.grid(nHetero=allN,nCon=allN,country=c('N','B'),
                          temp=allTemp,clone=NA)

vr <- buildVitalrates(treatments=treatments,
                      allsizes,
                      ndraws=ndraws)


## Build IPM for each combination
for (i in 1:nrow(treatments)) {

    incF <- vr$df$nHetero == treatments$nH[i] & vr$df$nCon == treatments$nC[i] &
        vr$df$country == treatments$country[i] &
        vr$df$temp == treatments$temp[i] &
        vr$df$sex == 'f'

    ## Which rows to select for other country
    incF2 <- vr$df$nHetero == treatments$nH[i] & vr$df$nCon == treatments$nC[i] &
        vr$df$country != treatments$country[i] & ## select other country
        vr$df$temp == treatments$temp[i] &
        vr$df$sex == 'f'

    for (k in 1:length(changeVR)) {

        for (j in 1:ndraws) {
            
            ipm <- createIPM(allsizes,
                             P=vr$P[j,get(ind[[k]][1])],
                             growth=vr$growth[j,get(ind[[k]][2])],
                             R1=vr$R1[j,get(ind[[k]][3])],
                             R2=vr$R2[j,get(ind[[k]][4])],
                             R3=vr$R3[j,get(ind[[k]][5])],
                             R4=vr$R4[j,get(ind[[k]][6])],
                             offsize=vr$offsize[j,get(ind[[k]][7])],
                             nc=nc)
            
            ## Save in dN array to ease plotting
            if (treatments$country[i] == 'B') {
                dN[which(allN == treatments$nC[i]),
                   which(allN == treatments$nH[i]),1,
                   which(allTemp == treatments$temp[i]),k,j] <- getLambda(ipm$ipm)
                ## If Belgian
                
            } else {
                dN[which(allN == treatments$nH[i]),
                   which(allN == treatments$nC[i]),2,
                   which(allTemp == treatments$temp[i]),k,j] <- getLambda(ipm$ipm)
                ## If Norwegian
            }

            cat('\r',i,'/',nrow(treatments),'-',k,'/',length(changeVR),'-',j,'/',ndraws)

        }  

    }
}

saveRDS(dN,file='Results/decomposition.rds')
