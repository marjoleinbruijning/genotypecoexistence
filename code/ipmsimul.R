
##########################################################################################
## Simulate population trajectories of competing genotypes
##########################################################################################

times <- 100 # number of time steps
NB <- NN <- matrix(NA,nrow=length(allTemp),ncol=times)

# start with observed neonate body size distribution
sizeoffcounts <- table(cut(dat$sizeoff,breaks=allsizes)) 
sizeoffcounts <- 6* (sizeoffcounts / sum(sizeoffcounts))

for (j in 1:length(allTemp)) {

    popdens <- array(NA,dim=c(nc,times,2))
    popdens[,1,] <- 0
    popdens[-1,1,1] <- popdens[-1,1,2] <- sizeoffcounts

    for (i in 1:(times-1)) {

        ## Create Southern IPM
        treatments <- expand.grid(nHetero=sum(popdens[,i,2]),
                                  nCon=sum(popdens[,i,1]),
                                  country='B',
                                  temp=allTemp[j],clone=NA)

        vr <- buildVitalrates(allsizes,
                              treatments=treatments,
                              ndraws=ndraws)

        ipm <- createIPM(allsizes,
                         P=colMeans(vr$P),
                         growth=colMeans(vr$growth),
                         R1=colMeans(vr$R1),
                         R2=colMeans(vr$R2),
                         R3=colMeans(vr$R3),
                         R4=colMeans(vr$R4),
                         offsize=colMeans(vr$offsize),
                         nc=nc)

        popdens[,i+1,1] <- ipm$ipm %*% popdens[,i,1]


        ## Create Northern IPM
        treatments <- expand.grid(nHetero=sum(popdens[,i,1]),
                                  nCon=sum(popdens[,i,2]),
                                  country='N',
                                  temp=allTemp[j],clone=NA)

        vr <- buildVitalrates(allsizes,
                              treatments=treatments,
                              ndraws=ndraws)

        ipm <- createIPM(allsizes,
                         P=colMeans(vr$P),
                         growth=colMeans(vr$growth),
                         R1=colMeans(vr$R1),
                         R2=colMeans(vr$R2),
                         R3=colMeans(vr$R3),
                         R4=colMeans(vr$R4),
                         offsize=colMeans(vr$offsize),
                         nc=nc)

        popdens[,i+1,2] <- ipm$ipm %*% popdens[,i,2]

        cat('\r',i,'-',times,'\n')
    }

    NB[j,] <- apply(popdens[,,1],2,sum)
    NN[j,] <- apply(popdens[,,2],2,sum)
 
}

save(NB,NN,file='Results/ipmsimultrajec.RData')
