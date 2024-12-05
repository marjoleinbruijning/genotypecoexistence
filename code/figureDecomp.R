

####################################################################
## Figure 4
####################################################################


dN <- readRDS('Results/decomposition.rds')
changeVR <- list(NA,'surv','growth','eggs','repr','clutch','sexoff','offsize')

propb <- array(NA,dim=c(length(allTemp),length(changeVR),ndraws))
pars <- array(NA,dim=c(length(allTemp),4,ndraws),
              dimnames=list(allTemp,c('a12','a11','a22','a21'),1:ndraws))

for (k in 1:length(allTemp)) {

    for (z in 1:length(changeVR)) {
        
        for (j in 1:ndraws) {

            mod <- list()

            for (i in 1:2) { # First B, then N
                df <- contourLines(allN,allN,log(dN[,,i,k,z,j]),levels=0) ## get values where log lambda=0
                if (length(df) > 0) {
                    df <- data.frame(n2=df[[1]]$y,n1=df[[1]]$x)
                    mod[[i]] <- lm(n2~n1,data=df)
                }
            }

            if (all(unlist(lapply(mod,length)) > 0)) {
                ## Predict eq proportions (where do nullclines intersect?)
                beq <- (coef(mod[[1]])[1] - coef(mod[[2]])[1]) / ( coef(mod[[2]])[2] - coef(mod[[1]])[2] )
                neq <- predict(mod[[1]],newdata=data.frame(n1=beq))

                if (beq < 0) beq <- 0
                if (neq < 0) neq <- 0
                propb[k,z,j] <- beq / (beq + neq)

              }
        }
    }
}



pdf('Results/figure4.pdf',width=12,height=8)

ly <- matrix(0,nrow=7,ncol=16)
ly[1:3,1:3] <- 1
ly[1:3,4:6] <- 2
ly[1:3,7:9] <- 3
ly[1:3,10:12] <- 4
ly[1:3,13:15] <- 5

ly[5:7,2:4] <- 6
ly[5:7,5:7] <- 7
ly[5:7,8:10] <- 8
ly[5:7,11:13] <- 9

layout(ly)

cols <-  wesanderson::wes_palette("GrandBudapest1", 7, type = "continuous")
namess <- c('Original','Survival','Growth','Carrying eggs','Egg development','Clutch size',
            'Female probability','Neonate size')


par(mar=c(1,2,1,1),oma=c(4,1,2,1))

plot(100,100,xlim=c(1,3),ylim=c(.5,7.5),bty='n',xaxt='n',yaxt='n',xlab='',ylab='',
     xaxs='i',yaxs='i')
text(rep(3,7),1:7,labels=namess[-1],pos=2,cex=1.5)
mtext('Exchanged vital rate',2,outer=FALSE,line=-3)

for (i in 1:length(allTemp)) {

    df <- data.frame(vr=rep(unlist(changeVR),ndraws),prop=c(propb[i,,]))
    df$vr[is.na(df$vr)] <- 'full'
    df$vr <- factor(df$vr,levels=c('full','surv','growth','eggs','repr','clutch','sexoff','offsize'))

    plot(100,100,xlim=c(0,1),ylim=c(.5,7.5),xlab='',yaxt='n',ylab='',
         bty='n',xaxs='i',yaxs='i',xaxt='n')
    abline(v=seq(0,1,.25),lty=2,col='grey')
    abline(v=.5,col='grey')
    axis(1,at=seq(0,1,.25))
    mtext(paste0(round(allTemp[i]*scTemp[2]+scTemp[1]),'°C'),3)

    abline(v=median(df$prop[df$vr == 'full'],na.rm=TRUE),lwd=3)
    abline(v=c(0,1))
    abline(h=7.5)

    for (j in 1:7) {
        inc <- df$vr == levels(df$vr)[j+1]
        quant <- quantile(df$prop[inc],prob=c(.05,.25,.5,.75,.95),na.rm=TRUE)
        points(quant[3],j,col=cols[j],cex=2.5,pch=16)
        segments(x0=quant[2],x1=quant[4],y0=j,col=cols[j],lwd=4)
        segments(x0=quant[1],x1=quant[5],y0=j,col=cols[j],lty=1,lwd=1)
    }
}



## Violin plot
temps <- (c(18,26) - scTemp[1]) / scTemp[2]

highdens <- quantile(dat$n,prob=.9)

df <- expand.grid(size=1,size2=1,country=c('B','N'),nCon=c(1,highdens),nHetero=c(1,highdens),temp=temps)
df <- df[!(df$nCon == highdens & df$nHetero == highdens),]
pred <- plogis(posterior_linpred(modSexoff,newdata=df,re_formula=NA))

namess <- c('Low density','High sympatric\ndensity','High allopatric\ndensity')

for (j in 1:2) { ## over temperature
    for (i in 1:2) { ## over country

        inc <- df$temp == temps[j] & df$country == names(colsCountry)[i]

        vioplot::vioplot(pred[,inc],col=colsCountry[i],ylim=c(0,1),yaxt='n',bty='l',
                         xaxt='n')
        
        if (i == 1 & j == 1) {
            axis(2)
            mtext('Neonate female probability',2,line=2)
        } else {        
            axis(2,labels=NA)
        }

        if (i == 1 & j == 1) text(1,0.05,col=colsCountry[1],labels='Southern\n genotypes',cex=1.2)
        if (i == 2 & j == 1) text(1,0.05,col=colsCountry[2],labels='Northern\n genotypes',cex=1.2)

   
        text(1:3, -.1,#par("usr")[3]-0.25, 
             srt = 45, adj = 1, xpd = NA,
             labels = namess, cex = 1)
    }
}



text(-10,2.7,labels='A) Decomposition analysis',xpd=NA,cex=2)
text(-8.8,1.2,labels='B) Density-dependent neonate sex ratios',xpd=NA,cex=2)
text(-1,1.35,labels='Equilibrium proportion Southern clones',xpd=NA,cex=1.5)

text(-8,1.1,labels=paste0(' ',scTemp[2]*temps[1]+scTemp[1],' °C'),cex=2,xpd=NA)
text(-.3,1.1,labels=paste0(' ',scTemp[2]*temps[2]+scTemp[1],' °C'),cex=2,xpd=NA)


dev.off()
