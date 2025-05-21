

####################################################################
## Figure 4
####################################################################


## Figure decomposition
dN <- readRDS('Results/decomposition.rds')

propb <- array(NA,dim=c(length(allTemp),length(changeVR),ndraws))
pars <- array(NA,dim=c(length(allTemp),4,ndraws,length(changeVR)),
              dimnames=list(allTemp,c('a12','a11','a22','a21'),
                            1:ndraws,changeVR))

for (k in 1:length(allTemp)) {

    for (z in 1:length(changeVR)) {
        
        for (j in 1:ndraws) {

            mod <- list()

            for (i in 1:2) { # First B, then N
                df <- contourLines(allN,allN,log(dN[,,i,k,z,j]),
                                   levels=0)
                ## get values where log lambda=0
                if (length(df) > 0) {
                    df <- data.frame(n2=df[[1]]$y,n1=df[[1]]$x)
                    mod[[i]] <- lm(n2~n1,data=df)
                }
            }

            if (all(unlist(lapply(mod,length)) > 0)) {
                ## Predict eq proportions (where do nullclines intersect?)
                beq <- (coef(mod[[1]])[1] - coef(mod[[2]])[1]) /
                    ( coef(mod[[2]])[2] - coef(mod[[1]])[2] )
                
                neq <- predict(mod[[1]],newdata=data.frame(n1=beq))

                if (beq < 0) beq <- 0
                if (neq < 0) neq <- 0
                propb[k,z,j] <- beq / (beq + neq)

                ## Get niche overlap and fitness difference
                pars[k,'a12',j,z] <- 1/coef(mod[[1]])[1]
                pars[k,'a11',j,z] <- 1/ (-(coef(mod[[1]])[1]) /
                                       (coef(mod[[1]])[2]))

                pars[k,'a22',j,z] <- 1/coef(mod[[2]])[1]
                pars[k,'a21',j,z] <- 1/ (-(coef(mod[[2]])[1]) /
                                       (coef(mod[[2]])[2]))

              }
        }
    }
}

nicheoverlap <- apply(pars,c(1,3,4),
                      function(x) sqrt((x['a21']*x['a12'])/(x['a11']*x['a22'])))
fitnessdiff <- apply(pars,c(1,3,4),
                     function(x) sqrt((x['a21']*x['a22'])/(x['a11']*x['a12'])))


pdf('Results/Figure4.pdf',width=12,height=8)

cols <-  wesanderson::wes_palette("GrandBudapest1", 7, type = "continuous")
namess <- c('Original','Survival','Growth','Carrying eggs',
            'Egg development','Clutch size',
            'Neonate sex ratio','Neonate size')


ly <- matrix(0,nrow=8,ncol=20)
ly[1:4,1:4] <- 2
ly[1:4,5:8] <- 3
ly[1:4,9:12] <- 4
ly[1:4,13:16] <- 5
ly[1:4,17:20] <- 1

ly[5:8,1:4] <- 6
ly[5:8,5:8] <- 7
ly[5:8,9:12] <- 8
ly[5:8,13:16] <- 9

ly[5:6,17:18] <- 10
ly[5:6,19:20] <- 11
ly[7:8,17:18] <- 12
ly[7:8,19:20] <- 13

layout(ly)

par(mar=c(3,3,3,1),oma=c(4,2,2,2))


## Equilibrium proportions
plot(100,100,xlim=c(1,3),ylim=c(.5,7.5),bty='n',xaxt='n',yaxt='n',xlab='',ylab='',
     xaxs='i',yaxs='i')
points(rep(1.1,7),1:7,pch=16,col=cols,cex=2)
text(rep(1.2,7),1:7,labels=namess[-1],pos=4,cex=1.5)
mtext('Exchanged vital rate',4,outer=FALSE,line=-1)


for (i in 1:length(allTemp)) {

    df <- data.frame(vr=rep(unlist(changeVR),ndraws),prop=c(propb[i,,]))
    df$vr[is.na(df$vr)] <- 'full'
    df$vr <- factor(df$vr,levels=c('full','surv','growth','eggs',
                                   'repr','clutch','sexoff','offsize'))

    plot(100,100,xlim=c(0,1),ylim=c(.5,7.5),xlab='',yaxt='n',ylab='',
         bty='n',xaxs='i',yaxs='i',xaxt='n')
    abline(v=seq(0,1,.25),lty=2,col='grey')
    abline(v=.5,col='grey')
    axis(1,at=seq(0,1,.25))

    if (i %in% 1:3) {
        mtext(paste0(round(allTemp[i]*scTemp[2]+scTemp[1]),'°C'),3)
    } else if (i == 4) {
        mtext(paste0('Extrapolate to: ',round(allTemp[i]*scTemp[2]+scTemp[1]),'°C'),3)
    }
    
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


## Coexistence planes
allniches <- seq(0,1,length.out=1001)
allk <-  seq(0,3.5,length.out=500)
mat <- outer(allniches,allk,
             function(x,y) x < y & y < 1/x)

for(k in 1:4) {
    plot(100,100,xlim=c(0,1),ylim=c(.4,2.1),xaxt='n',
         yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
    contour(allniches,allk,mat,add=TRUE,labels='',lwd=2)
    image(allniches,allk,mat,col=c(NA,'lightgrey'),add=TRUE)

    if (k == 1) mtext('Competitive ability difference',2,cex=0.8,line=2)
    
    axis(1,at=seq(0,1,.25))
    if (k == 1) {
        axis(2)
    } else {
        axis(2,labels=NA)
    }
    
    
    for (z in 2:8) {

        quant <- quantile(nicheoverlap[k,,z],prob=c(.25,.75),na.rm=TRUE)
        segments(x0=quant[1],x1=quant[2],
                 y0=median(fitnessdiff[k,,z],na.rm=TRUE),col=cols[z-1])

        quant <- quantile(fitnessdiff[k,,z],prob=c(.25,.75),na.rm=TRUE)
        segments(y0=quant[1],y1=quant[2],
                 x0=median(nicheoverlap[k,,z],na.rm=TRUE),col=cols[z-1])
        
        points(median(nicheoverlap[k,,z],na.rm=TRUE),
               median(fitnessdiff[k,,z],na.rm=TRUE),col=cols[z-1],pch=16,cex=2.5)        
        
    }

     points(median(nicheoverlap[k,,1],na.rm=TRUE),
        median(fitnessdiff[k,,1],na.rm=TRUE),cex=2.5,pch=15)

}



## Violin plot
temps <- (c(18,26) - scTemp[1]) / scTemp[2] ##unique(dat$temp)[c(2,4)]

highdens <- quantile(dat$n,prob=.9)

df <- expand.grid(size=1,size2=1,country=c('B','N'),nCon=c(1,highdens),nHetero=c(1,highdens),temp=temps)
df <- df[!(df$nCon == highdens & df$nHetero == highdens),]
pred <- plogis(posterior_linpred(modSexoff,newdata=df,re_formula=NA))

namess <- c('Low density','High sympatric density','High allopatric density')

for (j in 1:2) { ## over temperature
    for (i in 1:2) { ## over country

        inc <- df$temp == temps[j] & df$country == names(colsCountry)[i]

        vioplot::vioplot(pred[,inc],col=colsCountry[i],ylim=c(0,1),yaxt='n',bty='l',
                         xaxt='n')
        
        if (i == 1) {
            axis(2)
        } else {        
            axis(2,labels=NA)
        }

        if (j == 2) {
            text(1:3, -.1,
                 srt = 45, adj = 1, xpd = NA,
                 labels = namess, cex = 1)
        }
        
    }
}


text(-56,7.3,labels='A) Vital rate decomposition: Equilibrium proportions',
     xpd=NA,cex=1.5,pos=4)

text(-56,3.1,
     labels='B) Vital rate decomposition: Niche overlap and competitive ability difference',
     xpd=NA,cex=1.5,pos=4)

text(-8,3.3,labels='C) Density-dependent neonate\nfemale probabilities',
     xpd=NA,cex=1.5,
     pos=4)


text(-25,3.4,labels='Equilibrium proportion Southern clones',xpd=NA,cex=1.2)

text(-25,-.5,labels='Niche overlap',xpd=NA,cex=1.2)

text(-7.5,1.5,labels='Neonate female probability',xpd=NA,srt=90,cex=1.2)


text(4.5,2.2,labels=paste0(' ',scTemp[2]*temps[1]+scTemp[1],' °C'),cex=1.2,xpd=NA)
text(4.5,.5,labels=paste0(' ',scTemp[2]*temps[2]+scTemp[1],' °C'),cex=1.2,xpd=NA)


text(-3.9,1.45,col=colsCountry[1],labels='Southern\ngenotypes',xpd=NA,cex=1.2)
text(2.3,1.45,col=colsCountry[2],labels='Northern\ngenotypes',xpd=NA,cex=1.2)

dev.off()
