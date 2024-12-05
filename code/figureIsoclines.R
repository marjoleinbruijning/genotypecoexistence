
####################################################################
## Figure 2
####################################################################

pdf('Results/figure2.pdf',width=12,height=6)

## Figure with null isoclines and equilibrium proportions
ly <- matrix(0,nrow=6,ncol=12)
ly[1:3,1:3] <- 1
ly[1:3,4:6] <- 2
ly[4:6,1:3] <- 3
ly[4:6,4:6] <- 4

ly[1:6,7:12] <- 5

layout(ly)
par(mar=c(4,4,2,2))


## Load results
dNmean <- readRDS('Results/coexistenceMean.rds')

## Plot output
propb <- rep(NA,length(allTemp))
for (k in 1:length(allTemp)) {

    plot(-100,-100,xlim=c(0,1000),ylim=c(0,1000),xlab='',
         ylab='')
    fig_label(paste0(toupper(letters)[k],')'),pos='topright',cex=1.5)
    for (j in 1:n) { # for each Belgian pop size
        arrows(y0=allN, # all Norwegian sizes
               x0=rep(allN[j],n), # focal belgian
               y1=(allN*dNmean[j,,2,k]), # how much does Norwegian pop size change?
               x1=(allN[j]*dNmean[j,,1,k]), # how much does Belgian pop size change?
               length=.05,col='grey')
    }
    mtext(paste0(round(allTemp[k]*scTemp[2]+scTemp[1]),'°C'),3)

    ## Fit linear lines to obtain alpha coefficients
    mod <- list()
    for (i in 1:2) { ## First B, then N
         df <- contourLines(allN,allN,log(dNmean[,,i,k]),levels=0)
        ## get values where log lambda=0
        df <- data.frame(n2=df[[1]]$y,n1=df[[1]]$x)

        mod[[i]] <- lm(n2~n1,data=df)
        pred <- predict(mod[[i]],newdata=data.frame(n1=seq(-100,2000,1)),se.fit=TRUE)
        lines(-100:2000,pred$fit,lwd=4,lty=1,col=colsCountry[i])
        polygon(x=c(-100:2000,rev(-100:2000),-100),
                y=c(pred$fit+1.96*pred$se.fit,rev(pred$fit-1.96*pred$se.fit),
                    pred$fit[1]+1.96*pred$se.fit[1]),
                border=NA,col=paste0(colsCountry[i],'30'))
    }

    ## Predict eq proportions (where do nullclines intersect?)
    beq <- (coef(mod[[1]])[1] - coef(mod[[2]])[1]) /
        ( coef(mod[[2]])[2] - coef(mod[[1]])[2] )
    neq <- predict(mod[[1]],newdata=data.frame(n1=beq))

    if (!is.na(beq) & beq < 0) beq <- 0
    if (!is.na(neq) & neq < 0) neq <- 0
    propb[k] <- beq / (beq + neq)
    
    points(beq,neq,cex=1.5)

}
legend('topright',col=colsCountry,lwd=2,legend=c('South','North'),bg='white',cex=1.2)

text(grconvertX(1/4, "ndc", "user"), grconvertY(.02, "ndc", "user"),
     "Southern population size", cex=1.5, col="black", xpd=NA)

text(grconvertX(.01, "ndc", "user"), grconvertY(3/6, "ndc", "user"),
     "Northern population size", cex=1.5, col="black", xpd=NA,srt=90)



### Comparison genotype frequencies and IPM predictions
plot(100,100,
     ylim=c(0,1),
     xlim=c(18,34),
     xlab='',ylab='',
     bty='l',cex.axis=1.2)
fig_label('E)',pos='topright',cex=1.5)
mtext('Equilibrium proportion Southern genotypes',2,line=2.2)
mtext('Temperature',1,line=2.2)

pred <- posterior_epred(modFreq,newdata=data.frame(temp=seq(-2,3,.1),day=29),re_formula=NA)
lines(seq(-2,3,.1)*scTemp[2]+scTemp[1],
      colMeans(pred),lwd=4)
lines(seq(-2,3,.1)*scTemp[2]+scTemp[1],
      apply(pred,2,quantile,prob=.025),lty=2)
lines(seq(-2,3,.1)*scTemp[2]+scTemp[1],
      apply(pred,2,quantile,prob=.975),lty=2)


## Add IPM predictions
points(allTemp*scTemp[2]+scTemp[1],propb,cex=2.5,col='red',lwd=2,type='b',pch=16)
points(allTemp*scTemp[2]+scTemp[1],propb,cex=2.5,col='black',pch=1)
legend('bottomright',col=c('black','#f03b20'),bty='n',lwd=2,
       legend=c('Genotype frequency model','IPM predictions'),cex=1.2)


dev.off()
