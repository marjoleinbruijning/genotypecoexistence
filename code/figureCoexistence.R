

####################################################################
## Figure 3
####################################################################


dN <- readRDS('Results/coexistence.rds')

nicheoverlap <- fitnessdiff <- matrix(NA,nrow=length(allTemp),ncol=ndraws)

pars <- array(NA,dim=c(length(allTemp),4,ndraws),
              dimnames=list(allTemp,c('a12','a11','a22','a21'),1:ndraws))

for (k in 1:length(allTemp)) {

    for (j in 1:ncol(nicheoverlap)) {

        mod <- list()

        for (i in 1:2) { # First B, then N
            df <- contourLines(allN,allN,log(dN[,,i,k,j]),levels=0)
            ## get values where log lambda=0
            if (length(df) > 0) {
                df <- data.frame(n2=df[[1]]$y,n1=df[[1]]$x)
                mod[[i]] <- lm(n2~n1,data=df)
            }

        }
        if (length(mod) == 2 & all(unlist(lapply(mod,length)) > 0)) {

            ##  Get coefficients
            pars[k,'a12',j] <- 1/coef(mod[[1]])[1]
            pars[k,'a11',j] <- 1/ (-(coef(mod[[1]])[1]) / (coef(mod[[1]])[2]))

            pars[k,'a22',j] <- 1/coef(mod[[2]])[1]
            pars[k,'a21',j] <- 1/ (-(coef(mod[[2]])[1]) / (coef(mod[[2]])[2]))

        }
    }
}

ranges <- quantile(pars,na.rm=T,prob=c(.05,.95))
pars[pars < ranges[1] | pars > ranges[2]] <- NA

nicheoverlap <- apply(pars,c(1,3),function(x) sqrt((x['a21']*x['a12'])/(x['a11']*x['a22'])))
fitnessdiff <- apply(pars,c(1,3),function(x) sqrt((x['a21']*x['a22'])/(x['a11']*x['a12'])))


allniches <- seq(0,1,length.out=1001)
allk <-  seq(0,3.5,length.out=500)
mat <- outer(allniches,allk,
             function(x,y) x < y & y < 1/x)



pdf('Results/figure3.pdf')

par(mfrow=c(2,2),mar=c(1,1,1,1),oma=c(4,4,0,0))

for (i in 1:length(allTemp)) {

    exl <- which(is.na(nicheoverlap[i,]) | is.na(fitnessdiff[i,]))
    
    k <- with(df,MASS:::kde2d(nicheoverlap[i,-exl],fitnessdiff[i,-exl]))
    df <- contourLines(k$x,k$y,k$z,nlevels=10)

    plot(100,100,xlim=c(0,1),ylim=c(0,3),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
    contour(allniches,allk,mat,add=TRUE,labels='',lwd=2)
    image(allniches,allk,mat,col=c(NA,'lightgrey'),add=TRUE)
    
    for (j in 1:length(df)) {
        polygon(df[[j]]$x,df[[j]]$y,
                col=colorRampPalette(c('white',colsTemp[i]))(length(df)-1)[-1][j],
                border=NA)
    }
    contour(allniches,allk,mat,add=TRUE,drawlabels=FALSE,lwd=2,lty=1,col='darkgrey')

    if (i %in% 3:4) {
        axis(1)
    } else {
        axis(1,labels=NA)
    }

    if (i %in% c(1,3)) {
        axis(2)
    } else {
        axis(2,labels=NA)
    }
    
    
    legend('topright',bty='n',legend=paste0(round(allTemp[i]*scTemp[2]+scTemp[1]),'°C'),
           cex=1.2)

    if (i == 1) {
        text(c(.42,.42,.53),c(2.65,2.05,2.35),
             labels=c('Coexistence','Northern dominance','Southern dominance'),
             col=c('black',rev(colsCountry)),cex=1)

        arrows(x0=.78,x1=.85,y0=2.35,y1=2.25,col=colsCountry[1],length=.1)
        
        arrows(x0=.2,x1=.25,y0=1.9,y1=.1,col=colsCountry[2],length=.1)

    } 
}

mtext('Niche overlap',1,outer=TRUE,line=1)
mtext('Competitive ability difference',2,outer=TRUE,line=1)

for (i in 1:length(allTemp)) {

    par(fig = c(c(0,.5,0,0.5)[i], #left
                c(.2,.7,.2,.7)[i], # right
                c(0.8,.8,0.3,0.3)[i],# bottom
                c(1,1,0.5,0.5)[i]), # top
        new = TRUE)

    ## Coexistence outcomes
    exl <- which(is.na(nicheoverlap[i,]) | is.na(fitnessdiff[i,]))
    
    coexist <- (nicheoverlap[i,-exl] < fitnessdiff[i,-exl] &
                fitnessdiff[i,-exl]  < 1/nicheoverlap[i,-exl])

    out <- c(sum(coexist),table(factor(fitnessdiff[i,-exl][!coexist] > 1,levels=c(FALSE,TRUE))))
    pie(out,col=c('grey',rev(colsCountry)),labels='')
    
}
    
dev.off()
