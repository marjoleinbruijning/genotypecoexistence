
####################################################################
## Figure 5
####################################################################


dN <- readRDS(paste0('Results/clonalvariation.rds'))

## Panel a) coexistence outcomes
combis <- expand.grid(pop1=7:12,pop2=1:6)
nicheoverlap <- fitnessdiff <- propb <- matrix(NA,nrow=length(allTemp),ncol=nrow(combis))

pars <- array(NA,dim=c(length(allTemp),4,nrow(combis)),
              dimnames=list(allTemp,c('a12','a11','a22','a21'),1:nrow(combis)))

for (k in 1:length(allTemp)) {

    for (j in 1:ncol(nicheoverlap)) {

        mod <- list()
        for (i in 1:2) { # First B, then N
            df <- contourLines(allN,allN,log(dN[,,i,k,combis[j,i]]),levels=0)
            ## get values where log lambda=0
            if (length(df) > 0) {
                df <- data.frame(n2=df[[1]]$y,n1=df[[1]]$x)
                mod[[i]] <- lm(n2~n1,data=df)
            }
        }
        
        if (all(unlist(lapply(mod,length)) > 0)) {

            ##  Get coefficients
            pars[k,'a12',j] <- 1/coef(mod[[1]])[1]
            pars[k,'a11',j] <- 1/ (-(coef(mod[[1]])[1]) / (coef(mod[[1]])[2]))

            pars[k,'a22',j] <- 1/coef(mod[[2]])[1]
            pars[k,'a21',j] <- 1/ (-(coef(mod[[2]])[1]) / (coef(mod[[2]])[2]))

            ## Predict eq proportions (where do nullclines intersect?)
            beq <- (coef(mod[[1]])[1] - coef(mod[[2]])[1]) / ( coef(mod[[2]])[2] - coef(mod[[1]])[2] )
            neq <- predict(mod[[1]],newdata=data.frame(n1=beq))

            if (beq < 0) beq <- 0
            if (neq < 0) neq <- 0
            propb[k,j] <- beq / (beq + neq)            
        }
        
    }
}

nicheoverlap <- apply(pars,c(1,3),function(x) sqrt((x['a21']*x['a12'])/(x['a11']*x['a22'])))
fitnessdiff <- apply(pars,c(1,3),function(x) sqrt((x['a21']*x['a22'])/(x['a11']*x['a12'])))


## Prepare dataframe
scenarios <- do.call("rbind", replicate(4, combis, simplify = FALSE)) ## for each clone
scenarios$temp <- rep(allTemp,each=36)

scenarios$code <- paste(scenarios$pop2,scenarios$pop1,sep='-')
popinfo$code <- paste(popinfo$pop1,popinfo$pop2,sep='-')
scenarios$id  <- 1:nrow(scenarios)
scenarios <- merge(scenarios,popinfo[,c('code','popNo')],by='code')
scenarios <- scenarios[order(scenarios$id), ]

predmod <- posterior_epred(modFreq,newdata=cbind(scenarios,data.frame(day=29)))
scenarios$pred <- apply(predmod,2,quantile,prob=.5)
scenarios$ipmpred <- c(t(propb))
scenarios$nicheoverlap <- c(t(nicheoverlap))
scenarios$fitnessdiff <- c(t(fitnessdiff))


## Create plot
pdf('Results/figure5.pdf',width=12,height=9)

ly <- matrix(0,nrow=12,ncol=22)
ly[2:7,1:7] <- 1
ly[2:7,10:15] <- 2

## Density kernels
ly[1,10:15] <- 3
ly[2:7,16] <- 4

## Barplot
ly[2:7,18:22] <- 5

## All vital rate competitivneess correlations
ly[9:10,2:3] <- 6
ly[9:10,5:6] <- 7
ly[9:10,8:9] <- 8
ly[9:10,11:12] <- 9
ly[9:10,14:15] <- 10
ly[9:10,17:18] <- 11
ly[9:10,20:21] <- 12

ly[11:12,2:3] <- 13
ly[11:12,5:6] <- 14
ly[11:12,8:9] <- 15
ly[11:12,11:12] <- 16
ly[11:12,14:15] <- 17
ly[11:12,17:18] <- 18
ly[11:12,20:21] <- 19


layout(ly)

par(mar=c(0,0,0,0),oma=c(1,4,1,1))

## Colors
cols <- rev(c('#ca0020','#f4a582','#92c5de','#0571b0'))
names(cols) <- allTemp
colsclone <- c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02')
pchclone <- 1:6


allniches <- seq(0,1,length.out=1001)
allk <-  seq(0,3.5,length.out=500)
mat <- outer(allniches,allk,
             function(x,y) x < y & y < 1/x)
plot(100,100,xlim=c(0,1),ylim=c(0,3),xlab='',ylab='',
     xaxs='i',yaxs='i',cex.axis=1.2)
fig_label('A)',pos='topleft',cex=1.5)
contour(allniches,allk,mat,add=TRUE,labels='',lwd=2)
image(allniches,allk,mat,col=c(NA,'lightgrey'),add=TRUE)

mtext('Niche overlap',1,line=2.3,cex=1)
mtext('Fitness difference',2,line=2.3,cex=1)

legend('bottomright',col=cols,pch=1,legend=as.numeric(names(cols))*scTemp[2] + scTemp[1],
       title='Temperature',cex=1.2)
legend('topleft',pch=pchclone,title='Northern clones',legend=paste0('N',1:6),bg='white',cex=1.2)
legend('topright',pch=16,title='Southern clones',legend=paste0('S',1:6),col=colsclone,
       bg='white',cex=1.2)


points(scenarios$nicheoverlap,scenarios$fitnessdiff,
       col=cols[match(scenarios$temp,names(cols))],
       cex=2,lwd=2)
points(scenarios$nicheoverlap,scenarios$fitnessdiff,
       col=colsclone[scenarios$pop1-6],
       pch=pchclone[scenarios$pop2],lwd=1)
       


## Plot IPM predictions vs modFreq predictions on competitive outcome
plot(scenarios$ipmpred,scenarios$pred,ylim=c(0,1),xlim=c(0,1),pch=1,
     cex=2,
     xlab='',
     ylab='',
     col=cols[match(scenarios$temp,names(cols))],cex.axis=1.2,lwd=2)
fig_label('B)',pos='topleft',cex=1.5)

points(scenarios$ipmpred,scenarios$pred,
       col=colsclone[scenarios$pop1-6],
       lwd=1,cex=1,
       pch=pchclone[scenarios$pop2])
lines(-1:2,-1:2,lty=2,lwd=2)

mtext('Proportion Southern clones (IPM)',1,line=2.3,cex=1)
mtext('Proportion Southern clones (genotype frequency model)',2,line=2.3,cex=1)

legend('topleft',bty='n',
       legend=paste0('r= ',round(sqrt(summary(lm(pred~ipmpred,scenarios))$r.squared),2)),
       cex=1.2)


## Density kernels
plot(100,100,xlim=c(0,1),ylim=c(0,5),xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
for (i in 1:4) {
    inc <- scenarios$temp == allTemp[i]
    d <- density(scenarios$pred[inc],bw=.1)
    polygon(d, col = paste0(cols[i],50),border=NA)

}

plot(100,100,ylim=c(0,1),xlim=c(0,5),xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
for (i in 1:4) {
    inc <- scenarios$temp == allTemp[i]
    d <- density(scenarios$ipmpred[inc],bw=.1)
    polygon(d$y,d$x, col = paste0(cols[i],50),border=NA)

}


## Average proportions
a <- c(1-tapply(scenarios$pred,scenarios$pop2,mean),tapply(scenarios$pred,scenarios$pop1,mean))
b <- c(1-tapply(scenarios$ipmpred,scenarios$pop2,mean),tapply(scenarios$ipmpred,scenarios$pop1,mean))

x <- barplot(b,ylim=c(0,1),col=c(rep('grey',6),colsclone),names.arg=rep(1:6,2),
             cex.axis=1.2)
fig_label('C)',pos='topleft',cex=1.5)

abline(v=mean(x[6:7]),lwd=2)
abline(h=0)

points(x[1:6],rep(.1,6),pch=pchclone,cex=2.5)
         
mtext('Northern clones     Southern clones',1,line=2)
mtext('Competitiveness',2,line=2)


pchclone <- c(pchclone,rep(16,6))
colsclone <- c(rep('black',6),colsclone)

mods <- c('modSurv','modGrowth','modEggs','modRepr','modClutch','modSexoff','modOffsize')
namess <- c('Surival','Growth',
            'Carrying eggs','Egg development','Clutch size',
            'Female probability','Offspring size')


par(mar=c(1,0,1,0))

## Add panel plot with competitivness
for (j in 1:2) {

    if (j == 1) inc <- 1:6
    if (j == 2) inc <- 7:12
    
    for (i in 1:7) {
        y <- b[inc]
        x <- ranef(get(mods[i]))$clone[,,1][,1][inc]
        mod <- lm(y ~ x)
        p <- summary(mod)$coefficients[2,4]
        plot(x,y,bty='l',cex=ifelse(j==1,1.5,2),xaxt='n',yaxt='n',
             pch=pchclone[inc],col=colsclone[inc],xlab='',ylab='')
        if (j == 1) mtext(namess[i],3,line=.5)

        if (j == 1 & i == 1) fig_label('D)',pos='topleft',cex=1.5)

        if (p < .01) {
            lines(-10:10,predict(mod,data.frame(x=-10:10)),lwd=2)        
        } else if (p < .05) {
            lines(-10:10,predict(mod,data.frame(x=-10:10)),lwd=2,lty=2)
        } 
        
        if (i == 1) mtext('Competitiveness',2,outer=FALSE,line=1)

    }
}

mtext('Random clone effect',1,outer=TRUE,line=0)


dev.off()
