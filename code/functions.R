
################################################################################
## Functions
################################################################################

################################################################################
## Vital rate functions
################################################################################
## Model selection for each vital rate, and return only best model
selectmodel <- function (allmods) {

    modname <- deparse(substitute(allmods))

    npar <- unlist(lapply(allmods,function(x) nrow(fixef(x))))
    loo <- eval(parse(text=c('loo_compare(',c(paste0(modname,'[[',1:(length(allmods)-1),']],'),
                                              paste0(modname,'[[',(length(allmods)),']]')),
                             ')')))
    names(npar) <- paste0(modname,'[[',1:(length(allmods)),']]')

    namess <- rownames(as.matrix(loo))[[1]]
    return(allmods[[readr::parse_number(namess)]]) # return best model
}

################################################################################
## IPM functions
################################################################################
pxy <- function (mod=modSurv,df,ndraws=1) {

    reform <- ifelse(is.na(df$clone[1]),NA,'~(1|clone)')
    pred <- (posterior_linpred(mod,transform=TRUE,df,ndraws=ndraws,re_formula=reform))
    return(pred)

}

gxy <- function (mod=modGrowth,df,ndraws=1) {

    reform <- ifelse(is.na(df$clone[1]),NA,'~(1|clone)')
    pred <- (posterior_linpred(mod,transform=TRUE,df,ndraws=ndraws,re_formula=reform))
    return(pred)
    
}


fxy1 <- function (mod=modEggs,df,ndraws=1) {

    reform <- ifelse(is.na(df$clone[1]),NA,'~(1|clone)')
    R1 <- (posterior_linpred(mod,transform=TRUE,df,ndraws=ndraws,re_formula=reform))
    return(R1)

}

fxy2 <- function (mod=modRepr,df,ndraws=1) {

    reform <- ifelse(is.na(df$clone[1]),NA,'~(1|clone)')
    R2 <- (posterior_linpred(mod,transform=TRUE,df,ndraws=ndraws,re_formula=reform))
    return(R2)

}

fxy3 <- function (mod=modClutch,df,ndraws=1) {

    reform <- ifelse(is.na(df$clone[1]),NA,'~(1|clone)')
    R3 <- (posterior_linpred(mod,transform=TRUE,df,ndraws=ndraws,re_formula=reform) + 1)
    return(R3)

}

fxy4 <- function (mod=modSexoff,df,ndraws=1) {

    reform <- ifelse(is.na(df$clone[1]),NA,'~(1|clone)')
    R4 <- (posterior_linpred(mod,transform=TRUE,df,ndraws=ndraws,re_formula=reform))
    return(R4)
    
}

dxy <- function (mod=modOffsize,df,ndraws=1) {

    reform <- ifelse(is.na(df$clone[1]),NA,'~(1|clone)')
    pred <- (posterior_linpred(mod,transform=TRUE,df,ndraws=ndraws,re_formula=reform))
    return(pred)

}

buildVitalrates <- function (treatments,allsizes,
                             ndraws=1) {

    ## All body sizes
    df <- data.frame(size=allsizes,size2=allsizes^2,sex='f')

    ## Combine with treatments
    tmp <- as.data.frame(lapply(df, rep, each = nrow(treatments)))
    df <- cbind(tmp,treatments)
    df$tempSq <- df$temp

    ## Survival
    cat('\r','Vital rate 1/7')
    P <- pxy(df=df,ndraws=ndraws)

    ## Growth
    cat('\r','Vital rate 2/7')
    growth <- gxy(df=df,ndraws=ndraws)

    ## Reproduction
    cat('\r','Vital rate 3/7')
    R1 <- fxy1(df=df,ndraws=ndraws)
    
    cat('\r','Vital rate 4/7')
    R2 <- fxy2(df=df,ndraws=ndraws)
    
    cat('\r','Vital rate 5/7')
    R3 <- fxy3(df=df,ndraws=ndraws)
    
    cat('\r','Vital rate 6/7')
    R4 <- fxy4(df=df,ndraws=ndraws)
    
    ## Offspring size
    cat('\r','Vital rate 7/7')
    offsize <- dxy(df=df,ndraws=ndraws)

    return(list(P=P,
                growth=growth,
                R1=R1,
                R2=R2,R3=R3,R4=R4,
                offsize=offsize,df=df))
}

createIPM <- function (allsizes,P,growth,R1,R2,R3,R4,offsize,nc) {

    ## Survival
    P <- matrix(P,ncol=nc,nrow=nc,byrow=TRUE)
    
    ## Growth
    allsizesUnsc <- (allsizes * scSize[2]) + scSize[1]
    newsize <- allsizesUnsc + growth # predicted new size (in mm; unscaled)
    G <- sapply(1:length(newsize),function(x) {
        p <- dnorm(allsizesUnsc,newsize[x],resGrowth) * diff(allsizesUnsc[1:2])
        if (x < (length(newsize)/2)) {
            p[1] <- p[1] + (1-sum(p))
        } else {
            p[length(p)] <- p[length(p)] + (1-sum(p))
        }
        return(p)
    })

    ## Reproduction
    R <- matrix(R1*R2*R3*R4,ncol=nc,nrow=nc,byrow=TRUE)

    ## Offspring size
    D <- sapply(1:length(offsize),function(x) {
        p <- dnorm(allsizes,offsize[x],resOff) * cw
        p[1] <- p[1] + (1-sum(p))
        return(p)
    })


    # # Convert to daily IPMs
    P <- P^(1/3.5)
    R <- R / 3.5

    ## Construct full IPM
    ipm <- P*G + R*D

    return(list(ipm=ipm,G=G,P=P,D=D,R=R,R1=R1,R2=R2,R3=R3,R4=R4))
}


getLambda <- function (A) {
    ev <- eigen(A)
    lmax <- which.max(Re(ev$values))
    lambda <- Re(ev$values[lmax])
    return(lambda)
}


################################################################################
## PLOT LABELS
## from: https://www.r-bloggers.com/2017/03/adding-figure-labels-a-b-c-in-the-top-left-corner-of-the-plotting-region/
################################################################################

fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {

  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))

  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")

    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }

  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }

  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100

  x1 <- switch(pos,
    topleft     =x[1] + sw, 
    left        =x[1] + sw,
    bottomleft  =x[1] + sw,
    top         =(x[1] + x[2])/2,
    center      =(x[1] + x[2])/2,
    bottom      =(x[1] + x[2])/2,
    topright    =x[2] - sw,
    right       =x[2] - sw,
    bottomright =x[2] - sw)

  y1 <- switch(pos,
    topleft     =y[2] - sh,
    top         =y[2] - sh,
    topright    =y[2] - sh,
    left        =(y[1] + y[2])/2,
    center      =(y[1] + y[2])/2,
    right       =(y[1] + y[2])/2,
    bottomleft  =y[1] + sh,
    bottom      =y[1] + sh,
    bottomright =y[1] + sh)

  old.par <- par(xpd=NA)
  on.exit(par(old.par))

  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}
