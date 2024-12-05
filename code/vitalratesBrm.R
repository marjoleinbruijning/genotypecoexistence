

########################################################################
## Structure of all models to be tested
########################################################################
## List containing all models to be tested
forms <- list("size + size2 + temp + country + nCon + nHetero + (1|clone)", 
    "size + size2 + temp + country + sex + nCon + nHetero + (1|clone)", 
    "size + size2 + temp + country + nCon + nCon * temp + nHetero + nHetero * temp + (1|clone)", 
    "size + size2 + temp + country + sex + nCon + nCon * temp + nHetero + nHetero * temp + (1|clone)", 
    "size + size2 + temp + country + nCon + nCon * country + nHetero + nHetero * country + (1|clone)", 
    "size + size2 + temp + country + sex + nCon + nCon * country + nHetero + nHetero * country + (1|clone)", 
    "size + size2 + temp + country + nCon + nCon * temp + nCon * country + nHetero + nHetero * temp + nHetero * country + (1|clone)", 
    "size + size2 + temp + country + sex + nCon + nCon * temp + nCon * country + nHetero + nHetero * temp + nHetero * country + (1|clone)", 
    "size + size2 + temp + country + nCon + nCon * temp + nCon * country + nCon * temp * country + nHetero + nHetero * temp + nHetero * country + nHetero * temp * country + (1|clone)", 
    "size + size2 + temp + country + sex + nCon + nCon * temp + nCon * country + nCon * temp * country + nHetero + nHetero * temp + nHetero * country + nHetero * temp * country + (1|clone)")

incM <- c(2,4,6,8,10) ## models that include sex (for survival and growth)


#################################################################
## Fit models
#################################################################

if (runall) {
    
    modSexoff <- lapply(forms[-incM],function(x) {
        mod <- brm(as.formula(paste('sexoff ~ ',x)),data=dat,
                   family='bernoulli',
                   iter=niter,cores=ncores,chains=nchains,thin=thin,
                   silent=2,refresh=0,open_progress=FALSE,
                   save_pars = save_pars(all = TRUE)
                   )
        mod
    })
    saveRDS(modSexoff,file='Results/modSexoff-full.rds')
    rm(modSexoff); gc()


    modOffsize <- lapply(forms[-incM],function(x) {
        mod <- brm(as.formula(paste('sizeoff ~ ',x)),data=dat,
                   iter=niter,cores=ncores,chains=nchains,thin=thin,silent=2,
                   refresh=0,
                   open_progress=FALSE,
                   save_pars = save_pars(all = TRUE)
                   )
        mod
    })
    saveRDS(modOffsize,file='Results/modOffsize-full.rds')
    rm(modOffsize); gc()

    
    modSurv <- lapply(forms[incM],function(x) {
        mod <- brm(as.formula(paste('surv ~ ',x)),data=dat,
                   family='bernoulli',
                   iter=niter,cores=ncores,chains=nchains,thin=thin,
                   silent=2,refresh=0,open_progress=FALSE,
                   save_pars = save_pars(all = TRUE))
        mod
    })
    saveRDS(modSurv,file='Results/modSurv-full.rds')
    rm(modSurv); gc()


    modEggs <- lapply(forms[-incM],function(x) {
        mod <- brm(as.formula(paste('carryingeggs ~ ',x)),data=dat,
                   family='bernoulli',
                   iter=niter,cores=ncores,chains=nchains,thin=thin,
                   silent=2,refresh=0,open_progress=FALSE,
                   save_pars = save_pars(all = TRUE))
        mod
    })
    saveRDS(modEggs,file='Results/modEggs-full.rds')
    rm(modEggs); gc()

    
    modRepr <- lapply(forms[-incM],function(x) {
        mod <- brm(as.formula(paste('offspringproduced ~ ',x)),data=dat,
                   family='bernoulli',
                   iter=niter,cores=ncores,chains=nchains,
                   thin=thin,silent=2,refresh=0,open_progress=FALSE,
                   save_pars = save_pars(all = TRUE))
        mod
    })
    saveRDS(modRepr,file='Results/modRepr-full.rds')
    rm(modRepr); gc()

    
    modClutch <- lapply(forms[-incM],function(x) {
        mod <- brm(as.formula(paste('off0 ~ ',x)),data=dat,
                   family='poisson',
                   iter=niter,cores=ncores,chains=nchains,thin=thin,
                   silent=2,refresh=0,open_progress=FALSE,
                   save_pars = save_pars(all = TRUE))
        mod
    })
    saveRDS(modClutch,file='Results/modClutch-full.rds')
    rm(modClutch); gc()

    
    modGrowth <- lapply(forms[incM],function(x) {
        mod <- brm(as.formula(paste('growth ~ ',x)),data=dat,
                   iter=niter,cores=ncores,chains=nchains,thin=thin,
                   silent=2,refresh=0,open_progress=FALSE,
                   save_pars = save_pars(all = TRUE))
        mod
    })
    saveRDS(modGrowth,file='Results/modGrowth-full.rds')
    rm(modGrowth); gc()

    
    ##################################################################################
    ## Perform model selection
    ##################################################################################
    modnames <- c('modSurv','modGrowth','modEggs','modRepr','modClutch','modSexoff','modOffsize')
    
    for (i in 1:length(modnames)) {

        mod <- readRDS(paste0('Results/',modnames[i],'-full.rds'))
        
        for (j in 1:length(mod)) {
            mod[[j]] <- add_criterion(mod[[j]],'loo'))
            cat(all(rhat(mod[[j]]) < 1.01),'\n')
        }
        
        mod <- selectmodel(mod)
        saveRDS(mod,file=paste0('Results/',modnames[i],'.rds'))
    }
}


##################################################################################
## Load previously saved models for subsequent analysis
##################################################################################
modSurv <- readRDS('Results/modSurv.rds')
modGrowth <- readRDS('Results/modGrowth.rds')
modEggs <- readRDS('Results/modEggs.rds')
modRepr <- readRDS('Results/modRepr.rds')
modClutch <- readRDS('Results/modClutch.rds')
modSexoff <- readRDS('Results/modSexoff.rds')
modOffsize <- readRDS('Results/modOffsize.rds')

