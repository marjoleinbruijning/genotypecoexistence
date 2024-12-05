

###########################################################################
## All code for the analysis and figures presented in:
## Towards a mechanistic understanding of genotype coexistence under warming
## Bruijning et al., in prep.
###########################################################################

## All .R files should be in directory called 'code' (in working directory)
## .csv file named 'dataDaphnia20241205' containing raw data should be in directory called data
## directory named 'Results' should exist to store all output


###########################################################################
## Prepare everything
###########################################################################

require(brms); require(wesanderson) ## Load packages
source('code/functions.R') ## Load functions
source('code/preparedata.R') ## Load and prepare data, fit genotype freq model


#########################################################
## Fit vital rate models
#########################################################

## General settings
niter <- 50000
ncores <- 3
nchains <- 3
thin <- 10

runall <- TRUE # fit and save all vital rate models? (takes long, run only once)
source('code/vitalratesBrm.R') ## saves rds files for each vital rates

## Prepare everything for building IPM
nc <- 100 # number of size classes
allsizes <- seq(-2.5,2.5,length.out=nc) # size range, corresponds to: range((allsizes * scSize[2]) + scSize[1])
cw <- allsizes[2]-allsizes[1] # class width
resGrowth <- sd((modGrowth$data$growth-fitted(modGrowth)[,1]))
resOff <- sd((modOffsize$data$sizeoff-fitted(modOffsize)[,1]))


#########################################################
## Perform all IPM analysis and create figures
#########################################################

## Settings
n <- 6  # number of sym- and allopatric densities
allN <- seq(1,1000,length.out=n) # density range
allTemp <-  (c(18,22,26,32)-scTemp[1]) / scTemp[2] # evaluated temperatures (standardized)
ndraws <- 1000 # number of posterior draws
colsCountry <- c('#7570b3','#1b9e77')
names(colsCountry) <- c('B','N')
colsTemp <- wesanderson::wes_palette("Zissou1", 4, type = "continuous")


## Coexistence outcomes
source('code/vectorplot.R') # creates coexistence.rds and coexistenceMean.rds
source('code/figureIsoclines.R') # creates figure2.pdf
source('code/figureCoexistence.R') # creates figure3.pdf


## Decomposition
source('code/decomposition.R') # creates decomposition.rds
source('code/figureDecomp.R') # creates figure4.pdf


## Clonal variation in competition outcomes
source('code/clonalvar.R') # creates clonalvariation.rds
source('code/figureClonalvar.R') # creates figure5.pdf

#########################################################
