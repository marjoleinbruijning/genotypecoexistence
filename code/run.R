

###########################################################################
## All code for the analysis and figures presented in:
## Linking individual performance to density-dependent population dynamics
## to understand temperature-mediated genotype coexistence
## Bruijning et al., in prep.
###########################################################################


## Step-by-step instructions can be found in the README file


###########################################################################
## 1) Getting started
###########################################################################

## Create folders Results and Figures if they do not exist yet to store output
if (!file.exists('Results')) {dir.create('Results')}
if (!file.exists('Figures')) {dir.create('Figures')}


## General settings for the Bayesian regression models
niter <- 50000
ncores <- 3
nchains <- 4
thin <- 10

## Load required packages
require(brms); require(wesanderson) 

## Load functions
source('code/functions.R')

## Load and prepare data, and fit genotype freq model (takes a few minutes)
source('code/preparedata.R')


#########################################################
## 2) Fit vital rate models
#########################################################

runall <- FALSE
## If TRUE, fit and save all vital rate models (takes up to a day)
## If FALSE, existing Rmd files will be loaded from folder 'Results'

source('code/vitalratesBrm.R') 


#########################################################
## 3) Perform IPM analyses
#########################################################

## Requires fitted vital rate objects

## Settings across all IPM analyses
nc <- 100 # number of size classes
allsizes <- seq(-2.5,2.5,length.out=nc) # size range, corresponds to: range((allsizes * scSize[2]) + scSize[1])
cw <- allsizes[2]-allsizes[1] # class width
resGrowth <- sd((modGrowth$data$growth-fitted(modGrowth)[,1])) # residual growth
resOff <- sd((modOffsize$data$sizeoff-fitted(modOffsize)[,1])) # residual offspring body size

n <- 6  # number of sym- and allopatric densities
allN <- seq(1,1000,length.out=n) # density range
allTemp <-  (c(18,22,26,32)-scTemp[1]) / scTemp[2] # evaluated temperatures (standardized)
ndraws <- 1000 # number of posterior draws


## Temperature-dependent coexistence
## Output: Results/coexistence.rds and Results/coexistenceMean.rds
source('code/vectorplot.R')

## IPM simulations for each temperature
## Output: Results/ipmsimultrajec.RData
source('code/ipmsimul.R') 

## Decomposition analysis
## Output: Results/decomposition.rds
source('code/decomposition.R') 


#########################################################
## 4) Create figures
#########################################################

## Uses IPM settings set above

## Color settings
colsCountry <- c('#7570b3','#1b9e77')
names(colsCountry) <- c('B','N')
colsTemp <- wesanderson::wes_palette("Zissou1", 4, type = "continuous")

## Figure 2
## Requires Results/coexistenceMean.rds and Results/ipmsimultrajec.RData
## Output: Figures/figure2.pdf
source('code/figureIsoclines.R') 

## Figure 3
## Requires Results/coexistence.rds
## Output: Figures/figure3.pdf
source('code/figureCoexistence.R') 

## Figure 4
## Requires: Results/decomposition.rds
## Output: Figures/figure4.pdf
source('code/figureDecomp.R') 


#########################################################
