

dat <- read.csv('data/dataDaphnia20241205.csv') ## load data

## Scale variables
dat$temp <- scale(dat$temp)
dat$size <- scale(dat$size)

scTemp <- unlist(attributes(dat$temp)[2:3])
scSize <- unlist(attributes(dat$size)[2:3])

dat$sizeoff <- (dat$sizeoff-scSize[1]) / scSize[2]
dat$temp <- c(dat$temp)
dat$size <- c(dat$size)

dat$size2 <- dat$size^2
dat$off0 <- dat$clutchsize - 1
dat$sexoff <- ifelse(dat$sexoff == 'f',1,0)


## Create data frame that contains all treatments
popinfo <- data.frame(popNo=1:114,
                     temps=tapply(dat$temp,dat$popNo,function(x) x[1]),
                     pop1=tapply(dat$pop1,dat$popNo,function(x) x[1]),
                     pop2=tapply(dat$pop2,dat$popNo,function(x) x[1]))


## Within-latitude combinations and single genotypes
singlepop <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 
18, 19, 20, 21, 22, 23, 24, 25, 35, 42, 46, 47, 48, 49, 50, 51, 
52, 53, 54, 55, 56, 57, 58, 59, 69, 76, 80, 81, 82, 83, 84, 85, 
86, 87, 88, 89, 90, 91, 92, 96, 105, 111, 114, 26, 29, 32, 43, 
44, 45, 60, 63, 66, 77, 78, 79, 93, 97, 100, 110, 112, 113)

## Across latitude combinations
mixpop <- c(1:114)[!1:114 %in% singlepop]


## Fit model for genotype frequency as function of temp and day
inc <- dat$popNo %in% mixpop & !is.na(dat$countryBin)
modFreq <- brm(countryBin ~ day + day:temp-1 + (0+day|popNo),data=dat[inc,],
               family='bernoulli',iter=5000,cores=3,chains=3,thin=10)


## Predict proportion B and N based on fitted model to fill in hetero- and con densities
inc <- dat$popNo %in% mixpop
tmp <- dat[inc,]
tmp$popNo <- factor(tmp$popNo) # reset factors to only include mixed pops
predB <- fitted(modFreq,newdata=tmp,type='response')[,1] # predicted proportion Belgians
dat$predB <- NA # add new column to data
dat$predB[inc] <- predB

dat$predB[dat$pop1 %in% 7:12 & (is.na(dat$pop2)|dat$pop2%in%7:12) ] <- 1 # single country treatments
dat$predB[dat$pop1 %in% 1:6 & (is.na(dat$pop2)|dat$pop2%in%1:6) ] <- 0

dat$propCon <- NA # add new column
dat$propCon[dat$country == 'B' & !is.na(dat$country)] <- dat$predB[dat$country == 'B' & !is.na(dat$country)]
dat$propCon[dat$country == 'N' & !is.na(dat$country)] <- 1-dat$predB[dat$country == 'N' & !is.na(dat$country)]

dat$propHetero <- 1-dat$propCon
dat$nCon <- dat$n * dat$propCon
dat$nHetero <- dat$n * dat$propHetero
