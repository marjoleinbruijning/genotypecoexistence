# Data and code for: Towards a mechanistic understanding of genotype coexistence under warming (in prep.)
<i>Marjolein Bruijning, Luc De Meester, Marco D. Visser, Erlend I. F. Fossen, Helene Vanvelk, Joost A. M. Raeymaekers, Lynn Govaert, Kristien I. Brans, Sigurd Einum, Eelke Jongejans</i>

## Get started
In order to run all code in R, all .R files should be in directory called 'code' (placed in the in working directory); the data file named 'dataDaphnia20241205.csv' containing raw data should be in directory called `data`; a directory named 'Results' should exist to store all output.

## R-Code in code/…
Code to run all analyses and create all figures presented in the manuscript. All code runs from file run.R.

## Datafile data/dataDaphnia20241205.csv:
Every row in the dataset represents data on one individual. The dataset contains the following columns:

* day: day of the experiment when first measuring the individual. Individual measurements started at day 7 of the experiment, and continued until day 29.
* interval: time intervals between measurements (3 or 4 days)
* popNo: aquarium number (ranging between 1-114).
* temp: experimental temperature (in C).
* pop1: clonal treatment of the aquarium (N1 and N2 refer to clones originating from Norwegian ponds, B1 and B2 refer to clones originating from Belgian ponds).
* pop2: in the case of a mixed-population treatment, this column gives the second clonal treatment. In the case of a single-population treatment, set to NA.
* gen: location of origin of the measured individual, either identified using genetic markers or known if the aquarium contained only individuals from one location.
* n: population density as estimated using R-package trackdem.
* size: body size when first measuring the individual (mm).
* surv: survival status after one time interval (0 = dead, 1 = alive).
* growth: body size increment (mm/day).
* carryingeggs: binary variable indicating whether the individual carried eggs when first measuring (0 = no, 1 = yes).
* offspringproduced: binary variable indicating whether the individual produced any neonates when remeasuring, conditional on carrying eggs (0 = no, 1 = yes).
* clutchsize: number of produced offspring, conditional on producing offspring.
* sizeoff: body size of one of the neonates.
* clone: clonal identity of the measured individual, either identified using genetic markers or known if the aquarium contained only individuals from one latitude.
* sexoff: offspring sex (m/f).
* sex: sex (m/f) of the measured individual.
* country: location of origin of the measured individual at the latitude level (B=Belgium, N=Norway).
* countryBin: location of origin of the measured individual at the latitude level denoted as binary variable (0=Norway, 1=Belgium).
