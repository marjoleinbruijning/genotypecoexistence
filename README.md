# Data and code for: Linking individual performance to density-dependent population dynamics to understand temperature-mediated genotype coexistence
<i>Marjolein Bruijning, Luc De Meester, Marco D. Visser, Erlend I. F. Fossen, Helene Vanvelk, Joost A. M. Raeymaekers, Lynn Govaert, Kristien I. Brans, Sigurd Einum, Eelke Jongejans</i>

## Dependencies
* `R` (version 4.5.0)
* Package `brms` (version 2.22.0)
* Package `wesanderson` (0.3.7)
* Package `vioplot` (0.4.0)
* Package `MASS` (7.3.55)
* Package `readr` (2.1.5)

## Get started
In order to run all code in R, start by getting all required files and folders in place:

* All .R files should be in directory called 'code'.
* The .csv file named 'dataDaphnia20241205' containing raw data should be in directory called 'data'.
* Open an R workspace and set its working directory to the folder that contains the folders 'code' and 'data'.


## Run the analysis
All code can be run from the file '/code/run.R', which consists of four parts:

### 1) Getting started
Load all dependencies, functions and data, and create folders to store output (if they do not exist yet). Here, also the genotype frequency model is fitted (can take a few minutes).

### 2) Fit vital rate models
Fit each of the seven vital rate models. Set `runall <- TRUE` to fit and save all vital rate models. This can take up to a day, and should be run only once. Running this analysis produces two .rds files per vital rate, stored in /Results:
* `modXXX-full.rds`: each of the tested models
* `modXXX.rds`: the selected model based on leave-one-out cross-validation

(here, `XXX` is a placeholder for the name of each vital rate)

Set `runall <- FALSE` to load the pre-fitted vital rate models that are stored in /Results. (Note that these files are provided on Github, so that all subsequent analyses can be performed without fitting the vital rate models)

### 3) Perform IPM analyses
This part runs the analyses that are presented in the manuscript: a) Assessing the temperature-dependent coexistence outcomes by combining Integral Projection Models (IPMs) with Modern Coexistence Theory (MCT); b) Running IPMs to simulate the dynamics of competing clones; c) Decomposition analysis to quantify the contribution of latitude-specific vital rates in determining coexistence outcomes. Each of these scripts stores the output in the /Results folder. These analyses can take a few hours to complete and require the fitted and loaded vital rate objects obtained under 2).

### 4) Create figures
Finally, the figures 2-4 presented in the manuscript can be created, using the objects created under 3).


## R-Code in code/…
This folder contains all .R files to run the analyses and create the figures presented in the manuscript:
* decomposition.R: Performs decomposition analysis. *Output: Results/decomposition.rds*
* figureCoexistence.R: Reproduces figure 3. Requires: Results/coexistence.rds. *Output: Figures/figure3.pdf*
* figureDecomp.R: Reproduces figure 4. Requires: Results/decomposition.rds. *Output: Figures/figure4.pdf*
* figureIsoclines.R: Reproduces figure 2. Requires: Results/coexistenceMean.rds; Results/ipmsimultrajec.RData. *Output: Figures/figure2.pdf*
* functions.R: All functions related to building the IPMs
* ipmsimul.R: Simulate population trajectories of competing genotypes. *Output: Results/ipmsimultrajec.RData*
* preparedata.R: load the raw data and prepare data for all subsequent analyses, fit Genotype Frequency model.
* run.R: Script with the full workflow.
* vectorplot.R: Apply MCT to IPM to assess temperature-mediated genotype coexistence. *Output: Results/coexistence.rds; Results/coexistenceMean.rds*
* vitalratesBrm.R: Run or load all vital rate models.



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
