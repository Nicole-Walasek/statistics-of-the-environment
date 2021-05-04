# load packages -----------------------------------------------------------


if (!require("pacman"))
  install.packages("pacman")

pacman::p_load(
  ggplot2,
  ggthemes,
  boostmtree,
  MargCond,
  simstudy,
  lme4,
  MASS,
  sjPlot,
  effects,
  dplyr,
  reshape,
  scales,
  Hmisc,
  gsubfn,
  ggpubr
)

# set working directory

wdPath = file.path(
  "C:",
  "Users",
  "nicwa",
  "Dropbox",
  "PhD",
  "Publications",
  "statistics of the environment",
  "statistics_of_the_environment_code"
)
setwd(wdPath)

# load my own R functions

source('simRMData.R')


# generate data -----------------------------------------------------------


nSample = 100 # number of subjects

# between participant variances in intercept and slope, and covariance between the two
SigmaModel = rbind(c(1, 0.2),
                   c(0.2, 0.3))

popInt = 2
popSlope = 1.5

nObs = 50 # number of repeated measures

# changepoints in the mean
# number of change points; variance in the timing of change points;
# mean change in slope at a change point; variance in the change of the slope
changePoints = as.integer(3) # should be 'None if there is None'
changePointsVar = 0

slopeChangeMean = c(-0.5, 0.7, -0.1)
slopeChangeVar = 0.1


# change points in the variance
# change in variance; number of variance change points  

# ppVar can increasing, decreasing or a specified function 
varFun <- function(x) {
  (x ** 2)
}

ppVar = c(4, -2.5) #c(2, -2.5)#"increasing" 

ppVarChangePoints = as.integer(2)

noiseVar = 1 # noise around each participants mean trend

rangeMin = 1
rangeMax = 10

missingness = 'None'



list[data, populationChangePoints, populationPPVarChangePoints] <-
  simData(
    nSample,
    popInt,
    popSlope,
    sigmaModel,
    nObs,
    changePoints,
    changePointsVar,
    slopeChangeMean,
    slopeChangeVar,
    noiseVar,
    rangeMin,
    rangeMax,
    missingness,
    ppVar,
    ppVarChangePoints
  )

head(data, 70)

print(populationChangePoints)
print(populationPPVarChangePoints)
# plot raw data and descriptives ------------------------------------------


numParticipants = 30

p <-
  plotSimData(data, numParticipants, nSample, nObs,
              populationChangePoints,
              populationPPVarChangePoints)

print(p)

# these settings look good just need to change the path
ggsave(
  device = 'pdf',
  dpi = 600 ,
  width = 12,
  height = 8,
  filename = "figureStatsEnvVarChangePointsandSlopeChangePoints_50RM_Extreme.pdf"
)



# saving the dataframe if zou like it ----------------------------------------------------

write.csv(data, "simData_VarChangePointsandSlopeChangePoints_50RM_Extreme.csv", row.names=TRUE) 
  
# fit model ---------------------------------------------------------------


# fitting a mixed effects model

modelFitted <- lme4::lmer(y ~ time + (time | ppID), data)

# let's look at the mode results

summary(modelFitted)
sjPlot::tab_model(modelFitted)

plot_model(modelFitted, type = "pred", terms = "time")











# PANEL PLOTS FOR DIFFERENT PARAMETER COMBINATIONS ------------------

# syntax for figures is L or H for the autocorrelation, D(ecreasing),I(ncreasing),N(Null) for the main effect

nSample = 100 # number of subjects

# between participant variances in intercept and slope,
# and covariance between the two

# ((varIntercept, covar),
# (coVar,varSlope))
SigmaModel = rbind(c(5, 0.2),
                   c(0.2, 0.05))


nObs   = 20 # number of repeated measures


changePoints = as.integer(2) # should be 'None if there is None'
changePointsVar = 1
slopeChangeMean = c(-1, 1)
slopeChangeVar = 0.05


rangeMin = 1
rangeMax = 10

numParticipants = 20


# low autocorr, increasing trend


popInt = 5
popSlope = 1.5
noiseVar = 6

ppVar = -5
ppVarChangePoints = 17

list[data, populationChangePoints, populationPPVarChangePoints] <-
  simData(
    nSample,
    popInt,
    popSlope,
    sigmaModel,
    nObs,
    changePoints,
    changePointsVar,
    slopeChangeMean,
    slopeChangeVar,
    noiseVar,
    rangeMin,
    rangeMax,
    missingness,
    ppVar,
    ppVarChangePoints
  )

LI <-
  plotSimData(data,numParticipants,nSample, nObs, 
              populationChangePoints,
              populationPPVarChangePoints)


# high autocorr, increasing trend


popInt = 5
popSlope = 1.5
noiseVar = 1

ppVar = 5
ppVarChangePoints = 17

list[data, populationChangePoints, populationPPVarChangePoints] <-
  simData(
    nSample,
    popInt,
    popSlope,
    sigmaModel,
    nObs,
    changePoints,
    changePointsVar,
    slopeChangeMean,
    slopeChangeVar,
    noiseVar,
    rangeMin,
    rangeMax,
    missingness,
    ppVar,
    ppVarChangePoints
  )

HI <-
  plotSimData(data,numParticipants, nSample, nObs,
              populationChangePoints,
              populationPPVarChangePoints)


# low autocorr, decreasing trend


popInt = 5
popSlope = -1.5
noiseVar = 6

ppVar = -5
ppVarChangePoints = 17

list[data, populationChangePoints, populationPPVarChangePoints] <-
  simData(
    nSample,
    popInt,
    popSlope,
    sigmaModel,
    nObs,
    changePoints,
    changePointsVar,
    slopeChangeMean,
    slopeChangeVar,
    noiseVar,
    rangeMin,
    rangeMax,
    missingness,
    ppVar,
    ppVarChangePoints
  )

LD <-
  plotSimData(data,numParticipants, nSample, nObs, 
              populationChangePoints,
              populationPPVarChangePoints)


# high autocorr, Dncreasing trend


popInt = 5
popSlope = -1.5
noiseVar = 1

ppVar = 5
ppVarChangePoints = 17

list[data, populationChangePoints, populationPPVarChangePoints] <-
  simData(
    nSample,
    popInt,
    popSlope,
    sigmaModel,
    nObs,
    changePoints,
    changePointsVar,
    slopeChangeMean,
    slopeChangeVar,
    noiseVar,
    rangeMin,
    rangeMax,
    missingness,
    ppVar,
    ppVarChangePoints
  )

HD <-
  plotSimData(data,numParticipants,nSample, nObs, 
              populationChangePoints,
              populationPPVarChangePoints)


# low autocorr, no trend


popInt = 5
popSlope = 0
noiseVar = 6

ppVar = -5
ppVarChangePoints = 17

list[data, populationChangePoints, populationPPVarChangePoints] <-
  simData(
    nSample,
    popInt,
    popSlope,
    sigmaModel,
    nObs,
    changePoints,
    changePointsVar,
    slopeChangeMean,
    slopeChangeVar,
    noiseVar,
    rangeMin,
    rangeMax,
    missingness,
    ppVar,
    ppVarChangePoints
  )

LN <-
  plotSimData(data,numParticipants,nSample, nObs, 
              populationChangePoints,
              populationPPVarChangePoints)


# high autocorr, no trend


popInt = 5
popSlope = 0
noiseVar = 1

ppVar = 5
ppVarChangePoints = 17

list[data, populationChangePoints, populationPPVarChangePoints] <-
  simData(
    nSample,
    popInt,
    popSlope,
    sigmaModel,
    nObs,
    changePoints,
    changePointsVar,
    slopeChangeMean,
    slopeChangeVar,
    noiseVar,
    rangeMin,
    rangeMax,
    missingness,
    ppVar,
    ppVarChangePoints
  )

HN <-
  plotSimData(data,numParticipants,nSample, nObs, 
              populationChangePoints,
              populationPPVarChangePoints)





figure <- ggarrange(
  LI,
  LD,
  LN,
  HI,
  HD,
  HN,
  labels = c("LI", "LD", "LN", "HI", "HD", "HN"),
  ncol = 3,
  nrow = 2,
  font.label = list(color = "black", size = 10)
)

print(figure)

# these settings look good just need to change the path
ggsave(
  device = 'pdf',
  dpi = 1200 ,
  width = 12,
  height = 8,
  filename = "figureTest.pdf"
)



