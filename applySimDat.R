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

# load my own R functions

source('simRMData.R')

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



# generate data -----------------------------------------------------------


nSample = 100 # number of subjects

# between participant variances in intercept and slope, and covariance between the two
SigmaModel = rbind(c(1, 0.2),
                   c(0.2, 0.05))

popInt = 2
popSlope = 1.5

nObs   = 20 # number of repeated measures

changePoints = 10 # should be 'None if there is None'
changePointsVar = 1
slopeChangeMean = c(-1)
slopeChangeVar = 0.1

noiseVar = 1

rangeMin = 1
rangeMax = 10

missingness = 'None'


varFun <- function(x) {
  (x ** 2)
}

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
    c(2, -3, 2),
    c(5, 10, 15)
  )

head(data, 50)

print(populationChangePoints)
print(populationPPVarChangePoints)
# plot raw data and descriptives ------------------------------------------


numParticipants = 30
p <-
  plotSimData(data,
              numParticipants,
              populationChangePoints,
              populationPPVarChangePoints)

print(p)

# these settings look good just need to change the path
ggsave(
  device = 'pdf',
  dpi = 600 ,
  width = 12,
  height = 8,
  filename = "figureStatsEnvChangePointsInVariance.pdf"
)

# make a panel plot for different parameter combinations ------------------

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
  plotSimData(data,
              numParticipants,
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
  plotSimData(data,
              numParticipants,
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
  plotSimData(data,
              numParticipants,
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
  plotSimData(data,
              numParticipants,
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
  plotSimData(data,
              numParticipants,
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
  plotSimData(data,
              numParticipants,
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


# fit model ---------------------------------------------------------------


# fitting a mixed effects model

modelFitted <- lme4::lmer(y ~ time + (time | ppID), data)

# let's look at the mode results

summary(modelFitted)
sjPlot::tab_model(modelFitted)

plot_model(modelFitted, type = "pred", terms = "time")