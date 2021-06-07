# load packages -----------------------------------------------------------
if (!require("pacman"))
  install.packages("pacman")

pacman::p_load(
  ggplot2,
  ggthemes,
  dplyr,
  reshape,
  scales,
  Hmisc,
  forecast,
  nlme,
  pdp,
  piecewiseSEM,
  lmtest,
  changepoint,
  changepoint.np,
  tidyverse,
  RColorBrewer,
  hrbrthemes,
  broom
)


# set working directory

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()

# load other functions; this only works if the sourced file is
# in the same directory as this R script
source("simRMData.R")
source("explore.R")
source("preprocess.R")
source("framework.R")


# loading the different simulated data sets --------------------------------------------------------


# 30 repeated measures; increasing variance; switch points in slope
data30RM_IncreasingVar <-
  read.csv("for analysis//simData_VarIncreasingandSlopeChangePoints.csv")

head(data30RM_IncreasingVar)


# 30 repeated measures; more extreme increases variance; switch points in slope
data30RM_IncreasingVarExtr <-
  read.csv("for analysis//simData_VarIncreasingExtremeandSlopeChangePoints_30RM.csv")

head(data30RM_IncreasingVarExtr)

# 30 repeated measures; switch points in variance; switch points in slope
data30RM_VarChangePoints <-
  read.csv("for analysis//simData_VarChangePointsandSlopeChangePoints_30RM.csv")

head(data30RM_VarChangePoints)


# 50 repeated measures; increasing variance; switch points in slope
data50RM_IncreasingVar <-
  read.csv("for analysis//simData_VarIncreasingandSlopeChangePoints_50RM.csv")

head(data50RM_IncreasingVar)

# 50 repeated measures; more extreme increasing variance; switch points in slope
data50RM_IncreasingVarExtr <-
  read.csv("for analysis//simData_VarIncreasingExtremeandSlopeChangePoints_50RM.csv")

head(data50RM_IncreasingVarExtr)

# 50 repeated measures; switch points in variance; switch points in slope
data50RM_VarChangePoints <-
  read.csv("for analysis//simData_VarChangePointsandSlopeChangePoints_50RM.csv")

head(data50RM_VarChangePoints)

# 50 repeated measures; switch points in variance (more extreme); switch points in slope
data50RM_VarChangePointsExtr <-
  read.csv("for analysis//simData_VarChangePointsandSlopeChangePoints_50RM_Extreme.csv")

head(data50RM_VarChangePointsExtr)



# SPECIFY DATASET ---------------------------------------------------------

nSample = 100
dataSim <- data50RM_VarChangePointsExtr #data50RM_IncreasingVarExtr
nObs <- 51
numParticipants <- 30 # number of participants to show in the plot

head(dataSim)
#only keep the variables we need 
data <- data.frame(dataSim$ppID, dataSim$time, dataSim$yNorm)
colnames(data) <- c("ppID", "time", "yNorm")
head(data)

# user needs to specify their variable names
data_ppID <- "ppID"
data_x <- "time"
data_y <- "yNorm"

# assign the variable names that are neceesary for the subsequent functions
names(data)[names(data) == data_ppID] <- "ppID"
names(data)[names(data) == data_x] <- "time"
names(data)[names(data) == data_y] <- "y"


# EXPLORE -----------------------------------------------------------------

# PART ONE
# plot the original data
# need to have a list slope and variance changepoints for plotting

head(data)

list[populationChangePoints, populationPPVarChangePoints] <-
  detect_ChangePoints(dataSim)

yLabel = "morbidity-mortality"

pOrg <-
  explore.plotTimeSeries(
    data,
    numParticipants,
    nSample,
    nObs,
    populationChangePoints,
    populationPPVarChangePoints, yLabel
  )

pOrg


# if the data are not simulated with my function, both changepoint arguments
# should be set to "None"
# alternatively the user can specify the location of the change points if they like

pOrg <-
  explore.plotTimeSeries(data,
              numParticipants,
              nSample,
              nObs,
              "None",
              "None", yLabel)

pOrg

# PART TWO
# time series decomposition
# create one time series per participant 

# what is the number of observations per time unit
# this information is necessary to remove/quantify seasonal effects
# it could be for example x observations per day, week, or year. only the x
# needs to be specified and not the unit
# the unit matters for the interpretation by the user
# freq needs to be at least 2 to be able to estimate seasonal effects
freq = 5
freqUnit = "week"

# If the seasonality and residual components are independent of the trend,
# then you have an additive series. If the seasonality and residual components
# are in fact dependent, meaning they fluctuate on trend, then you have a multiplicative series

type <- 'additive'
numParticipantsDecompose <- 5 #how many lines to plot

# this plots the different components of a time series for the
# specified number of participants
# TS_decompose contains one time series object per participant 
list[pTS, dataTS] <-
  explore.decompose_TS(data , numParticipantsDecompose, type, freq, freqUnit)
pTS


# PART 3
# autocorrelation and partial autocorrelation at different lags

# partial autocorrelation displays the autocorrelation of each lag after controlling for all preceding lags
# the PACF is used to determine the number of autoregressive terms to include in an ARIMA model
# in this case we would include one autoregressive term

lagMax = 20

# create an ACF matrix for all participants with a specified number of lags

list[acPlot,pacPlot] <- explore.visualize_autocorr(data = data, lagMax, nSample)



# explore change points
p_cp <- explore.changePoints(data = data)
p_cp

# PREPROCESS --------------------------------------------------------------

# various options are offered to the user

#1. seasonally adjust the time series

type <- 'additive'
data_SeasonalaAdj <- preprocess.adjust_season(dataTS, type)
#inspect the adjusted data
list[pTS, dataTSNew] <-
  explore.decompose_TS(data_SeasonalaAdj , numParticipantsDecompose, type, freq,freqUnit)
pTS

#2. difference the time series to achieve stationarity in mean or apply 
#   a logarithmic transformation to achieve stationarity in variance

# difference once
data_diffAdj <- preprocess.difference(data, degree =1, FALSE)
head(data_diffAdj)
list[pTS,dataTSNew] <-
  explore.decompose_TS(data_diffAdj , numParticipantsDecompose, type, freq,freqUnit)
pTS


# difference twice
data_diffAdj <- preprocess.difference(data, degree =2, FALSE)
head(data_diffAdj)
list[pTS,dataTSNew] <-
  explore.decompose_TS(data_diffAdj , numParticipantsDecompose, type, freq,freqUnit)
pTS

# logarithmic transformation and no differencing
data_diffAdj <- preprocess.difference(data, degree =0, TRUE)
head(data_diffAdj)
list[pTS,dataTSNew] <-
  explore.decompose_TS(data_diffAdj , numParticipantsDecompose, type, freq,freqUnit)
pTS

# logarithmic transformation and one round of differencing
data_diffAdj <- preprocess.difference(data, degree =1, TRUE)
head(data_diffAdj)
list[pTS,dataTSNew] <-
  explore.decompose_TS(data_diffAdj , numParticipantsDecompose, type, freq, freqUnit)
pTS
# add seasonal adjustment
data_SeasonalaAdj <- preprocess.adjust_season(dataTSNew, type)
#inspect the adjusted data
list[pTS, xx] <-
  explore.decompose_TS(data_SeasonalaAdj , numParticipantsDecompose, type, freq, freqUnit)
pTS

# split the data set into different subsets based on the specified time indices  

splitValues <- c(18,36)
dataSplitted <- preprocess.splitData(data = data,splitValues)
# these individual data sets may be used as input to the framework function
dataEarly <- dataSplitted[[1]]
dataMiddle <- dataSplitted[[2]]
dataLate <- dataSplitted[[3]]



# FRAMEWORK ---------------------------------------------------------------

# if interested in polynomial effects 
data$time_c <- scale(data$time, center = TRUE, scale = FALSE)
data$time2 <- data$time_c^2
data$time3 <- data$time_c^3


meanModel <- function(data){
  #specify the model of the mean here 
  lm(y ~ time, data = data)
}

varModel <- function(data){
  #specify the model of the variance here 
  lm(yVar ~ time_c +time2, data = data)
}

acLAG = 4
freq = 5
freqUnit = "week"
yLabel = "morbidity-mortality"

modelResults <-
  framework.applyModel(data, acLAG, meanModel, varModel)


head(modelResults)
view(modelResults)

# to request the "statistics of the environment"
modelResults$Statistics_Mean
modelResults$Statistics_Variance
modelResults$Statistics_ChangePoints

# plotting the models of the mean and variance
pModelResults <- framework.plotModel(modelResults, 2, "mean", "mean", freq, freqUnit, yLabel)

pModelResults <- framework.plotModel(modelResults, 1, "variance", "variance", freq, freqUnit, yLabel)

#plot color of noise distribution 

pColorOfNoise <- framework.plotColorOfNoise(modelResults, binwidth = 0.05)
pColorOfNoise

