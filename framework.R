# load packages -----------------------------------------------------------


if (!require("pacman"))
  install.packages("pacman")

pacman::p_load(
  ggplot2,
  ggthemes,
  lme4,
  sjPlot,
  effects,
  dplyr,
  reshape,
  scales,
  Hmisc,
  forecast,
  nlme,
  pdp
)

# load my own R functions
simRMDataPath <- file.path(
  "C:",
  "Users",
  "nicwa",
  "Dropbox",
  "PhD",
  "Publications",
  "statistics of the environment",
  "statistics_of_the_environment_code",
  "simRMData.R"
)
source(simRMDataPath)

# helper functions

# detect the change points
# has to be a loop across participants


detect_Slope_And_Varaince_ChangePoints <- function(data) {
  # detect slope change points
  outcomeSlope <-
    data.frame(sapply(split(data$RE.s, data$ppID), function(vec) {
      return(which(abs(diff(vec)) > 0))
    }))
  
  # next average across participants
  resultSlope <- round(rowMeans(outcomeSlope))
  
  if (length(resultSlope) == (nrow(data[data$ppID == 1,]) - 1)) {
    resultSlope = 'None'
  } else if ((is.integer(resultSlope) ||
              is.numeric(resultSlope)) &&
             length(resultSlope) == 0L) {
    resultSlope = 'None'
  }
  
  
  
  # detect variance change points
  outcomeVar <-
    data.frame(sapply(split(data$RE.var, data$ppID), function(vec) {
      return(which(abs(diff(vec)) > 0))
    }))
  
  # next average across participants
  resultVar <- round(rowMeans(outcomeVar))
  
  
  if (length(resultVar) == (nrow(data[data$ppID == 1,]) - 1)) {
    resultVar = 'None'
  } else if ((is.integer(resultVar) ||
              is.numeric(resultVar)) && length(resultVar) == 0L) {
    resultVar = 'None'
  }
  
  
  return(
    list(
      "populationChangePoints" = resultSlope,
      "populationPPVarChangePoints" = resultVar
    )
  )
  
}


TS_decompose <- function(data_TS, numParticipants, nSample, type){
  
  data_TS_subet <- list()
  counter = 1
  for (ppID in sample(c(1:nSample), numParticipants)) {
    data_TS_subet[counter] = data_TS[ppID]
    counter <- counter + 1
  }
  
  data_TS_decomposed <- sapply(data_TS_subet, function(x) {decompose(x, type = type)})
  
  data_TS_decomposed_x <-
    data.frame((sapply(data_TS_decomposed["x",], c)))
  data_TS_decomposed_seasonal <-
    data.frame((sapply(data_TS_decomposed["seasonal",], c)))
  data_TS_decomposed_trend <-
    data.frame((sapply(data_TS_decomposed["trend",], c)))
  data_TS_decomposed_random <-
    data.frame((sapply(data_TS_decomposed["random",], c)))
  
  # next tarnsofrm each into long format and then merge columns
  # reshape into long format
  
  data_TS_decomposed_x_Long <-
    reshape(
      data_TS_decomposed_x,
      direction = 'long',
      varying = list(1:numParticipants),
      v.names = "originalData",
      timevar = "ppID",
      idvar = "time"
    )
  
  
  # next merge the dataframes together
  data_TS_decomposed_DF <-
    cbind(
      data_TS_decomposed_x_Long,
      as.numeric(unlist(data_TS_decomposed_seasonal)),
      as.numeric(unlist(data_TS_decomposed_trend)),
      as.numeric(unlist(data_TS_decomposed_random))
    )
  
  names(data_TS_decomposed_DF) <-
    c("ppID", "originalData", "time", "seasonal", "trend", "random")
  data_TS_decomposed_DF$ppID <- as.factor(data_TS_decomposed_DF$ppID)
  head(data_TS_decomposed_DF) #looks good
  
  # reshape into long format
  data_TS_decomposed_DF_Long <-
    reshape(
      data_TS_decomposed_DF,
      direction = 'long',
      varying = c("originalData", "seasonal", "trend", "random"),
      v.names = "value",
      timevar = 'dataType',
    )
  
  data_TS_decomposed_DF_Long$dataType <-
    as.factor(data_TS_decomposed_DF_Long$dataType)
  levels(data_TS_decomposed_DF_Long$dataType) <-
    c("originalData", "seasonal", "trend", "random")
  
  
  p <-
    ggplot(data = data_TS_decomposed_DF_Long, aes(x = time, y = value, group = ppID)) +  geom_line(color = "grey", size = (1/numParticipants)) + geom_point(size = (1.5/numParticipants)) + facet_wrap(
      ~
        dataType,
      ncol = 1,
      strip.position = "left",
      scales = "free_y"
    )
  p <-
    p + labs(y = "value", x = "week") + scale_x_continuous(breaks =
                                                             seq(0, nObs, 5), labels = c(0:(nObs / 5)))
  p <- p +       theme(legend.position = "none") + theme_bw() + theme(
    plot.title = element_text(face = "bold", size = 8) ,
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(size = 1, colour = 'grey70'),
    axis.line.y = element_line(size = 1, colour = 'grey70'),
    legend.position = c(0.06, 0.75)
  )
  
  return(list("p" = p,"data_TS_decomposed" = data_TS_decomposed))
  
  
}







# loading the data --------------------------------------------------------
# set working directory

wdPath = file.path(
  "C:",
  "Users",
  "nicwa",
  "Dropbox",
  "PhD",
  "Publications",
  "statistics of the environment",
  "statistics_of_the_environment_code",
  "for analysis"
)
setwd(wdPath)


# 30 repeated measures; increasing variance; switch points in slope
data30RM_IncreasingVar <-
  read.csv("simData_VarIncreasingandSlopeChangePoints.csv")

head(data30RM_IncreasingVar)


# 30 repeated measures; more extreme increases variance; switch points in slope
data30RM_IncreasingVarExtr <-
  read.csv("simData_VarIncreasingExtremeandSlopeChangePoints_30RM.csv")

head(data30RM_IncreasingVarExtr)

# 30 repeated measures; switch points in variance; switch points in slope
data30RM_VarChangePoints <-
  read.csv("simData_VarChangePointsandSlopeChangePoints_30RM.csv")

head(data30RM_VarChangePoints)


# 50 repeated measures; increasing variance; switch points in slope
data50RM_IncreasingVar <-
  read.csv("simData_VarIncreasingandSlopeChangePoints_50RM.csv")

head(data50RM_IncreasingVar)

# 50 repeated measures; increasing variance; switch points in slope
data50RM_IncreasingVarExtr <-
  read.csv("simData_VarIncreasingExtremeandSlopeChangePoints_50RM.csv")

head(data50RM_IncreasingVarExtr)

# 50 repeated measures; switch points in variance; switch points in slope
data50RM_VarChangePoints <-
  read.csv("simData_VarChangePointsandSlopeChangePoints_50RM.csv")

head(data50RM_VarChangePoints)


# Jeb et al. 2015 ----------------------------------------------------

# I am consciously not explicitly showing example code for seasonal effects
# I don't think there are many theoretical justfications for the questions we ask
# something like yearly flu season might be interesting; however, this would require us
# to have multiple years of data each ideally with at 12 observations (1 per month)
# I think it is sufficient for this purpose to refer the reader to the jebb paper and other resources

# additionally I don't think that we want to build an arima model per paerticipant as
# we are not interested in forecasting but rather descriptive modeling (see table in the jebb et al paper)


### preprocessing steps ###

# specify here which data set you would like to use
nSample = 100
data <- data50RM_IncreasingVarExtr
nObs <- 50
numParticipants <- 30 # number of participants shown in the plot



# plot the original data
# need to have a list slope and variance changepoints for plotting
list[populationChangePoints, populationPPVarChangePoints] <-
  detect_Slope_And_Varaince_ChangePoints(data)


p <-
  plotSimData(
    data,
    numParticipants,
    nSample,
    nObs,
    populationChangePoints,
    populationPPVarChangePoints
  )

print(p)



# this will create one time series object per participant
data_TS <- list()

# assuming data have been collected across 6 weeks with 5 measurements per week
# frequency specifies the number of observations per unit of time
# this is necessary for identifying cycles and seasons
for (idx in 1:nSample) {
  data_TS[[idx]] <-
    ts(data[data$ppID == idx,]$yNorm[-1], frequency = 5)
}


### time series decomposition ###

# used to identify seasonality/cycles, trend and irregular components

# additive: most appropriate when the magnitude of the trend-cycle and seasonal
# components remain constant across time

# multiplicative: most appropriate when magnitude varies but still appears
# proportional over time (change by a multiplicative factor)

# If the seasonality and residual components are independent of the trend,
# then you have an additive series. If the seasonality and residual components
# are in fact dependent, meaning they fluctuate on trend, then you have a multiplicative series.


# create a function for this

numParticipantsDecompose <- 1
nSample <- 100
type <- 'additive'

list[p,TS_decompose_subset] <- TS_decompose(data_TS , numParticipantsDecompose, nSample,type)

print(p)

currPP <-
  sample(c(1:numParticipantsDecompose), 1) # specifies a random participant; use could also specify a participant



tsdecomp <- TS_decompose_subset[,currPP]
class(tsdecomp) <- "decomposed.ts"

### seasonal adjustment of the time series ###

# as variance is dependent on time we use a multiplicative decomposition
tsadj = seasadj(tsdecomp)
combined <- data.frame(as.numeric(data_EXAMPLE), tsadj)
names(combined) <- c("originalData", "adjustedData")
combined$time <- c(1:nObs)

# reshape into long format
combinedLong <-
  reshape(
    combined,
    direction = 'long',
    varying = c("originalData", "adjustedData"),
    v.names = "morbidity_mortality",
    timevar = 'dataType',
    idvar = "time"
  )
combinedLong$dataType <- as.factor(combinedLong$dataType)
levels(combinedLong$dataType) <- c('original data', 'adjusted data')

p <-
  ggplot(data = combinedLong, aes(x = time, y = morbidity_mortality)) + geom_line() + geom_point() + facet_wrap( ~
                                                                                                                   dataType, ncol = 1, strip.position = "left")
p <-
  p + labs(y = "morbidity-mortality", x = "week") + scale_x_continuous(breaks =
                                                                         seq(0, nObs, 5), labels = c(0:(nObs / 5)))
p <- p +       theme(legend.position = "none") + theme_bw() + theme(
  plot.title = element_text(face = "bold", size = 8) ,
  axis.ticks = element_line(colour = "grey70", size = 0.2),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line.x = element_line(size = 1, colour = 'grey70'),
  axis.line.y = element_line(size = 1, colour = 'grey70'),
  legend.position = c(0.06, 0.75)
)
print(p)


### stationarity ###
# I don't think we want the data to be stationary as we are explicitly interested in the trend

### autocorrelation ###
# plotting the ACF and PACF
# partial autocorrelation displays the autocorrelation of each lag after controlling for all preceding lags
# the PACF is used to determine the number of autoregressive terms to inlcude in an ARIMA model
# in this case we would include one autoregressive term

plot(acf(as.numeric(tsadj), plot = FALSE), main = 'ACF of the seasonally adjusted time series')
plot(pacf(as.numeric(tsadj), plot = FALSE), main = 'PACF of the seasonally adjusted time series')

### modeling the mean trend ###

# because we are not interested in seasonal effects, we use the seasonally adjusted data to model the trend
# we use generalized least squares to account for the autocorrelation in the data

time = c(1:nObs)
timeModel <- gls(tsadj ~ time)
summary(timeModel)
fittedTimeModel <- fitted(timeModel)
residualsTimeModel <- as.numeric(residuals(timeModel))
absResidualsTimeModel <- abs(residualsTimeModel)
coefTimeModel <-
  coef(timeModel) #this could be the descriptive population intercept and slope; is this better than a mixed effects model?
originalData <- as.numeric(tsadj)

results <-
  data.frame(originalData,
             fittedTimeModel,
             residualsTimeModel,
             absResidualsTimeModel,
             time)


# residual time model to assess varaince across time
timeModelRes <- gls(absResidualsTimeModel ~ time)
summary(timeModelRes)
fittedTimeModelRes <- fitted(timeModelRes)
coefTimeModelRes <- coef(timeModelRes)

# plotting the fitted line and residuals

pModel <-
  ggplot(data = results) + geom_line(aes(x = time, y = originalData), color = "grey") + geom_point(aes(x =
                                                                                         time, y = originalData)) + geom_line(aes(x = time, y = fittedTimeModel), color = "red")
pModel <-
  pModel + labs(y = "morbidity-mortality", x = "week") + scale_x_continuous(breaks =
                                                                              seq(0, nObs, 5), labels = c(0:(nObs / 5)))
pModel <- pModel + labs(title = "linear model")
pModel <-
  pModel +       theme(legend.position = "none") + theme_bw() + theme(
    plot.title = element_text(face = "bold", size = 8) ,
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(size = 1, colour = 'grey70'),
    axis.line.y = element_line(size = 1, colour = 'grey70'),
    legend.position = c(0.06, 0.75)
  )


pRes <-
  ggplot(data = results, aes(x = time, y = residualsTimeModel)) + geom_line(color = "grey") + geom_point() 
pRes <-
  pRes + labs(y = "residual", x = "week") + scale_x_continuous(breaks =
                                                                 seq(0, nObs, 5), labels = c(0:(nObs / 5)))
pRes <- pRes + labs(title = "residual error series")
pRes <-
  pRes +       theme(legend.position = "none") + theme_bw() + theme(
    plot.title = element_text(face = "bold", size = 8) ,
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(size = 1, colour = 'grey70'),
    axis.line.y = element_line(size = 1, colour = 'grey70'),
    legend.position = c(0.06, 0.75)
  )


pAbsRes <-
  ggplot(data = results) + geom_line(aes(x = time, y = absResidualsTimeModel), color = "grey") + geom_point(aes(x = time, y = absResidualsTimeModel)) + geom_line(aes(x = time, y = fittedTimeModelRes), color = "red")

pAbsRes <-
  pAbsRes + labs(y = "residual", x = "week") + scale_x_continuous(breaks =
                                                                    seq(0, nObs, 5), labels = c(0:(nObs / 5)))
pAbsRes <- pAbsRes + labs(title = "absolute residual error series")
pAbsRes <-
  pAbsRes +       theme(legend.position = "none") + theme_bw() + theme(
    plot.title = element_text(face = "bold", size = 8) ,
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(size = 1, colour = 'grey70'),
    axis.line.y = element_line(size = 1, colour = 'grey70'),
    legend.position = c(0.06, 0.75)
  )

# combining the plots

pModelEval <- grid.arrange(pModel, pRes, pAbsRes, nrow = 1)

# TODO: create a data set containing per pp statistics
# TODO search for change point analysis
# TODO read the jebb 2017 paper
# TODO continue with Gerstorf
# create framework functions once we (Ethan, Willen and I) agree on all the components








### framework component ###
# a.	From spreadsheet / papers:
#   i.	Trend in mean of morbidity-mortality (mean, slope, intercept) (done; is better to do this via mixed effects)
# ii.	Autocorrelation at different lags and perhaps at different life-stages (done)
# 1.	Color of noise (Willem would like to have that) (still to do)
# iii.	Identification of change points in mean and variance (still to do)
# iv.	Trend in variance of morbidity-mortality (mean, slope, intercept) (done)
# b.	Extra:
#   i.	Identify seasonality and cycles (in principle done; although I don't see the theoretical justification)




# fitting a mixed model [IGNORE THIS FOR NOW] ---------------------------------------------------


# fitting a mixed effects model


modelFitted <- lme4::lmer(y ~ time + (time | ppID), data)

# let's look at the mode results

summary(modelFitted)
sjPlot::tab_model(modelFitted)

plot_model(modelFitted, type = "pred", terms = "time")


