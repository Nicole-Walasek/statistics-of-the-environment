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
  pdp, piecewiseSEM,lmtest, changepoint, changepoint.np
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
  
  data_TS_subet <- data_TS
  # counter = 1
  # for (ppID in sample(c(1:nSample), numParticipants)) {
  #   data_TS_subet[counter] = data_TS[ppID]
  #   counter <- counter + 1
  # }
  
  data_TS_decomposed <- sapply(data_TS_subet, function(x) {decompose(x, type = type)})
  
  data_TS_decomposed_x <-
    data.frame((sapply(data_TS_decomposed["x",], c)))
  data_TS_decomposed_seasonal <-
    data.frame((sapply(data_TS_decomposed["seasonal",], c)))
  data_TS_decomposed_trend <-
    data.frame((sapply(data_TS_decomposed["trend",], c)))
  data_TS_decomposed_random <-
    data.frame((sapply(data_TS_decomposed["random",], c)))
  

  # perform seasonal adjustment 
  data_TS_decomposed_seasonalAdj <- list()
  if(type == 'additive'){
    data_TS_decomposed_seasonalAdj <- data_TS_decomposed_x - data_TS_decomposed_seasonal   
  
    } else{
    
      data_TS_decomposed_seasonalAdj <- data_TS_decomposed_x/data_TS_decomposed_seasonal   
  }

  
  # next transform each into long format and then merge columns
  # reshape into long format
  
  data_TS_decomposed_x_Long <-
    reshape(
      data_TS_decomposed_x,
      direction = 'long',
      varying = list(1:nSample),
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
      as.numeric(unlist(data_TS_decomposed_random)),
      as.numeric(unlist(data_TS_decomposed_seasonalAdj))
    )
  
  names(data_TS_decomposed_DF) <-
    c("ppID", "originalData", "time", "seasonal", "trend", "random","adjusted")
  data_TS_decomposed_DF$ppID <- as.factor(data_TS_decomposed_DF$ppID)
  head(data_TS_decomposed_DF) #looks good
  
  # reshape into long format
  data_TS_decomposed_DF_Long <-
    reshape(
      data_TS_decomposed_DF,
      direction = 'long',
      varying = c("originalData", "seasonal", "trend", "random","adjusted"),
      v.names = "value",
      timevar = 'dataType',
    )
  
  data_TS_decomposed_DF_Long$dataType <-
    as.factor(data_TS_decomposed_DF_Long$dataType)
  levels(data_TS_decomposed_DF_Long$dataType) <-
    c("originalData", "seasonal", "trend", "random", "adjusted")
  
  scaleFactor <- 1
  
  if(numParticipants > 1) {
    scaleFactor <- numParticipants/2
  } 
  
  
  data_TS_decomposed_DF_LongSubset <- filter(data_TS_decomposed_DF_Long, ppID %in% sample(data_TS_decomposed_DF_Long$ppID, numParticipants)) #select random participants
  p <-
    ggplot(data = data_TS_decomposed_DF_LongSubset, aes(x = time, y = value, group = ppID)) +  geom_line(color = "grey", size = (1/scaleFactor)) + geom_point(size = (1.5/scaleFactor)) + facet_wrap(
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
  
  # for convenience, we will transpose the matrix before returning 
  data_TS_decomposed_seasonalAdj <- as.data.frame(t(as.matrix(data_TS_decomposed_seasonalAdj)))
  
  return(list("p" = p,"data_TS_decomposed" = data_TS_decomposed, "data_TS_decomposed_seasonalAdj" =data_TS_decomposed_seasonalAdj))
  
  
}


# framework functions

meanModel <- function(data, nObs ,degree =1,type){
  # data: seasonally adjusted time series data. assumes one row per pp 
  # degree: degree of time predictor
  # 1 for linear, 2 for quadratic, 3 for polynomial
  
  time <- c(1:nObs)
  time_c <- time - mean(time)
  time2 <- time_c^2
  time3 <- time_c^3
  
  alpha <- 0.05 # to check whether residuals are correlated; could be changed to an absolute criterion  
  
  # framework components that we will compute and store
  interceptRes <- c()
  time1Res <- c()
  time2Res <- c()
  time3Res <- c()
  rsquareRes <- c()
  meanRes <- c()
  stdRes <- c()
  acRes <- c()
  
  
  # for plotting DF
  fittedAll <- c()
  residualsAll <- c()
  
  #counter
  current <-1
  # iterate through participants
  for(currentPP in unique(data$ppID)){
    
    originalData <- c()
    if(type == 'mean'){
      originalData <- with(data, data[ppID == currentPP,]$y)
    }else{
      originalData <- with(data, data[ppID == currentPP,]$yVar)
    }
    
    # not sure how to specify the best ARMA model for the residuals; the tutorial is not helpful
    # there also doesn't seem to be a way to specify differencing for the residuals 
    # @Ethan
    
    
    if(degree == 1) {
      timeModel <- lm(originalData ~ time) 
      
      if(dwtest(timeModel)$p.value < alpha){ # in this case there is autocorrelation present and will model this autocorrelation using GLS
        timeModel <- gls(originalData ~ time, correlation = corARMA(p=1, q=0, form = ~ 1)) 
      }
      
    }else if(degree == 2){
      timeModel <- lm(originalData ~ time_c+time2) 
      if(dwtest(timeModel)$p.value < alpha){ 
        timeModel <- gls(originalData ~ time_c+time2, correlation = corARMA(p=1, q=0, form = ~ 1)) 
      } 
      
    }else{ # degree corresponds to 3
      # first model without gls 
      timeModel <- lm(originalData ~ time_c+time2+time3) 
      if(dwtest(timeModel)$p.value < alpha){ 
        timeModel <- gls(originalData ~ time_c+time2+time3, correlation = corARMA(p=1, q=0, form = ~ 1)) 
      }
      
    }
    
    # for the framework components 
    
    #this could be the descriptive population intercept and slope; is this better than a mixed effects model
    coefTimeModel <-coef(timeModel) 
    interceptRes[current] <- as.numeric(coefTimeModel[1])
    
    
    if (degree == 1){
      time1Res[current] <- as.numeric(coefTimeModel[2])  
    } else if(degree == 2) {
      time1Res[current] <- as.numeric(coefTimeModel[2])
      time2Res[current] <- as.numeric(coefTimeModel[3])
    } else if(degree == 3){
      time1Res[current] <- as.numeric(coefTimeModel[2])
      time2Res[current] <- as.numeric(coefTimeModel[3])
      time3Res[current] <- as.numeric(coefTimeModel[4])
    }
    
    # crude r-squared estimate 
    rsquareRes[current] <- cor(originalData, fitted(timeModel))^2 
    meanRes[current] <- mean(originalData)
    stdRes[current] <- sd(originalData)
    acRes[current] <- pacf(as.numeric(originalData), plot = FALSE)$acf[1]
    
    
    # for plotting
    fittedAll <- c(fittedAll,as.numeric(fitted(timeModel)))
    residualsAll <- c(residualsAll,as.numeric(residuals(timeModel)))
    
    #advance the counter 
    current <- current + 1
  }
  meanResults <- data.frame()
  
  if(degree ==1){
    meanResults <- data.frame(unique(data$ppID),interceptRes, time1Res,rsquareRes, meanRes, stdRes, acRes)
    names(meanResults) <-  c("ppID", "intercept", "time1", "rsquare", "mean", "sd", "ac")
    
  }else if(degree  ==2){
    meanResults <- data.frame(unique(data$ppID),interceptRes, time1Res,time2Res,rsquareRes, meanRes, stdRes, acRes)  
    names(meanResults) <-  c("ppID", "intercept", "time1", "time2","rsquare", "mean", "sd", "ac")
    
  }else{
    meanResults <- data.frame(unique(data$ppID),interceptRes, time1Res,time2Res,time3Res,rsquareRes, meanRes, stdRes, acRes)
    names(meanResults) <-  c("ppID", "intercept", "time1", "time2", "time3", "rsquare", "mean", "sd", "ac")
    
  }
  data$fittedAll <- fittedAll
  data$residualsAll <- residualsAll
  
  return(list(meanResults, data))
  
}


# plotting the fitted line and residuals
plotModel <- function(data, numPP, modelTitle, type){
  
  results <-
    filter(data, ppID %in% sample(data$ppID, numPP)) #select random participants for plotting
  results$ppID <- as.factor(results$ppID)
  names(results)[names(results) == 'fittedAll'] <- 'fittedTimeModel'
  names(results)[names(results) == 'residualsAll'] <- 'residualsTimeModel'
  
  if(type =="mean"){
    names(results)[names(results) == 'y'] <- 'originalData'
    yLab <- "morbidity-mortality"
  }else{
    names(results)[names(results) == 'yVar'] <- 'originalData'
    yLab <- "morbidity-mortality variance"
  }
  
  pModel <-
    ggplot(data = results) + geom_line(aes(x = time, y = originalData, group = ppID), color = "grey") + geom_point(aes(x =
                                                                                                                         time, y = originalData, group = ppID), size = (1.5/numPP)) + geom_line(aes(x = time, y = fittedTimeModel, group = ppID), color = "firebrick")
  pModel <-
    pModel + labs(y =yLab, x = "week") + scale_x_continuous(breaks =
                                                                                seq(0, nObs, 5), labels = c(0:(nObs / 5)))
  pModel <- pModel + labs(title = paste(modelTitle, " model"))
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
    ggplot(data = results, aes(x = time, y = residualsTimeModel, group = ppID)) + geom_line(color = "grey", aes(group = ppID)) + geom_point(aes(group = ppID), size = (1.5/numPP)) 
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
  
  # combining the plots
  pModelEval <- grid.arrange(pModel, pRes, nrow = 1)
  
  return(pModelEval)
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

# 50 repeated measures; more extreme increasing variance; switch points in slope
data50RM_IncreasingVarExtr <-
  read.csv("simData_VarIncreasingExtremeandSlopeChangePoints_50RM.csv")

head(data50RM_IncreasingVarExtr)

# 50 repeated measures; switch points in variance; switch points in slope
data50RM_VarChangePoints <-
  read.csv("simData_VarChangePointsandSlopeChangePoints_50RM.csv")

head(data50RM_VarChangePoints)

# 50 repeated measures; switch points in variance (more extreme); switch points in slope
data50RM_VarChangePointsExtr <-
  read.csv("simData_VarChangePointsandSlopeChangePoints_50RM_Extreme.csv")

head(data50RM_VarChangePointsExtr)



# Jeb et al. 2015 + change point analysis ----------------------------------------------------

# I am consciously not explicitly showing example code for seasonal effects
# I don't think there are many theoretical justfications for the questions we ask
# something like yearly flu season might be interesting; however, this would require us
# to have multiple years of data each ideally with at 12 observations (1 per month)
# I think it is sufficient for this purpose to refer the reader to the jebb paper and other resources

# additionally I don't think that we want to build an arima model per paerticipant as
# we are not interested in forecasting but rather descriptive modeling (see table in the jebb et al paper)


### preprocessing steps + visual inspection of the data ###

# specify here which data set you would like to use
nSample = 100
data <- data50RM_VarChangePointsExtr #data50RM_IncreasingVarExtr
nObs <- 50
numParticipants <- 30 # number of participants to show in the plot



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


### time series decomposition + seasonal adjustment ###

# used to identify seasonality/cycles, trend and irregular components

# additive: most appropriate when the magnitude of the trend-cycle and seasonal
# components remain constant across time

# multiplicative: most appropriate when magnitude varies but still appears
# proportional over time (change by a multiplicative factor)

# If the seasonality and residual components are independent of the trend,
# then you have an additive series. If the seasonality and residual components
# are in fact dependent, meaning they fluctuate on trend, then you have a multiplicative series.
# this is a preprocessing step to identify the individual framework components


numParticipantsDecompose <- 3 #how many lines to plot 
nSample <- 100
type <- 'additive'

list[p,TS_decompose_subset, seasonallyAdjData] <- TS_decompose(data_TS , numParticipantsDecompose, nSample,type)

print(p)

### stationarity ###
# I don't think we want the data to be stationary as we are explicitly interested in the trend

### autocorrelation ###
# plotting the ACF and PACF
# partial autocorrelation displays the autocorrelation of each lag after controlling for all preceding lags
# the PACF is used to determine the number of autoregressive terms to inlcude in an ARIMA model
# in this case we would include one autoregressive term

# one example pp
currPP <-
  sample(c(1:numParticipantsDecompose), 1) # specifies a random participant; use could also specify a participant
tsadj <- seasonallyAdjData[currPP,]


plot(acf(as.numeric(tsadj), plot = FALSE), main = 'ACF of the seasonally adjusted time series')
plot(pacf(as.numeric(tsadj), plot = FALSE), main = 'PACF of the seasonally adjusted time series')

### modeling the mean trend ###

# @Ethan: can we always fit polynomial coefficients to use them as per-individual predictors?  
# because we are not interested in seasonal effects, we use the seasonally adjusted data to model the trend
# we use generalized least squares to account for the autocorrelation in the data
# without specifying a correlation object GLS() is identical to lm() and GLM()


#reshape seasonally adjusted data into long format for analysis and plotting
seasonallyAdjDataLong <- reshape(seasonallyAdjData, idvar = "ppID", varying = list(1:nObs), v.names = 'y', timevar = 'time' ,direction = "long")
seasonallyAdjDataLong <- seasonallyAdjDataLong[
  with(seasonallyAdjDataLong, order(ppID)),]
head(seasonallyAdjDataLong,55)



# run mean-trend analysis for whole dataset
list[meanResultsFramework,dataMeanModel] <- meanModel(seasonallyAdjDataLong, nObs ,degree =3,"mean")
# plot one example particpant 
pMeanModel <- plotModel(dataMeanModel,5, "mean", "mean")





### model changes in variance (intercept, mean, slope)
# @Ethan: should we treat the variance also as a time series?; should we always fit a polynomial model


# run mean-trend analysis for whole dataset
head(seasonallyAdjDataLong)
seasonallyAdjDataLong$yVar <-as.numeric(sapply(split(seasonallyAdjDataLong$y, seasonallyAdjDataLong$ppID), function(x) {as.numeric(scale(x, scale = FALSE, center = TRUE)^2)}))
head(seasonallyAdjDataLong)

list[varResultsFramework,dataVarModel] <- meanModel(seasonallyAdjDataLong, nObs ,degree =3, "variance")
# plot one example particpant 
pMeanModel <- plotModel(dataVarModel,4, "variance", "variance")



### change point analysis (both mean and variance) ###

# see also https://www.marinedatascience.co/blog/2019/09/28/comparison-of-change-point-detection-methods/
# see also https://lindeloev.github.io/mcp/articles/packages.html
# see also https://onlinelibrary.wiley.com/doi/full/10.1002/env.2576 (publication explaining the methodology of the non-parametric approach)
# this https://lindeloev.github.io/mcp/ looks really cool but also non-trivial to use; does require user input on the number 
# of segments and therefore maybe not exactly what we need; however, may be good to reference in the MS

# @Ethan

# change points in slope and variance in the simulated data
list[populationChangePoints, populationPPVarChangePoints] <-
  detect_Slope_And_Varaince_ChangePoints(data)

populationChangePoints
populationPPVarChangePoints

results <-
  filter(dataMeanModel, ppID %in% sample(data$ppID, 1))

# parametric version
# will be using the seasonally adjusted data 

# PELT stands for pruned linear exact time and is computationally very efficient 
# looking at the plot the result is actually pretty decent 
# always ends with n

# identify changes in mean
head(results)
changePointsMean <- cpt.mean(results$y, method = 'PELT', penalty = "AIC")
print(changePointsMean@cpts) 
plot(changePointsMean)

# identify changes in variance
changePointsVar <- cpt.var(results$y, method = 'PELT', penalty = "AIC")
print(changePointsVar@cpts) 
plot(changePointsVar)


# identify changes in both mean and variance (in the scale parameter)
# a bit too much 
changePointsMeanVar <- cpt.meanvar(results$y, method = 'PELT', penalty = "MBIC")
print(changePointsMeanVar@cpts) 
plot(changePointsMeanVar)

# non-parametric version
# non-parametric changepoint detection
# very good recovery of slope changes
changePointsNP <- cpt.np(results$y, method = 'PELT', penalty = "MBIC")
print(changePointsNP@cpts) 
plot(changePointsNP)

out <- cpt.var(results$y,pen.value=c(log(length(results$y)),100*log(length(results$y))),penalty="CROPS",method="PELT")
cpts.full(out)
plot(out,diagnostic=TRUE)
plot(out,ncpts=2)

# TODO: create a data set containing per pp statistics; nearly done only need to add changepopints

# TODO read the jebb 2017 paper
# TODO continue with Gerstorf
# jaab publications 








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


