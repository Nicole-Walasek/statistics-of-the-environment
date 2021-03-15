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

options(digits = 5)
# functions for plotting and simulating data --------------------------------




plotSimData <-
  function(data,
           numParticipants,
           populationChangePoints) {
    # compute the sample autocorrelation and its standard deviation
    
    dataWide <- cast(data, ID ~
                       time, value = "yNorm")
    
    results = apply(dataWide[,-1], 1, acf, lag.max = 1, plot = FALSE)
    sampleAutoCorrelations <- c()
    for (i in 1:nSample) {
      sampleAutoCorrelations[i] <- results[i][[1]]$acf[2]
    }
    autoCorrMean <- mean(sampleAutoCorrelations)
    autoCorrSD <- sd(sampleAutoCorrelations)
    
    
    
    dataNew <-
      filter(data, ppID %in% sample(data$ppID, numParticipants)) #select random participants
    
    
    firstPart <-
      ggplot(dataNew, aes(x = time, y = yNorm, group = ppID)) + geom_line(colour = 'steelblue3',
                                                                          size = 0.6,
                                                                          alpha = 0.2) +
      theme(legend.position = "none") + theme_bw() + theme(
        plot.title = element_text(face = "bold", size = 8) ,
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1, colour = 'grey70'),
        axis.line.y = element_line(size = 1, colour = 'grey70'),
        legend.position = c(0.06, 0.75)
      )
    
    finalPlot <-
      firstPart + stat_summary(
        data = data,
        fun = mean,
        geom = "line",
        size = 0.8,
        linetype = "solid",
        colour = 'darkred' ,
        aes(x = time, y = yNorm),
        inherit.aes = FALSE
      ) +
      stat_summary(
        data = data,
        fun = mean,
        fun.min = function(x)
          mean(x) - sd(x),
        fun.max = function(x)
          mean(x) + sd(x),
        geom = "ribbon",
        aes(x = time, y = yNorm),
        alpha = 0.3,
        fill = "grey70",
        colour = NA,
        inherit.aes = FALSE
      ) + stat_summary(
        data = data,
        fun.data = mean_se,
        geom = "ribbon",
        aes(x = time, y = yNorm),
        inherit.aes = FALSE,
        alpha = 0.3
      )
    
    finalPlotLabeled <- finalPlot + labs(
      x = "time",
      y = "harshness",
      title = paste(
        'autocorrelation = ',
        round(autoCorrMean, 3),
        "(",
        round(autoCorrSD, 3),
        ")"
      )
    ) + scale_x_continuous(
      expand = c(0, 0),
      limits = c(0, nObs),
      breaks = round(seq(0, nObs, by = 1), 1)
    ) + scale_y_continuous(
      expand = c(0, 0),
      limits = c(1, max(data$yNorm)),
      breaks = round(seq(1, max(data$yNorm), by = 1), 1)
    )
    
    
    
    # + annotate(
    #                                        'text',
    #                                        x = 10,
    #                                        y = 9,
    #                                        size = 4.5,
    #                                        color = 'black'
    #                                      )
    #
    
    # check whether there were change points
    
    
    
    if (!is.character(populationChangePoints)) {
      yMin = c()
      yMax = c()
      for (idx in 1:length(populationChangePoints)) {
        yMin[idx] <-
          min(dataNew$yNorm[dataNew$time == populationChangePoints[idx]])
        yMax[idx] <-
          max(dataNew$yNorm[dataNew$time == populationChangePoints[idx]])
      }
      
      
      #marks <- data.frame(x = rep(populationChangePoints, each = 2), y = rep(c(2,8), times = length(populationChangePoints)))
      for (idx in 1:length(populationChangePoints)) {
        currPoint <- populationChangePoints[idx]
        finalPlotLabeled <-
          finalPlotLabeled + geom_segment(
            x = currPoint,
            y = yMin[idx],
            xend = currPoint,
            yend = yMax[idx],
            alpha = 0.2,
            linetype = 'longdash'
          )
        
      }
      
    }
    
    return(finalPlotLabeled)
    
    
  }


simData <-
  function(nSample,
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
           rangeMax) {
    ### parameters
    
    # nSample: desired sample size
    # popInt: population intercept
    # popSlope: population slope
    # sigmaModel: is the variance-covaraince matrix
    # nObs: number of observations/repeated measures
    # changePoints: number of change points in the population trend; can also be use specified individual changePoints       or 'None' if the user does not want any
    
    # slopeChangeMean: mean value of the slope changing at each change point
    # slopeChangeVar: variance around the mean slope change
    # noiseVar: varaince around the per-participant regression line; also controls the autocorrelation
    
    RE = mvrnorm(nSample, mu = c(popInt, popSlope), Sigma = SigmaModel)
    
    
    colnames(RE) = c("ints", "slopes")
    t(round(RE, 2))
    
    populationChangePoints = changePoints
    
    
    particpantSlopesAcrossTime = rep(RE[, 2], each = (nObs + 1))
    
    if (is.character(changePoints)) {
      # if we are not interested in random change points
      particpantSlopesAcrossTime = rep(RE[, 2], each = (nObs + 1))
      
    } else {
      if (!is.numeric(changePoints)) {
        print('change points must be a single numeric value or a vector of numeric vaues')
        
        return()
        
      }
      
      
      # we have to deal with changepoints
      
      if (length(changePoints) == 1) {
        # if this is true we know that the user only specified the number of change points
        # we sample the use specified number of change points
        
        populationChangePoints = sort(sample(2:nObs-1, changePoints))
        
        
      } else {
        # the user specified the change points him/her self
        
        populationChangePoints = changePoints
        
      }
      
      if (length(slopeChangeMean) < length(populationChangePoints)) {
        slopeChangeMeanList = rep(mean(slopeChangeMean),
                                  each = length(populationChangePoints))
      } else {
        slopeChangeMeanList = slopeChangeMean
      }
      
      if (length(slopeChangeVar) < length(populationChangePoints)) {
        slopeChangeVarList = rep(mean(slopeChangeVar), each = length(populationChangePoints))
      } else {
        slopeChangeVarList = slopeChangeVar
      }
      
      
      if ((length(slopeChangeMean) > length(populationChangePoints)) |
          (length(slopeChangeVar) > length(populationChangePoints))) {
        print('cannot have more slope changes than changepoints')
        return()
        
      }
      
      # let's derive the participant slope vectors
      
      
      
      particpantSlopesAcrossTimeDF = data.frame(
        ID   = rep(1:nSample,   each = (nObs + 1)),
        time = rep(0:nObs, times = nSample),
        RE.s = particpantSlopesAcrossTime
      )
      
      
      particpantSlopesAcrossTimeDFWide = cast(particpantSlopesAcrossTimeDF, ID ~
                                                time, value = "RE.s")
      
      # now apply the manipulations
      
      loopDF = data.frame(x = populationChangePoints, y = slopeChangeMeanList, z = slopeChangeVarList)
      
      if (changePointsVar == 'None') {
        for (currRowIDX in (1:nrow(loopDF))) {
          # x refers to the change point in units of time
          # y refers to the mean slope change
          # z refers to the variance around the mean slope change
          currRow <- loopDF[currRowIDX, ]
          
          # TODO consider adding per participant variance in the change point
          # woud be easiest to accomplish with another loop through the entire dataframe
          
          
          particpantSlopesAcrossTimeDFWide[, (currRow$x + 1):(nObs + 2)] <-
            particpantSlopesAcrossTimeDFWide[, (currRow$x + 1):(nObs + 2)] + rnorm(n = nSample ,
                                                                                   mean = currRow$y,
                                                                                   sd =
                                                                                     currRow$z)
        }
        
      } else {
        for (currRowIDX in (1:nrow(loopDF))) {
          # x refers to the change point in units of time
          # y refers to the mean slope change
          # z refers to the variance around the mean slope change
          currRow <- loopDF[currRowIDX, ]
          
          # TODO consider adding per participant variance in the change point
          # woud be easiest to accomplish with another loop through the entire dataframe
          
          for (ppIDX in (1:nSample)) {
            currCol <-  abs(round(rnorm(1, currRow$x, changePointsVar)))
            if (currCol >= nObs+2){ 
              currCol <- nObs+1}
            particpantSlopesAcrossTimeDFWide[particpantSlopesAcrossTimeDFWide$ID == ppIDX, (currCol + 1):(nObs + 2)] <-
              particpantSlopesAcrossTimeDFWide[particpantSlopesAcrossTimeDFWide$ID == ppIDX, (currCol + 1):(nObs + 2)] +
              rnorm(1, mean = currRow$y, sd = currRow$z)
            
            
          }
          
        }
        
      }
      
      
      
      particpantSlopesAcrossTimeDFLong = melt(particpantSlopesAcrossTimeDFWide, id =
                                                c("ID"))
      particpantSlopesAcrossTimeDFLongSorted <-
        particpantSlopesAcrossTimeDFLong[order(particpantSlopesAcrossTimeDFLong$ID,
                                               particpantSlopesAcrossTimeDFLong$time), ]
      
      
      particpantSlopesAcrossTime = particpantSlopesAcrossTimeDFLongSorted$value
      
    }
    
    data = data.frame(
      ID   = rep(1:nSample,   each = (nObs + 1)),
      time = rep(0:nObs,   times = nSample),
      RE.i = rep(RE[, 1], each = (nObs + 1)),
      RE.s = particpantSlopesAcrossTime,
      y    = NA
    )
    
    
    #simulate the response variable
    
    # add random noise around the linear per participant line
    # the sd controls the autocorrelation
    
    # an alternative may be to add noise around the slope for each time step
    # y = with(data, (0 + RE.i) + (0 + RE.s + rnorm(n=nSample*nObs, mean=0, sd=0.2))*time)
    
    
    y = with(data,
             (0 + RE.i) + (0 + RE.s) * time + rnorm(
               n = nSample * (nObs + 1),
               mean = 0,
               sd = noiseVar
               
             ))
    
    
    
    # do you want to add missing data at random?
    
    #m  = rbinom(n=ni*nj, size=1, prob=.1)
    #y[m==1] = NA
    
    data$y  = y
    
    #create a factor ppID
    data$ppID <- as.factor(data$ID)
    
    
    # add a rescaled outcome variable to the dataframe
    
    data$yNorm <- rescale(data$y, c(rangeMin, rangeMax))
    
    
    return(list("data" = data, "populationChangePoints" = populationChangePoints))
    
    
    
    
  }



# generate data -----------------------------------------------------------


nSample = 100 # number of subjects

# between participant variances in intercept and slope, and covariance between the two
SigmaModel = rbind(c(1, 0.2),
                   c(0.2, 0.1))

popInt = 2
popSlope = 1.5

nObs   = 20 # number of repeated measures

changePoints = 2 # should be 'None if there is None'
changePointsVar = 1
slopeChangeMean = c(-0.2, 0.2)
slopeChangeVar = 0.1

noiseVar = 5

rangeMin = 1
rangeMax = 10

list[data, populationChangePoints] <- simData(
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
  rangeMax
)


head(data, 50)

print(populationChangePoints)

# plot raw data and descriptives ------------------------------------------



numParticipants = 30
p <- plotSimData(data, numParticipants, populationChangePoints)

print(p)



# make a panel plot for different parameter combinations ------------------

# syntax for figures is L or H for the autocorrelation, D(ecreasing),I(ncreasing),N(Null) for the main effect

nSample = 100 # number of subjects
SigmaModel = rbind(c(1, 0.3),
                   c(0.3, 0.5))


nObs   = 20 # number of repeated measures

changePoints = 2 # should be 'None if there is None'
changePointsVar = 1
slopeChangeMean = c(-0.6, 1)
slopeChangeVar = 0.1

rangeMin = 1
rangeMax = 10




## @ Ethan: new parameters that create less variance:


nSample = 100 # number of subjects

# between participant variances in intercept and slope, and covariance between the two

# ((varIntercept, covar), (coVar,varSlope))
SigmaModel = rbind(c(5, 0.2),
                   c(0.2, 0.05))


nObs   = 20 # number of repeated measures


changePoints = 2 # should be 'None if there is None'
changePointsVar = 1
slopeChangeMean = c(-0.05,0.05)
slopeChangeVar = 0.05


rangeMin = 1
rangeMax = 10






# low autocorr, increasing trend


popInt = 5
popSlope = 1.5
noiseVar = 7

list[data, populationChangePoints] <- simData(
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
  rangeMax
)

numParticipants = 30
LI <- plotSimData(data, numParticipants, populationChangePoints)


# high autocorr, increasing trend


popInt = 5
popSlope = 1.5
noiseVar = 1

list[data, populationChangePoints] <- simData(
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
  rangeMax
)

numParticipants = 30
HI <- plotSimData(data, numParticipants, populationChangePoints)


# low autocorr, decreasing trend


popInt = 5
popSlope = -1.5
noiseVar = 6

list[data, populationChangePoints] <- simData(
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
  rangeMax
)

numParticipants = 30
LD <- plotSimData(data, numParticipants, populationChangePoints)


# high autocorr, Dncreasing trend


popInt = 5
popSlope = -1.5
noiseVar = 1

list[data, populationChangePoints] <- simData(
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
  rangeMax
)

numParticipants = 30
HD <- plotSimData(data, numParticipants, populationChangePoints)


# low autocorr, no trend


popInt = 5
popSlope = 0
noiseVar = 6

list[data, populationChangePoints] <- simData(
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
  rangeMax
)

numParticipants = 30
LN <- plotSimData(data, numParticipants, populationChangePoints)


# high autocorr, no trend


popInt = 5
popSlope = 0
noiseVar = 1

list[data, populationChangePoints] <- simData(
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
  rangeMax
)

numParticipants = 30
HN <- plotSimData(data, numParticipants, populationChangePoints)





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
  font.label =list(color="black",size=10)
)
figure

# these settings look good just need to change the path 
ggsave(device='pdf',dpi=600 , width = 12, height = 8, filename="figureStatsEnvLessVar2.pdf")


# fit model ---------------------------------------------------------------


# fitting a mixed effects model

modelFitted <- lme4::lmer(y ~ time + (time | ppID), data)

# let's look at the mode results

summary(modelFitted)
sjPlot::tab_model(modelFitted)

plot_model(modelFitted, type = "pred", terms = "time")
