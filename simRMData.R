# load packages -----------------------------------------------------------

pacman::p_load(
  ggplot2,
  ggthemes,
  MASS,
  dplyr,
  reshape,
  scales,
  Hmisc,
  gsubfn,
  ggpubr
)




options(digits = 5)
# functions for plotting and simulating data --------------------------------
detect_ChangePoints <- function(data) {
  # detect slope change points
  outcomeSlope <-
    data.frame(sapply(split(data$RE.s, data$ppID), function(vec) {
      return(which(abs(diff(vec)) > 0))
    }))
  
  # next average across participants
  resultSlope <- round(rowMeans(outcomeSlope))
  
  if (length(resultSlope) == (nrow(data[data$ppID == 1, ]) - 1)) {
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
  
  
  if (length(resultVar) == (nrow(data[data$ppID == 1, ]) - 1)) {
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

plotSimData <-
  function(data,
           numParticipants, nSample, nObs, 
           populationChangePoints,
           populationPPVarChangePoints) {
    
    # compute the sample autocorrelation and its standard deviation
    
    dataWide <- cast(data, ppID ~
                       time, value = "y")
    
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
      ggplot(dataNew, aes(x = time, y = y, group = ppID)) + geom_line(colour = 'steelblue3',
                                                                          size = 0.4,
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
        aes(x = time, y = y),
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
        aes(x = time, y = y),
        alpha = 0.3,
        fill = "grey70",
        colour = NA,
        inherit.aes = FALSE
      ) + stat_summary(
        data = data,
        fun.data = mean_se,
        geom = "ribbon",
        aes(x = time, y = y),
        inherit.aes = FALSE,
        alpha = 0.3
      )

    finalPlotLabeled <- finalPlot + labs(
      x = "measurement point",
      y = "morbidity-mortality",
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
      limits = c(1, max(data$y)),
      breaks = round(seq(1, max(data$y), by = 1), 1)
    )
    
    ## marke change points in per participant variance
    
    if (!is.character(populationPPVarChangePoints)) {
      yMin = c()
      yMax = c()
      
      for (idx in 1:length(populationPPVarChangePoints)) {
        currMean = (mean(data$y[data$time == populationPPVarChangePoints[idx]]))
        currSD = (sd(data$y[data$time == populationPPVarChangePoints[idx]]))
        
        
        yMin[idx] <- currMean + currSD + 0.5
        yMax[idx] <- currMean - currSD - 0.5
      }
      
      for (idx in 1:length(populationPPVarChangePoints)) {
        currPoint <- populationPPVarChangePoints[idx]
        finalPlotLabeled <-
          finalPlotLabeled + geom_segment(
            x = currPoint,
            y = yMin[idx],
            xend = currPoint,
            yend = yMax[idx],
            alpha = 0.2,
            linetype = 'longdash',
            color = 'darkred',
            size = 0.8
          )
        
      }
      
    }
    
    
    if (!is.character(populationChangePoints)) {
      yMin = c()
      yMax = c()
      
      
      for (idx in 1:length(populationChangePoints)) {
        currMean = (mean(data$y[data$time == populationChangePoints[idx]]))
        currSD = (sd(data$y[data$time == populationChangePoints[idx]]))
        
        
        yMin[idx] <- currMean + currSD + 0.15
        yMax[idx] <- currMean - currSD - 0.15
      }
      
      
      
      for (idx in 1:length(populationChangePoints)) {
        currPoint <- populationChangePoints[idx]
        finalPlotLabeled <-
          finalPlotLabeled + geom_segment(
            x = currPoint,
            y = yMin[idx],
            xend = currPoint,
            yend = yMax[idx],
            alpha = 0.2,
            linetype = 'longdash',
            color = 'seagreen',
            size = 0.8
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
           rangeMax,
           missingness,
           ppVar,
           ppVarChangePoints) {
    ### parameters
    
    # nSample: desired sample size
    # popInt: population intercept
    # popSlope: population slope
    # sigmaModel: is the variance-covaraince matrix
    # nObs: number of observations/repeated measures
    # changePoints: number of change points in the population trend; can also be use specified individual changePoints
    # or 'None' if the user does not want any
    
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
      
      #TODO does it make sense to add time dependent changes in slope here? 
      
      if (length(changePoints) == 1 & is.integer(changePoints)) {
        # if this is true we know that the user only specified the number of change points
        # we sample the use specified number of change points
        
        populationChangePoints = sort(sample(2:(nObs - 1), changePoints))
        
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
            if (currCol >= nObs + 2) {
              currCol <- nObs + 1
            }
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
    
    
    # you can manipulate the variance within participants here
    
    populationPPVarChangePoints = "None"
    if (is.function(ppVar)) {
      x <- 0:nObs
      x <- x + 0.1
      
      if (min(ppVar(x)) < 0) {
        noiseFunctionIB <- ppVar(x) + (abs(min(ppVar(x))) + 0.1)
        noiseFunction <-
          (noiseFunctionIB) / abs(sum(noiseFunctionIB)) * (noiseVar * max(x))
        
      } else {
        noiseFunction <-
          ((ppVar(x)) / abs(sum(ppVar(x)))) * (noiseVar * max(x))
        
      }
      plot(noiseFunction)
      
      RE.var = rep(noiseFunction, times = nSample)
      
    } else if (is.character(ppVar))  {
      if (ppVar == 'increasing') {
        # implement an increasing patterns such that the mean corresponds to
        # noiseVar
        
        x <- 0:nObs
        noiseIncreasing <- x + 0.1
        noiseIncreasing <-
          ((noiseIncreasing) / abs(sum(noiseIncreasing))) * (noiseVar * max(x))
        RE.var = rep(noiseIncreasing, times = nSample)
        
      } else if (ppVar == 'decreasing') {
        x <- 0:nObs
        noiseDecreasing <- sort(x + 0.1, decreasing = T)
        noiseDecreasing <-
          ((noiseDecreasing) / abs(sum(noiseDecreasing))) * (noiseVar * max(x))
        RE.var = rep(noiseDecreasing, times = nSample)
        
      }else{
        # assume that the user put in None, so variance is the same across time
        RE.var = rep(noiseVar, each = (nObs + 1) * nSample)
        
      }
      
    } else if (is.numeric(ppVar)) {
      # change points in varaince
      if (length(ppVarChangePoints) == 1 &
          is.integer(ppVarChangePoints)) {
        populationPPVarChangePoints = sort(sample(2:(nObs - 1), ppVarChangePoints))
        
        
        if (length(ppVar) == ppVarChangePoints){
          ppVarArray = ppVar
        } else {
          ppVarArray = rep(ppVar, each = length(populationPPVarChangePoints))
          
        }
        
        
      } else {
        populationPPVarChangePoints = ppVarChangePoints
        
        if (length(ppVar) < length(populationPPVarChangePoints)) {
          ppVarArray = rep(mean(ppVar),
                           each = length(populationPPVarChangePoints))
        } else {
          ppVarArray = ppVar
        }
        
        
        
      }
      
      ppVarAcrossTime <- rep(noiseVar, each = (nObs + 1))
      
      for (idx in 1:length(populationPPVarChangePoints)) {
        ppVarIDX = populationPPVarChangePoints[idx]
        ppVarVal = ppVarArray[idx]
        ppVarAcrossTime[ppVarIDX:length(ppVarAcrossTime)] <-
          ppVarAcrossTime[ppVarIDX:length(ppVarAcrossTime)] + ppVarVal
        
      }
      RE.var = rep(ppVarAcrossTime, times = nSample)
      
      
      
    } 
    
    
    data = data.frame(
      ID   = rep(1:nSample,   each = (nObs + 1)),
      time = rep(0:nObs,   times = nSample),
      RE.i = rep(RE[, 1], each = (nObs + 1)),
      RE.s = particpantSlopesAcrossTime,
      RE.var = RE.var ,
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
               sd = abs(RE.var)
               
             ))
    
    
    # do you want to add missing data at random?
    
    if (missingness != 'None') {
      m  = rbinom(n = sample * nObs,
                  size = 1,
                  prob = missingness)
      y[m == 1] = NA
      
      
    }
    
    data$y  = y
    
    #create a factor ppID
    data$ppID <- as.factor(data$ID)
    
    
    # add a rescaled outcome variable to the dataframe
    
    data$yNorm <- rescale(data$y, c(rangeMin, rangeMax))
    
    
    
    
    # for vizualisation purposes
    if(is.numeric(populationChangePoints)){
      populationChangePoints = populationChangePoints -1 
    }
    
    if (is.numeric(populationPPVarChangePoints)){
      populationPPVarChangePoints = populationPPVarChangePoints -1
      
    }
    
    return(
      list(
        "data" = data,
        "populationChangePoints" = populationChangePoints,
        "populationPPVarChangePoints
" = populationPPVarChangePoints
      )
    )
    
    
  }
