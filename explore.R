# load packages -----------------------------------------------------------
if (!require("pacman"))
  install.packages("pacman")

pacman::p_load(
  ggplot2,
  ggthemes,
  gsubfn,
  ggpubr
)



# plot original data 
explore.plotTimeSeries <-
  function(data,
           numParticipants, nSample, nObs, 
           populationChangePoints,
           populationPPVarChangePoints, ylab) {
    
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
      y = ylab,
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



# decompose time series 
explore.decompose_TS <- function(data, numParticipants,type, freq, freqUnit) {
  nSample <- max(data$ppID)
  data_TS <- list()
  
  for (idx in 1:nSample) {
    data_TS[[idx]] <-
      ts(data[data$ppID == idx, ]$y[-1], frequency = freq)
  }
  

  data_TS_decomposed <-
    sapply(data_TS, function(x) {
      decompose(x, type = type)
    })
  
  data_TS_decomposed_x <-
    data.frame((sapply(data_TS_decomposed["x", ], c)))
  data_TS_decomposed_seasonal <-
    data.frame((sapply(data_TS_decomposed["seasonal", ], c)))
  data_TS_decomposed_trend <-
    data.frame((sapply(data_TS_decomposed["trend", ], c)))
  data_TS_decomposed_random <-
    data.frame((sapply(data_TS_decomposed["random", ], c)))

  
  # next transform each into long format and then merge columns
  # reshape into long format
  
  data_TS_decomposed_x_Long <-
    reshape(
      data_TS_decomposed_x,
      direction = 'long',
      varying = list(1:nSample),
      v.names = "data",
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
    c("ppID",
      "data",
      "time",
      "seasonal",
      "trend",
      "random")
  data_TS_decomposed_DF$ppID <-
    as.factor(data_TS_decomposed_DF$ppID)
  head(data_TS_decomposed_DF) #looks good
  
  # reshape into long format
  data_TS_decomposed_DF_Long <-
    reshape(
      data_TS_decomposed_DF,
      direction = 'long',
      varying = c("data", "seasonal", "trend", "random"),
      v.names = "value",
      timevar = 'dataType',
    )
  
  data_TS_decomposed_DF_Long$dataType <-
    as.factor(data_TS_decomposed_DF_Long$dataType)
  levels(data_TS_decomposed_DF_Long$dataType) <-
    c("data", "seasonal", "trend", "random")
  
  scaleFactor <- 1
  
  if (numParticipants > 1) {
    scaleFactor <- numParticipants / 2
  }
  
  
  data_TS_decomposed_DF_LongSubset <-
    filter(
      data_TS_decomposed_DF_Long,
      ppID %in% sample(data_TS_decomposed_DF_Long$ppID, numParticipants)
    ) #select random participants
  p <-
    ggplot(data = data_TS_decomposed_DF_LongSubset, aes(x = time, y = value, group = ppID)) +  geom_line(color = "grey", size = (1 /
                                                                                                                                   scaleFactor)) + geom_point(size = (1.5 / scaleFactor)) + facet_wrap(
                                                                                                                                     ~
                                                                                                                                       dataType,
                                                                                                                                     ncol = 1,
                                                                                                                                     strip.position = "left",
                                                                                                                                     scales = "free_y"
                                                                                                                                   )
  p <-
    p + labs(y = "value", x = freqUnit) + scale_x_continuous(breaks =
                                                             seq(0, nObs, freq), labels = c(0:(nObs / freq)))
  p <-
    p +       theme(legend.position = "none") + theme_bw() + theme(
      plot.title = element_text(face = "bold", size = 8) ,
      axis.ticks = element_line(colour = "grey70", size = 0.2),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(size = 1, colour = 'grey70'),
      axis.line.y = element_line(size = 1, colour = 'grey70'),
      legend.position = c(0.06, 0.75)
    )

  return(
    list(
      "p" = p,
      "data_TS_decomposed" = data_TS_decomposed
    )
  )
  
  
}

# visualize autocorrelation and partial autocorrelation patterns 
explore.visualize_autocorr <- function(data, lagMax, nSample) {
  acfMat <-
    data.frame(sapply(split(data$y, data$ppID), function(vec) {
      return(as.numeric(acf(
        as.numeric(vec), plot = FALSE, lag.max = lagMax
      )$acf)[2:(lagMax + 1)])
    }))
  
  
  pacfMat <-
    data.frame(sapply(split(data$y, data$ppID), function(vec) {
      return(as.numeric(pacf(
        as.numeric(vec), plot = FALSE, lag.max = lagMax
      )$acf))
    }))
  
  
  acfMatLong <- reshape(
    acfMat,
    idvar = "lag",
    varying = list(1:nSample),
    v.names = 'value',
    timevar = 'participant' ,
    direction = "long"
  )
  acfMatLong$lag <- as.factor(acfMatLong$lag)
  acfMatLong$autocorrelation <- round(acfMatLong$value, 2)
  
  p_acf <- ggplot(acfMatLong, aes(lag, participant, fill = value)) +
    geom_tile() + scale_fill_gradient2(low = '#0752A2', high = "#cc0000",
                                       limits = round(c(min(acfMat), max(acfMat)), 2),
                                       breaks = round(c(min(acfMat), 0, max(acfMat)), 2)
    ) + theme(legend.position = "none") + theme_bw() + theme(
      plot.title = element_text(face = "bold", size = 8) ,
      axis.ticks = element_line(colour = "grey70", size = 0.2),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(size = 1, colour = 'grey70'),
      axis.line.y = element_line(size = 1, colour = 'grey70'),
      panel.border = element_blank()
    ) + ggtitle("autocorrelation")
  
  
  pacfMatLong <- reshape(
    pacfMat,
    idvar = "lag",
    varying = list(1:nSample),
    v.names = 'value',
    timevar = 'participant' ,
    direction = "long"
  )
  pacfMatLong$lag <- as.factor(pacfMatLong$lag)
  pacfMatLong$autocorrelation <- round(pacfMatLong$value, 2)
  
  
  p_pacf <-
    ggplot(pacfMatLong, aes(lag, participant, fill = value)) +
    geom_tile() + scale_fill_gradient2(low = '#0752A2', high = "#cc0000",
                                       limits = round(c(min(pacfMat), max(pacfMat)), 2),
                                       breaks = round(c(min(pacfMat), 0, max(pacfMat)), 2)
    ) + theme(legend.position = "none") + theme_bw() + theme(
      plot.title = element_text(face = "bold", size = 8) ,
      axis.ticks = element_line(colour = "grey70", size = 0.2),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(size = 1, colour = 'grey70'),
      axis.line.y = element_line(size = 1, colour = 'grey70'),
      panel.border = element_blank()
    ) + ggtitle("partial autocorrelation")
  
  # combining the acf and pacf plots
  pACF_PLOTS <- grid.arrange(p_acf, p_pacf, nrow = 1)
  return(list("p_acf" = p_acf, "p_pacf" = p_pacf))
  
}



identifyChangePoints <- function(data, penalty) {
  cpMEAN <-cpt.mean(data$y, method = 'PELT', penalty = penalty)
  cpVAR <- cpt.var(data$y, method = 'PELT', penalty = penalty)
  cpNP_ALL <- cpt.np(data$y, method = 'PELT', penalty = penalty)
  result <- list("cpMEAN" = cpMEAN, "cpVAR" = cpVAR, "cpNP_ALL" = cpNP_ALL)
  return(result)
  
}



createVector <- function(x,type,len){
  outcome <- c(rep(0,len))
  
  if(type == 'var'){
    currvec <- x$cpVAR@cpts
    currvec <- currvec[-length(x$cpVAR@cpts)]
    
    
  }else if(type == 'mean'){
    currvec <- x$cpMEAN@cpts
    currvec <- currvec[-length(x$cpMEAN@cpts)]
    
  }else{
    currvec <- x$cpNP_ALL@cpts
    currvec <- currvec[-length(x$cpNP_ALL@cpts)]
    
  }
  outcome[currvec] <- 1
  return(outcome)
}


explore.changePoints <- function(data){
  
  data_nested <- data %>%
    group_by(ppID) %>%
    nest()
  
  len <- nrow(data_nested$data[[1]])
  
  data_nested <- data_nested %>%
    group_by(ppID) %>%
    mutate(
      changePoints_MBIC = map2(data,"MBIC", identifyChangePoints),
      changePoints_AIC  = map2(data, "AIC", identifyChangePoints),
      meanCP_MBIC      = pmap(list(changePoints_MBIC,'mean' ,len),createVector),
      varCP_MBIC      = pmap(list(changePoints_MBIC,'var' ,len),createVector),
      npCP_MBIC      = pmap(list(changePoints_MBIC,'np',len),createVector),
      meanCP_AIC      = pmap(list(changePoints_AIC,'mean',len),createVector),
      varCP_AIC      = pmap(list(changePoints_AIC,'var',len),createVector),
      npCP_AIC      = pmap(list(changePoints_AIC,'np',len),createVector),
    )
  
  
  #let's plot this 
  plottingDF <- data_nested %>%
    unnest(cols = c(data,meanCP_MBIC, varCP_MBIC, 
                    npCP_MBIC, meanCP_AIC, varCP_AIC, npCP_AIC)) %>%
    dplyr::select(c(ppID,time,meanCP_MBIC, varCP_MBIC, 
                    npCP_MBIC, meanCP_AIC, varCP_AIC, 
                    npCP_AIC))
  plottingDF <- data.frame(plottingDF)
  names(plottingDF) <- c("participant", "time","meanMBIC", "varianceMBIC", "bothMBIC", "meanAIC", "varianceAIC", "bothAIC")
  plottingDF #transform into long format 
  
  
  
  plottingDFLong <- reshape(
    plottingDF,
    idvar = c("time","participant"),
    varying = list(names(plottingDF)[3:8]),
    v.names = 'value',
    timevar = 'method' ,
    direction = "long",
    times = c("mean (MBIC)", "variance (MBIC)", "both (MBIC)", "mean (AIC)", "variance (AIC)", "both (AIC)")
  )
  
  
  newOrder <- c("mean (MBIC)", "variance (MBIC)", "both (MBIC)", "mean (AIC)", "variance (AIC)", "both (AIC)")
  plottingDFLong <- transform(plottingDFLong, method = factor(method, levels = newOrder))
  
  p<- ggplot(plottingDFLong, aes(time, participant, fill = value)) +
    geom_tile(show.legend = FALSE) + scale_fill_gradient2(low = 'white', high = "black") +
    facet_wrap("method",nrow = 2, ncol =3) +
    theme_bw() + theme(
      plot.title = element_text(face = "bold", size = 8) ,
      axis.ticks = element_line(colour = "grey70", size = 0.2),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = 'grey70', size = 1)
    ) + ggtitle("change points")
  
  return(p)
  
}










