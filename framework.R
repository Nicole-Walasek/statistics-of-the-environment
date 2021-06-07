#this script contains the framework functions 


# load packages -----------------------------------------------------------
if (!require("pacman"))
  install.packages("pacman")

pacman::p_load(
  ggplot2,
  ggthemes,
  gsubfn,
  ggpubr,
  tidyverse
)

# creates the framework

# helper function to compute to change points
# penalties can be AIC or MBIC, where MBIC applies a strciter change point criterion
# resulting in fewer changepoints 

identifyChangePoints <- function(data, penalty) {
  cpMEAN <-cpt.mean(data$y, method = 'PELT', penalty = penalty)
  cpVAR <- cpt.var(data$y, method = 'PELT', penalty = penalty)
  cpNP_ALL <- cpt.np(data$y, method = 'PELT', penalty = penalty)
  result <- list("cpMEAN" = cpMEAN, "cpVAR" = cpVAR, "cpNP_ALL" = cpNP_ALL)
  return(result)
}

spectral_analysis <- function(data){
  data.spect <- spectrum(data$y, log= "no", plot = FALSE)
  data.spect$freq_log <- log(data.spect$freq)
  data.spect$spec_log <- -log(data.spect$spec)
  return(data.spect)
}

spectralModel <- function(data){
  #specify the model of the mean here 
  lm(data$spec_log ~ data$freq_log)
}

extraxSpectralCoef <- function(data){
  return(data$coefficients[2])
}


framework.applyModel <- function(data, acLAG, meanModel = function(data){lm(y~time, data = data)}, varModel =function(data){lm(yVar~time, data = data)}) {
  
  # add variance in the outcome variable 
  data$yVar <-
    as.numeric(sapply(split(data$y, data$ppID), function(x) {
      as.numeric(scale(x, scale = FALSE, center = TRUE) ^ 2)
    }))
  
  data_nested <- data %>% 
    group_by(ppID) %>% 
    nest()
  
  data_models <- data_nested %>% 
    group_by(ppID) %>% 
    mutate(
      #mean stats
      data_lm      = map(data,meanModel),
      data_ac      = map(data,function(data){
        
        result <- data.frame(acf(as.numeric(data$y), lag.max = nrow(data)-1, plot = FALSE)$acf[2:(nrow(data))], 
                             pacf(as.numeric(data$y), plot = FALSE)$acf[1:(nrow(data)-1)])
        names(result) <- c("acf", "pacf") 
        result$lag <- as.factor(c(1:(nrow(data)-1)))
        result <- as_tibble(result)
      }),
      lm_statsMean = map(data_lm,glance),
      lm_coefsMean = map(data_lm,tidy),
      lm_valsMean  = map(data_lm,augment),
      
      # variance statistics
      dataVar_lm      = map(data,varModel),
      dataVar_ac      = map(data,function(data){
        
        result <- data.frame(acf(as.numeric(data$yVar), lag.max = nrow(data)-1, plot = FALSE)$acf[2:(nrow(data))], 
                             pacf(as.numeric(data$yVar), plot = FALSE)$acf[1:(nrow(data)-1)])
        names(result) <- c("acf", "pacf") 
        result$lag <- as.factor(c(1:(nrow(data)-1)))
        result <- as_tibble(result)
      }),
      lm_statsVar = map(dataVar_lm,glance),
      lm_coefsVar = map(dataVar_lm,tidy),
      lm_valsVar  = map(dataVar_lm,augment),
      
      # add changepoints
      changePoints_MBIC = map2(data,"MBIC", identifyChangePoints),
      changePoints_AIC  = map2(data, "AIC", identifyChangePoints),
      
      # add spectral analysis 
      dataSpectral = map(data,spectral_analysis),
      dataSpectralLM = map(dataSpectral, spectralModel),
      dataSpectralCoef = map_dbl(dataSpectralLM, extraxSpectralCoef)
    )
  
  #adding all the other model coefficients for the mean 
  
  # add model coefficients  
  #add ppID
  ppID <- data_models$ppID
  df <- data.frame(ppID)
  df$m_intercept <- map_dbl(data_models$lm_coefsMean, function(x) unlist(x[1,"estimate"]))
  degree = nrow(data_models$lm_coefsMean[[1]])-1
  for(idx in 1:degree){
    current <-map_dbl(data_models$lm_coefsMean, function(x) unlist(x[idx+1,"estimate"]))
    df <- cbind(df,current)
    
    names(df)[names(df) == "current"] <-  sprintf("m_slope%s",idx)
  }
  
  #add rsquared
  df$m_rsqrd <- data_models$lm_statsMean %>% map_dbl("r.squared")
  
  # add mean and sd
  df$m_mean <- map_dbl(data_models$data, function(data){mean(data$y)})
  df$m_sd <- map_dbl(data_models$data, function(data){sd(data$y)})
  
  
  dfTibble <- as_tibble(df)
  
  # adding all the ac and pac values
  dfAC <- data.frame(data.frame(matrix(ncol = 0, nrow = nrow(data_models))))
  dfPAC <- data.frame(data.frame(matrix(ncol = 0, nrow = nrow(data_models))))
  for(idx in 1:acLAG){
    current <-map_dbl(data_models$data_ac, function(x) unlist(x[idx,"acf"]))
    dfAC <- cbind(dfAC,current)
    names(dfAC)[names(dfAC) == "current"] <-  sprintf("m_ac%s",idx)
    
    current <-map_dbl(data_models$data_ac, function(x) unlist(x[idx,"pacf"]))
    dfPAC <- cbind(dfPAC,current)
    names(dfPAC)[names(dfPAC) == "current"] <-  sprintf("m_pac%s",idx)
  }
  dfACTibble <- as_tibble(dfAC)
  dfPACTibble <- as_tibble(dfPAC)
  meanStats = bind_cols("model_coef" = dfTibble, "AC"= dfACTibble, "PAC" = dfPACTibble)
  
  
  #adding all the other model coefficients for the variance 
  
  ppID <- data_models$ppID
  df <- data.frame(ppID)
  df$v_intercept <- map_dbl(data_models$lm_coefsVar, function(x) unlist(x[1,"estimate"]))
  degree = nrow(data_models$lm_coefsVar[[1]])-1
  for(idx in 1:degree){
    current <-map_dbl(data_models$lm_coefsVar, function(x) unlist(x[idx+1,"estimate"]))
    df <- cbind(df,current)
    
    names(df)[names(df) == "current"] <-  sprintf("v_slope%s",idx)
  }
  
  #add rsquared
  df$v_rsqrd <- data_models$lm_statsVar %>% map_dbl("r.squared")
  
  # add mean and sd
  df$mean <- map_dbl(data_models$data, function(data){mean(data$yVar)})
  df$sd <- map_dbl(data_models$data, function(data){sd(data$yVar)})
  
  #add ppID
  df$ppID <- data_models$ppID
  
  
  dfTibble <- as_tibble(df)
  
  # adding all the ac values
  dfAC <- data.frame(data.frame(matrix(ncol = 0, nrow = nrow(data_models))))
  dfPAC <- data.frame(data.frame(matrix(ncol = 0, nrow = nrow(data_models))))
  for(idx in 1:acLAG){
    current <-map_dbl(data_models$dataVar_ac, function(x) unlist(x[idx,"acf"]))
    dfAC <- cbind(dfAC,current)
    names(dfAC)[names(dfAC) == "current"] <-  sprintf("v_ac%s",idx)
    
    current <-map_dbl(data_models$dataVar_ac, function(x) unlist(x[idx,"pacf"]))
    dfPAC <- cbind(dfPAC,current)
    names(dfPAC)[names(dfPAC) == "current"] <-  sprintf("v_pac%s",idx)
  }
  dfACTibble <- as_tibble(dfAC)
  dfPACTibble <- as_tibble(dfPAC)
  
  varStats = bind_cols("model_coef" = dfTibble, "AC"= dfACTibble, "PAC" = dfPACTibble)
  
  # create change point statistics 
  meanNumberMBIC <- map_dbl(data_models$changePoints_MBIC, function(x){length(x$cpMEAN@cpts)-1})
  df <- data.frame(meanNumberMBIC)
  df$meanMedianMBIC <- map_dbl(data_models$changePoints_MBIC, function(x){median(x$cpMEAN@cpts[-length(x$cpMEAN@cpts)])})
  df$varNumberMBIC <- map_dbl(data_models$changePoints_MBIC, function(x){length(x$cpVAR@cpts)-1})
  df$varMedianMBIC <- map_dbl(data_models$changePoints_MBIC, function(x){median(x$cpVAR@cpts[-length(x$cpVAR@cpts)])})
  df$npNumberMBIC <- map_dbl(data_models$changePoints_MBIC, function(x){length(x$cpNP_ALL@cpts)-1})
  df$npMedianMBIC <- map_dbl(data_models$changePoints_MBIC, function(x){median(x$cpNP_ALL@cpts[-length(x$cpNP_ALL@cpts)])})
  
  df$meanNumberAIC <- map_dbl(data_models$changePoints_AIC, function(x){length(x$cpMEAN@cpts)-1})
  df$meanMedianAIC <- map_dbl(data_models$changePoints_AIC, function(x){median(x$cpMEAN@cpts[-length(x$cpMEAN@cpts)])})
  df$varNumberAIC <- map_dbl(data_models$changePoints_AIC, function(x){length(x$cpVAR@cpts)-1})
  df$varMedianAIC <- map_dbl(data_models$changePoints_AIC, function(x){median(x$cpVAR@cpts[-length(x$cpVAR@cpts)])})
  df$npNumberAIC <- map_dbl(data_models$changePoints_AIC, function(x){length(x$cpNP_ALL@cpts)-1})
  df$npMedianAIC <- map_dbl(data_models$changePoints_AIC, function(x){median(x$cpNP_ALL@cpts[-length(x$cpNP_ALL@cpts)])})
  df <- as_tibble(df)
  
  
  
  data_models <- data_models %>%
    add_column("Statistics_Mean" = meanStats)
  
  data_models <- data_models %>%
    add_column("Statistics_Variance" = varStats)
  
  data_models <- data_models %>%
    add_column("Statistics_ChangePoints" = df)
  
  
  
  return(data_models)
  
}


# plots the framework models

framework.plotModel <- function(dataAll, numPP, modelTitle, type, freq = 1, xlab = "time", ylab) {
  timeDF <- dataAll %>%
    dplyr::select(data, ppID) %>%
    unnest(cols = c(data))%>%
    dplyr::select(time, ppID)
  
  timePlot <- timeDF$time
  
  if(type == "mean"){
    
    # plot mean
    data <- dataAll %>%
      dplyr::select(lm_valsMean,ppID) %>% 
      unnest(cols = c(lm_valsMean)) %>%
      add_column(timePlot = timePlot) 
    
    data <- as.data.frame(data)
    
  }else{
    
    # plot variance
    data <- dataAll %>%
      dplyr::select(lm_valsVar,ppID) %>%  
      unnest(cols = c(lm_valsVar))%>%
      add_column(timePlot = timePlot) 
    
    data <- as.data.frame(data)
    
    
  }
  
  results <-
    filter(data, ppID %in% sample(data$ppID, numPP)) #select random participants for plotting
  results$ppID <- as.factor(results$ppID)
  names(results)[names(results) == '.fitted'] <- 'fittedTimeModel'
  names(results)[names(results) == '.resid'] <-
    'residualsTimeModel'
  
  if (type == "mean") {
    names(results)[names(results) == 'y'] <- 'originalData'
    yLab <- ylab
  } else{
    names(results)[names(results) == 'yVar'] <- 'originalData'
    yLab <- sprintf("%s variance", ylab)
  }
  
  
  pModel <-
    ggplot(data = results) + geom_segment(aes(x = timePlot, xend = timePlot, y = originalData, yend= originalData - residualsTimeModel), color = "grey")+
    geom_line(aes(x = timePlot, y = originalData, group = ppID), color = "black", size = 0.6) + 
    geom_point(aes(x =
                     
                     timePlot, y = originalData, group = ppID), size = (1.5 /(0.8*numPP))) + geom_line(aes(x = timePlot, y = fittedTimeModel, group = ppID),size =0.8, color = "firebrick")
  pModel <-
    pModel + labs(y = yLab, x = xlab) + scale_x_continuous(breaks =
                                                             seq(0, nObs, freq), labels = c(0:(nObs / freq)))
  
  pModel <- pModel + labs(title = paste(modelTitle, " model"))
  pModel <-
    pModel +       theme(legend.position = "none") + theme_bw() + theme(
      plot.title = element_text(face = "bold", size = 8) ,
      axis.ticks = element_line(colour = "grey70", size = 0.2),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = 'grey70', size = 1),
      legend.position = c(0.06, 0.75)
    )
  
  pRes <-
    ggplot(data = results, aes(x = timePlot, y = residualsTimeModel, group = ppID)) + geom_line(color = "#535754", aes(group = ppID), size = 0.6) + geom_point(aes(group = ppID), size = (1.5 /
                                                                                                                                                                                            (0.8*numPP)))
  pRes <-
    pRes + labs(y = "residual", x = xlab) + scale_x_continuous(breaks =
                                                                 seq(0, nObs, freq), labels = c(0:(nObs / freq)))
  pRes <- pRes + labs(title = "residual error series")
  pRes <-
    pRes +       theme(legend.position = "none") + theme_bw() + theme(
      plot.title = element_text(face = "bold", size = 8) ,
      axis.ticks = element_line(colour = "grey70", size = 0.2),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = 'grey70', size = 1),
      legend.position = c(0.06, 0.75)
    )
  
  # combining the plots
  pModelEval <- grid.arrange(pModel, pRes, nrow = 1)
  
  return(pModelEval)
}


framework.plotColorOfNoise <- function(data, binwidth) {
  df_spectralCoef <- data.frame(data$dataSpectralCoef)
  names(df_spectralCoef) <- c("beta")
  color <- rep("", nrow(df_spectralCoef))
  color[df_spectralCoef < 0] <- "blue"
  color[df_spectralCoef > -0.5 &df_spectralCoef < 0.5] <- "white"
  color[df_spectralCoef >= 0.5 &df_spectralCoef < 1.5] <- "pink"
  color[df_spectralCoef >= 1.5 ] <- "brown"
  df_spectralCoef$color <- color
  
  group.colors <- c(blue = "#00008B", white = "white", pink = "#E30B5C", brown = "#A52A2A" )
  
  p_spectralCoef <- ggplot(df_spectralCoef, aes(x = beta, fill = color)) + 
    geom_histogram(binwidth = binwidth,alpha=0.5, position="identity", color = 'black')+
    scale_fill_manual(values=group.colors) +
    
    geom_vline(aes(xintercept=0),
               color="black", linetype="dashed", size=1) +
    geom_vline(aes(xintercept=1),
               color="black", linetype="dashed", size=1)+
    geom_vline(aes(xintercept=2),
               color="black", linetype="dashed", size=1)
  p_spectralCoef  <- p_spectralCoef + theme(legend.position = "none") + theme_bw() + theme(
    plot.title = element_text(face = "bold", size = 8) ,
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(size = 1, colour = 'grey70'),
    axis.line.y = element_line(size = 1, colour = 'grey70'),
    panel.border = element_blank(),
    legend.position = "bottom"
  )
  
  return(p_spectralCoef)
  
  
  
  
}
