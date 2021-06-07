

preprocess.adjust_season <- function(dataTS, type){
  
  data_TS_x <-
    data.frame((sapply(dataTS["x", ], c)))
  data_TS_seasonal <-
    data.frame((sapply(dataTS["seasonal", ], c)))

 
  # perform seasonal adjustment
  dataAdj <- list()
  if (type == 'additive') {
    dataAdj <-
      data_TS_x - data_TS_seasonal

  } else{
    dataAdj <-
      data_TS_x / data_TS_seasonal
  }
  dataAdj <- as.data.frame(t(as.matrix(dataAdj)))
  #reshape to long format for further processing
  adjDataLong <-
    reshape(
      dataAdj,
      idvar = "ppID",
      varying = list(1:nrow(data_TS_x)),
      v.names = 'y',
      timevar = 'time' ,
      direction = "long"
    )
  adjDataLong <- adjDataLong[with(adjDataLong, order(ppID)), ]

  return(adjDataLong)
}


# want to difference per participant 

# helper function
myDiff <- function(y, degree){
  for(idx in 1:degree){
    y <- diff(y)
  }
  
  return(y)
}

preprocess.difference <- function(data, degree = 1 ,logBool = FALSE){
  
  # first check for logarithmic transformations
  # if we first difference the negative values will lead to NaN's
  if(logBool) {
    new_y <- data %>%
      group_by(ppID) %>%
      dplyr::select(y,ppID)%>%
      group_modify(~log(.x))
    
    data$y <-new_y$y
    
  }
  
  
  if(degree >0){
    dataNested <- data %>% 
      group_by(ppID) %>%
      nest()
    
    dataNew <- as_tibble(dataNested$ppID)
    dataNew<- dplyr::rename(dataNew, ppID = value)
    
    y <- map2(dataNested$data,degree,function(data,degree){myDiff(data$y, degree)})
    
    
    time <- map(dataNested$data,function(data){data$time[-c(length(data$time):(length(data$time)-degree+1))]})
    
    # create the new data set 
    data <- dataNew %>%
      mutate(
        y = y,
        time = time
      ) %>%
      unnest(cols = c(y, time))

  }
    

  return(data)
}


preprocess.splitData <- function(data, splitValues){
  splitVector <- data$time
  
  for(idx in 1:length(splitValues)){
    splitVal <- splitValues[idx]
    
    if(idx == 1) {
      splitVector[data$time <= splitVal] <- idx 
    }else{
      splitVector[data$time > splitValues[idx-1] & data$time <= splitVal] <- idx
    }
    
    if(idx == length(splitValues)){ # we have reached the last value
      splitVector[data$time > splitVal] <- idx+1
    }
  }
  data$split <- as.factor(splitVector)
  
  dataSplitted <- data %>%
    group_split(split)
  
  for(idx in 1:length(splitValues)){
    splitVal <- splitValues[idx]
    dataSplitted[[idx+1]]$time <- dataSplitted[[idx+1]]$time - (splitVal+1) 
    
  }
  
  
  return(dataSplitted)
  
}


