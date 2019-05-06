
# Model fiting function
fit.model <- function(df, mymodel, initpar, low, up){
  # Fits a model to each participant of the dataframe using the optim function.
  
  # @param df Dataframe with all participants.
  # @param model Model to fit. Model has to return cuadratic loss between predictions and ratings.
  # @param initpar Initial guess for parameter values.
  # @param low Lower bound for parameter search.
  # @param up Upper bound for parameter search.
  # @return average parameter values aross participants.
  
  # Get unique participant values
  participants <- unique(df$participant)
  
  # Initialize matrix to save results
  all.vals <- matrix(NA, nrow=length(participants), ncol=length(initpar))
  
  # Iterate over participants
  for (i in participants){
    # Get data for fminsearch
    i.data <- subset(df, participant == i)
    i.x <- i.data$A  # chech if you need this as.vector
    i.y <- i.data$B
    i.z <- i.data$C
    i.d <-i.data$D
    i.r <- i.data$rating
    
    # Fit data
    # I tried fminsearch (pracma) but it was unstable
    i.vals = result <- optim(par=initpar,
                             fn=mymodel,
                             x=i.x,
                             y=i.y,
                             z=i.z,
                             d=i.d,
                             ratings=i.r,
                             method="L-BFGS-B",
                             lower=low,
                             upper=up,
                             control = list(maxit = 1000, factr=10^-10))
    
    # Append to list
    all.vals[i,] <- unlist(i.vals[1])
  }
  
  # Get average
  all.means <- colMeans(all.vals, na.rm = TRUE, dims = 1)
  
  return(list(means=all.means, data=all.vals))
}

# Model-data correlation function
avg.model.data.corr <- function(df, pred.model, part.pars){
  # Calculates average correlation over particioants between model preditions and ratings.
  
  # @param df Dataframe with all participants.
  # @param pred.model Model to fit. Model has to return the predicted rating, not the loss.
  # @param part.pars Estimated parameter values per participant (matrix).
  # @return average model-ratings correlation across participants.
  
  # Get unique participant values
  participants <- unique(df$participant)
  
  # Initialize vector to save results
  all.corr <- rep(NA, length=length(participants))
  
  # Iterate over participants
  for (i in participants){
    # Get data for each participant
    i.pars <- part.pars[i,]
    i.data <- subset(df, participant == i)
    i.x <- i.data$A
    i.y <- i.data$B
    i.z <- i.data$C
    i.d <- i.data$D
    i.r <- i.data$rating
    
    # Predict
    i.pred <- pred.model(i.pars, i.x, i.y, i.z, i.d)
    
    # Correlate
    i.corr <- cor(i.pred, i.r)

    # Append to vector
    all.corr[i] <- i.corr
  }
  
  # Get average
  corr.mean <- mean(all.corr, na.rm = TRUE)
  
  return(corr.mean)
}

# Akaike's information criterion
akaike.ic <- function(model.pred, ratings, n.pars){
  # Calculates the Akaike's information criterion for a single participant.
  
  # @param model.pre Model predictions (numeric vector).
  # @param ratings Criteria to calculate error against (numeric vector).
  # @param n.pars Number of model parameters
  # @return Scalar
  
  n.obs = length(ratings)
  sse = sum((model.pred - ratings)^2)
  # avoid division by 0, akaike.ic will resturn NA in this case
  # and this value will not be used to calculate the mean (see avg.aic)
  sse[sse==0] <- NA  
  aic = log(sse / n.obs) + 2 * (n.pars + 1)
  return(aic)
}

# AIC over participants
avg.aic <- function(df, pred.model, part.pars){
  
  # @param df Dataframe with all participants.
  # @param pred.model Model to fit. Model has to return the predicted rating, not the loss.
  # @param part.pars Estimated parameter values per participant (matrix).
  
  
  # Get unique participant values
  participants <- unique(df$participant)
  
  # Initialize vector to save results
  all.aic <- rep(NA, length=length(participants))
  
  
  # Iterate over participants
  for (i in participants){
    # Get data for each participant
    i.pars <- part.pars[i,]
    i.data <- subset(df, participant == i)
    i.x <- i.data$A
    i.y <- i.data$B
    i.z <- i.data$C
    i.d <- i.data$D
    i.r <- i.data$rating
    
    # Predict
    i.pred <- pred.model(i.pars, i.x, i.y, i.z, i.d)
    
    # Get AIC
    i.aic <- akaike.ic(model.pred=i.pred, ratings=i.r, n.pars=length(i.pars))
    
    # Append to vector
    all.aic[i] <- i.aic
  }
  
  # Get average
  aic.mean <- mean(all.aic, na.rm = TRUE)
  
  return(aic.mean)
}
