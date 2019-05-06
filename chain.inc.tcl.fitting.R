# Chain model, complete information, tradictional category labels

# Load utility functions
source("utilities.R")

# Load data
chain.inc.data <- read.csv("exp4_chain_incomplete.csv")
# str(chain.inc.data)

# Pure inference model chain functions 
im.chain.inc.loss <- function(pars, x, y, z, d, ratings){
  # Calculates quadratic loss between im.chain predictions and ratings...
  # for all exemplars for a single participant in the incomplete info condition.
  # chain model: X -> Y -> Z? and D; D is isolated and Z's state is unknown.
  
  # @param pars Numeric vector, includes the following:
  # bZ Base rate of the effect attribute Z.
  # mYZ Strength of the causal link from Y to Z.
  # gamma Power factor that provides a non-linear mapping of the model's predictions...
  # onto the response scales.
  # @param x Numeric vector, indicates presense (1) absence (0) of feature X.
  # @param y Numeric vector, indicates presense (1) absence (0) of feature Y.
  # @param z Numeric vector, indicates presense (1) absence (0) of feature Z.
  # @param d Numeric vector, indicates presense (1) absence (0) of feature D.
  # @return squared sum of the differences between exemplar's probabilities and ratings.
  
  # Extract parameters, gm_pars is a vector
  bZ <- pars[1]
  mYZ <- pars[2]
  gamma <- pars[3]
  
  # Calculate feature probabilities
  p_z = 1 - (1 - bZ) * (1 - mYZ * y)
  
  # In this model feature Z's p is equal to exemplar p
  exemplars_p = p_z
  
  # Ratings are in scale 0 to 6
  result <- sum((6 * exemplars_p ^ gamma - ratings)^2)
  
  return(result)
}

im.chain.inc.pred <- function(pars, x, y, z, d){
  # Extract parameters, gm_pars is a vector
  bZ <- pars[1]
  mYZ <- pars[2]
  gamma <- pars[3]
  
  # Calculate feature probabilities
  p_z = 1 - (1 - bZ) * (1 - mYZ * y)
  
  # In this model feature Z's p is equal to exemplar p
  exemplars_p = p_z
  
  # Ratings are in scale 0 to 6
  result <- 6 * exemplars_p ^ gamma
  return(result)
  
}

# Generative model chain functions 
gm.chain.inc.loss <- function(pars, x, y, z, d, ratings) {
  # Calculates quadratic loss between im.chain predictions and ratings...
  # for all exemplars for a single participant in the incomplete info condition.
  # chain model: X -> Y -> Z? and D; D is isolated and Z's state is unknown.
  # This version infers p(Z) using the generative model equation and then...
  # follows the same procedure as the standard generative model.
  
  # @param gm_pars Numeric vector, includes the following:
  # cX Base rate of the cause attribute X.
  # bY Base rate of the effect attribute Y.
  # bZ Base rate of the effect attribute Z.
  # bD Base rate of isolated feature D.
  # mXY Strength of the causal link from X to Y.
  # mYZ Strength of the causal link from Y to Z.
  # gamma Power factor that provides a non-linear mapping of the model's predictions onto...
  # the response scales.
  # @param x Numeric vector, indicates presense (1) absence (0) of feature X.
  # @param y Numeric vector, indicates presense (1) absence (0) of feature Y.
  # @param d Numeric vector, indicates presense (1) absence (0) of feature D.
  # @return squared sum of the differences between exemplar's probabilities and ratings.
  
  # Extract parameters, gm_pars is a vector
  cX <- pars[1]
  bY <- pars[2]
  bZ <- pars[3]
  bD <- pars[4]
  mXY <- pars[5]
  mYZ <- pars[6]
  gamma <- pars[7]
  
  # Calculate feature probabilities
  # if X is present... if X is absent
  p_x <- cX ^ x * (1 - cX) ^ (1 - x)
  # if Y is present... if Y is absent
  p_y_x <- (1 - (1 - mXY) ^ x * (1 - bY)) ^ y * ((1 - mXY) ^ x * (1 - bY)) ^ (1 - y)
  # is only pk(Z=1) not pk(Z=0)
  p_z = 1 - (1 - bZ) * (1 - mYZ * y)  
  # if D is present... if D is absent
  p_d <- bD ^ d * (1 - bD) ^ (1 - d)
  
  # Exemplar probability is the product of all (conditional) feature probabilities
  exemplars_p <- p_x * p_y_x * p_z * p_d
  
  # Ratings are in scale 0 to 6
  result <- sum((6 * exemplars_p ^ gamma - ratings)^2)
  return(result)
}

gm.chain.inc.pred <- function(pars, x, y, z, d){
  # Extract parameters, gm_pars is a vector
  cX <- pars[1]
  bY <- pars[2]
  bZ <- pars[3]
  bD <- pars[4]
  mXY <- pars[5]
  mYZ <- pars[6]
  gamma <- pars[7]
  
  # Calculate feature probabilities
  # if X is present... if X is absent
  p_x <- cX ^ x * (1 - cX) ^ (1 - x)
  # if Y is present... if Y is absent
  p_y_x <- (1 - (1 - mXY) ^ x * (1 - bY)) ^ y * ((1 - mXY) ^ x * (1 - bY)) ^ (1 - y)
  # is only pk(Z=1) not pk(Z=0)
  p_z = 1 - (1 - bZ) * (1 - mYZ * y)  
  # if D is present... if D is absent
  p_d <- bD ^ d * (1 - bD) ^ (1 - d)
  
  # Exemplar probability is the product of all (conditional) feature probabilities
  exemplars_p <- p_x * p_y_x * p_z * p_d
  
  # Ratings are in scale 0 to 6
  result <- 6 * exemplars_p ^ gamma
  return(result)
}

# Piecemeal model chain functions
pm.chain.inc.loss <- function(pars, x, y, z, d, ratings){
  # Calculates quadratic loss between pm.chain predictions and ratings...
  # for all exemplars for a single participant in the incomplete info condition.
  # chain model: X -> Y -> Z? and D; D is isolated and Z's state is unknown.
  # Assumes that all attributes are present in the prototype (X=1, Y=1, Z=1, D=1).
  
  # This version infers p(Z) using the generative model equations...
  # and uses it as the state of Z. Otherwise it follows the same...
  # procedure as the standard PS model. 
  
  # Extract parameters, gm_pars is a vector
  p_x <- pars[1]
  p_y <- pars[2]
  p_d <- pars[3]
  bZ <- pars[4]
  mYZ <- pars[5]
  beta <- pars[6]
  
  # Correct conceptual weight of attributes used twice in the model
  p_y <- p_y/2.0
  
  # Infer p_c from the state of its cause b
  # Note that the causal strength parameter is multiplied by p(B), which in this case is known (so 0 or 1)
  p_z = 1 - (1 - bZ) * (1 - mYZ * y)  # is only pk(Z=1) not pk(Z=0)
  
  # Discount inferences flowing from effect to cause (Eq. 1)
  # Not necessary here as I'm dealing with the chain model
  
  # Compute distances for each directly connected pair (Eq. 2)
  dist_xy <- (p_x / (p_x + p_y)) * abs(x - 1) + (p_y / (p_x + p_y)) * abs(y - 1)
  dist_yz <- (p_y / (p_y + p_z)) * abs(y - 1) + (p_z / (p_y + p_z)) * abs(p_z - 1) # p_z takes the place of c
  
  # Compute distances for each isolated feature (Eq. 3)
  dist_d <- p_d * abs(d - 1)
  
  # Transform distances into similarities (Ep. 4)
  sim_xy <- exp(-dist_xy) # b parameter is set to 1
  sim_yz <- exp(-dist_yz)
  sim_d <- exp(-dist_d)
  
  # Compute overall similarity (eq. 5)
  sim <- rowMeans(cbind(sim_xy, sim_yz, sim_d))
  
  # Take the minimum sim
  min_sim = min(sim)
  
  # Ratings are in scale 0 to 6
  result <- sum(((sim - min_sim) / beta - ratings)^2)
  return(result)
}

pm.chain.inc.pred <- function(pars, x, y, z, d){
  # Extract parameters, gm_pars is a vector
  p_x <- pars[1]
  p_y <- pars[2]
  p_d <- pars[3]
  bZ <- pars[4]
  mYZ <- pars[5]
  beta <- pars[6]
  
  # Correct conceptual weight of attributes used twice in the model
  p_y <- p_y/2.0
  
  # Infer p_c from the state of its cause b
  # Note that the causal strength parameter is multiplied by p(B), which in this case is known (so 0 or 1)
  p_z = 1 - (1 - bZ) * (1 - mYZ * y)  # is only pk(Z=1) not pk(Z=0)
  
  # Discount inferences flowing from effect to cause (Eq. 1)
  # Not necessary here as I'm dealing with the chain model
  
  # Compute distances for each directly connected pair (Eq. 2)
  dist_xy <- (p_x / (p_x + p_y)) * abs(x - 1) + (p_y / (p_x + p_y)) * abs(y - 1)
  dist_yz <- (p_y / (p_y + p_z)) * abs(y - 1) + (p_z / (p_y + p_z)) * abs(p_z - 1) # p_z takes the place of c
  
  # Compute distances for each isolated feature (Eq. 3)
  dist_d <- p_d * abs(d - 1)
  
  # Transform distances into similarities (Ep. 4)
  sim_xy <- exp(-dist_xy) # b parameter is set to 1
  sim_yz <- exp(-dist_yz)
  sim_d <- exp(-dist_d)
  
  # Compute overall similarity (eq. 5)
  sim <- rowMeans(cbind(sim_xy, sim_yz, sim_d))
  
  # Take the minimum sim
  min_sim = min(sim)
  
  # Ratings are in scale 0 to 6
  result <- (sim - min_sim) / beta
  return(result)
}


# Fit generative model
gm.pars <- fit.model(df=chain.inc.data,
                     mymodel=gm.chain.inc.loss,
                     initpar=c(0.4, 0.4, 0.9, 0.3, 0.1, 0.1, 1.0),
                     low=c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                     up=c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0))
print("mean gm paramaters:")
print(round(gm.pars$means, 3))

# Correlate gm predictions to data
gm.corr = avg.model.data.corr(df=chain.inc.data,
                              pred.model=gm.chain.inc.pred,
                              part.pars=gm.pars$data)
print("mean gm corr:")
print(round(gm.corr, 2))

# Get AIC
gm.aic <- avg.aic(df=chain.inc.data,
                  pred.model=gm.chain.inc.pred,
                  part.pars=gm.pars$data) 
print("mean gm aic:")
print(round(gm.aic, 1))
writeLines("\n")


# Fit piecemeal model
pm.pars <- fit.model(df=chain.inc.data,
                     mymodel=pm.chain.inc.loss,
                     initpar=c(0.5, 0.5, 0.5, 0.5, 0.5, 1.0),
                     low=c(10^-10, 10^-10, 10^-10, 10^-10, 10^-10, 10^-10),
                     up=c(1.0, 1.0, 1.0, 1.0, 1.0, 10.0))
print("mean pm paramaters:")
print(round(pm.pars$means, 3))

# Correlate pm predictions to data
pm.corr = avg.model.data.corr(df=chain.inc.data, 
                              pred.model=pm.chain.inc.pred,
                              part.pars=pm.pars$data)
print("mean pm corr:")
print(round(pm.corr, 2))

# Get AIC
pm.aic <- avg.aic(df=chain.inc.data, 
                  pred.model=pm.chain.inc.pred,
                  part.pars=pm.pars$data) 
print("mean pm aic:")
print(round(pm.aic, 1))
writeLines("\n")
