# Al model functions go here
# .loss functions are for model fitting
# .pred functions are for geting the correlations with the actual ratings


# Chain
# Generative model chain (complete)
gm.chain.compl.loss <- function(pars, x, y, z, d, ratings) {
  # Calculates quadratic loss between gm.chain predictions and ratings...
  # for all exemplars for a single participant in the complete info condition.
  # chain model: X -> Y -> Z and D; D is isolated and Z's state is known.
  # Instead of modeling each individual exemplar, this function calculates p(E|C)...
  # and then multiplies the conditional probabilities to get each exemplar p.
  # This function uses bernoulli distributions to calculate in a single equation each effect's p.
  
  # @param gm_pars Numeric vector, includes the following:
    # cX Base rate of the cause attribute X.
    # bY Base rate of the effect attribute Y.
    # bZ Base rate of the effect attribute Z.
    # bD Base rate of isolated feature D.
    # mXY Strength of the causal link from X to Y.
    # mYZ Strength of the causal link from Y to Z.
    # gamma Power factor that provides a non-linear mapping of the model's predictions
      # onto the response scales.
  # @param x Numeric vector, indicates presense (1) absence (0) of feature X.
  # @param y Numeric vector, indicates presense (1) absence (0) of feature Y.
  # @param z Numeric vector, indicates presense (1) absence (0) of feature Z.
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
  # if Z is present... if Z is absent
  p_z_y <- (1 - (1 - mYZ) ^ y * (1 - bZ)) ^ z * ((1 - mYZ) ^ y * (1 - bZ)) ^ (1 - z)
  # if D is present... if D is absent
  p_d <- bD ^ d * (1 - bD) ^ (1 - d)
  
  # Exemplar probability is the product of all (conditional) feature probabilities
  exemplars_p <- p_x * p_y_x * p_z_y * p_d
  
  # Ratings are in scale 0 to 6
  result <- sum((6 * exemplars_p ^ gamma - ratings)^2)
  return(result)
}

gm.chain.compl.pred <- function(pars, x, y, z, d){
  # Extract parameters, gm_pars is a vector
  cX <- pars[1]
  bY <- pars[2]
  bZ <- pars[3]
  bD <- pars[4]
  mXY <- pars[5]
  mYZ <- pars[6]
  gamma <- pars[7]
  
  # Calculate attribute probabilities
  
  # if X is present... if X is absent
  p_x <- cX ^ x * (1 - cX) ^ (1 - x)
  # if Y is present... if Y is absent
  p_y_x <- (1 - (1 - mXY) ^ x * (1 - bY)) ^ y * ((1 - mXY) ^ x * (1 - bY)) ^ (1 - y)
  # if Z is present... if Z is absent
  p_z_y <- (1 - (1 - mYZ) ^ y * (1 - bZ)) ^ z * ((1 - mYZ) ^ y * (1 - bZ)) ^ (1 - z)
  # if D is present... if D is absent
  p_d <- bD ^ d * (1 - bD) ^ (1 - d)
  
  # Exemplar probability is the product of all (conditional) feature probabilities
  exemplars_p <- p_x * p_y_x * p_z_y * p_d
  
  # Ratings are in scale 0 to 6
  result <- 6 * exemplars_p^ gamma
  return(result)
}

# Piecemeal model chain (complete)
pm.chain.compl.loss <- function(pars, x, y, z, d, ratings){
  # Calculates quadratic loss between pm.chain predictions and ratings...
  # for all exemplars for a single participant in the complete info condition.
  # chain model: X -> Y -> Z and D; D is isolated and Z's state is known.
  # Assumes that all attributes are present in the prototype (A=1, B=1, C=1, D=1).
  
  # @param gm_pars Numeric vector, includes the following:
    # p_a Base probability of the cause A.
    # p_b: p(effect given cause) of B.
    # p_c: p(effect given cause) of C.
    # p_d: conceptual weight of D.
    # beta: parameter for linear transformation od similarities into ratings.
  # @param x Numeric vector, indicates presense (1) absence (0) of feature X.
  # @param y Numeric vector, indicates presense (1) absence (0) of feature Y.
  # @param z Numeric vector, indicates presense (1) absence (0) of feature Z.
  # @param d Numeric vector, indicates presense (1) absence (0) of feature D.
  # @return squared sum of the differences between exemplar's sim to prototype and ratings.
  
  # Extract parameters, gm_pars is a vector
  p_x <- pars[1]
  p_y <- pars[2]
  p_z <- pars[3]
  p_d <- pars[4]
  beta <- pars[5]
  
  # Correct conceptual weight of attributes used twice in the model
  p_y <- p_y/2.0
  
  # Discount inferences flowing from effect to cause (Eq. 1)
  # Not necessary here as I'm dealing with the chain model
  
  # Compute distances for each directly connected pair (Eq. 2)
  dist_xy <- (p_x / (p_x + p_y)) * abs(x - 1) + (p_y / (p_x + p_y)) * abs(y - 1)
  dist_yz <- (p_y / (p_y + p_z)) * abs(y - 1) + (p_z / (p_y + p_z)) * abs(z - 1)
  
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

pm.chain.compl.pred <- function(pars, x, y, z, d){
  # Extract parameters, gm_pars is a vector
  p_x <- pars[1]
  p_y <- pars[2]
  p_z <- pars[3]
  p_d <- pars[4]
  beta <- pars[5]
  
  # Correct conceptual weight of attributes used twice in the model
  p_y <- p_y/2.0
  
  # Discount inferences flowing from effect to cause (Eq. 1)
  # Not necessary here as I'm dealing with the chain model
  
  # Compute distances for each directly connected pair (Eq. 2)
  dist_xy <- (p_x / (p_x + p_y)) * abs(x - 1) + (p_y / (p_x + p_y)) * abs(y - 1)
  dist_yz <- (p_y / (p_y + p_z)) * abs(y - 1) + (p_z / (p_y + p_z)) * abs(z - 1)
  
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

# Generative model chain (incomplete) 
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

# Piecemeal model chain (incomplete)
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
  
  # Infer p_z from the state of its cause b
  # Note that the causal strength parameter is multiplied by p(Y), which in this case is known (so 0 or 1)
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
