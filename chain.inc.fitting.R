# Chain model, complete information, defining category labels

# Load utility functions
source("utilities.R")
source("models.R")

# Load data
chain.inc.data <- read.csv("exp1_chain_incomplete.csv")
# str(chain.inc.data)

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

