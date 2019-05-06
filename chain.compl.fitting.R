# Chain model, complete information, defining category labels

# Load utility functions
source("utilities.R")
source("models.R")
        
# Load data
exp1.data <- read.csv("exp1_chain_complete.csv")
str(exp1.data)

# Fit generative model
gm.pars <- fit.model(df=exp1.data,
                     mymodel=gm.chain.compl.loss,
                     initpar=c(0.4, 0.4, 0.9, 0.3, 0.1, 0.1, 0.3),
                     low=c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                     up=c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0))
print("mean gm paramaters:")
print(round(gm.pars$means, 3))

# Correlate gm predictions to data
gm.corr <- avg.model.data.corr(df=exp1.data, 
                              pred.model=gm.chain.compl.pred,
                              part.pars=gm.pars$data)
print("mean gm corr:")
print(round(gm.corr, 2))

# Get AIC
gm.aic <- avg.aic(df=exp1.data, 
                  pred.model=gm.chain.compl.pred,
                  part.pars=gm.pars$data) 
print("mean gm aic:")
print(round(gm.aic, 1))
writeLines("\n")


# Fit piecemeal model
pm.pars <- fit.model(df=exp1.data,
                          mymodel=pm.chain.compl.loss,
                          initpar=c(0.4, 0.4, 0.9, 0.3, 0.1),
                          low=c(10^-10, 10^-10, 10^-10, 10^-10, 10^-10),
                          up=c(1.0, 1.0, 1.0, 1.0, 10.0))
print("mean pm paramaters:")
print(round(pm.pars$means, 3))

# Correlate pm predictions to data
pm.corr = avg.model.data.corr(df=exp1.data, 
                              pred.model=pm.chain.compl.pred,
                              part.pars=pm.pars$data)
print("mean pm corr:")
print(round(pm.corr, 2))

# Get AIC
pm.aic <- avg.aic(df=exp1.data, 
                  pred.model=pm.chain.compl.pred,
                  part.pars=pm.pars$data) 
print("mean pm aic:")
print(round(pm.aic, 1))
writeLines("\n")










