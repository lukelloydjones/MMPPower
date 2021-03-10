# =========================================
# PowerBoot functions
# Author: Luke Lloyd-Jones adapted from 
#         Petra Kuhnert's MMP code.
# Date started: 12/10/2020
# Date updated: 24/11/2020
# =========================================



powerBoot2 <- function(fit, nSim, delta,
                       trend.nm   = "Date",
                       seas.1.nm  = "SEAS_C1",
                       seas.2.nm  = "SEAS_C2",
                       seas.trm   = 1,
                       mod.res    = 1)
{
  print(paste0("Preforming ", nSim, " simulations for % year-on-year change of ", delta))
  e    <- resid(fit)
  pval <- rep(0, nSim)
  y0   <- fit$coefficients["(Intercept)"]
  mod.ress <- e / sqrt(1 - hatvalues(fit))
  mod.res.cent <- mod.ress - mean(mod.ress)
  if (mod.res == 1)
  {
    e <- mod.res.cent
  }
  if (seas.trm == 1)
  {
    beta.seas.1 <- fit$coefficients[seas.1.nm]
    beta.seas.2 <- fit$coefficients[seas.2.nm]
  }
  for (i in 1:nSim) 
  {
    if ((i %% 100) == 0) {print(paste0("Finished iteration ", i))}
    # There is an ordering in time assumption here.  
    # Resample the residuals
    resid     <- sample(e, size = length(e), replace = TRUE)
    ystar     <- vector(mode = "numeric", length = length(e))
    m         <- fit$model
    trend     <- fit$model[, trend.nm]                         
    if (seas.trm == 1)
    {
      seas.c1   <- fit$model[, seas.1.nm]
      seas.c2   <- fit$model[, seas.2.nm]    
      epx.ystar <- y0 + trend * log(1 - delta) / 365.25 + beta.seas.1 * seas.c1 + beta.seas.2 * seas.c2
    } else {
      epx.ystar <- y0 + trend * log(1 - delta) / 365.25
    }

    # Is this imputation just for 0s or for missing data
    # Can distort things if pass lots of NAs
    newY    <- epx.ystar + resid 
    if (seas.trm == 1)
    {
      bootdat <- data.frame(e = resid, y = newY, trend = trend, seas.c1 = seas.c1, seas.c2 = seas.c2)
      mod     <- lm(y ~ trend + seas.c1 + seas.c2, data = bootdat)
    } else {
      bootdat <- data.frame(e = resid, y = newY, trend = trend)
      mod     <- lm(y ~ trend, data = bootdat)
    }
    tval    <- summary(mod)$coefficients["trend", 3]
    rdf     <- summary(mod)$df[2]
    pval[i] <- 2 * pt(-1 * abs(tval), rdf)
  }
  pval
}

powerBoot2Regional <- function(fit, nSim, delta, trend.nm   = "Date", mod.res    = 1)
{
  # Inputs
  #     - fit   - fitted linear model from lm
  #     - nSim  - number of Bootstrap resamples
  #     - delta - fractional year on year change
  #     - trend.nm - name of linear trend model component 
  #     - mod.res - logical for whether the modified residuals
  #                  should be used. Default is true.
  # Outputs
  #     - A set of pvals evaluating the significance of the
  #       two sided test for the linear trends coefficient for 
  #       the  nSim Bootstrap replicates
  print(paste0("Preforming ", nSim, " simulations for % year-on-year change of ", delta))
  e    <- resid(fit)
  pval <- rep(0, nSim)
  mod.ress     <- e / sqrt(1 - hatvalues(fit))
  mod.res.cent <- mod.ress - mean(mod.ress)
  if (mod.res == 1)
  {
    e <- mod.res.cent
  }
  # ----------------------
  # Using the model matrix
  # ----------------------
  #beta.est <- rep(0, nSim)
  for (i in 1:nSim) 
  {
    if ((i %% 100) == 0) {print(paste0("Finished iteration ", i))}
    # There is an ordering in time assumption here.  
    # Resample the residuals
    resid.new <- sample(e, size = length(e), replace = TRUE)
    
    # Replace the trend coefficients
    coef.new           <- fit$coefficients
    coef.new[trend.nm] <- log(1 - delta) / 365.25
    epx.ystar <- model.matrix(fit) %*% matrix(coef.new, nrow = length(coef.new), ncol = 1)
    # Is this imputation just for 0s or for missing data
    # Can distort things if pass lots of NAs
    newY       <- epx.ystar + resid.new 
    bootdat    <- data.frame(y = newY, fit$model[, -1])
    model.form <- paste(colnames(bootdat[,-1]), collapse="+")
    mod        <- lm(paste0("y ~ ", model.form), data = bootdat)
    tval       <- summary(mod)$coefficients[trend.nm, 3]
    rdf        <- summary(mod)$df[2]
    #beta.est[i] <- summary(mod)$coefficients[trend.nm, 1]
    pval[i]    <- 2 * pt(-1 * abs(tval), rdf)
  }
  return(pval)
  #return(list(pval = pval, beta.est = beta.est))
}

powerBoot2Regional2 <- function(fit, nSim, delta, trend.nm   = "Date", mod.res    = 1)
{
  print(paste0("Preforming ", nSim, " simulations for % year-on-year change of ", delta))
  e    <- resid(fit)
  pval <- rep(0, nSim)
  mod.ress     <- e / sqrt(1 - hatvalues(fit))
  mod.res.cent <- mod.ress - mean(mod.ress)
  if (mod.res == 1)
  {
    e <- mod.res.cent
  }
  # ----------------------
  # Using the model matrix
  # ----------------------
  #beta.est <- rep(0, nSim)
  for (i in 1:nSim) 
  {
    if ((i %% 100) == 0) {print(paste0("Finished iteration ", i))}
    # There is an ordering in time assumption here.  
    # Resample the residuals
    if (length(fit$model$SHORT_NAME) > 0)
    {
      e2 <- data.frame(resids = e, site = fit$model$SHORT_NAME)
      site.unq <- as.character(unique(fit$model$SHORT_NAME))
      for (k in seq(1, length(site.unq)))
      {
        sub.e       <- e2[which(e2$site == site.unq[k]), ]
        resid.new.k <- base::sample(sub.e$resids, size = length(sub.e$resids), replace = TRUE)
        if (k == 1)
        {
          resid.new <- resid.new.k
        } else {
          resid.new <- c(resid.new, resid.new.k)
        }
      }
    } else {
      resid.new <- base::sample(e, size = length(e), replace = TRUE)
    }
    # Replace the trend coefficients
    coef.new           <- fit$coefficients
    coef.new[trend.nm] <- log(1 - delta) / 365.25
    epx.ystar <- model.matrix(fit) %*% matrix(coef.new, nrow = length(coef.new), ncol = 1)
    # Is this imputation just for 0s or for missing data
    # Can distort things if pass lots of NAs
    newY       <- epx.ystar + resid.new 
    bootdat    <- data.frame(y = newY, fit$model[, -1])
    model.form <- paste(colnames(bootdat[,-1]), collapse="+")
    mod        <- lm(paste0("y ~ ", model.form), data = bootdat)
    tval       <- summary(mod)$coefficients[trend.nm, 3]
    rdf        <- summary(mod)$df[2]
    #beta.est[i] <- summary(mod)$coefficients[trend.nm, 1]
    pval[i]    <- 2 * pt(-1 * abs(tval), rdf)
  }
  return(pval)
  #return(list(pval = pval, beta.est = beta.est))
}


