# Select Build, Build and reload to build and lode into the R-session.
#library(car)
#data(SLID,package = "car")
#SLID = SLID[complete.cases(SLID), ]
#ds = SLID
#colnames(ds)
#dim(ds)
#summary(ds)
#levels(ds$sex)
#levels(ds$language)
#ggpairs(ds)
#lm = lm(formula = wages ~ education+age+sex+language,data=SLID)
#summary(lm)

mylm <- function(formula, data = list(), contrasts = NULL, ...){
  # Extract model matrix & responses
  mf <- model.frame(formula = formula, data = data)
  X  <- model.matrix(attr(mf, "terms"), data = mf, contrasts.arg = contrasts)
  y  <- model.response(mf)
  terms <- attr(mf, "terms")


  # Add code here to calculate coefficients, residuals, fitted values, etc...
  # and store the results in the list est

  #Coefficients, residuals
  coefficients = solve(t(X)%*%X)%*%t(X)%*%y
  Y_hat = X%*%coefficients
  residuals = y-Y_hat
  n = nrow(mf)
  p = ncol(mf)
  k = p-1
  residual_quantiles = quantile(residuals)
  # min max quantiles of residuals

  #SSE, SST, R squared, residual standard error
  SSE = t(y-X%*%coefficients)%*%(y-X%*%coefficients)
  I = matrix(0,nrow=n,ncol=n)
  I[row(I)==col(I)] = 1
  J = matrix(1,nrow=n,ncol=n)
  SST = t(y)%*%(I-J/n)%*%y
  R_squared = 1-SSE/SST
  R_adjusted = 1 - (n-1)*(1-R_squared)/(n-p)
  residual_standard_error = sqrt(SSE/(n-p))
  linear_corr_coeff = cor(X,method="pearson")

  #z-test test of significance
  sigmasq = SSE/(n-p)
  sigmasq = sigmasq[1,1]
  covmatrix = sigmasq*solve(t(X)%*%X)
  std_coefficients = sqrt(diag(covmatrix))
  z = coefficients/std_coefficients
  pvalue_z = 2*pnorm(z,lower.tail = FALSE)

  #chisq test significance of the regression
  chisq = (SST-SSE)/(SSE/(n-p))
  pvalue_chisq = pchisq(chisq, df=k, lower.tail=FALSE)

  #ciritcal values of the tests
  alpha = 0.05
  critical_z = c(qnorm(alpha/2),qnorm(1-alpha/2))
  critical_chi = c(qchisq(alpha/2,df=k),qchisq(1-alpha/2,df=k))

  est <- list(terms = terms, model = mf)

  # Store call and formula used
  est$call <- match.call()
  est$formula <- formula
  est$X = X
  est$coefficients = coefficients
  est$residuals = residuals
  est$sums_sq_errors = SSE
  est$sums_sq_total = SST
  est$R_squared = R_squared
  est$R_adjusted = R_adjusted
  est$n = n
  est$p = p
  est$sigma_sq = sigmasq
  est$covmatrix = covmatrix
  est$std_error_coeff = std_coefficients
  est$z_stat = z
  est$pvalue_z = pvalue_z
  est$chisq_stat = chisq
  est$pvalue_chisq = pvalue_chisq
  est$residual_std_err = residual_standard_error
  est$residual_quantiles = residual_quantiles
  est$critical_z = critical_z
  est$critical_chi = critical_chi
  est$linear_corr_coeff = linear_corr_coeff




  # Set class name. This is very important!
  class(est) <- 'mylm'

  # Return the object with all results
  return(est)
}

print.mylm <- function(object, ...){
  # Code here is used when print(object) is used on objects of class "mylm"
  # Useful functions include cat, print.default and format
  cat('Info about object\n')
  cat('Call:\n')
  print(object$call)
  cat('\n')
  cat('Coefficients:\n')
  print(t(object$coefficients))
}

summary.mylm <- function(object, ...){
  # Code here is used when summary(object) is used on objects of class "mylm"
  # Useful functions include cat, print.default and format
  cat('Summary of object\n')
  cat('Call\n')
  print(object$call)
  cat('\n')
  cat('Residuals:\n')
  print.default(object$residual_quantiles,digits=3,quote=FALSE)
}

plot.mylm <- function(object, ...){
  # Code here is used when plot(object) is used on objects of class "mylm"

  library(ggplot2)
  # ggplot requires that the data is in a data.frame, this must be done here
  ggplot() + geom_point()

  # if you want the plot to look nice, you can e.g. use "labs" to add labels, and add colors in the geom_point-function

}

anova.mylm <- function(object, ...){
  # Code here is used when anova(object) is used on objects of class "mylm"

  # Components to test
  comp <- attr(object$terms, "term.labels")

  # Name of response
  response <- deparse(object$terms[[2]])

  # Fit the sequence of models
  txtFormula <- paste(response, "~", sep = "")
  model <- list()
  for(numComp in 1:length(comp)){
    if(numComp == 1){
      txtFormula <- paste(txtFormula, comp[numComp])
    }
    else{
      txtFormula <- paste(txtFormula, comp[numComp], sep = "+")
    }
    formula <- formula(txtFormula)
    model[[numComp]] <- lm(formula = formula, data = object$model)
  }

  # Print Analysis of Variance Table
  cat('Analysis of Variance Table\n')
  cat(c('Response: ', response, '\n'), sep = '')
  cat('          Df  Sum sq X2 value Pr(>X2)\n')
  for(numComp in 1:length(comp)){
    # Add code to print the line for each model tested
  }

  return(model)

}

