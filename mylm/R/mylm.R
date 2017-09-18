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

  #Coefficients, residuals, df, residual quantiles
  coefficients <- solve(t(X)%*%X)%*%t(X)%*%y
  Y_hat <- X%*%coefficients
  residuals <- y-Y_hat
  n <- nrow(Y_hat)
  p <- nrow(coefficients)
  k <- p-1
  df <- n-p
  residual_quantiles <- quantile(residuals)

  #SSE, SST, R squared, residual standard error, Pearsons linear corr coefficient
  SSE <- t(residuals)%*%residuals
  I <- diag(n)
  J <- matrix(1,nrow=n,ncol=n)
  SST <- t(y)%*%(I-J/n)%*%y
  R_squared <- 1-SSE/SST
  residual_standard_error <- sqrt(SSE/(n-p))
  linear_corr_coeff <- cor(X,method="pearson")
  if (p <=2){
    linear_corr_coeff <- linear_corr_coeff[1,1]
  } else {
    linear_corr_coeff <- linear_corr_coeff[2:p,2:p]
  }

  #sigma^2, covariance matrix, z-test test of significance
  sigmasq <- SSE/(n-p)
  sigmasq <- sigmasq[1,1]
  covmatrix <- sigmasq*solve(t(X)%*%X)
  std_coefficients <- sqrt(diag(covmatrix))
  z <- coefficients/std_coefficients
  pvalue_z <- 2*pnorm(abs(z),lower.tail = FALSE)

  #chisq test significance of the regression
  chisq <- (SST-SSE)/(SSE/(n-p))
  pvalue_chisq <- pchisq(chisq, df=k, lower.tail = FALSE)

  #ciritcal values of the tests
  alpha <- 0.05
  critical_z <- c(qnorm(alpha/2),qnorm(1-alpha/2))
  critical_chi <- c(qchisq(alpha/2,df=k), qchisq(1-alpha/2, df=k))

  est <- list(terms = terms, model = mf)

  # Store call and formula used
  est$call <- match.call()
  est$formula <- formula
  est$X <- X
  est$y <- y
  est$Y_hat <- Y_hat
  est$coefficients <- coefficients
  est$residuals <- residuals
  est$n <- n
  est$p <- p
  est$k <- k
  est$df <- df
  est$residual_quantiles <- residual_quantiles

  est$SSE <- SSE
  est$SST <- SST
  est$R_squared <- R_squared
  est$residual_std_err <- residual_standard_error

  est$sigma_sq <- sigmasq
  est$covmatrix <- covmatrix
  est$std_error_coeff <- std_coefficients
  est$z_stat <- z
  est$pvalue_z <- pvalue_z
  est$chisq_stat <- chisq
  est$pvalue_chisq <- pvalue_chisq

  est$critical_z <- critical_z
  est$critical_chi <- critical_chi
  est$linear_corr_coeff <- linear_corr_coeff

  # Set class name. This is very important!
  class(est) <- 'mylm'

  # Return the object with all results
  return(est)
}

print.mylm <- function(object, ...){
  # Code here is used when print(object) is used on objects of class "mylm"
  cat('Call:\n')
  print(object$call)
  cat('\n')
  cat('Coefficients:\n')
  print(t(object$coefficients))
}

summary.mylm <- function(object, ...){
  # Code here is used when summary(object) is used on objects of class "mylm"
  cat('Call\n')
  print(object$call)
  cat('\nResiduals:\n')
  quantiles<-t(round(object$residual_quantiles,2))
  colnames(quantiles)=c("Min","1Q","Median","3Q","Max")
  print(quantiles)
  cat('\nCoefficients:\n')
  signif <- vector(mode="character",length=nrow(object$coefficients))

  for (i in 1:nrow(object$coefficients)){
    if (object$pvalue_z[i] < 0.001){
      signif[i] = "***"
    } else if (object$pvalue_z[i] < 0.01){
      signif[i] = "**"
    } else if (object$pvalue_z[i] < 0.05){
      signif[i] = "*"
    } else if (object$pvalue_z[i] < 0.1){
      signif[i] = "."
    } else if (object$pvalue_z[i] <= 1)
      signif[i] = " "
  }

  coeffmatrix <- matrix(c(round(object$coefficients,4),round(data.matrix(object$std_error_coeff),4),round(object$z_stat,2),object$pvalue_z),nrow<-nrow(object$coefficients))
  coeffmatrix <- data.frame(coeffmatrix)
  coeffmatrix <- cbind(coeffmatrix,signif)
  rownames(coeffmatrix)<-rownames(object$coefficients)
  colnames(coeffmatrix)<-c("Estimate","Std. Error","z value","Pr(>|z|)","   ")
  print(coeffmatrix)

  cat('\nSignif. codes: 0 ´***´ 0.001 ´**´ 0.01 ´*´ 0.05 ´.´ 0.1 ´ ´ 1\n')
  cat('\nResidual standard error: ', object$residual_std_err, 'on',object$df,'degrees of freedom\n')
  cat('Multiple R-squared: ', object$R_squared, '\n') #    Adjusted R-squared: ',object$R_adjusted,'\n')
  cat('Chisquared-statistics: ', object$chisq_stat, 'on', object$k, 'and', object$df, 'DF, ', '  p-value: ', object$pvalue_chisq)

}

plot.mylm <- function(object, ...){
  # Code here is used when plot(object) is used on objects of class "mylm"
  plotfit <- data.frame(fitvals = object$Y_hat, residuals = object$residuals,
                       obsvals = object$y)
  library(ggplot2)

  fittedvals_plot <- ggplot(plotfit,aes(fitvals, residuals)) + geom_point(pch = 21) + geom_hline(yintercept = 0, linetype ="dashed") + labs(x = "Fitted values", y = "Residuals", title = "Residuals vs Fitted values",subtitle = deparse(object$call))
  #+ geom_hline(yintercept = 0,linetype = "dashed") + geom_smooth(se = FALSE, col = "red", size = 0.5, method = "loess")
  obsvals_plot <- ggplot(plotfit,aes(obsvals, residuals)) + geom_point(pch = 21) + geom_hline(yintercept = 0, linetype ="dashed") + labs(x = "Observed values", y = "Residuals", title = "Residuals vs Observed values",subtitle = deparse(object$call))
  # if you want the plot to look nice, you can e.g. use "labs" to add labels, and add colors in the geom_point-function
  list_plot <- list(fittedvals_plot,obsvals_plot)

  return(list_plot)
}

anova.mylm <- function(object, ...){
  # Code here is used when anova(object) is used on objects of class "mylm"

  # Components to test
  comp <- attr(object$terms, "term.labels")

  # Name of response
  response <- deparse(object$terms[[2]])

  # Fit the sequence of models
  txtFormula <- paste(response, "~", sep = "")

  # Fit model with only intercept
  no = txtFormula
  no = paste(no,1)
  formulano = formula(no)
  model_nocoeff <- lm(formula = formulano,data=object$model)

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

  SSE <- vector()
  SSEdiff <- vector()
  SSE_nocoeff <- t(model_nocoeff$residuals)%*%model_nocoeff$residuals
  SSE[1] <- SSE_nocoeff
  Res.Df <- vector()
  Res.Df[1] <- model_nocoeff$df.residual
  Df <-vector()
  MeanSSE <- vector()
  X2_value <- vector()
  pvalue_chisqX2 <- vector()
  # Print Analysis of Variance Table
  cat('Analysis of Variance Table\n')
  cat(c('Response: ', response, '\n'), sep = '')
  #cat('          Df  sq X2 value Pr(>X2)\n')

  for(numComp in 1:length(comp)){
    SSE[numComp+1] <- t(model[[numComp]]$residuals)%*%model[[numComp]]$residuals
    SSEdiff[numComp] <- SSE[numComp]-SSE[numComp+1]
    Res.Df[numComp+1] <- model[[numComp]]$df.residual
    Df[numComp] <- Res.Df[numComp] - Res.Df[numComp+1]
    MeanSSE[numComp] <- SSEdiff[numComp]/Df[numComp]
  }
  for(numComp2 in 1:length(comp)){
    X2_value[numComp2] <- SSEdiff[numComp2]/(SSE[length(comp)]/Res.Df[length(comp)])
    pvalue_chisqX2[numComp2] <- pchisq(X2_value[numComp2], df=Df[numComp2], lower.tail = FALSE)
  }
  signif2 <- vector(mode="character",length=length(comp))
  for (i in 1:length(pvalue_chisqX2)){
    if (pvalue_chisqX2[i] < 0.001){
      signif2[i] = "***"
    } else if (pvalue_chisqX2[i] < 0.01){
      signif2[i] = "**"
    } else if (pvalue_chisqX2[i] < 0.05){
      signif2[i] = "*"
    } else if (pvalue_chisqX2[i] < 0.1){
      signif2[i] = "."
    } else if (pvalue_chisqX2[i] <= 1)
      signif2[i] = " "
  }
  anovamatrix <- cbind(Df,round(SSEdiff),round(X2_value,2),pvalue_chisqX2)
  anovamatrix <- data.frame(anovamatrix)
  anovamatrix <- cbind(anovamatrix,signif2)
  rownames(anovamatrix) <- comp
  colnames(anovamatrix) <- c("Df","Sum sq","X2 value","Pr(>X2)","   ")
  print(anovamatrix)
  cat("Residuals ", tail(Res.Df,n=1), tail(SSE,n=1), round(tail(SSE,n=1)/tail(Res.Df,n=1)))

  #return(model)

}

