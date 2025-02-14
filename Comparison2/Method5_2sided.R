# 2-sided version of power calculation based on method 5: conjunctive IU test

calc_pwr_conj_test_2sided <- function(K,            # Number of clusters in treatment arm
                                      m,            # Individuals per cluster
                                      alpha = 0.05, # Significance level
                                      beta1,        # Effect for outcome 1
                                      beta2,        # Effect for outcome 2
                                      varY1,        # Variance for outcome 1
                                      varY2,        # Variance for outcome 2
                                      rho01,        # ICC for outcome 1
                                      rho02,        # ICC for outcome 2
                                      rho1,         # Inter-subject between-endpoint ICC
                                      rho2,         # Intra-subject between-endpoint ICC
                                      r = 1,        # Treatment allocation ratio
                                      cv = 0,       # If equal cluster size, cv=0
                                      deltas = c(0,0),
                                      dist = "MVN"   # Distribution to be used,
){

  # Check that input values are valid
  if(!is.numeric(c(K, m, alpha, beta1, beta2, varY1, varY2, rho01, rho02, rho1, rho2, r, cv))){
    stop("All input parameters must be numeric values.")
  }
  if(r <= 0){
    stop("Treatment allocation ratio should be a number greater than 0.")
  }
  if(K < 1 | K != round(K)){
    stop("'K' must be a positive whole number.")
  }
  if(m < 1 | m != round(m)){
    stop("'m' must be a positive whole number.")
  }

  # Helper functions requires ratio be defined as K1/K rather than K2/K1,
  # so define new ratio variable based on the one that was inputted by user
  r_alt <- 1/(r + 1)
  K_total <- ceiling(K/r_alt)
  betas = c(beta1, beta2)
  Q = 2

  # Variance of trt assignment
  sigmaz.square <- r_alt*(1 - r_alt)

  # Helper Function 1: Construct covariance matrix Sigma_E for Y_k -------------
  constrRiE <- function(rho01, rho2, Q, vars){
    rho0q <- diag(rho01)
    SigmaE_Matrix <- diag((1-rho0q)*vars)
    for(row in 1:Q){
      for(col in 1:Q){
        if(row != col){
          SigmaE_Matrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*(rho2[row,col]-rho01[row,col])
        }
      }
    }
    # Check for matrix positive definite
    if(min(eigen(SigmaE_Matrix)$values) <= 1e-08){
      warning("The resulting covariance matrix Sigma_E is not positive definite. Check the inputs for the correlation values.")
    }
    return(SigmaE_Matrix)
  }

  # Helper Function 2: Construct covariance matrix Sigma_phi for Y_k -----------
  constrRiP <- function(rho01, Q, vars){
    rho0q <- diag(rho01)
    SigmaP_Matrix <- diag(rho0q*vars)
    for(row in 1:Q){
      for(col in 1:Q){
        if(row != col){
          SigmaP_Matrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*rho01[row,col]
        }
      }
    }
    # Check for matrix positive definite
    if(min(eigen(SigmaP_Matrix)$values) <= 1e-08){
      warning("The resulting covariance matrix Sigma_phi is not positive definite. Check the input of rho01 and rho2.")
    }
    return(SigmaP_Matrix)
  }

  # Helper Function 3: Calculate covariance between betas ----------------------
  calCovbetas <- function(vars, rho01, rho2, cv, sigmaz.square, m, Q){
    sigmaE <- constrRiE(rho01, rho2, Q, vars)
    sigmaP <- constrRiP(rho01, Q, vars)
    tmp <- solve(diag(1,Q) - cv^2*(m*sigmaP %*% solve(sigmaE + m*sigmaP) %*% sigmaE %*% solve(sigmaE + m*sigmaP)))
    covMatrix <- 1/(m*sigmaz.square)*(sigmaE + m*sigmaP)%*%tmp
    covMatrix <- (covMatrix + t(covMatrix))/2  # symmerize the off-diagonal
    return(covMatrix)
  }

  # Helper Function 4: Calculate correlation between test statistics -----------
  calCorWks <- function(vars, rho01, rho2, sigmaz.square, cv, m, Q){
    top <- calCovbetas(vars, rho01, rho2, cv, sigmaz.square, m, Q)
    wCor <- diag(Q)
    for(row in 1:Q){
      for(col in 1:Q){
        if(row != col){
          wCor[row,col] <- top[row,col]/sqrt(top[row,row]*top[col,col])
        }
      }
    }
    return(wCor)
  }

  # Define necessary parameters
  sigmaks.sq <- diag(calCovbetas(vars = c(varY1, varY2),
                                 rho01 = matrix(c(rho01, rho1,
                                                  rho1, rho02),
                                                2, 2),
                                 rho2 = matrix(c(1, rho2,
                                                 rho2, 1),
                                               2, 2),
                                 cv = cv,
                                 sigmaz.square = sigmaz.square,
                                 m = m,
                                 Q = Q))
  meanVector <- sqrt(K_total)*(betas - deltas)/sqrt(sigmaks.sq)
  wCor <- calCorWks(vars = c(varY1, varY2),
                    rho01 = matrix(c(rho01, rho1,
                                     rho1, rho02),
                                   2, 2),
                    rho2 = matrix(c(1, rho2,
                                    rho2, 1),
                                  2, 2),
                    sigmaz.square = sigmaz.square,
                    cv = cv,
                    m = m,
                    Q = Q)

  if(dist == "MVN"){ # Using multivariate normal distribution
    # Calculate critical value and power
    # Here, we use the 2-sided critical value
    criticalValue <- qnorm(p = 1 - alpha/2,
                           mean = 0,
                           sd = 1) # lower.tail = TRUE is default

    # We want (|Z1| > z_{\alpha/2}) AND (|Z2| > z_{\alpha/2})

    # Pr(reject) = Pr((Z1, Z2) in Region1 u Region2 u Region3 u Region 4)
    #            = \sum_{i=1^4} Pr((Z_1, Z_2) in Region_i)

    # Region 1. Probability both outcomes > +criticalValue
    p_upper_upper <- pmvnorm(lower = c(criticalValue, criticalValue),
                             upper = c(Inf, Inf),
                             mean = meanVector,
                             corr = wCor)[1]

    # Region 2. Probability outcome1 > +criticalValue and outcome2 < -criticalValue
    p_upper_lower <- pmvnorm(lower = c(criticalValue, -Inf),
                             upper = c(Inf, -criticalValue),
                             mean = meanVector,
                             corr = wCor)[1]

    # Region 3. Probability outcome1 < -criticalValue and outcome2 > +criticalValue
    p_lower_upper <- pmvnorm(lower = c(-Inf, criticalValue),
                             upper = c(-criticalValue, Inf),
                             mean = meanVector,
                             corr = wCor)[1]

    # Region 4. Probability both outcomes < -criticalValue
    p_lower_lower <- pmvnorm(lower = c(-Inf, -Inf),
                             upper = c(-criticalValue, -criticalValue),
                             mean = meanVector,
                             corr = wCor)[1]

    # Sum them up to get the total probability that both endpoints
    # are outside of the interval [-criticalValue, +criticalValue]
    # i.e. 2-sided power
    power <- p_upper_upper + p_upper_lower + p_lower_upper + p_lower_lower

  } else{
    stop("Please choose a valid input parameter for 'dist', i.e. 'MVN' for Multivariate Normal Distribution.")
  }
  return(round(power, 4))
} # End calc_pwr_conj_test_2sided()
