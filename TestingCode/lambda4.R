calc_lambda_disj_2dftest <- function(dist = "Chi2",# Distribution to base calculation from
                                  K,            # Number of clusters in treatment arm
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
                                  r = 1         # Treatment allocation ratio
){

  # Check that input values are valid
  if(!is.numeric(c(K, m, alpha, beta1, beta2, varY1, varY2, rho01, rho02, rho1, rho2, r))){
    stop("All input parameters must be numeric values (with the exception of 'dist').")
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

  # Define small dependent functions -------------------------------------------
  # Dependent Function 1: Calculates covariance between betas
  calCovbetas <- function(vars, rho01, rho2, cv, sigmaz.square, m, Q){
    sigmaE <- constrRiE(rho01, rho2, Q, vars)
    sigmaP <- constrRiP(rho01, Q, vars)
    tmp <- solve(diag(1,Q) - cv^2*(m*sigmaP %*% solve(sigmaE + m*sigmaP) %*% sigmaE %*% solve(sigmaE + m*sigmaP) ))
    covMatrix <- 1/(m*sigmaz.square)*(sigmaE + m*sigmaP)%*%tmp
    covMatrix <- (covMatrix + t(covMatrix))/2  # symmerize the off-diagonal
    return(covMatrix)
  }

  # Dependent Function 2: Constructs covariance matrix Sigma_E for Y_i
  constrRiE <- function(rho01, rho2, Q, vars){
    rho0k <- diag(rho01)
    SigmaE_Matrix <- diag((1 - rho0k)*vars)
    for(row in 1:Q) {
      for(col in 1:Q) {
        if(row != col){
          SigmaE_Matrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*(rho2[row,col] - rho01[row,col])
        }
      }
    }
    return(SigmaE_Matrix)
  }

  # Dependent Function 3: Constructs covariance matrix Sigma_phi for Y_i
  constrRiP <- function(rho01, Q, vars) {
    rho0k <- diag(rho01)
    SigmaP_Matrix <- diag(rho0k*vars)
    for(row in 1:Q) {
      for(col in 1:Q) {
        if(row != col){
          SigmaP_Matrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*rho01[row,col]
        }
      }
    }
    return(SigmaP_Matrix)
  }

  # Defining necessary parameters based on input values
  betas <- c(beta1, beta2) # Beta vector
  r_alt <- 1/(r + 1)
  Q <- 2 # Number of outcomes, could extend this to more than 2 in the future
  K_total <- ceiling(K/r_alt) # Total number of clusters
  vars <- c(varY1, varY2) # Vector of variances
  rho01_mat <- matrix(c(rho01, rho1, rho1, rho02), 2, 2)
  rho2_mat <- matrix(c(1, rho2, rho2, 1), 2, 2)
  clus_var <- 0 # cluster variation, placeholder for future extensions

  # Start calculating power
  sigmaz.square <- r_alt*(1 - r_alt) # Variance of treatment assignment
  omega <- calCovbetas(vars, rho01_mat, rho2_mat, clus_var, sigmaz.square, m, Q)
  tau <- K_total*t(betas) %*% solve(omega) %*% betas

  if(dist == "Chi2"){ # Using Chi2
    cv <- qchisq(1 - alpha, df = 2, ncp = 0, lower.tail = TRUE, log.p = FALSE)
    power <- round(1 - pchisq(cv, df = 2, ncp = tau, lower.tail = TRUE), 4)
  } else if(dist == "F"){ # Using F
    Fscore <- qf(1 - alpha, df1 = Q, df2 = K_total - 2*Q, ncp = 0,
                 lower.tail = TRUE, log.p = FALSE)
    power <- round(1 - pf(Fscore, df1 = Q, df2 = K_total - 2*Q, tau,
                          lower.tail = TRUE, log.p = FALSE), 4)
  } else{
    stop("Please choose a valid input parameter for 'dist', either 'Chi2' for Chi-Square or 'F' for F-distribution.")
  }

  return(tau)
} # End calc_pwr_disj_2dftest()
