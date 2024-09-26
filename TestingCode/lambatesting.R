source("./RequiredPackages.R")

# Table of all Parameters ------------------------------------------------------
numParameters <- expand.grid(K = c(6, 10),
                             m = 70,
                             betas = c(paste("0.1 0.4"),
                                       paste("0.25 0.4"),
                                       paste("0.4 0.4")),
                             vars = c(paste("0.5 10"),
                                      paste("10 0.5")),
                             rho0 = c(paste("0.05 0.1"),
                                      paste("0.1 0.1"),
                                      paste("0.1 0.05")),
                             rho1 = c(0.005, 0.01, 0.02, 0.05),
                             rho2 = c(0.1, 0.3, 0.5, 0.7),
                             alpha = c(0.05),
                             r = c(1)
) %>%
  separate(betas, sep = " ", into = c("beta1", "beta2")) %>%
  separate(rho0, sep = " ", into = c("rho01", "rho02")) %>%
  separate(vars, sep = " ", into = c("varY1", "varY2")) %>%
  mutate_if(is.character, as.numeric) %>%
  rowid_to_column(., "Scenario")

totalScenarios <- nrow(numParameters)

# Run power calculations on all true parameters --------------------------------
powerTable <- numParameters %>%
  rowwise() %>%
  mutate('method1_bonf' = calc_pwr_pval_adj(K = K, m = m, alpha = alpha,
                                            beta1 = beta1, beta2 = beta2,
                                            varY1 = varY1, varY2 = varY2,
                                            rho01 = rho01, rho02 = rho02,
                                            rho2  = rho2, r = r,
                                            dist = "Chi2")$'Final Power'[[1]],
         'method1_sidak' = calc_pwr_pval_adj(K = K, m = m, alpha = alpha,
                                             beta1 = beta1, beta2 = beta2,
                                             varY1 = varY1, varY2 = varY2,
                                             rho01 = rho01, rho02 = rho02,
                                             rho2  = rho2, r = r,
                                             dist = "Chi2")$'Final Power'[[2]],
         'method1_dap' = calc_pwr_pval_adj(K = K, m = m, alpha = alpha,
                                           beta1 = beta1, beta2 = beta2,
                                           varY1 = varY1, varY2 = varY2,
                                           rho01 = rho01, rho02 = rho02,
                                           rho2  = rho2, r = r,
                                           dist = "Chi2")$'Final Power'[[3]],
         'method2' = calc_pwr_comb_outcome(K = K, m = m, alpha = alpha,
                                           beta1 = beta1, beta2 = beta2,
                                           varY1 = varY1, varY2 = varY2,
                                           rho01 = rho01, rho02 = rho02,
                                           rho1 = rho1, rho2  = rho2, r = r,
                                           dist = "Chi2"),
         'method3' = calc_pwr_single_1dftest(K = K, m = m, alpha = alpha,
                                             beta1 = beta1, beta2 = beta2,
                                             varY1 = varY1, varY2 = varY2,
                                             rho01 = rho01, rho02 = rho02,
                                             rho1 = rho1, rho2  = rho2, r = r,
                                             dist = "Chi2"),
         'method4_Chi2' = calc_pwr_disj_2dftest(K = K, m = m, alpha = alpha,
                                                beta1 = beta1, beta2 = beta2,
                                                varY1 = varY1, varY2 = varY2,
                                                rho01 = rho01, rho02 = rho02,
                                                rho1 = rho1, rho2  = rho2,
                                                r = r, dist = "Chi2"),
         'method5_T' = calc_pwr_conj_test(K = K, m = m, alpha = alpha,
                                          beta1 = beta1, beta2 = beta2,
                                          varY1 = varY1, varY2 = varY2,
                                          rho01 = rho01, rho02 = rho02,
                                          rho1 = rho1, rho2  = rho2,
                                          r = r, dist = "T"),
         'method5_MVN' = calc_pwr_conj_test(K = K, m = m, alpha = alpha,
                                            beta1 = beta1, beta2 = beta2,
                                            varY1 = varY1, varY2 = varY2,
                                            rho01 = rho01, rho02 = rho02,
                                            rho1 = rho1, rho2  = rho2,
                                            r = r, dist = "MVN")) %>%
  mutate_at(vars(contains('method')), funs(.*100))

powerTable2 <- powerTable %>%
  mutate(VIF1 = 1 + (m-1)*rho01,
         VIF2 = 1 + (m-1)*rho02,
         VIF12 = rho2 + (m-1)*rho1) %>%
  mutate(method1_lambda1 = (beta1^2)/((1 + 1/r)*(varY1/(K*m))*(1 + (m-1)*rho01)),
         method1_lambda2 = (beta2^2)/((1 + 1/r)*(varY2/(K*m))*(1 + (m-1)*rho02)),
         method2_lambda = ((beta1 + beta2)^2)/((1 + 1/r)*((round(varY1 + varY2 + 2*rho2*sqrt(varY1)*sqrt(varY2), 4))/(K*m))*(1 + (m - 1)*((rho01*varY1 + rho02*varY2 + 2*rho1*sqrt(varY1*varY2))/
                                                                                                                                            (varY1 + varY2 + 2*rho2*sqrt(varY1*varY2))))),
         method3_lambda = ((sqrt((beta1^2)/((((1 + 1/r)*varY1)/(K*m))*(1 + (m - 1)*rho01))) + sqrt((beta2^2)/((((1 + 1/r)*varY2)/(K*m))*(1 + (m - 1)*rho02))))^2)/(2*(1 + (rho2 + (m - 1)*rho1)/sqrt((1 + (m - 1)*rho01)*(1 + (m - 1)*rho02)))),
         method4_lambda = calc_lambda_disj_2dftest(K = K, m = m, alpha = alpha,
                                                beta1 = beta1, beta2 = beta2,
                                                varY1 = varY1, varY2 = varY2,
                                                rho01 = rho01, rho02 = rho02,
                                                rho1 = rho1, rho2  = rho2,
                                                r = r, dist = "Chi2")) %>%
  mutate(lambda4GTlamba3 = ifelse('method4_lambda[, 1]' > method3_lambda, "Yes", "No"))
 # rename(method4_lambda)


View(powerTable2)
View(filter(powerTable2, method4_lambda < 0))




# mutate(pwr_method4 = round(1 - pchisq(qchisq(1 - alpha, df = 2, ncp = 0, lower.tail = TRUE, log.p = FALSE), df = 2, ncp = method4_lambda, lower.tail = TRUE), 4)*100)
#   View(dplyr::select(powerTable2, method4_Chi2, pwr_method4) %>%
#          filter(method4_Chi2 != pwr_method4))

  #   mutate(pwr_method3 = round(1 - pchisq(qchisq(p = alpha, df = 1, lower.tail = FALSE), ncp = method3_lambda, df = 1, lower.tail = TRUE), 4)*100)
  #
  # View(dplyr::select(powerTable2, method3, pwr_method3) %>%
  #        filter(method3 != pwr_method3))

#   mutate(pwr_method2 = round(1 - pchisq(qchisq(p = alpha, df = 1, lower.tail = FALSE), 1, ncp = method2_lambda, lower.tail = TRUE), 4)*100)
#
# View(dplyr::select(powerTable2, method2, pwr_method2) %>%
#        filter(method2 != pwr_method2))

  # mutate(pwr_method1_dap = round(min(c(1 - pchisq(qchisq(1 - (1 - (1 - alpha)^(1/(2^(1 - rho2)))), df = 1, ncp = 0), 1, ncp = method1_lambda1, lower.tail = TRUE),
  #                                  1 - pchisq(qchisq(1 - (1 - (1 - alpha)^(1/(2^(1 - rho2)))), df = 1, ncp = 0), 1, ncp = method1_lambda2, lower.tail = TRUE)))*100, 2))
# View(dplyr::select(powerTable2, method1_dap, pwr_method1_dap) %>%
#        filter(method1_dap != pwr_method1_dap))


counterExamples <- powerTable2 %>%
  mutate(p_greaterthan_2 = if_else(any(across(c(method1_bonf, method1_sidak, method1_dap), ~ . > method2)), "Yes", "No"),
         p_greaterthan_3 = if_else(any(across(c(method1_bonf, method1_sidak, method1_dap), ~ . > method3)), "Yes", "No"),
         p_greaterthan_4 = if_else(any(across(c(method1_bonf, method1_sidak, method1_dap), ~ . > method4_Chi2)), "Yes", "No")) %>%
  mutate(p_lambda_greaterthan_2 = if_else(min(method1_lambda1, method1_lambda2) >  method2_lambda, "Yes", "No"),
         p_lambda_greaterthan_3 = if_else(min(method1_lambda1, method1_lambda2) > method3_lambda, "Yes", "No"),
         p_lambda_greaterthan_4 = if_else(min(method1_lambda1, method1_lambda2) > method4_lambda, "Yes", "No"))

filter(counterExamples, VIF12 > VIF1)
filter(counterExamples, VIF12 > VIF2)

View(counterExamples)
