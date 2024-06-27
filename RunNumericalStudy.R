source("./RequiredPackages.R")

# Table of all Parameters ------------------------------------------------------
simParameters <- expand.grid(K = c(8, 10, 12),
                             m = c(60, 80),
                             betas = c(paste("0.1 0.4"),
                                       paste("0.25 0.4"),
                                       paste("0.4 0.4")),
                             vars = c(paste("0.5 1"),
                                      paste("1 1"),
                                      paste("1 0.5")),
                             rho0 = c(paste("0.05 0.1"),
                                      paste("0.1 0.1"),
                                      paste("0.1 0.05")),
                             rho1 = c(0.005, 0.01, 0.05, 0.1, 0.2),
                             rho2 = c(0.1, 0.2, 0.5, 0.7, 0.9),
                             alpha = c(0.05),
                             r = c(1)
) %>%
  separate(betas, sep = " ", into = c("beta1", "beta2")) %>%
  separate(rho0, sep = " ", into = c("rho01", "rho02")) %>%
  separate(vars, sep = " ", into = c("varY1", "varY2")) %>%
  mutate_if(is.character, as.numeric) %>%
  rowid_to_column(., "Scenario")

# Run power calculations on all true parameters --------------------------------
powerTable <- simParameters %>%
  rowwise() %>%
  mutate('method1_bonf' = calc_pwr_pval_adj(K = K, m = m, alpha = alpha,
                                            beta1 = beta1, beta2 = beta2,
                                            varY1 = varY1, varY2 = varY2,
                                            rho01 = rho01, rho02 = rho02,
                                            rho2  = rho2, r = r)$'Final Power'[[1]],
         'method1_sidak' = calc_pwr_pval_adj(K = K, m = m, alpha = alpha,
                                             beta1 = beta1, beta2 = beta2,
                                             varY1 = varY1, varY2 = varY2,
                                             rho01 = rho01, rho02 = rho02,
                                             rho2  = rho2, r = r)$'Final Power'[[2]],
         'method1_dap' = calc_pwr_pval_adj(K = K, m = m, alpha = alpha,
                                           beta1 = beta1, beta2 = beta2,
                                           varY1 = varY1, varY2 = varY2,
                                           rho01 = rho01, rho02 = rho02,
                                           rho2  = rho2, r = r)$'Final Power'[[3]],
         'method2' = calc_pwr_comb_outcome(K = K, m = m, alpha = alpha,
                                           beta1 = beta1, beta2 = beta2,
                                           varY1 = varY1, varY2 = varY2,
                                           rho01 = rho01, rho02 = rho02,
                                           rho1 = rho1, rho2  = rho2, r = r),
         'method3' = calc_pwr_single_1dftest(K = K, m = m, alpha = alpha,
                                             beta1 = beta1, beta2 = beta2,
                                             varY1 = varY1, varY2 = varY2,
                                             rho01 = rho01, rho02 = rho02,
                                             rho1 = rho1, rho2  = rho2, r = r),
         'method4_chi2' = calc_pwr_disj_2dftest(K = K, m = m, alpha = alpha,
                                                beta1 = beta1, beta2 = beta2,
                                                varY1 = varY1, varY2 = varY2,
                                                rho01 = rho01, rho02 = rho02,
                                                rho1 = rho1, rho2  = rho2,
                                                r = r, dist = "Chi2"),
         'method4_F' = calc_pwr_disj_2dftest(K = K, m = m, alpha = alpha,
                                             beta1 = beta1, beta2 = beta2,
                                             varY1 = varY1, varY2 = varY2,
                                             rho01 = rho01, rho02 = rho02,
                                             rho1 = rho1, rho2  = rho2,
                                             r = r, dist = "F"),
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

View(powerTable)

# Visualizations ---------------------------------------------------------------
# Frequency of how many times a method is most powerful
methodList <- c("method1_bonf", "method1_sidak", "method1_dap",
                "method2", "method3", "method4_chi2", "method4_F",
                "method5_T", "method5_MVN")
mostPowerful <- powerTable %>%
  pivot_longer(cols = c("method1_bonf", "method1_sidak", "method1_dap",
                        "method2", "method3", "method4_chi2", "method4_F",
                        "method5_T", "method5_MVN"), names_to = "Method",
               values_to = "Power") %>%
  group_by(Scenario) %>%
  slice(which.max(Power)) %>%
  ungroup() %>%
  group_by(Method) %>%
  summarize(n = n()) %>%
  mutate(Method = factor(Method, levels = methodList)) %>%
  complete(Method = levels(Method), fill = list(n = 0))

# Histogram of power results for methods
powerLong <- powerTable %>%
  pivot_longer(cols = c("method1_bonf", "method1_sidak", "method1_dap",
                        "method2", "method3", "method4_chi2", "method4_F",
                        "method5_T", "method5_MVN"), names_to = "Method",
               values_to = "Power")

ggplot(data = powerLong, aes(Power)) + geom_histogram(bins = 35) +
  facet_wrap(~Method)

# Table of frequencies separated by rho1 and rho2
rho1rho2MostPowerful <- powerTable %>%
  pivot_longer(cols = c("method1_bonf", "method1_sidak", "method1_dap",
                        "method2", "method3", "method4_chi2", "method4_F",
                        "method5_T", "method5_MVN"), names_to = "Method",
               values_to = "Power") %>%
  group_by(Scenario, rho1, rho2) %>%
  slice(which.max(Power)) %>%
  ungroup() %>%
  group_by(Method, rho1, rho2) %>%
  dplyr::summarize(n = n()) %>%
  mutate(rho1 = as.factor(rho1), rho2 = paste("rho2 =", rho2))

ggplot(rho1rho2MostPowerful, aes(x = rho1, y = n, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") + facet_wrap(~rho2) +
  ylab("Frequency") + ggtitle("Number of input scenarios a \nmethod is found to have highest power, \nseparated by rho1 and rho2") + geom_text(aes(label = n), vjust = -1, position = position_dodge(width = 0.9))


#
