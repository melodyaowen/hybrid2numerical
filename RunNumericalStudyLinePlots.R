source("./RequiredPackages.R")

# Varying Rho1 -----------------------------------------------------------------
# Table of all Parameters
numParameters_Rho1 <- expand.grid(K = 8,
                             m = 50,
                             betas = c(paste("0.4 0.4")),
                             vars = c(paste("1 1")),
                             rho0 = c(paste("0.1 0.1")),
                             rho1 = seq(0.001, 0.1, by = 0.001),
                             rho2 = c(0.1),
                             alpha = c(0.05),
                             r = c(1)
) %>%
  separate(betas, sep = " ", into = c("beta1", "beta2")) %>%
  separate(rho0, sep = " ", into = c("rho01", "rho02")) %>%
  separate(vars, sep = " ", into = c("varY1", "varY2")) %>%
  mutate_if(is.character, as.numeric) %>%
  rowid_to_column(., "Scenario")

# Run power calculations on all true parameters
powerTable_Rho1 <- numParameters_Rho1 %>%
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

View(powerTable_Rho1)

plotData_Rho1 <- powerTable_Rho1 %>%
  pivot_longer(cols = c("method1_bonf", "method1_sidak", "method1_dap",
                        "method2", "method3", "method4_Chi2",
                        "method5_T", "method5_MVN"), names_to = "Method",
               values_to = "Power") %>%
  mutate(Method = fct_recode(Method,
                             "1. P-Value Adjustment (Bonferroni)" = "method1_bonf",
                             "1. P-Value Adjustment (Sidak)" = "method1_sidak",
                             "1. P-Value Adjustment (D/AP)" = "method1_dap",
                             "2. Combined Outcomes" = "method2",
                             "3. Single 1-DF Test" = "method3",
                             "4. Disjunctive 2-DF" = "method4_Chi2",
                             "5. Conjunctive IU Test (t-Dist)" = "method5_T",
                             "5. Conjunctive IU Test (MVN-Dist)" = "method5_MVN"))

# Line graphs of Rho1
rho1_graph <- ggplot(plotData_Rho1, aes(x = rho1, y = Power,
                                        group = Method, color = Method)) +
  geom_line(linewidth = 1) +
  theme(text = element_text(size = 15)) + xlab(expression(rho[1]))
rho1_graph

ggsave(filename = "~/Desktop/2. Hybrid Software and Simulation/Word Drafts/Figures/Line Graphs/rho1_plot.png",
       plot = rho1_graph,
       width = 10, height = 6, dpi = 300)

# Varying Rho2 -----------------------------------------------------------------
# Table of all Parameters
numParameters_Rho2 <- expand.grid(K = 8,
                                  m = 50,
                                  betas = c(paste("0.4 0.4")),
                                  vars = c(paste("1 1")),
                                  rho0 = c(paste("0.1 0.1")),
                                  rho1 = c(0.01),
                                  rho2 = seq(0.1, 0.9, by = 0.01),
                                  alpha = c(0.05),
                                  r = c(1)
) %>%
  separate(betas, sep = " ", into = c("beta1", "beta2")) %>%
  separate(rho0, sep = " ", into = c("rho01", "rho02")) %>%
  separate(vars, sep = " ", into = c("varY1", "varY2")) %>%
  mutate_if(is.character, as.numeric) %>%
  rowid_to_column(., "Scenario")

# Run power calculations on all true parameters
powerTable_Rho2 <- numParameters_Rho2 %>%
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

View(powerTable_Rho2)

plotData_Rho2 <- powerTable_Rho2 %>%
  pivot_longer(cols = c("method1_bonf", "method1_sidak", "method1_dap",
                        "method2", "method3", "method4_Chi2",
                        "method5_T", "method5_MVN"), names_to = "Method",
               values_to = "Power") %>%
  mutate(Method = fct_recode(Method,
                             "1. P-Value Adjustment (Bonferroni)" = "method1_bonf",
                             "1. P-Value Adjustment (Sidak)" = "method1_sidak",
                             "1. P-Value Adjustment (D/AP)" = "method1_dap",
                             "2. Combined Outcomes" = "method2",
                             "3. Single 1-DF Test" = "method3",
                             "4. Disjunctive 2-DF" = "method4_Chi2",
                             "5. Conjunctive IU Test (t-Dist)" = "method5_T",
                             "5. Conjunctive IU Test (MVN-Dist)" = "method5_MVN"))

# Line graphs of Rho2
rho2_graph <- ggplot(plotData_Rho2, aes(x = rho2, y = Power,
                                        group = Method, color = Method)) +
  geom_line(linewidth = 1) +
  theme(text = element_text(size = 15)) + xlab(expression(rho[2]))

rho2_graph

ggsave(filename = "~/Desktop/2. Hybrid Software and Simulation/Word Drafts/Figures/Line Graphs/rho2_plot.png",
       plot = rho2_graph,
       width = 10, height = 6, dpi = 300)

# Varying Betas -----------------------------------------------------------------
# Table of all Parameters
numParameters_Betas <- expand.grid(K = 8,
                                  m = 50,
                                  betas = paste("0.4", seq(0.1, 0.8, by = 0.01)),
                                  vars = c(paste("1 1")),
                                  rho0 = c(paste("0.1 0.1")),
                                  rho1 = c(0.01),
                                  rho2 = c(0.1),
                                  alpha = c(0.05),
                                  r = c(1)
) %>%
  separate(betas, sep = " ", into = c("beta1", "beta2")) %>%
  separate(rho0, sep = " ", into = c("rho01", "rho02")) %>%
  separate(vars, sep = " ", into = c("varY1", "varY2")) %>%
  mutate_if(is.character, as.numeric) %>%
  rowid_to_column(., "Scenario")

# Run power calculations on all true parameters
powerTable_Betas <- numParameters_Betas %>%
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

View(powerTable_Betas)

plotData_Betas <- powerTable_Betas %>%
  mutate(BetaRatio = beta2/beta1) %>%
  pivot_longer(cols = c("method1_bonf", "method1_sidak", "method1_dap",
                        "method2", "method3", "method4_Chi2",
                        "method5_T", "method5_MVN"), names_to = "Method",
               values_to = "Power") %>%
  mutate(Method = fct_recode(Method,
                             "1. P-Value Adjustment (Bonferroni)" = "method1_bonf",
                             "1. P-Value Adjustment (Sidak)" = "method1_sidak",
                             "1. P-Value Adjustment (D/AP)" = "method1_dap",
                             "2. Combined Outcomes" = "method2",
                             "3. Single 1-DF Test" = "method3",
                             "4. Disjunctive 2-DF" = "method4_Chi2",
                             "5. Conjunctive IU Test (t-Dist)" = "method5_T",
                             "5. Conjunctive IU Test (MVN-Dist)" = "method5_MVN"))

# Line graphs of Beta ratios
betas_graph <- ggplot(plotData_Betas, aes(x = BetaRatio, y = Power,
                          group = Method, color = Method)) +
  geom_line(linewidth = 1) +
  theme(text = element_text(size = 15)) + xlab(expression(beta[2]/beta[1]))

betas_graph

ggsave(filename = "~/Desktop/2. Hybrid Software and Simulation/Word Drafts/Figures/Line Graphs/betas_plot.png",
       plot = betas_graph,
       width = 10, height = 6, dpi = 300)

# Varying ICCs -----------------------------------------------------------------
# Table of all Parameters
numParameters_ICCs <- expand.grid(K = 8,
                                   m = 50,
                                   betas = c(paste("0.4 0.4")),
                                   vars = c(paste("1 1")),
                                   rho0 = paste("0.1", seq(0.01, 0.5, by = 0.01)),
                                   rho1 = c(0.01),
                                   rho2 = c(0.1),
                                   alpha = c(0.05),
                                   r = c(1)
) %>%
  separate(betas, sep = " ", into = c("beta1", "beta2")) %>%
  separate(rho0, sep = " ", into = c("rho01", "rho02")) %>%
  separate(vars, sep = " ", into = c("varY1", "varY2")) %>%
  mutate_if(is.character, as.numeric) %>%
  rowid_to_column(., "Scenario")

# Run power calculations on all true parameters
powerTable_ICCs <- numParameters_ICCs %>%
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

View(powerTable_ICCs)

plotData_ICCs <- powerTable_ICCs %>%
  mutate(ICCsRatio = rho02/rho01) %>%
  pivot_longer(cols = c("method1_bonf", "method1_sidak", "method1_dap",
                        "method2", "method3", "method4_Chi2",
                        "method5_T", "method5_MVN"), names_to = "Method",
               values_to = "Power") %>%
  mutate(Method = fct_recode(Method,
                             "1. P-Value Adjustment (Bonferroni)" = "method1_bonf",
                             "1. P-Value Adjustment (Sidak)" = "method1_sidak",
                             "1. P-Value Adjustment (D/AP)" = "method1_dap",
                             "2. Combined Outcomes" = "method2",
                             "3. Single 1-DF Test" = "method3",
                             "4. Disjunctive 2-DF" = "method4_Chi2",
                             "5. Conjunctive IU Test (t-Dist)" = "method5_T",
                             "5. Conjunctive IU Test (MVN-Dist)" = "method5_MVN"))

# Line graphs of ICC ratios
ICCs_graph <- ggplot(plotData_ICCs, aes(x = ICCsRatio, y = Power,
                           group = Method, color = Method)) +
  geom_line(linewidth = 1) +
  theme(text = element_text(size = 15)) +
  xlab(expression(rho[0]^{(2)} / rho[0]^{(1)}))

ICCs_graph

ggsave(filename = "~/Desktop/2. Hybrid Software and Simulation/Word Drafts/Figures/Line Graphs/ICCs_plot.png",
       plot = ICCs_graph,
       width = 10, height = 6, dpi = 300)

# Varying Vars -----------------------------------------------------------------
# Table of all Parameters
numParameters_Vars <- expand.grid(K = 8,
                                  m = 50,
                                  betas = c(paste("0.4 0.4")),
                                  vars = paste("1", seq(0.01, 2, by = 0.01)),
                                  rho0 = paste("0.1 0.1"),
                                  rho1 = c(0.01),
                                  rho2 = c(0.1),
                                  alpha = c(0.05),
                                  r = c(1)
) %>%
  separate(betas, sep = " ", into = c("beta1", "beta2")) %>%
  separate(rho0, sep = " ", into = c("rho01", "rho02")) %>%
  separate(vars, sep = " ", into = c("varY1", "varY2")) %>%
  mutate_if(is.character, as.numeric) %>%
  rowid_to_column(., "Scenario")

# Run power calculations on all true parameters
powerTable_Vars <- numParameters_Vars %>%
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

View(powerTable_Vars)

plotData_Vars <- powerTable_Vars %>%
  mutate(VarsRatio = varY2/varY1) %>%
  pivot_longer(cols = c("method1_bonf", "method1_sidak", "method1_dap",
                        "method2", "method3", "method4_Chi2",
                        "method5_T", "method5_MVN"), names_to = "Method",
               values_to = "Power") %>%
  mutate(Method = fct_recode(Method,
                             "1. P-Value Adjustment (Bonferroni)" = "method1_bonf",
                             "1. P-Value Adjustment (Sidak)" = "method1_sidak",
                             "1. P-Value Adjustment (D/AP)" = "method1_dap",
                             "2. Combined Outcomes" = "method2",
                             "3. Single 1-DF Test" = "method3",
                             "4. Disjunctive 2-DF" = "method4_Chi2",
                             "5. Conjunctive IU Test (t-Dist)" = "method5_T",
                             "5. Conjunctive IU Test (MVN-Dist)" = "method5_MVN"))

# Line graphs of Variance ratios
vars_graph <- ggplot(plotData_Vars, aes(x = VarsRatio, y = Power,
                          group = Method, color = Method)) +
  geom_line(linewidth = 1) +
  theme(text = element_text(size = 15)) +
  xlab(expression(sigma[2]^{2} / sigma[1]^{2}))

vars_graph

ggsave(filename = "~/Desktop/2. Hybrid Software and Simulation/Word Drafts/Figures/Line Graphs/vars_plot.png",
       plot = vars_graph,
       width = 10, height = 6, dpi = 300)

# Varying K -----------------------------------------------------------------
# Table of all Parameters
numParameters_K <- expand.grid(K = c(4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24),
                                  m = 50,
                                  betas = c(paste("0.4 0.4")),
                                  vars = paste("1 1"),
                                  rho0 = paste("0.1 0.1"),
                                  rho1 = c(0.01),
                                  rho2 = c(0.1),
                                  alpha = c(0.05),
                                  r = c(1)
) %>%
  separate(betas, sep = " ", into = c("beta1", "beta2")) %>%
  separate(rho0, sep = " ", into = c("rho01", "rho02")) %>%
  separate(vars, sep = " ", into = c("varY1", "varY2")) %>%
  mutate_if(is.character, as.numeric) %>%
  rowid_to_column(., "Scenario")

# Run power calculations on all true parameters
powerTable_K <- numParameters_K %>%
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

View(powerTable_K)

plotData_K <- powerTable_K %>%
  pivot_longer(cols = c("method1_bonf", "method1_sidak", "method1_dap",
                        "method2", "method3", "method4_Chi2",
                        "method5_T", "method5_MVN"), names_to = "Method",
               values_to = "Power") %>%
  mutate(Method = fct_recode(Method,
                             "1. P-Value Adjustment (Bonferroni)" = "method1_bonf",
                             "1. P-Value Adjustment (Sidak)" = "method1_sidak",
                             "1. P-Value Adjustment (D/AP)" = "method1_dap",
                             "2. Combined Outcomes" = "method2",
                             "3. Single 1-DF Test" = "method3",
                             "4. Disjunctive 2-DF" = "method4_Chi2",
                             "5. Conjunctive IU Test (t-Dist)" = "method5_T",
                             "5. Conjunctive IU Test (MVN-Dist)" = "method5_MVN"))

# Line graphs of K
K_graph <- ggplot(plotData_K, aes(x = K, y = Power,
                          group = Method, color = Method)) +
  geom_line(linewidth = 1) +
  theme(text = element_text(size = 15))

K_graph

ggsave(filename = "~/Desktop/2. Hybrid Software and Simulation/Word Drafts/Figures/Line Graphs/K_plot.png",
       plot = K_graph,
       width = 10, height = 6, dpi = 300)


# Varying m -----------------------------------------------------------------
# Table of all Parameters
numParameters_m <- expand.grid(K = 8,
                               m = seq(20, 500, by = 5),
                               betas = c(paste("0.4 0.4")),
                               vars = paste("1 1"),
                               rho0 = paste("0.1 0.1"),
                               rho1 = c(0.01),
                               rho2 = c(0.1),
                               alpha = c(0.05),
                               r = c(1)
) %>%
  separate(betas, sep = " ", into = c("beta1", "beta2")) %>%
  separate(rho0, sep = " ", into = c("rho01", "rho02")) %>%
  separate(vars, sep = " ", into = c("varY1", "varY2")) %>%
  mutate_if(is.character, as.numeric) %>%
  rowid_to_column(., "Scenario")

# Run power calculations on all true parameters
powerTable_m <- numParameters_m %>%
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

View(powerTable_m)

plotData_m <- powerTable_m %>%
  pivot_longer(cols = c("method1_bonf", "method1_sidak", "method1_dap",
                        "method2", "method3", "method4_Chi2",
                        "method5_T", "method5_MVN"), names_to = "Method",
               values_to = "Power") %>%
  mutate(Method = fct_recode(Method,
                             "1. P-Value Adjustment (Bonferroni)" = "method1_bonf",
                             "1. P-Value Adjustment (Sidak)" = "method1_sidak",
                             "1. P-Value Adjustment (D/AP)" = "method1_dap",
                             "2. Combined Outcomes" = "method2",
                             "3. Single 1-DF Test" = "method3",
                             "4. Disjunctive 2-DF" = "method4_Chi2",
                             "5. Conjunctive IU Test (t-Dist)" = "method5_T",
                             "5. Conjunctive IU Test (MVN-Dist)" = "method5_MVN"))

# Line graphs of m
m_graph <- ggplot(plotData_m, aes(x = m, y = Power,
                       group = Method, color = Method)) +
  geom_line(linewidth = 1) +
  theme(text = element_text(size = 15))

m_graph

ggsave(filename = "~/Desktop/2. Hybrid Software and Simulation/Word Drafts/Figures/Line Graphs/m_plot.png",
       plot = m_graph,
       width = 10, height = 6, dpi = 300)


