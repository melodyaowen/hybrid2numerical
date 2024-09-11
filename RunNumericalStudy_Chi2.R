source("./RequiredPackages.R")

# Table of all Parameters ------------------------------------------------------
numParameters <- expand.grid(K = c(6, 8, 10),
                             m = c(50, 70),
                             betas = c(paste("0.1 0.4"),
                                       paste("0.25 0.4"),
                                       paste("0.4 0.4")),
                             vars = c(paste("0.5 1"),
                                      paste("1 1"),
                                      paste("1 0.5")),
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

View(powerTable)
View(numParameters)
nrow(numParameters)

# Check for cases where power is 100
scenarios100 <- powerTable %>%
  dplyr::filter(if_any(c(method1_bonf, method1_sidak, method1_dap,
                         method2, method3,
                         method4_Chi2, method5_T, method5_MVN), ~ . == 100))

nrow(scenarios100)
View(scenarios100)

# Visualizations ---------------------------------------------------------------
# Frequency of how many times a method is most powerful
methodList <- c("method1_bonf", "method1_sidak", "method1_dap",
                "method2", "method3", "method4_Chi2",
                "method5_T", "method5_MVN")
# mostPowerful <- powerTable %>%
#   pivot_longer(cols = c("method1_bonf", "method1_sidak", "method1_dap",
#                         "method2", "method3", "method4_Chi2",
#                         "method5_T", "method5_MVN"), names_to = "Method",
#                values_to = "Power") %>%
#   group_by(Scenario) %>%
#   slice(which.max(Power)) %>%
#   ungroup() %>%
#   group_by(Method) %>%
#   summarize(n = n()) %>%
#   mutate(Method = factor(Method, levels = methodList)) %>%
#   complete(Method = levels(Method), fill = list(n = 0)) %>%
#   mutate(Percent = paste0(round((n/nrow(numParameters))*100, 2), "%")) %>%
#   mutate("Method" = fct_recode(Method,
#                                "1. P-Value Adjustment (Bonferroni)" = "method1_bonf",
#                                "1. P-Value Adjustment (Sidak)" = "method1_sidak",
#                                "1. P-Value Adjustment (D/AP)" = "method1_dap",
#                                "2. Combined Outcomes" = "method2",
#                                "3. Single 1-DF Test" = "method3",
#                                "4. Disjunctive 2-DF" = "method4_Chi2",
#                                "5. Conjunctive IU Test (t-Dist)" = "method5_T",
#                                "5. Conjunctive IU Test (MVN-Dist)" = "method5_MVN"))
# View(mostPowerful)

mostPowerful <- powerTable %>%
  pivot_longer(cols = c("method1_bonf", "method1_sidak", "method1_dap",
                        "method2", "method3", "method4_Chi2",
                        "method5_T", "method5_MVN"), names_to = "Method",
               values_to = "Power") %>%
  group_by(Scenario) %>%
  filter(Power == max(Power)) %>%
  mutate(unique_id = row_number()) %>%
  pivot_wider(names_from = unique_id, values_from = c("Method", "Power")) %>%
  group_by(Method_1, Method_2) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(Method_1 = factor(Method_1, levels = methodList)) %>%
  complete(Method_1 = levels(Method_1), fill = list(n = 0)) %>%
  mutate(Percent = paste0(round((n/nrow(numParameters))*100, 2), "%"))

View(mostPowerful)

write.csv(mostPowerful, file = "./Results/Chi2/MostPowerful.csv")

# Frequency of how many times a method is least powerful
leastPowerful <- powerTable %>%
  pivot_longer(cols = c("method1_bonf", "method1_sidak", "method1_dap",
                        "method2", "method3", "method4_Chi2",
                        "method5_T", "method5_MVN"), names_to = "Method",
               values_to = "Power") %>%
  group_by(Scenario) %>%
  filter(Power == min(Power)) %>%
  mutate(unique_id = row_number()) %>%
  pivot_wider(names_from = unique_id, values_from = c("Method", "Power")) %>%
  group_by(Method_1, Method_2) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(Method_1 = factor(Method_1, levels = methodList)) %>%
  complete(Method_1 = levels(Method_1), fill = list(n = 0)) %>%
  mutate(Percent = paste0(round((n/nrow(numParameters))*100, 2), "%"))

View(leastPowerful)
write.csv(leastPowerful, file = "./Results/Chi2/LeastPowerful.csv")

# Histogram of power results for methods
powerLong <- powerTable %>%
  pivot_longer(cols = c("method1_bonf", "method1_sidak", "method1_dap",
                        "method2", "method3", "method4_Chi2",
                        "method5_T", "method5_MVN"), names_to = "Method",
               values_to = "Power") %>%
  mutate("Method Label" = fct_recode(Method,
                                     "1. P-Value Adjustment (Bonferroni)" = "method1_bonf",
                                     "1. P-Value Adjustment (Sidak)" = "method1_sidak",
                                     "1. P-Value Adjustment (D/AP)" = "method1_dap",
                                     "2. Combined Outcomes" = "method2",
                                     "3. Single 1-DF Test" = "method3",
                                     "4. Disjunctive 2-DF" = "method4_Chi2",
                                     "5. Conjunctive IU Test (t-Dist)" = "method5_T",
                                     "5. Conjunctive IU Test (MVN-Dist)" = "method5_MVN"))

summaryStats <- powerLong %>%
  dplyr::select(`Method Label`, Power) %>%
  group_by(`Method Label`) %>%
  summarize(SD = round(sd(Power), 2),
            Mean = round(mean(Power), 2),
            Median = round(median(Power), 2),
            Min = round(min(Power), 2),
            Max = round(max(Power), 2))

ggplot(data = powerLong, aes(Power)) +
  geom_histogram(bins = 30, fill = 'blue') +
  facet_wrap(~`Method Label`) +
  ylab("Count") +
  xlab("Statistical Power") +
  geom_text(data = summaryStats,
            aes(x = -Inf, y = Inf,
                label = paste0("Mean: ", Mean,
                               "   Min: ", Min,
                               "   Max: ", Max
                )
            ),
            hjust = -0.26, vjust = 1.5,
            size = 4) +
  theme(text = element_text(size = 20))

# Average ranking among the methods
rankData <- powerTable %>%
  dplyr::select(Scenario, starts_with("method")) %>%
  pivot_longer(cols = starts_with("Method"),
               names_to = "Method",
               values_to = "Power") %>%
  arrange(Scenario, Power) %>%
  group_by(Scenario) %>%
  mutate(Rank = rank(-Power, ties.method = "min")) %>%
  ungroup()

rankDataSummary <- rankData %>%
  group_by(Method) %>%
  summarize(Mean = round(mean(Rank), 2),
            SD = round(sd(Rank), 2),
            Median = round(median(Rank), 2),
            Min = round(min(Rank), 2),
            Max = round(max(Rank), 2)) %>%
  arrange(Mean)

# Create a summary table
rank_summary_table <- dcast(rankData, Scenario ~ Method, value.var = "Rank")

# Melt the table for ggplot
rank_summary_melted <- melt(rank_summary_table, id.vars = "Scenario") %>%
  mutate(Method = variable) %>%
  dplyr::select(-variable) %>%
  mutate(Method = recode(Method,
                         "method3" = "3. Single 1-DF Weighted",
                         "method2" = "2. Combined Outcomes",
                         "method4_Chi2" = "4. Disj. 2-DF (Chi2)",
                         "method5_MVN" = "5. Conj. IU (MVN)",
                         "method5_T" = "5. Conj. IU (T)",
                         "method1_dap" = "1. P-Val Adj. (D/AP)",
                         "method1_sidak" = "1. P-Val Adj. (Sidak)",
                         "method1_bonf" = "1. P-Val Adj. (Bonf.)"))

mean_ranks <- rankData %>%
  group_by(Method) %>%
  summarize(`Mean Ranking` = round(mean(Rank, na.rm = TRUE), 2)) %>%
  arrange(`Mean Ranking`) %>%
  mutate(Method = recode(Method,
                         "method3" = "3. Single 1-DF Weighted",
                         "method2" = "2. Combined Outcomes",
                         "method4_Chi2" = "4. Disj. 2-DF (Chi2)",
                         "method5_MVN" = "5. Conj. IU (MVN)",
                         "method5_T" = "5. Conj. IU (T)",
                         "method1_dap" = "1. P-Val Adj. (D/AP)",
                         "method1_sidak" = "1. P-Val Adj. (Sidak)",
                         "method1_bonf" = "1. P-Val Adj. (Bonf.)"))

table_grob <- tableGrob(mean_ranks)

# Plot with heatmap
heatmap_plot <- ggplot(rank_summary_melted, aes(x = Method, y = Scenario, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue", name = "Rank of Power\n(Lower is More Powerful)") +
  scale_y_continuous(breaks = seq(0, 3000, by = 500)) +
  xlab("Design Method") +
  ylab("Scenario Index") +
  labs(fill = "Rank of Power\n(Lower is More Powerful)") +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 25, hjust = 1),
        plot.margin = unit(c(1, 4, 1, 1), "lines")) +  # Adjust margins to make space for table
  annotation_custom(grob = table_grob,
                    xmin = 7, xmax = 13.1, ymin = -400, ymax = 600)  # Adjust these values to place the table
print(heatmap_plot)


# Boxplots for differences in betas, sigmas, rhos
diffData <- powerTable %>%
  mutate(`rho02 - rho01` = rho02 - rho01,
         `Var(Y2) - Var(Y1)` = varY2 - varY1,
         `Beta2 - Beta1` = beta2 - beta1) %>%
  pivot_longer(cols = starts_with("method"),
               names_to = "Method", values_to = "Power") %>%
  pivot_longer(cols = c("rho02 - rho01", "Var(Y2) - Var(Y1)", "Beta2 - Beta1"),
               names_to = "Variable", values_to = "Difference") %>%
  mutate(Variable = paste0(Variable, " = ", Difference)) %>%
  dplyr::select(-Difference) %>%
  arrange(Variable)

ggplot(diffData, aes(x = Power, y = Method, group = Method, fill = Method)) +
  geom_boxplot() + facet_wrap(~Variable)

# Boxplots for K and m
sizeData <- powerTable %>%
  arrange(K, m) %>%
  mutate(`Sample Size` = paste0("K = ", K, ", m = ", m),
         `Sample Size` = factor(`Sample Size`, levels = c("K = 6, m = 50",
                                                          "K = 6, m = 70",
                                                          "K = 8, m = 50",
                                                          "K = 8, m = 70",
                                                          "K = 10, m = 50",
                                                          "K = 10, m = 70"))) %>%
  pivot_longer(cols = starts_with("method"),
               names_to = "Method", values_to = "Power")

ggplot(sizeData, aes(x = Power, y = Method, group = Method, fill = Method)) +
  geom_boxplot() + facet_wrap(~`Sample Size`, ncol = 2)

# Most Powerful Methods --------------------------------------------------------

methods234 <- powerTable %>%
  dplyr::select(-method1_bonf, -method1_sidak, -method1_dap,
                -method5_T, -method5_MVN) %>%
  mutate(rho02minus01 = rho02 - rho01,
         var2minus1 = varY2 - varY1,
         beta2minus1 = beta2 - beta1,
         mostPower = pmax(method2, method3, method4_Chi2)) %>%
  mutate(best = pmap_chr(list(method2, method3, method4_Chi2), ~ {
    values <- c(method2 = ..1, method3 = ..2, method4_Chi2 = ..3)
    tied_methods <- names(values)[values == max(values)]
    paste(tied_methods, collapse = " = ")
  })) %>%
  ungroup() %>%
  arrange(beta2minus1, var2minus1, rho02minus01) %>%
  mutate(newID = paste(rho02minus01, var2minus1, beta2minus1)) %>%
  group_by(newID) %>%
  mutate(group_id = cur_group_id()) %>%
  ungroup()

# Getting counts for groups
scenarioGroups <- methods234 %>%
  dplyr::select(Scenario, group_id) %>%
  group_by(group_id) %>%
  summarize(n = n())

numCases <- if(all(scenarioGroups$n == 96)) {96} else {NA}

mosaic::tally(best ~ group_id, data = methods234)

bestSummary <- methods234 %>%
  dplyr::select(beta1, beta2, varY1, varY2, rho01, rho02, best,
                best, newID, group_id) %>%
  group_by(best, group_id,
           beta1, beta2, varY1, varY2, rho01, rho02
  ) %>%
  mutate(n = n()) %>%
  distinct() %>%
  spread(best, n) %>%
  mutate(across(contains("method"), ~ replace_na(.x, 0))) %>%
  arrange(group_id)

View(bestSummary)

allBest <- bestSummary %>%
  mutate(beta2minus1 = beta2 - beta1,
         var2minus1 = varY2 - varY1,
         rho02minus01 = rho02 - rho01) %>%
  mutate(beta1beta2 = paste0("(", beta1, ", ", beta2, ")"),
         varY1varY2 = paste0("(", varY1, ", ", varY2, ")"),
         rho01rho02 = paste0("(", rho01, ", ", rho02, ")")) %>%
  ungroup() %>%
  mutate(BetaCase = ifelse(beta1 < beta2, paste0("beta1 < beta2"),
                           ifelse(beta1 > beta2, paste0("beta1 > beta2"),
                                  ifelse(beta1 == beta2, paste0("beta1 = beta2"), NA))),
         VarCase = ifelse(varY1 < varY2, paste0("varY1 < varY2"),
                          ifelse(varY1 > varY2, paste0("varY1 > varY2"),
                                 ifelse(varY1 == varY2, paste0("varY1 = varY2"), NA))),
         RhoCase = ifelse(rho01 < rho02, paste0("rho01 < rho02"),
                          ifelse(rho01 > rho02, paste0("rho01 > rho02"),
                                 ifelse(rho01 == rho02, paste0("rho01 = rho02"), NA)))) %>%
  dplyr::select(group_id, beta1beta2, beta2minus1,
                varY1varY2, var2minus1,
                rho01rho02, rho02minus01,
                contains("="), contains("method"),
                BetaCase, VarCase, RhoCase)

allBestRaw <- allBest %>%
  dplyr::select(-BetaCase, -VarCase, -RhoCase) %>%
  arrange(beta2minus1, var2minus1, rho02minus01) %>%
  mutate(across(contains("method"), ~ paste0(round(.x/numCases*100, 2),
                                             "%", " (n = ", .x, ")"))) %>%
  mutate(n = numCases)

allBestCases <- allBest %>%
  dplyr::select(BetaCase, VarCase, RhoCase, starts_with("method")) %>%
  mutate(n = numCases) %>%
  group_by(BetaCase, VarCase, RhoCase) %>%
  summarise_all(sum) %>%
  mutate(across(contains("method"), ~ paste0(round(.x/n*100, 2),
                                             "%", " (n = ", .x, ")")))

View(allBest)

View(allBestRaw)
View(allBestCases)

write.csv(allBestRaw, file = "./Results/Chi2/BestAll_Raw.csv")
write.csv(allBestCases, file = "./Results/Chi2/BestAll_Cases.csv")
