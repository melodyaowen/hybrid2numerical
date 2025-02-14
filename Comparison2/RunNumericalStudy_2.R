source("./RequiredPackages.R")
source("./Comparison1/Method5_2sided.R")

# "Comparison 1"
# "2-sided" comparison using the Chi^2 distribution and MVN distribution
# Method 5 is two 2-sided tests, we use a new function defined for this purpose
# and not the package function

# Table of all Parameters ------------------------------------------------------
numParameters <- expand.grid(K = c(4, 6, 8, 10),
                             m = c(50, 70, 100),
                             betas = c(paste("0.1 0.4"),
                                       paste("0.2 0.4"),
                                       paste("0.3 0.4"),
                                       paste("0.4 0.4")),
                             vars = c(paste("0.5 1.5"),
                                      paste("0.5 1"),
                                      paste("1 1"),
                                      paste("1 0.5"),
                                      paste("1.5 0.5")),
                             rho0 = c(paste("0.05 0.1"),
                                      paste("0.07 0.1"),
                                      paste("0.1 0.1"),
                                      paste("0.1 0.07"),
                                      paste("0.1 0.05")),
                             rho1 = c(0.005, 0.01, 0.02, 0.05, 0.07),
                             rho2 = c(0.1, 0.3, 0.5, 0.7, 0.9),
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
         'method5_MVN' = calc_pwr_conj_test_2sided(K = K, m = m, alpha = alpha,
                                                   beta1 = beta1, beta2 = beta2,
                                                   varY1 = varY1, varY2 = varY2,
                                                   rho01 = rho01, rho02 = rho02,
                                                   rho1 = rho1, rho2  = rho2,
                                                   r = r, dist = "MVN")) %>%
  mutate_at(vars(contains('method')), funs(.*100))

View(head(powerTable))
View(head(numParameters))
nrow(numParameters)

# Check for cases where power is 100
scenarios100 <- powerTable %>%
  dplyr::filter(if_any(c(method1_bonf, method1_sidak, method1_dap,
                         method2, method3,
                         method4_Chi2, method5_MVN), ~ . == 100))

nrow(scenarios100)
View(scenarios100)

# Most and Least Powerful Tables -----------------------------------------------
# Frequency of how many times a method is the most powerful
methodList <- c("method1_bonf", "method1_sidak", "method1_dap",
                "method2", "method3", "method4_Chi2", "method5_MVN")

mostPowerful <- powerTable %>%
  pivot_longer(cols = c("method1_bonf", "method1_sidak", "method1_dap",
                        "method2", "method3", "method4_Chi2", "method5_MVN"),
               names_to = "Method",
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

write.csv(mostPowerful, file = "./Comparison1/Results1/MostPowerful_1.csv")

# Frequency of how many times a method is least powerful
leastPowerful <- powerTable %>%
  pivot_longer(cols = c("method1_bonf", "method1_sidak", "method1_dap",
                        "method2", "method3", "method4_Chi2", "method5_MVN"),
               names_to = "Method",
               values_to = "Power") %>%
  group_by(Scenario) %>%
  filter(Power == min(Power)) %>%
  mutate(unique_id = row_number()) %>%
  pivot_wider(names_from = unique_id, values_from = c("Method", "Power")) %>%
  group_by(Method_1) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(Method_1 = factor(Method_1, levels = methodList)) %>%
  complete(Method_1 = levels(Method_1), fill = list(n = 0)) %>%
  mutate(Percent = paste0(round((n/nrow(numParameters))*100, 2), "%"))

View(leastPowerful)
write.csv(leastPowerful, file = "./Comparison1/Results1/LeastPowerful_1.csv")

# Power Histogram --------------------------------------------------------------

# Histogram of power results for methods
powerLong <- powerTable %>%
  pivot_longer(cols = c("method1_bonf", "method1_sidak", "method1_dap",
                        "method2", "method3", "method4_Chi2", "method5_MVN"),
               names_to = "Method",
               values_to = "Power") %>%
  mutate("Method Label" = fct_recode(Method,
                                     "1. P-Value Adjustment (Bonferroni)" = "method1_bonf",
                                     "1. P-Value Adjustment (Sidak)" = "method1_sidak",
                                     "1. P-Value Adjustment (D/AP)" = "method1_dap",
                                     "2. Combined Outcomes" = "method2",
                                     "3. Single Weighted 1-DF Test" = "method3",
                                     "4. Disjunctive 2-DF" = "method4_Chi2",
                                     "5. Conjunctive IU Test (MVN-Dist)" = "method5_MVN"))

summaryStats <- powerLong %>%
  dplyr::select(`Method Label`, Power) %>%
  group_by(`Method Label`) %>%
  summarize(SD = round(sd(Power), 2),
            Mean = round(mean(Power), 2),
            Median = round(median(Power), 2),
            Min = round(min(Power), 2),
            Max = round(max(Power), 2))

powerHistogram <- ggplot(data = powerLong, aes(Power)) +
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
            hjust = -0.5, vjust = 1.5,
            size = 4) +
  theme(text = element_text(size = 20))

ggsave(filename = "./Comparison1/Results1/PowerHistogram_1.png",
       plot = powerHistogram,
       width  = 5000, height = 3000, units  = "px")

# Ranking Heatmap --------------------------------------------------------------
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
                         "method3" = "Single Weighted 1-DF",
                         "method2" = "Combined Outcomes",
                         "method4_Chi2" = "Disjunctive 2-DF",
                         "method5_MVN" = "Conjunctive IU (MVN)",
                         "method1_dap" = "P-Val Adj. (D/AP)",
                         "method1_sidak" = "P-Val Adj. (Sidak)",
                         "method1_bonf" = "P-Val Adj. (Bonf.)"))

mean_ranks <- rankData %>%
  group_by(Method) %>%
  summarize(`Mean Ranking` = round(mean(Rank, na.rm = TRUE), 2)) %>%
  arrange(`Mean Ranking`) %>%
  mutate(Method = recode(Method,
                         "method3" = "Single Weighted 1-DF",
                         "method2" = "Combined Outcomes",
                         "method4_Chi2" = "Disjunctive 2-DF",
                         "method5_MVN" = "Conjunctive IU (MVN)",
                         "method1_dap" = "P-Val Adj. (D/AP)",
                         "method1_sidak" = "P-Val Adj. (Sidak)",
                         "method1_bonf" = "P-Val Adj. (Bonf.)"))

table_grob <- tableGrob(mean_ranks, rows = NULL)

# Plot with heatmap
rankHeatmap <- ggplot(rank_summary_melted, aes(x = Method, y = Scenario, fill = value)) +
  geom_tile() +
  scale_y_continuous(breaks = seq(0, 30000, by = 5000)) +
  xlab("Design Method") +
  ylab("Scenario Index") +
  scale_fill_gradient(
    low = "white", high = "blue",
    name = "Rank of Power\n(Smaller number = higher power)",
    guide = guide_colorbar(reverse = TRUE)
  ) +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 25, hjust = 1),
        plot.margin = unit(c(1, 4, 1, 1), "lines")) +
  annotation_custom(grob = table_grob,
                    xmin = 7.5, xmax = 10, ymin = -400, ymax = 600)

ggsave(filename = "./Comparison1/Results1/RankHeatmap_1.png",
       plot = rankHeatmap,
       width  = 5000, height = 2500, units  = "px")

# Most Powerful Method Tables --------------------------------------------------

# Start by grabbing the names of all methods that have at least
# one scenario of being the most powerful among all the methods
mostPowerfulMethodNames <- powerTable %>%
  pivot_longer(cols = c("method1_bonf", "method1_sidak", "method1_dap",
                        "method2", "method3", "method4_Chi2", "method5_MVN"),
               names_to = "Method",
               values_to = "Power") %>%
  group_by(Scenario) %>%
  filter(Power == max(Power)) %>%
  ungroup() %>%
  dplyr::select(Method) %>%
  distinct(Method)

# Vector of most powerful method names
mostPowerfulNames <- mostPowerfulMethodNames$Method

# Full power table with just most powerful methods
mostPowerfulFull <- powerTable %>%
  dplyr::select(Scenario, K, m, beta1, beta2, varY1, varY2, rho01, rho02,
                rho1, rho2, alpha, r, all_of(mostPowerfulNames)) %>%
  rowwise() %>%
  mutate(rho02minus01 = rho02 - rho01,
         var2minus1 = varY2 - varY1,
         beta2minus1 = beta2 - beta1,
         mostPower = max(c_across(all_of(mostPowerfulNames)))) %>%
  ungroup() %>%
  mutate(
    best = pmap_chr(
      dplyr::select(cur_data(), starts_with("method")),
      ~ {
        rowvals <- c(...)
        names(rowvals) <- names(dplyr::select(cur_data(), starts_with("method")))
        tied_methods <- names(rowvals)[rowvals == max(rowvals)]
        paste(tied_methods, collapse = ", ")
      }
    )
  ) %>%
  arrange(beta2minus1, var2minus1, rho02minus01) %>%
  mutate(newID = paste(rho02minus01, var2minus1, beta2minus1)) %>%
  group_by(newID) %>%
  mutate(group_id = cur_group_id()) %>%
  ungroup()

# Getting counts for groups
scenarioGroups <- mostPowerfulFull %>%
  dplyr::select(Scenario, group_id) %>%
  group_by(group_id) %>%
  summarize(n = n())

# Summary of methods and how many scenarios they're most powerful for
mosaic::tally(best ~ group_id, data = mostPowerfulFull)

# Summary data
bestSummary <- mostPowerfulFull %>%
  dplyr::select(beta1, beta2, varY1, varY2, rho01, rho02,
                best, newID, group_id) %>%
  group_by(best, group_id,
           beta1, beta2, varY1, varY2, rho01, rho02
  ) %>%
  mutate(n = n()) %>%
  distinct() %>%
  spread(best, n) %>%
  mutate(across(contains("method"), ~ replace_na(.x, 0))) %>%
  arrange(group_id)

#View(bestSummary)

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
  left_join(., scenarioGroups, by = "group_id") %>%
  rowwise() %>%
  mutate(across(contains("method"),
                ~ if_else(.x == 0, "0%",  # the special case: just "0%"
                          paste0(round(.x/n*100, 0), "% (n = ", .x, ")"))))

allBestCases <- allBest %>%
  dplyr::select(group_id, BetaCase, VarCase, RhoCase, starts_with("method")) %>%
  left_join(., scenarioGroups, by = "group_id") %>%
  dplyr::select(-group_id) %>%
  group_by(BetaCase, VarCase, RhoCase) %>%
  summarise_all(sum) %>%
  rowwise() %>%
  mutate(across(contains("method"),
                ~ if_else(.x == 0, "0%",  # the special case: just "0%"
                          paste0(round(.x/n*100, 0), "% (n = ", .x, ")"))))

#View(allBest)
#View(allBestRaw)
#View(allBestCases)

write.csv(allBestRaw, file = "./Comparison1/Results1/BestAll_1.csv")
write.csv(allBestCases, file = "./Comparison1/Results1/BestAllCases_1.csv")

# Result table based on standardized effect sizes ------------------------------

mostPowerfulFull_std <- powerTable %>%
  dplyr::select(Scenario, K, m, beta1, beta2, varY1, varY2, rho01, rho02,
                rho1, rho2, alpha, r, all_of(mostPowerfulNames)) %>%
  rowwise() %>%
  mutate(effect1std = round(beta1/sqrt(varY1), 2),
         effect2std = round(beta2/sqrt(varY2), 2),
         rho02minus01 = rho02 - rho01,
         mostPower = max(c_across(all_of(mostPowerfulNames)))) %>%
  ungroup() %>%
  mutate(
    best = pmap_chr(
      dplyr::select(cur_data(), starts_with("method")),
      ~ {
        rowvals <- c(...)
        names(rowvals) <- names(dplyr::select(cur_data(), starts_with("method")))
        tied_methods <- names(rowvals)[rowvals == max(rowvals)]
        paste(tied_methods, collapse = ", ")
      }
    )
  ) %>%
  arrange(effect1std, effect2std, rho02minus01) %>%
  mutate(newID = paste(effect1std, effect2std, rho02minus01)) %>%
  group_by(newID) %>%
  mutate(group_id = cur_group_id()) %>%
  ungroup()

# Getting counts for groups
scenarioGroups_std <- mostPowerfulFull_std %>%
  dplyr::select(Scenario, group_id) %>%
  group_by(group_id) %>%
  summarize(n = n())

# Summary of methods of how many scenarios they're most powerful for
mosaic::tally(best ~ group_id, data = mostPowerfulFull_std)

bestSummary_std <- mostPowerfulFull_std %>%
  dplyr::select(beta1, beta2, varY1, varY2, rho01, rho02,
                best, newID, group_id,
                effect1std, effect2std, rho02minus01) %>%
  group_by(best, group_id,
           beta1, beta2, varY1, varY2, rho01, rho02
  ) %>%
  mutate(n = n()) %>%
  distinct() %>%
  spread(best, n) %>%
  mutate(across(contains("method"), ~ replace_na(.x, 0))) %>%
  arrange(group_id)

#View(bestSummary_std)

allBest_std <- bestSummary_std %>%
  mutate(eff2minus1 = effect2std - effect1std,
         rho02minus01 = rho02 - rho01) %>%
  mutate(eff1eff2 = paste0("(", effect1std, ", ", effect2std, ")"),
         rho01rho02 = paste0("(", rho01, ", ", rho02, ")")) %>%
  ungroup() %>%
  mutate(EffectCase = ifelse(effect1std < effect2std, paste0("beta1/sigma1 < beta2/sigma2"),
                             ifelse(effect1std > effect2std, paste0("beta1/sigma1 > beta2/sigma2"),
                                    ifelse(effect1std == effect2std, paste0("beta1/sigma1 = beta2/sigma2"), NA))),
         RhoCase = ifelse(rho01 < rho02, paste0("rho01 < rho02"),
                          ifelse(rho01 > rho02, paste0("rho01 > rho02"),
                                 ifelse(rho01 == rho02, paste0("rho01 = rho02"), NA)))) %>%
  dplyr::select(group_id, eff1eff2, eff2minus1,
                rho01rho02, rho02minus01,
                contains("="), contains("method"),
                EffectCase, RhoCase)

allBestRaw_std <- allBest_std %>%
  dplyr::select(-EffectCase, -RhoCase) %>%
  arrange(eff2minus1, rho02minus01) %>%
  left_join(., scenarioGroups_std, by = "group_id") %>%
  mutate(across(contains("method"),
                ~ if_else(.x == 0, "0%",  # the special case: just "0%"
                          paste0(round(.x/n*100, 0), "% (n = ", .x, ")"))))

allBestCases_std <- allBest_std %>%
  dplyr::select(group_id, EffectCase, RhoCase, starts_with("method")) %>%
  left_join(., scenarioGroups_std, by = "group_id") %>%
  dplyr::select(-group_id) %>%
  group_by(EffectCase, RhoCase) %>%
  summarise_all(sum) %>%
  rowwise() %>%
  mutate(across(contains("method"),
                ~ if_else(.x == 0, "0%",  # the special case: just "0%"
                          paste0(round(.x/n*100, 0), "% (n = ", .x, ")"))))


write.csv(allBestRaw_std, file = "./Comparison1/Results1/BestAll_1_STD.csv")
write.csv(allBestCases_std, file = "./Comparison1/Results1/BestAllCases_1_STD.csv")

