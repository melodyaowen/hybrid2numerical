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

# Run power calculations on all true parameters --------------------------------
powerTable <- numParameters %>%
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
View(numParameters)
nrow(numParameters)

# Check for cases where power is 100
scenarios100 <- powerTable %>%
  dplyr::filter(if_any(c(method1_bonf, method1_sidak, method1_dap,
                         method2, method3,
                         method4_chi2, method4_F,
                         method5_T, method5_MVN), ~ . == 100))

View(scenarios100)

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
  complete(Method = levels(Method), fill = list(n = 0)) %>%
  mutate(Percent = round((n/nrow(numParameters))*100, 2))
View(mostPowerful)

# Frequency of how many times a method is least powerful
leastPowerful <- powerTable %>%
  pivot_longer(cols = c("method1_bonf", "method1_sidak", "method1_dap",
                        "method2", "method3", "method4_chi2", "method4_F",
                        "method5_T", "method5_MVN"), names_to = "Method",
               values_to = "Power") %>%
  group_by(Scenario) %>%
  slice(which.min(Power)) %>%
  ungroup() %>%
  group_by(Method) %>%
  summarize(n = n()) %>%
  mutate(Method = factor(Method, levels = methodList)) %>%
  complete(Method = levels(Method), fill = list(n = 0)) %>%
  mutate(Percent = round((n/nrow(numParameters))*100, 2))
View(leastPowerful)

# Histogram of power results for methods
powerLong <- powerTable %>%
  pivot_longer(cols = c("method1_bonf", "method1_sidak", "method1_dap",
                        "method2", "method3", "method4_chi2", "method4_F",
                        "method5_T", "method5_MVN"), names_to = "Method",
               values_to = "Power") %>%
  mutate("Method Label" = fct_recode(Method,
                                     "1. P-Value Adjustment (Bonferroni)" = "method1_bonf",
                                     "1. P-Value Adjustment (Sidak)" = "method1_sidak",
                                     "1. P-Value Adjustment (D/AP)" = "method1_dap",
                                     "2. Combined Outcomes" = "method2",
                                     "3. Single 1-DF Test" = "method3",
                                     "4. Disjunctive 2-DF (Chi2-Dist)" = "method4_chi2",
                                     "4. Disjunctive 2-DF (F-Dist)" = "method4_F",
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

# Boxplots of the methods

ggplot(data = powerLong, aes(x = Power, y = fct_rev(`Method Label`))) +
  geom_boxplot() +
  ylab("Design Method") +
  xlab("Statistical Power")



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
                         "method4_chi2" = "4. Disj. 2-DF (Chi2)",
                         "method2" = "2. Combined Outcomes",
                         "method4_F" = "4. Disj. 2-DF (F)",
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
                         "method4_chi2" = "4. Disj. 2-DF (Chi2)",
                         "method2" = "2. Combined Outcomes",
                         "method4_F" = "4. Disj. 2-DF (F)",
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
                    xmin = 9.1, xmax = 13.1, ymin = -400, ymax = 600)  # Adjust these values to place the table
print(heatmap_plot)

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


# Most Powerful Methods --------------------------------------------------------

methods234 <- powerTable %>%
  dplyr::select(-method1_bonf, -method1_sidak, -method1_dap,
                -method4_F, -method5_T, -method5_MVN) %>%
  mutate(rho02minus01 = rho02 - rho01,
         var2minus1 = varY2 - varY1,
         beta2minus1 = beta2 - beta1,
         mostPower = max(method2, method3, method4_chi2)) %>%
  mutate(best = pmap_chr(list(method2, method3, method4_chi2), ~ {
    values <- c(method2 = ..1, method3 = ..2, method4_chi2 = ..3)
    names(values)[which.max(values)]})) %>%
  ungroup %>%
  arrange(beta2minus1, var2minus1, rho02minus01) %>%
  mutate(newID = paste(rho02minus01, var2minus1, beta2minus1)) %>%
  group_by(newID) %>%
  mutate(group_id = cur_group_id())

mosaic::tally(best ~ group_id, data = methods234)

bestSummary <- methods234 %>%
  #dplyr::filter(best == "method2") %>%
  dplyr::select(beta1, beta2, varY1, varY2, rho01, rho02, best,
                #method2, method3, method4_chi2,
                best, newID, group_id) %>%
  group_by(best, group_id,
           beta1, beta2, varY1, varY2, rho01, rho02
           ) %>%
  mutate(n = n()) %>%
  distinct() %>%
  spread(best, n) %>%
  arrange(group_id) %>%
  mutate(best2 = ifelse(!is.na(method2), 1, 0),
         best3 = ifelse(!is.na(method3), 1, 0),
         best4 = ifelse(!is.na(method4_chi2), 1, 0),
         count = best2 + best3 + best4) %>%
  dplyr::select(-best2, -best3, -best4)

View(bestSummary)


allBest <- bestSummary %>%
  mutate(beta2minus1 = beta2 - beta1,
         var2minus1 = varY2 - varY1,
         rho02minus01 = rho02 - rho01) %>%
  mutate(beta1beta2 = paste0("(", beta1, ", ", beta2, ")"),
         varY1varY2 = paste0("(", varY1, ", ", varY2, ")"),
         rho01rho02 = paste0("(", rho01, ", ", rho02, ")")) %>%
  ungroup() %>%
  dplyr::select(group_id, beta1beta2, beta2minus1,
                varY1varY2, var2minus1,
                rho01rho02, rho02minus01,
                method2, method3, method4_chi2, count) %>%
  arrange(beta2minus1, var2minus1, rho02minus01) %>%
  mutate(method2 = ifelse(is.na(method2), 0, method2),
         method3 = ifelse(is.na(method3), 0, method3),
         method4_chi2 = ifelse(is.na(method4_chi2), 0, method4_chi2))

View(allBest)

write.csv(allBest, file = "./Results/BestAll.csv")

# Cases where only one best method given the scenario
bestOverall <- bestSummary %>%
  dplyr::filter(count == 1) %>%
  mutate(beta2minus1 = beta2 - beta1,
         var2minus1 = varY2 - varY1,
         rho02minus01 = rho02 - rho01) %>%
  mutate(beta1beta2 = paste0("(", beta1, ", ", beta2, ")"),
         varY1varY2 = paste0("(", varY1, ", ", varY2, ")"),
         rho01rho02 = paste0("(", rho01, ", ", rho02, ")")) %>%
  ungroup() %>%
  mutate(BestMethod = ifelse(!is.na(method2), "Method 2",
                             ifelse(!is.na(method3), "Method 3",
                                    ifelse(!is.na(method4_chi2), "Method 4 (Chi2)", "ERROR")))) %>%
  dplyr::select(group_id, beta1beta2, beta2minus1,
                varY1varY2, var2minus1,
                rho01rho02, rho02minus01, BestMethod,
                method2, method3, method4_chi2, count) %>%
  arrange(beta2minus1, var2minus1, rho02minus01)

write.csv(bestOverall, file = "./Results/BestOverall.csv")

# Cases where there are multiple best methods for a given scenario
bestSplit <- bestSummary %>%
  dplyr::filter(count != 1) %>%
  mutate(beta2minus1 = beta2 - beta1,
         var2minus1 = varY2 - varY1,
         rho02minus01 = rho02 - rho01) %>%
  mutate(beta1beta2 = paste0("(", beta1, ", ", beta2, ")"),
         varY1varY2 = paste0("(", varY1, ", ", varY2, ")"),
         rho01rho02 = paste0("(", rho01, ", ", rho02, ")")) %>%
  ungroup() %>%
dplyr::select(group_id, beta1beta2, beta2minus1,
              varY1varY2, var2minus1,
              rho01rho02, rho02minus01,
              method2, method3, method4_chi2, count) %>%
  arrange(beta2minus1, var2minus1, rho02minus01)

write.csv(bestSplit, file = "./Results/BestSplit.csv")


# ggplot(rho1rho2MostPowerful, aes(x = rho1, y = n, fill = Method)) +
#   geom_bar(stat = "identity", position = "dodge") + facet_wrap(~rho2) +
#   ylab("Frequency") + ggtitle("Number of input scenarios a \nmethod is found to have highest power, \nseparated by rho1 and rho2") + geom_text(aes(label = n), vjust = -1, position = position_dodge(width = 0.9))




# Least Powerful Methods -------------------------------------------------------

methods145 <- powerTable %>%
  dplyr::select(-method1_dap, -method1_sidak,
                -method2, -method3, -method4_chi2,
                -method5_MVN) %>%
  mutate(rho02minus01 = rho02 - rho01,
         var2minus1 = varY2 - varY1,
         beta2minus1 = beta2 - beta1,
         leastPower = min(method1_bonf, method4_F, method5_T)) %>%
  mutate(worst = pmap_chr(list(method1_bonf, method4_F, method5_T), ~ {
    values <- c(method1_bonf = ..1, method4_F = ..2, method5_T = ..3)
    names(values)[which.min(values)]})) %>%
  ungroup() %>%
  arrange(beta2minus1, var2minus1, rho02minus01) %>%
  mutate(newID = paste(rho02minus01, var2minus1, beta2minus1)) %>%
  group_by(newID) %>%
  mutate(group_id = cur_group_id())

mosaic::tally(worst ~ group_id, data = methods145)

worstSummary <- methods145 %>%
  #dplyr::filter(best == "method2") %>%
  dplyr::select(beta1, beta2, varY1, varY2, rho01, rho02,
                #method2, method3, method4_chi2,
                worst, newID, group_id) %>%
  group_by(worst, group_id,
           beta1, beta2, varY1, varY2, rho01, rho02
  ) %>%
  mutate(n = n()) %>%
  distinct() %>%
  spread(worst, n) %>%
  arrange(group_id) %>%
  mutate(worst1 = ifelse(!is.na(method1_bonf), 1, 0),
         worst4 = ifelse(!is.na(method4_F), 1, 0),
         worst5 = ifelse(!is.na(method5_T), 1, 0),
         count = worst1 + worst4 + worst5) %>%
  dplyr::select(-worst1, -worst4, -worst5)

View(worstSummary)


allWorst <- worstSummary %>%
  mutate(beta2minus1 = beta2 - beta1,
         var2minus1 = varY2 - varY1,
         rho02minus01 = rho02 - rho01) %>%
  mutate(beta1beta2 = paste0("(", beta1, ", ", beta2, ")"),
         varY1varY2 = paste0("(", varY1, ", ", varY2, ")"),
         rho01rho02 = paste0("(", rho01, ", ", rho02, ")")) %>%
  ungroup() %>%
  dplyr::select(group_id, beta1beta2, beta2minus1,
                varY1varY2, var2minus1,
                rho01rho02, rho02minus01,
                method1_bonf, method4_F, method5_T, count) %>%
  arrange(beta2minus1, var2minus1, rho02minus01) %>%
  mutate(method1_bonf = ifelse(is.na(method1_bonf), 0, method1_bonf),
         method4_F = ifelse(is.na(method4_F), 0, method4_F),
         method5_T = ifelse(is.na(method5_T), 0, method5_T))

View(allWorst)

write.csv(allWorst, file = "./Results/WorstAll.csv")


