source("./RequiredPackages.R")

# Figure 2. Comparison of methods 1 and 4 power (bonf and sidak)

method14dat <- expand.grid(lambda1 = seq(1, 30, by = 0.5),
                           lambda4 = seq(1, 30, by = 0.5),
                           rho2 = c(0.1, 0.3, 0.5, 0.7),
                           alpha = 0.05
                           ) %>%
  filter(lambda1 < lambda4) %>%
  mutate(lambda4minus1 = lambda4 - lambda1) %>%
  mutate(alpha_B = alpha/2, # Bonferroni
         alpha_S = 1 - (1 - alpha)^(1/2), # Sidak
         alpha_D = 1 - (1 - alpha)^(1/(2^(1 - rho2))) # D/AP
         ) %>%
  mutate(`Bonferonni` = 1 - pchisq(qchisq(1 - alpha_B, df = 1, ncp = 0),
                                        df = 1, ncp = lambda1, lower.tail = TRUE),
         `Sidak` = 1 - pchisq(qchisq(1 - alpha_S, df = 1, ncp = 0),
                                        df = 1, ncp = lambda1, lower.tail = TRUE),
         `DAP` = 1 - pchisq(qchisq(1 - alpha_D, df = 1, ncp = 0),
                                       df = 1, ncp = lambda1, lower.tail = TRUE)) %>%
  mutate(`Disjunctive` = 1 - pchisq(qchisq(1 - alpha, df = 2, ncp = 0,
                                                lower.tail = TRUE, log.p = FALSE),
                                         df = 2, ncp = lambda4, lower.tail = TRUE)) %>%
  mutate(`Bonferroni Correction` = Disjunctive - Bonferonni,
         `Sidak Correction` = Disjunctive - Sidak,
         `D/AP Correction` = Disjunctive - DAP) %>%
  pivot_longer(cols = c(`Bonferroni Correction`, `Sidak Correction`, `D/AP Correction`),
               names_to = "Methods", values_to = "PowerDiff") %>%
  mutate(PowerDiff = round(PowerDiff, 4)) %>%
  mutate(`HighestPower` = ifelse(PowerDiff > 0, "Disjunctive 2-DF",
                                     ifelse(PowerDiff < 0, "P-Value Adjustment",
                                            "Equal"))) %>%
  mutate(HighestPower = factor(HighestPower, levels = c("Disjunctive 2-DF", "P-Value Adjustment", "Equal")))

plot14data1 <- method14dat %>%
  filter(Methods != "D/AP Correction")

plot14data2 <- method14dat %>%
  filter(Methods == "D/AP Correction") %>%
  mutate(rho2 = factor(rho2, levels = c("0.1", "0.3", "0.5", "0.7")))

figure2 <- ggplot(plot14data1, aes(x = lambda4minus1, y = PowerDiff,
                                          group = HighestPower,
                                   color = HighestPower)) +
  geom_point() + theme(text = element_text(size = 15)) +
  facet_wrap(~Methods) +
  xlab(TeX("$\\lambda^{DIS2DF} - \\lambda^{PADJ}$")) +   # LaTeX for x-axis using TeX()
  ylab(TeX("$\\pi^{DIS2DF} - \\pi^{PADJ}$")) +
  labs(color = TeX("Method with highest power"))

# Saving figure 3 as png
ggsave(filename = "./MathematicalComparisons/Output/Figure2.png",
       plot = figure2,
       width  = 3200, height = 1500, units  = "px")




# Figure 3. Comparison of methods 1 and 4 power (d/ap)

latex_labels <- c(
  `0.1` = "rho[2] == 0.1",
  `0.3` = "rho[2] == 0.3",
  `0.5` = "rho[2] == 0.5",
  `0.7` = "rho[2] == 0.7"
)

figure3 <- ggplot(plot14data2, aes(x = lambda4minus1, y = PowerDiff,
                                          group = HighestPower,
                                   color = HighestPower)) +
  geom_point() + theme(text = element_text(size = 15)) +
  facet_wrap(~rho2, ncol = 2,
             labeller = as_labeller(latex_labels, label_parsed)) +
  theme(strip.text = element_text(size = 15)) +
  xlab(TeX("$\\lambda^{DIS2DF} - \\lambda^{PADJ}$")) +   # LaTeX for x-axis using TeX()
  ylab(TeX("$\\pi^{DIS2DF} - \\pi^{PADJ}$")) +
  labs(color = TeX("Method with highest power"))

# Saving figure 4 as png
ggsave(filename = "./MathematicalComparisons/Output/Figure3.png",
       plot = figure3,
       width  = 3200, height = 1800, units  = "px")
