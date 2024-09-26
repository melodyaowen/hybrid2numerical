# Comparison of power between Methods 1 and 4
source("./RequiredPackages.R")

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
  mutate(Method4Bigger = ifelse(PowerDiff > 0, "True", ifelse(PowerDiff < 0, "False", "Equal")))

plot14data1 <- method14dat %>%
  filter(Methods != "D/AP Correction")

plot14data2 <- method14dat %>%
  filter(Methods == "D/AP Correction") %>%
  mutate(rho2 = factor(rho2, levels = c("0.1", "0.3", "0.5", "0.7")))

method14_plot1 <- ggplot(plot14data1, aes(x = lambda4minus1, y = PowerDiff,
                                          group = Method4Bigger, color = Method4Bigger)) +
  geom_point() + theme(text = element_text(size = 15)) +
  facet_wrap(~Methods) +
  xlab(TeX("$\\lambda^{(Method 4)} - \\lambda^{(Method 1)}$")) +   # LaTeX for x-axis using TeX()
  ylab(TeX("$\\pi^{(Method 4)} - \\pi^{(Method 1)}$")) +
  labs(color = TeX("$\\pi^{(Method 1)} < \\pi^{(Method 4)}$"))

method14_plot1
ggsave(filename = "~/Desktop/2. Hybrid Software and Simulation/Word Drafts/Figures/methods14plot1.png",
       plot = method14_plot1,
       width = 10, height = 4, dpi = 300)

latex_labels <- c(
  `0.1` = "rho[2] == 0.1",
  `0.3` = "rho[2] == 0.3",
  `0.5` = "rho[2] == 0.5",
  `0.7` = "rho[2] == 0.7"
)

method14_plot2 <- ggplot(plot14data2, aes(x = lambda4minus1, y = PowerDiff,
                                          group = Method4Bigger, color = Method4Bigger)) +
  geom_point() + theme(text = element_text(size = 15)) +
  facet_wrap(~rho2, ncol = 2,
             labeller = as_labeller(latex_labels, label_parsed)) +
  theme(strip.text = element_text(size = 15)) +
  xlab(TeX("$\\lambda^{(Method 4)} - \\lambda^{(Method 1)}$")) +   # LaTeX for x-axis using TeX()
  ylab(TeX("$\\pi^{(Method 4)} - \\pi^{(Method 1)}$")) +
  labs(color = TeX("$\\pi^{(Method 1)} < \\pi^{(Method 4)}$"))

method14_plot2
ggsave(filename = "~/Desktop/2. Hybrid Software and Simulation/Word Drafts/Figures/methods14plot2.png",
       plot = method14_plot2,
       width = 10, height = 6, dpi = 300)



  # mutate(`P-Value (Bonf.)` = 1 - pchisq(qchisq(1 - alpha_B, df = 1, ncp = 0),
  #                                       df = 1, ncp = lambda1, lower.tail = TRUE),
  #        `P-Value (Sidak)` = 1 - pchisq(qchisq(1 - alpha_S, df = 1, ncp = 0),
  #                                       df = 1, ncp = lambda1, lower.tail = TRUE),
  #        `P-Value (D/AP)` = 1 - pchisq(qchisq(1 - alpha_D, df = 1, ncp = 0),
  #                                     df = 1, ncp = lambda1, lower.tail = TRUE)) %>%
  # mutate(`Disjunctive 2-DF` = 1 - pchisq(qchisq(1 - alpha, df = 2, ncp = 0,
  #                                               lower.tail = TRUE, log.p = FALSE),
  #                                        df = 2, ncp = lambda4, lower.tail = TRUE)) %>%
  # pivot_longer(cols = c(`P-Value (Bonf.)`, `P-Value (Sidak)`,
  #                       `P-Value (D/AP)`, `Disjunctive 2-DF`),
  #              names_to = "Method", values_to = "Power") %>%
  # mutate(Method = factor(Method,
  #                        levels = c("P-Value (Bonf.)", "P-Value (Sidak)",
  #                                   "P-Value (D/AP)", "Disjunctive 2-DF"))) %>%
  # mutate(rho2 = factor(rho2, levels = c("0.1", "0.3", "0.5", "0.7")))

# method14_plot <- ggplot(method14dat, aes(x = lambda4minus1, y = Power,
#                                          group = Method, color = Method)) +
#   geom_point() + theme(text = element_text(size = 15)) +
#   facet_wrap(~rho2)
#
# method14_plot

