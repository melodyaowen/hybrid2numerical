# rm(list = ls())
#
# # Testing out comparisons
#
# beta1 <- 0.2
# beta2 <- 0.1
# varY1 <- 0.25
# varY2 <- 0.4
# rho01 <- 0.05
# rho02 <- 0.02
# rho1 <- 0.01
# rho2 <- 0.05
# K <- 10
# m <- 50
# VIF1_new <- 1 + (m - 1)*rho01
# VIF2_new <- 1 + (m - 1)*rho02
# VIF12_new <- rho2 + (m - 1)*rho1
#
# # Method 3 NCP
# method3_NCP <- ((sqrt((beta1^2)/(((2*varY1)/(K*m))*VIF1_new)) +
#                    sqrt((beta2^2)/(((2*varY2)/(K*m))*VIF2_new)))/
#                   (sqrt(2*(1+(VIF12_new)/(sqrt(VIF1_new*VIF2_new))))))^2; method3_NCP
#
# # Method 4 NCP
# method4_NCP <- (K*m*((beta1^2)*varY2*VIF2_new -
#                        2*beta1*beta2*sqrt(varY1*varY2)*VIF12_new +
#                        (beta2^2)*varY1*VIF1_new))/
#   (2*varY1*varY2*(VIF1_new*VIF2_new - VIF12_new^2)); method4_NCP
#
# # Method 3 Power
# method3_crit <- qchisq(p = 0.05, df = 1, lower.tail = FALSE)
# method3_power   <- 1 - pchisq(method3_crit, ncp = method3_NCP,
#                               df = 1, lower.tail = TRUE); method3_power
#
# # Method 4 Power (Chi2-Distribution)
# method4_crit_chi2 <- qchisq(1 - 0.05, df = 2,
#                             ncp = 0, lower.tail = TRUE,
#                             log.p = FALSE)
# method4_power_chi2 <- 1 - pchisq(method4_crit_chi2, df = 2,
#                                  ncp = method4_NCP,
#                                  lower.tail = TRUE); method4_power_chi2
#
# # Method 4 Power (F-Distribution)
# method4_Fscore <- qf(1 - 0.05, df1 = 2, df2 = K*2 - 2*2, ncp = 0,
#                      lower.tail = TRUE, log.p = FALSE)
# method4_power_F <- 1 - pf(method4_Fscore, df1 = 2, df2 = K*2 - 2*2,
#                           method4_NCP, lower.tail = TRUE,
#                           log.p = FALSE); method4_power_F
#
# cat(" Method 3 NCP = ", method3_NCP, "\n",
#     "Method 4 NCP = ", method4_NCP, "\n",
#     "Method 3 Power (Chi2) = ", method3_power, "\n",
#     "Method 4 Power (Chi2) = ", method4_power_chi2, "\n",
#     "Method 4 Power (F)    = ", method4_power_F)

# Comparing Power as a function of NCP and Critical Value ----------------------

# Case when the NCP's are equal for Method 4 and 5
ncpEqualData <- expand.grid(alpha = c(0.01, 0.025, 0.05, 0.1),
                            NCP = seq(0, 30, by = 0.1)) %>%
  mutate(c3 = qchisq(p = alpha, df = 1,
                      ncp = 0, lower.tail = FALSE),
         c4 = qchisq(p = alpha, df = 2,
                      ncp = 0, lower.tail = FALSE)) %>%
  mutate(Method3 = 1 - pchisq(c3, ncp = NCP,
                                   df = 1, lower.tail = TRUE),
         Method4 = 1 - pchisq(c4, df = 2,
                                   ncp = NCP, lower.tail = TRUE))

ncpEqualDataLong <- ncpEqualData %>%
  pivot_longer(cols = c(Method3, Method4),
               values_to = "Power", names_to = "Method") %>%
  mutate(`Design Method` = str_replace(Method, "(\\D)(\\d)", "\\1 \\2")) %>%
  mutate(alpha_char = paste0(alpha))

ncpEqualDataLong$alpha_char <- factor(ncpEqualDataLong$alpha_char,
                                      labels = c('0.01' = parse(text = TeX('$\\alpha =$ 0.01')),
                                                 '0.025' = parse(text = TeX('$\\alpha =$ 0.025')),
                                                 '0.05' = parse(text = TeX('$\\alpha =$ 0.05')),
                                                 '0.1' = parse(text = TeX('$\\alpha =$ 0.1'))))

# Figure 3
ggplot(data = ncpEqualDataLong, aes(x = NCP,
                                    y = Power,
                                    color = `Design Method`)) +
  geom_point() + facet_wrap(~alpha_char, labeller = label_parsed) +
  xlab(TeX("Non-Centrality Parameter")) +
  ylab(TeX("Statistical Power")) +
  scale_color_manual(labels = c(TeX("Method 3"),
                                TeX("Method 4")),
                     values = c("blue", "violet")) +
  theme(text = element_text(size = 25))


MostPowerEqual <- mutate(ncpEqualData, `Highest Power` = ifelse(Method3 > Method4, "Method 3",
                                     ifelse(Method4 > Method3, "Method 4",
                                            ifelse(Method3 == Method4, "Equal", "ERROR")))) %>%
  group_by(`Highest Power`) %>%
  summarize(n = n())
MostPowerEqual

# Case when the NCP's are unequal for Method 3 and 4

# ncpUnequalData <- expand.grid(alpha = c(0.01, 0.025, 0.05, 0.1),
#                               NCP3 = seq(0, 25, by = 1),
#                               NCP4 = seq(0, 25, by = 1))

# NCP4Vector <- sort(seq(0, 25, by = 0.01), decreasing = TRUE)
# ncpUnequalData <- data.frame(alpha = c(rep(0.01, length(NCP4Vector)),
#                                         rep(0.025, length(NCP4Vector)),
#                                         rep(0.05, length(NCP4Vector)),
#                                         rep(0.1, length(NCP4Vector))),
#                               NCP3 = rep(seq(0, 25, by = 0.01), 4),
#                              NCP4 = rep(NCP4Vector, 4)) %>%

# ncpUnequalData <- expand.grid(alpha = c(0.01, 0.025, 0.05, 0.1),
#                                NCP3 = seq(0, 30, by = 0.5),
#                                NCP4 = seq(0, 30, by = 0.5)) %>%
#   dplyr::filter(NCP3 != NCP4) %>%
#   mutate(NCP_4minus3 = NCP4 - NCP3) %>%
#   mutate(c3 = qchisq(p = alpha, df = 1,
#                       ncp = 0, lower.tail = FALSE),
#          c4 = qchisq(p = alpha, df = 2,
#                       ncp = 0, lower.tail = FALSE)) %>%
#   mutate(Method3 = 1 - pchisq(c3, ncp = NCP3,
#                               df = 1, lower.tail = TRUE),
#          Method4 = 1 - pchisq(c4, df = 2,
#                               ncp = NCP4, lower.tail = TRUE)) %>%
#   mutate(alpha_char = paste0(alpha))
#
# ncpUnequalDataLong <- ncpUnequalData %>%
#   pivot_longer(cols = c(Method3, Method4),
#                values_to = "Power", names_to = "Method") %>%
#   mutate(`Design Method` = str_replace(Method, "(\\D)(\\d)", "\\1 \\2"))
#
# ncpUnequalDataLong$alpha_char <- factor(ncpUnequalDataLong$alpha_char,
#                                         labels = c('0.01' = parse(text = TeX('$\\alpha =$ 0.01')),
#                                                    '0.025' = parse(text = TeX('$\\alpha =$ 0.025')),
#                                                    '0.05' = parse(text = TeX('$\\alpha =$ 0.05')),
#                                                    '0.1' = parse(text = TeX('$\\alpha =$ 0.1'))))
#
# # Find where lines intersect
#
# # Finding the intersection points (where the difference is minimized)
# intersections <- ncpUnequalData %>%
#   mutate(PowerDifference = abs(Method3 - Method4)) %>%
#   group_by(alpha_char) %>%
#   slice(which.min(PowerDifference)) %>%
#   dplyr::select(alpha_char, NCP_4minus3, Method3, Method4, PowerDifference)
# intersections$alpha_char <- factor(intersections$alpha_char,
#                                         labels = c('0.01' = parse(text = TeX('$\\alpha =$ 0.01')),
#                                                    '0.025' = parse(text = TeX('$\\alpha =$ 0.025')),
#                                                    '0.05' = parse(text = TeX('$\\alpha =$ 0.05')),
#                                                    '0.1' = parse(text = TeX('$\\alpha =$ 0.1'))))
#
# intersections
#
# # Figure 4
# ggplot(data = ncpUnequalDataLong, aes(x = NCP_4minus3, y = Power, color = `Design Method`)) +
#   geom_point() +
#   facet_wrap(~alpha_char, labeller = label_parsed) +
#   geom_vline(data = intersections, aes(xintercept = NCP_4minus3),
#              linetype = "dashed", color = "black") +
#   geom_text(data = intersections, aes(x = NCP_4minus3,
#                                       y = max(Method3, Method4) + 0.02,
#                                       label = round(NCP_4minus3, 2)),
#             color = "black", angle = 0, vjust = 8, hjust = -0.5) +
#   ylab(TeX("Statistical Power")) +
#   xlab(TeX("$\\lambda^{(Method 4)} - \\lambda^{(Method 3)}$")) +
#   scale_color_manual(labels = c(TeX("Method 3"),
#                                 TeX("Method 4")),
#                      values = c("blue", "violet")) #+
#   #theme(text = element_text(size = 25))
#
#
#
# # ggplot(data = ncpUnequalDataLong, aes(x = NCP_4minus3, y = Power, color = Method)) +
# #   geom_point() + facet_wrap(~alpha)
#
# MostPowerUnequal <- mutate(ncpUnequalData, `Highest Power` = ifelse(Method3 > Method4, "Method 3",
#                                                            ifelse(Method4 > Method3, "Method 4",
#                                                                   ifelse(Method3 == Method4, "Equal", "ERROR")))) %>%
#   group_by(`Highest Power`) %>%
#   summarize(n = n())
#
# MostPowerUnequal
#
# filteredTest <- ncpUnequalData %>%
#   dplyr::filter(NCP4 > NCP3) %>%
#   mutate(`Highest Power` = ifelse(Method3 > Method4, "Method 3",
#                                   ifelse(Method4 > Method3, "Method 4",
#                                          ifelse(Method3 == Method4, "Equal", "ERROR"))))
#
#
# # Critical value for method 3, df = 1
# method3_crit <- qchisq(p = 0.05, df = 1, lower.tail = FALSE)
#
# # Method 3 Power
# method3_power   <- 1 - pchisq(method3_crit, ncp = method3_lambda,
#                               df = 1, lower.tail = TRUE)
#
# # Critical value for method 4, df = 2
# method4_crit_chi2 <- qchisq(1 - 0.05, df = 2,
#                             ncp = 0, lower.tail = TRUE, log.p = FALSE)
#
# # Method 4 Chi2 Power
# method4_power_chi2 <- 1 - pchisq(method4_crit_chi2, df = 2,
#                                  ncp = method4_lambda, lower.tail = TRUE)
#
#
#
#
# marcumq(a = sqrt(method3_NCP), b = sqrt(method3_crit), m = .5) # m is order
# marcumq(a = sqrt(method4_NCP), b = sqrt(method4_crit_chi2), m = 1)
#
#
#
#
#
# # -------
# # More DF's mean more power
# 1 - pchisq(2, ncp = 8,
#            df = 1, lower.tail = TRUE)
# 1 - pchisq(2, ncp = 8,
#            df = 5, lower.tail = TRUE)
#
# # Lower NCP means lower power
# 1 - pchisq(method3_crit, ncp = 20,
#            df = 1, lower.tail = TRUE)
#
# 1 - pchisq(method3_crit, ncp = 1,
#            df = 1, lower.tail = TRUE)
