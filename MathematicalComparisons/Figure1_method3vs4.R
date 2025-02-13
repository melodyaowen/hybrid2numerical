source("./RequiredPackages.R")

# Figure 2. Comparison of methods 3 and 4 when their non-centrality parameters
#           are equal

# Case when the NCP's are equal for Method 3 and 4
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

figure1 <- ggplot(data = ncpEqualDataLong,
                  aes(x = NCP, y = Power, color = `Design Method`)) +
  geom_point() + facet_wrap(~alpha_char, labeller = label_parsed) +
  xlab(TeX("Non-Centrality Parameter")) +
  ylab(TeX("Statistical Power")) +
  scale_color_manual(labels = c(TeX("Single Weighted 1-DF"),
                                TeX("Disjunctive 2-DF")),
                     values = c("blue", "violet")) +
  theme(text = element_text(size = 25))

# Saving figure 2 as png
ggsave(filename = "./MathematicalComparisons/Output/Figure1.png",
       plot = figure1,
       width  = 5000, height = 3000, units  = "px")

# Checking to ensure Method 3 is always more powerful (it is)
# Power method 3 = power method 4 only due to lack of significant digits in R
MostPowerEqual <- mutate(ncpEqualData, `Highest Power` = ifelse(Method3 > Method4, "Method 3",
                                     ifelse(Method4 > Method3, "Method 4",
                                            ifelse(Method3 == Method4, "Equal", "ERROR")))) %>%
  group_by(`Highest Power`) %>%
  summarize(n = n())
MostPowerEqual
