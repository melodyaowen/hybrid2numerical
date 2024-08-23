require(latex2exp)

critValueData34 <- data.frame(alpha = seq(from = 0.001, to = 0.999, by = 0.001)) %>%
  mutate(`c4` = qchisq(p = alpha, df = 1,
                     ncp = 0, lower.tail = FALSE),
         `c5` = qchisq(p = alpha, df = 2,
                     ncp = 0, lower.tail = FALSE))

critValueData34Long <- critValueData34 %>%
  pivot_longer(cols = c(c4, c5),
               values_to = "Value", names_to = "Distribution")

ggplot(data = critValueData34Long, aes(x = alpha, y = Value,
                                       color = Distribution)) +
  scale_color_manual(labels = c(TeX("Central $\\chi^2$ with 1-DF"),
                                TeX("Central $\\chi^2$ with 2-DF")),
                     values = c("violet", "blue")) +
  geom_line(linewidth = 1.5) +
  xlab(TeX("Overall False-Positive Rate ($\\alpha$)")) +
  theme(text=element_text(size=25))

View(TestData)

