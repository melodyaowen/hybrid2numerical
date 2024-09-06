# Figure 2: Comparison of the Marcum Q-functions 1/2 and 1

# Creating Marcum Q data to compare and find threshold
ncpEqualData <- expand.grid(c3 = seq(0, 15, by = 0.05),
                            c4 = seq(0, 15, by = 0.05),
                            NCP = seq(0, 15, by = 5)) %>%
  mutate(Method3 = 1 - pchisq(c3, ncp = NCP,
                              df = 1, lower.tail = TRUE),
         Method4 = 1 - pchisq(c4, df = 2,
                              ncp = NCP, lower.tail = TRUE)) %>%
  rowwise() %>%
  mutate(Inequality = ifelse(Method3 > Method4, "Method 3 > Method 4",
                            ifelse(Method4 > Method3, "Method 3 < Method 4",
                                   ifelse(Method3 == Method4, "Equal", "ERROR")))) %>%
  mutate(c4_minus_c3 = c4 - c3,
         power4_minus_power3 = Method4 - Method3) %>%
  dplyr::filter(c4 > c3) %>% # Only cases where c4 > c3
  mutate(NCP_Character = paste0("a = ", NCP))

ncpEqualData_filtered <- ncpEqualData %>%
  dplyr::filter(NCP %in% c(0, 5, 10, 15)) %>%
  dplyr::filter(c4_minus_c3 <= 4)

ncpEqualData_filtered$NCP_Character <- factor(ncpEqualData_filtered$NCP_Character,
                                              levels = c("a = 0", "a = 5", "a = 10", "a = 15"))

# Figure 2
ggplot(data = ncpEqualData_filtered, #dplyr::filter(ncpEqualData, NCP == 5),
       aes(x = c4_minus_c3,
           y = power4_minus_power3,
           color = Inequality)) +
  geom_point() + facet_wrap(~NCP_Character) +
  ylab(TeX("$Q_1(a, b_2) - Q_{1/2}(a, b_1)$")) +
  xlab(TeX("$b_2 - b_1$")) +
  scale_color_manual(labels = c(TeX("$Q_{1/2}(a, b_1) < Q_1(a, b_2)$"),
                                TeX("$Q_1(a, b_2) < Q_{1/2}(a, b_1)$")),
                     values = c("blue", "violet")) +
  theme(text = element_text(size = 25))





# Extra crap that I tried ---

ggplot(data = ncpEqualData, #dplyr::filter(ncpEqualData, NCP == 5),
       aes(x = c4_minus_c3,
           y = power4_minus_power3,
           color = NCP,
           group = NCP)) +
  geom_point()

ggplot(data = ncpEqualData,
       aes(x = c4_minus_c3)) +
  geom_histogram() + facet_wrap(~MostPower)

ncpEqualDataLong <- ncpEqualData %>%
  pivot_longer(cols = c(Method3, Method4),
               values_to = "Power", names_to = "Method")

ggplot(data = ncpEqualDataLong,
       aes(x = c4_minus_c3, y = NCP, group = Method)) +
  geom_point() + facet_wrap(~MostPower)

View(dplyr::filter(ncpEqualData, MostPower != "Equal"))

View(dplyr::filter(ncpEqualData, c4 > c3))




mosaic::favstats(c4_minus_c3 ~ MostPower, data = ncpEqualData)

# Making sure using MarcumQ function is the same
# ncpEqualData <- expand.grid(c3 = seq(0, 50, by = 1),
#                             c4 = seq(0, 50, by = 1),
#                             NCP = seq(0, 100, by = 1)) %>%
#   mutate(Method3 = 1 - pchisq(c3, ncp = NCP,
#                               df = 1, lower.tail = TRUE),
#          Method4 = 1 - pchisq(c4, df = 2,
#                               ncp = NCP, lower.tail = TRUE),
#          Method3_MQ = marcumq(a = sqrt(NCP), b = sqrt(c3), m = .5), # m is order
#          Method4_MQ = marcumq(a = sqrt(NCP), b = sqrt(c4), m = 1)) %>%
#   mutate(Method3_Check = ifelse(round(Method3, 5) == round(Method3_MQ, 5), "Same", "DIFFERENT"),
#          Method4_Check = ifelse(round(Method4, 5) == round(Method4_MQ, 5), "Same", "DIFFERENT")) %>%
#   dplyr::select(c3, c4, NCP,
#                 Method3, Method3_MQ, Method3_Check,
#                 Method4, Method4_MQ, Method4_Check)
