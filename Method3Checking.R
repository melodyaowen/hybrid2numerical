# Checking my derivation for Method 3
beta1 <- 0.1
beta2 <- 0.2
varY1 <- 0.25
varY2 <- 0.5
rho01 <- 0.08
rho02 <- 0.05
rho1 <- 0.01
rho2 <- 0.05
K <- 10
m <- 50
VIF1_new <- 1 + (m - 1)*rho01
VIF2_new <- 1 + (m - 1)*rho02
VIF12_new <- rho2 + (m - 1)*rho1

# Method 3 NCP
method3_NCP <- ((sqrt((beta1^2)/(((2*varY1)/(K*m))*VIF1_new)) +
                   sqrt((beta2^2)/(((2*varY2)/(K*m))*VIF2_new)))/
                  (sqrt(2*(1+(VIF12_new)/(sqrt(VIF1_new*VIF2_new))))))^2; method3_NCP

# Method 4 NCP
method4_NCP <- (K*m*((beta1^2)*varY2*VIF2_new -
                       2*beta1*beta2*sqrt(varY1*varY2)*VIF12_new +
                       (beta2^2)*varY1*VIF1_new))/
  (2*varY1*varY2*(VIF1_new*VIF2_new - VIF12_new^2)); method4_NCP

# Method 3 Power
method3_crit <- qchisq(p = 0.05, df = 1, lower.tail = FALSE)
method3_power   <- 1 - pchisq(method3_crit, ncp = method3_NCP,
                              df = 1, lower.tail = TRUE); method3_power


(1/2)*(erfc((sqrt(method3_crit)-sqrt(method3_NCP))/(sqrt(2))) + erfc((sqrt(method3_crit)+sqrt(method3_NCP))/(sqrt(2))))


# Method 4 Power
method4_crit_chi2 <- qchisq(1 - 0.05, df = 2,
                            ncp = 0, lower.tail = TRUE,
                            log.p = FALSE)
method4_power_chi2 <- 1 - pchisq(method4_crit_chi2, df = 2,
                                 ncp = method4_NCP,
                                 lower.tail = TRUE); method4_power_chi2


exp(-method4_NCP/2)*


marcumq(a = sqrt(method3_NCP), b = sqrt(method3_crit), m = .5) # m is order
marcumq(a = sqrt(method4_NCP), b = sqrt(method4_crit_chi2), m = 1)

marcumq(a = 9, b = 3.84, m = .5) # m is order
marcumq(a = 9, b = 5.99, m = 1)

marcumq(a = 4, b = 2, m = .5)
marcumq(a = 4, b = 2, m = 1)
marcumq(a = 4, b = 2, m = 2)
marcumq(a = 4, b = 2, m = 3)
marcumq(a = 4, b = 2, m = 4)

marcumq(a = 4, b = 1, m = 1)
marcumq(a = 4, b = 2, m = 1)
marcumq(a = 4, b = 3, m = 1)
marcumq(a = 4, b = 4, m = 1)

grid <- expand.grid(seq(0, 10, by = 0.5), seq(0, 10, by = 0.5), seq(0, 30, by = 1))
marcumDat <- data.frame(c3 = grid$Var1,
                        c4 = grid$Var2,
                        lambda = grid$Var3) %>%
  mutate(c4_minus_c3 = c4 - c3,
         Method3 = marcumq(a = sqrt(lambda),
                           b = sqrt(c3),
                           m = .5),
         Method4 = marcumq(a = sqrt(lambda),
                           b = sqrt(c4),
                           m = 1))

marcumDatLong <- marcumDat %>%
  pivot_longer(cols = c(Method3, Method4),
               values_to = "Power", names_to = "Method")

ggplot(data = marcumDatLong, aes(x = c4_minus_c3, y = Power, color = Method)) +
  geom_point()



testData <- data.frame(s = seq(0, 500, by = 1)) %>%
  rowwise() %>%
  mutate(oneVal = (1/(factorial(s)^2))*((method4_NCP/2)^s)*gammainc(s + 1, method4_crit_chi2/2)[[2]])
sum(testData$oneVal)

s <- 1
(1/(factorial(s)^2))*((method4_NCP/2)^s)*gammainc(s + 1, method4_crit_chi2/2)[[2]]
