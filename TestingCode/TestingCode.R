

testData <- expand.grid(b1 = seq(1, 10, by = 0.1),
                        b2 = seq(1, 10, by = 0.1),
                        v1 = seq(0, 1, by = 0.01),
                        v2 = seq(0, 1, by = 0.01)) %>%
  mutate(val1 = round((b1^2)*v2, 6),
         val2 = round((b2^2)*v1, 6),
         compare = ifelse(val1 < val2, "Yes", ifelse(val1 > val2, "No", "Equal"))) %>%
  mutate(val1sqrt = round(b1*sqrt(v2), 6),
         val2sqrt = round(b2*sqrt(v1), 6),
         comparesqrt = ifelse(val1sqrt < val2sqrt, "Yes", ifelse(val1sqrt > val2sqrt, "No", "Equal")))
View(head(testData))

%>%
  filter(compare != comparesqrt)

nrow(testData)

View(testData)

all(testData$compare == testData$comparesqrt)

nrow(testData)
View(head(testData))



alpha_B <- alpha/2
cv_B <- qchisq(1 - alpha_B, df = 1, ncp = 0)

alpha <- 0.05

# Smaller alpha, smaller critical value, smaller power

qchisq(1 - alpha/2, df = 1, ncp = 0)
qchisq(1 - alpha, df = 2, ncp = 0)

1 - pchisq(qchisq(1 - 0.05, df = 2, ncp = 0), df = 2, ncp = 8, lower.tail = TRUE)
1 - pchisq(qchisq(1 - 0.05, df = 1, ncp = 0), df = 1, ncp = 8, lower.tail = TRUE)


1 - pchisq(6, df = 1, ncp = 8, lower.tail = TRUE)
1 - pchisq(3, df = 1, ncp = 8, lower.tail = TRUE)

1 - pchisq(qchisq(1 - 0.05, df = 1, ncp = 0), 1, ncp = 8, lower.tail = TRUE)
1 - pchisq(qchisq(1 - 0.05/2, df = 1, ncp = 0), 1, ncp = 8, lower.tail = TRUE)
