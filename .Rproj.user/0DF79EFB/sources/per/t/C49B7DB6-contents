library(dplyr)
set.seed(111)

subject_random <- rnorm(500, mean = 0, sd = 5)
baseline_severity <- sample(1:5, 500, replace = TRUE, prob = c(0.1, 0.2, 0.4, 0.2, 0.1))
pgi_s <- data.frame(
  SSID = rep(1:500,5),
  timepoint = rep(NA,500)
)
pgi_s <- pgi_s%>%group_by(SSID)%>%
  mutate(timepoint = rep(1:5))
within_noise <- rnorm(nrow(pgi_s), mean = 0, sd = 1.2)

pgi_s$PGI_S <- round(pmin(pmax(baseline_severity + within_noise, 1), 5))

pgi_s$subject_intercept <- subject_random[pgi_s$SSID]
pgi_s$PROM <- round(
  100 - (pgi_s$PGI_S * 15) + 
    rnorm(nrow(pgi_s), mean = 0, sd = 10) + 
    pgi_s$subject_intercept,
  1
)
pgi_s <- pgi_s%>%select(-subject_intercept)
write.csv(pgi_s, "example_data/pgi_s.csv")
