library(tseries)
library(dplyr)
library(urca)
library(dynlm)
library(car)
library(ggplot2)

CPI <- read.csv("CPIAUCSL_W.csv")
CPI <- na.omit(CPI)
CPI$DATE <- as.Date(CPI$DATE)
CPI1970_2012 <- subset(CPI, DATE >= "1970-01-01" & DATE <= "2012-12-31")

inflation_rate <- (log(CPI1970_2012$CPIAUCSL) - log(lag(CPI1970_2012$CPIAUCSL, 1))) * 12 * 100

CPI1970_2012$infl <- inflation_rate

cpi_dec1969 <- CPI$CPIAUCSL[which(CPI$DATE == "1969-12-01")]
cpi_nov1969 <- CPI$CPIAUCSL[which(CPI$DATE == "1969-11-01")]
Infl_dec1969 <- (log(cpi_dec1969) - log(cpi_nov1969)) * 12 * 100
CPI1970_2012$infl[1] <- (log(CPI1970_2012$CPIAUCSL[1]) - log(cpi_dec1969)) * 12 * 100

plot(CPI1970_2012$DATE, CPI1970_2012$infl, 
     type = "l", xlab = "Date", ylab = "Inflation Rate", main = "Inflation Rate from 1970 to 2012")

acf(ts(CPI1970_2012$infl), lag.max = 12, plot = F, main = "Autocorrelation of Inflation Rate")
acf(ts(CPI1970_2012$infl), lag.max = 12, plot = T, main = "Autocorrelation of Inflation Rate")

odel1 <- lm(infl ~ lag(infl, 1), data = CPI1970_2012)
summary(Model1)

Model2 <- lm(infl ~ lag(infl, 1) + lag(infl, 2), data = CPI1970_2012)

summary(Model2)
IC <- function(model) {
  SSR <- sum(model$residuals^2)
  T <- length(model$residuals)
  p <- length(model$coef)  
  
  return(
    round(c("Lags" = p - 1,
            "AIC" = AIC(model),
            "BIC" = BIC(model),
            "R.Square" = summary(model)$r.squared), 4)
  )
}
max_lag <- 1:6

table <- sapply(max_lag, function(max_lag) {
  formula <- as.formula(paste("infl ~", paste("lag(infl,", 1:max_lag, ")", collapse = " + ")))
  
  model <- lm(formula, data = CPI1970_2012)
  IC(model)
})
IC_table <- as.data.frame(t(table))
print(IC_table)
best_by_AIC <- IC_table[which.min(IC_table$AIC),]
best_by_BIC <- IC_table[which.min(IC_table$BIC),]
best_by_R2 <- IC_table[which.max(IC_table$R.Square),]
cat("Best model based on AIC:\n")
print(best_by_AIC)

cat("Best model based on BIC:\n")
print(best_by_BIC)

cat("Best model based on R-Squared:\n")
print(best_by_R2)

inflation_nov_dec_2012 <- filter(CPI1970_2012, format(CPI1970_2012$DATE, "%Y-%m") %in% c("2012-11", "2012-12"))

infl_2012_12 <- inflation_nov_dec_2012 |>
  filter(format(DATE, "%m") == "12") |>
  pull(infl)  

infl_2012_11 <- inflation_nov_dec_2012 |>
  filter(format(DATE, "%m") == "11") |>
  pull(infl)

infl_2012_12
infl_2012_11

beta0 <- coef(Model2)[2]  
beta1 <- coef(Model2)[3]  

pred_infl_2013_01 <- beta0 * infl_2012_12 + beta1* infl_2012_11
pred_infl_2013_01


CPI1970_2012 <- CPI1970_2012 |>
  mutate(Delta_Infl = c(NA, diff(infl)))  
CPI1970_2012$Delta_Infl[1] <- CPI1970_2012$infl[1] - Infl_dec1969
CPI1970_2012 <- CPI1970_2012 %>%
  mutate(
    Delta_Infl_LAG1 = lag(Delta_Infl, 1),  
    Delta_Infl_LAG2 = lag(Delta_Infl, 2)   
  )

adf_test_without_t <- ur.df(CPI1970_2012 $Delta_Infl, type = "drift", lags = 2, selectlags = "Fixed")

summary(adf_test_without_t)

t <- 1:nrow(CPI1970_2012 )

eq02 <- lm(Delta_Infl ~ t + Delta_Infl_LAG1 + Delta_Infl_LAG2, data = CPI1970_2012 )

residuals_eq02 <- residuals(eq02)

adf_test_with_t <- ur.df(residuals_eq02, type = "drift", lags = 2, selectlags = "Fixed")
summary(adf_test_with_t)

# Ensure infl is a time series as dynlm deals with ts objects only
CPI1970_2012$infl <- ts(CPI1970_2012$infl, start = c(1970, 1), frequency = 12)

# Create a time vector from the time series
Annual_Infl_ts <- CPI1970_2012$infl
TV <- time(Annual_Infl_ts)

# 15% triming as required
w <- 0.15
NOB <- length(TV)
tau0 <- TV[round(NOB * w, 0)]
tau1 <- TV[round(NOB * (1 - w), 0)]
tau <- seq(tau0, tau1, 0.083) 

# Initialize vector for F-statistics
Fstats <- numeric(length(tau))

# Loop to compute F-statistics
for (i in 1 : length(tau) ) {
  # Create dummy variable D
  D <- ifelse(TV >= tau[i], 1, 0)
  # ADL2 model with all possible interactions
  eq <- dynlm(infl ~ L(infl)+L(infl,2) + 
                D*L(infl) , CPI1970_2012)
  coefs <- names(coef(eq))
  
  # compute and save the F-statistic
  Fstats[i] <- linearHypothesis(eq,c("D=0","L(infl):D"))$F[2]
  
}

QLR <- max(Fstats)
QLR

# identify the time period where the QLR-statistic is observed
as.yearqtr(tau[which.max(Fstats)])

dates <- CPI1970_2012$DATE   
forecast_start_date <- as.Date("2005-12-01")
forecast_end_date <- as.Date("2012-12-01")
forecast_period <- which(dates >= forecast_start_date & dates <= forecast_end_date)
estimation_period <- which(dates < forecast_start_date)
pseudo_forecasts <- numeric(length(forecast_period))
forecast_errors <- numeric(length(forecast_period))
for (i in 1:length(forecast_period)) {
  estimation_data <- CPI1970_2012[estimation_period, ]
  ar2_model <- dynlm(infl ~ L(infl, 1) + L(infl, 2), data = estimation_data)
  forecast_value <- predict(ar2_model, newdata = CPI1970_2012[forecast_period[i], , drop = FALSE])
  actual_value <- CPI1970_2012$infl[forecast_period[i]]
  pseudo_forecasts[i] <- forecast_value
  forecast_errors[i] <- actual_value - forecast_value
}
forecast_results <- data.frame(
  Date = dates[forecast_period],
  Actual = CPI1970_2012$infl[forecast_period],
  Forecast = pseudo_forecasts,
  Error = forecast_errors
)
print(forecast_results) #I was interested in seeing the results which were quite similar to the actual values
mean_ferror <- mean(forecast_errors)
t_test_result <- t.test(forecast_errors)
print(t_test_result) #
RMSFE <- sqrt(mean(forecast_errors^2))
cat("RMSFE:", RMSFE, "\n") 
