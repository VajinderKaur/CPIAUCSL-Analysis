#Part 3 
library(tseries)
library(dplyr)
library(urca)
library(dynlm)
library(car)
library(ggplot2)

BXP1 <- read.csv("BXP.csv")
AMT1 <- read.csv("AMT.csv")
BXP1$Date <- as.Date(BXP1$Date, format = "%d-%m-%Y")
AMT1$Date <- as.Date(AMT1$Date, format = "%d-%m-%Y")

BXP1$log_price <- log(BXP1$Price)
AMT1$log_price <- log(AMT1$Price)

BXP1$log_return <- c(NA, diff(log(BXP1$Price)))

AMT1$log_return <- c(NA, diff(log(AMT1$Price)))
BXP1 <- na.omit(BXP1)
AMT1 <- na.omit(AMT1)

merged_data1 <- merge(BXP1, AMT1, by = "Date", suffixes = c("_BXP", "_AMT"))

# Plot
ggplot(merged_data1, aes(x = Date)) +
  geom_line(aes(y = log_return_BXP, color = "BXP")) +
  geom_line(aes(y = log_return_AMT, color = "AMT")) +
  labs(title = "Log Returns of BXP and AMT", y = "Log Return", color = "Stock") +
  theme_minimal()

regression_model <- lm(log_return_BXP ~ log_return_AMT, data = merged_data1)
summary(regression_model)
# Extract residuals
merged_data1$residuals <- regression_model$residuals

plot(regression_model)

adf_test_residuals <- ur.df(merged_data1$residuals, type = "none", selectlags = "AIC")
summary(adf_test_residuals)

ecm_model <- lm(log_return_BXP ~ 
                  lag(log_return_BXP, 1) +   # Positive lag for spot return
                  lag(log_return_AMT, 1) +   # Positive lag for futures return
                  (log_return_BXP - log_return_AMT) +  # Cointegration term
                  residuals, 
                data = merged_data1)

# Summary of the ECM model
summary(ecm_model)

merged_data1$z_score <- (merged_data1$residuals - mean(merged_data1$residuals)) / sd(merged_data1$residuals)

ggplot(merged_data1, aes(x = Date, y = z_score)) +
  geom_line() +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  geom_hline(yintercept = -2, linetype = "dashed", color = "blue") +
  labs(title = "Z-Scores of Residuals for Pair-Trading", y = "Z-Score") +
  theme_minimal()

threshold <- 0.02  

merged_data1 <- merged_data1 %>%
  mutate(
    signal = case_when(
      residuals > threshold ~ "Sell BXP and Buy AMT",
      residuals < -threshold ~ "Buy BXP and Sell AMT",
      TRUE ~ "Hold"
    )
  )

ggplot(merged_data1, aes(x = Date)) +
  geom_line(aes(y = residuals), color = "blue") +
  geom_point(aes(y = residuals, color = signal)) +
  labs(title = "Pair-Trading Strategy Signals for BXP and AMT",
       y = "Residuals",
       color = "Signal") +
  theme_minimal()