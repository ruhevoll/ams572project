# AMS 572 Project Question 1


#IMPORTANT: Be sure to run getwd() and ensure that fred_gas_tx.csv, rbob_futures_hist.csv are present in your working directory before running this program.
library(readr)
library(dplyr)
library(lubridate)
library(tseries)

# Retail gas prices dataframe
df_retail <- read_csv("fred_gas_tx.csv") %>%
  rename(retail_price = APUS37A74714,
         date = observation_date) %>%
  mutate(date = ymd(date),
         YearMonth = floor_date(date, "month")) %>%
  select(YearMonth, retail_price)

# Futures data dataframe
df_futures <- read_csv("rbob_futures_hist.csv") %>%
  rename(rbob_price = Price,
         date = Date) %>%
  mutate(date = mdy(date), 
         YearMonth = floor_date(date, "month")) %>%
  select(YearMonth, rbob_price)

# Need to convert RBOB daily into RBOB monthly (rbob_futures_hist.csv is daily data)
# For this, we just take the monthly mean
df_futures_monthly <- df_futures %>%
  group_by(YearMonth) %>%
  summarise(rbob_price = mean(rbob_price, na.rm = TRUE))

# Combine data sets together and remove rows which have missing data
df_merged <- inner_join(df_retail, df_futures_monthly, by = "YearMonth") %>%
  na.omit() # Remove any rows with missing data

# To study the relationship between RBOB and the retail rprice, we model the linear regression
# retail_price = beta_0 + beta_1 * rbob_price

model <- lm(retail_price ~ rbob_price, data = df_merged)

# Print the summary to see R-squared, t-statistic, etc.
# We want to test H_0: beta_1 = 0 vs. H_1 : beta =/= 0
summary(model)

# Here's a plot of the regression line we can use in the report:
# Will edit the formatting later so the plot looks nicer for the report. 
plot(retail_price ~ rbob_price, data = df_merged)
abline(model, col = "red")

# --------------------------------------------------------------------------------------

# PART 2

# --------------------------------------------------------------------------------------

# To justify our OLS model, we want to test if RBOB and retail prices are cointegrated (I think this is part of assuming no-arbitrage). So, in this part
# we test to see if the linear combination is stationary. 

# We'll test if retail_price - beta_1 * rbob_price is a stationary process N(0, s)
# Testing H0: Residuals are non-stationary vs. H1: Residuals are stationary
# For this, we use the Augmented Dickey-Fuller test (to be explained in our report)

model_residuals <- resid(model)
adf_test_result <- adf.test(model_residuals)
print(adf_test_result)

# Result of the ADF test: Residuals are stationary (reject Ho w/ p-value 0.01)
# For the paper, I'm also generating a plot of the residuals to show visually the residuals are stationary

df_merged_with_resid <- df_merged
df_merged_with_resid$residuals <- model_residuals

plot(df_merged_with_resid$YearMonth, df_merged_with_resid$residuals,
     type = "l", # 'l' for line plot
     col = "blue",
     main = "Residuals",
     xlab = "Date",
     ylab = "Residual Value ($)")
# Add a horizontal line at 0
abline(h = 0, col = "red", lty = 2, lwd = 2) 
