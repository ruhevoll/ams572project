# AMS 572 Project Question 1


#IMPORTANT: Be sure to run getwd() and ensure that fred_gas.csv, rbob_futures_hist.csv are present in your working directory before running this program.
library(readr)
library(dplyr)
library(lubridate)

# Retail gas prices dataframe
df_retail <- read_csv("fred_gas.csv") %>%
  rename(retail_price = APU000074714,
         date = observation_date) %>%
  mutate(date = mdy(date), 
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