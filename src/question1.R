# AMS 572 Project Question 1

## !!!! Quantmod may not always work, we should consider downloading the data manually, storing as .csv (for our final submission) so there is no need for an internet connection
library(quantmod)
library(tidyverse)

# You can set the randomness to a certain "seed" so our sim results are reproducible 
set.seed(42)

# Download data
start_date <- "2021-01-01"
tickers <- c("SPY", "BTC-USD")
data_env <- new.env()
getSymbols(tickers, src = "yahoo", from = start_date, auto.assign = TRUE, env = data_env)


log_returns <- do.call(merge, lapply(data_env, function(x) dailyReturn(Ad(x), type = "log"))) # Compute log returns
colnames(log_returns) <- c("spy", "btc")

# Data cleaning below
clean_data <- as.data.frame(log_returns) %>% 
  na.omit() %>%
  filter(is.finite(spy) & is.finite(btc)) # Remove Inf values

print(paste("Baseline N:", nrow(clean_data)))

# We create a function so we can easily repeat the test 3 times
analyze_correlation <- function(df, label) {
  # Handling Method: Complete Case Analysis (na.omit)
  # We drop the rows where data is missing.
  df_clean <- na.omit(df)
  
  test <- cor.test(df_clean$spy, df_clean$btc, method = "pearson")
  
  return(data.frame(
    Scenario = label,
    N_Remaining = nrow(df_clean),
    Correlation = round(test$estimate, 4),
    P_Value = format.pval(test$p.value, digits = 3)
  ))
}


# Baseline residuals 
res_baseline <- analyze_correlation(clean_data, "1. Baseline (Full Data)")


# MCAR 
# Delete 15% of SPY data, I guess this is like a packet loss in downloading the data from the internet
# This implies no bias, just loss of sample size

mcar_data <- clean_data
n_rows <- nrow(mcar_data)
missing_indices <- sample(1:n_rows, size = 0.15 * n_rows) # Randomly select 15% of indices to set to NA
mcar_data$spy[missing_indices] <- NA

res_mcar <- analyze_correlation(mcar_data, "2. MCAR (Random Deletion)")


# MNAR
# For MNAR, we delete 80% of the bottom 10% of SPY returns

mnar_data <- clean_data
# Identify the bottom 10% of SPY returns (crashes)
threshold <- quantile(mnar_data$spy, 0.10)

# Delete 80% of the data ONLY if it is below that threshold
crash_indices <- which(mnar_data$spy < threshold)
delete_indices <- sample(crash_indices, size = 0.8 * length(crash_indices))
mnar_data$spy[delete_indices] <- NA

res_mnar <- analyze_correlation(mnar_data, "3. MNAR (Missing Low Values)")


# Results:
final_results <- rbind(res_baseline, res_mcar, res_mnar)
print(final_results)
