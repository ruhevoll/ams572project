# AMS 572 Project Question 2

library(quantmod)
library(tidyverse)
library(broom)

# Ensure reproducibility
set.seed(42)
options(scipen = 999)

# Extract interaction term after some data cleaning (making sure each category has enough data)
run_mlr <- function(df, label) {
  df_clean <- na.omit(df)

  if (length(levels(factor(df_clean$is_memecoin))) < 2) {
    warning(paste("Skipping", label, ": Not enough data in one category."))
    return(NULL)
  }
  
  model <- lm(log_return_coin ~ log_return_btc * is_memecoin, data = df_clean)
  stats <- tidy(model)
  
  target_row <- stats[grep(":", stats$term), ]
  
  return(data.frame(
    Scenario = label,
    N_Remaining = nrow(df_clean),
    Interaction_Beta = round(target_row$estimate, 4),
    Std_Error = round(target_row$std.error, 5),
    P_Value = format.pval(target_row$p.value, digits = 3)
  ))
}


# Download the data

start_date <- "2021-01-01"
native_coins <- c("ETH-USD", "SOL-USD", "XRP-USD", "ADA-USD", "AVAX-USD", 
                  "DOT-USD", "LINK-USD", "POL-USD", "LTC-USD", "ATOM-USD")
meme_coins <- c("DOGE-USD", "SHIB-USD", "PEPE-USD", "WIF-USD", "BONK-USD", 
                "FLOKI-USD", "MEME-USD", "ELON-USD", "BABYDOGE-USD", "SAFEMOON-USD")
btc_ticker <- "BTC-USD"
all_tickers <- c(btc_ticker, native_coins, meme_coins)

data_env <- new.env()
valid_tickers <- c()

print("Starting data download...")
for (ticker in all_tickers) {
  tryCatch({
    getSymbols(ticker, src = "yahoo", from = start_date, env = data_env, auto.assign = TRUE)
    if (nrow(get(ticker, envir = data_env)) > 0) valid_tickers <- c(valid_tickers, ticker)
  }, error = function(e) {})
}

# Calculate Returns
log_returns_list <- list()
for (ticker in valid_tickers) {
  d <- get(ticker, envir = data_env)
  d <- na.omit(d) 
  try({
    r <- dailyReturn(Ad(d), type = "log")
    colnames(r) <- ticker
    log_returns_list[[ticker]] <- r
  }, silent = TRUE)
}

# Merge 
all_log_returns <- do.call(merge, c(log_returns_list, all = TRUE))
all_df <- as.data.frame(all_log_returns) %>% rownames_to_column("date") %>% as_tibble()
colnames(all_df) <- gsub(".Adjusted", "", colnames(all_df))

# Extract BTC returns (Predictor)
btc_col_name <- grep("BTC", colnames(all_df), value = TRUE)
if (length(btc_col_name) == 0) stop("Bitcoin data is missing!")

# Rename the BTC column explicitly by name to preserve 'date' column
btc_data_clean <- all_df %>% 
  select(date, all_of(btc_col_name)) %>%
  rename(log_return_btc = all_of(btc_col_name))

# Stack data
stacked_data <- all_df %>% 
  select(-all_of(btc_col_name)) %>% 
  pivot_longer(cols = -date, names_to = "coin", values_to = "log_return_coin")

# Final Baseline Data Creation
final_data <- stacked_data %>%
  left_join(btc_data_clean, by = "date") %>%
  mutate(
    clean_ticker = gsub("\\.USD|-USD", "", coin),
    is_memecoin = ifelse(coin %in% meme_coins | clean_ticker %in% gsub("-USD","",meme_coins), 1, 0),
    is_memecoin = factor(is_memecoin, levels = c(0, 1), labels = c("Native", "Memecoin"))
  ) %>%
  na.omit() %>% # Remove NAs
  filter(is.finite(log_return_coin) & is.finite(log_return_btc)) # Remove Inf

print(paste("Baseline Data Loaded. N =", nrow(final_data)))


# Now we run the analysis
res_baseline <- run_mlr(final_data, "1. Baseline (Full Data)")


# MCAR
# Randomly lose 20% of Memecoin returns
mcar_data <- final_data
meme_indices <- which(mcar_data$is_memecoin == "Memecoin")
remove_indices <- sample(meme_indices, size = 0.2 * length(meme_indices))
mcar_data$log_return_coin[remove_indices] <- NA

res_mcar <- run_mlr(mcar_data, "2. MCAR (Random Memecoin Loss)")


# MNAR 
# Delete 90% of the "very" negative returns (0.05 percentile)
mnar_data <- final_data
threshold <- quantile(mnar_data$log_return_coin, 0.05)
crash_indices <- which(mnar_data$is_memecoin == "Memecoin" & mnar_data$log_return_coin < threshold)


delete_indices <- sample(crash_indices, size = 0.9 * length(crash_indices))
mnar_data$log_return_coin[delete_indices] <- NA

res_mnar <- run_mlr(mnar_data, "3. MNAR (Crash Data Missing)")

# Filter out any NULL results if a scenario failed the factor check
final_results <- list(res_baseline, res_mcar, res_mnar) %>%
  bind_rows()

# Final results:
print(final_results)