# AMS 572 Project Question 1

## !!!! Quantmod may not always work, we should consider downloading the data manually, storing as .csv (for our final submission) so there is no need for an internet connection
# (for our final submission)
library(quantmod)
library(tidyverse)
library(gridExtra) 
library(zoo) 

# You can set the randomness to a certain "seed" so our sim results are reproducible 
set.seed(42)

options(scipen = 999)

# Creates directory in getwd() directory to store plots
dir.create("project_plots", showWarnings = FALSE)

# Download data
start_date <- "2023-11-01" 
tickers <- c("SPY", "BTC-USD")
data_env <- new.env()

print("Downloading data.")
tryCatch({
  getSymbols(tickers, src = "yahoo", from = start_date, auto.assign = TRUE, env = data_env)
}, error = function(e) {
  stop("Error downloading data. Bad internet or typo in ticker name.")
})

log_returns <- do.call(merge, lapply(data_env, function(x) dailyReturn(Ad(x), type = "log")))
colnames(log_returns) <- c("spy", "btc")

# Data cleaning below

clean_data <- as.data.frame(log_returns) %>% 
  rownames_to_column(var = "date") %>%
  mutate(date = as.Date(date)) %>%
  na.omit() %>%
  filter(is.finite(spy) & is.finite(btc)) # Remove Inf values

print(paste("Baseline N:", nrow(clean_data)))


# Data Exploration 

graphics.off()

# Prepare Data
cumulative_data <- clean_data %>%
  mutate(Growth_SPY = exp(cumsum(spy)), Growth_BTC = exp(cumsum(btc))) %>%
  pivot_longer(cols = c("Growth_SPY", "Growth_BTC"), names_to = "Asset", values_to = "Value")

long_returns <- clean_data %>%
  pivot_longer(cols = c("spy", "btc"), names_to = "Asset", values_to = "Log_Return")

# Growth of $1 invested in S&P and BTC
p1 <- ggplot(cumulative_data, aes(x = date, y = Value, color = Asset)) +
  geom_line(size = 1) + scale_color_manual(values = c("Growth_BTC" = "orange", "Growth_SPY" = "blue")) +
  labs(title = "1. Growth of $1", y = "Value ($)") + theme_minimal() + theme(legend.position = "bottom")
ggsave("project_plots/q1_1_growth.png", plot = p1, width = 10, height = 6)

# Daily volatility plots
p2 <- ggplot(long_returns, aes(x = date, y = Log_Return, color = Asset)) +
  geom_line(alpha = 0.8) + scale_color_manual(values = c("btc" = "orange", "spy" = "blue")) +
  labs(title = "2. Daily Volatility", y = "Log Return") + theme_minimal()
ggsave("project_plots/q1_2_volatility.png", plot = p2, width = 10, height = 6)

# Correlation scatter plots
p3 <- ggplot(clean_data, aes(x = spy, y = btc)) +
  geom_point(alpha = 0.3, color = "darkblue") + geom_smooth(method = "lm", color = "red") +
  labs(title = "3. Correlation Scatter", x = "SPY", y = "BTC") + theme_minimal()
ggsave("project_plots/q1_3_scatter.png", plot = p3, width = 8, height = 6)

# Rolling correlation plots
clean_data$rolling_corr <- rollapply(clean_data[,c("spy", "btc")], width = 30,
                                     function(x) cor(x[,1], x[,2]), by.column = FALSE, fill = NA)
p4 <- ggplot(clean_data, aes(x = date, y = rolling_corr)) +
  geom_line(color = "purple", size = 1) + geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "4. Rolling 30-Day Correlation", y = "Correlation") + theme_minimal()
ggsave("project_plots/q1_4_rolling.png", plot = p4, width = 10, height = 6)

# Return distribution boxplot
p5 <- ggplot(long_returns, aes(x = Asset, y = Log_Return, fill = Asset)) +
  geom_boxplot() +
  labs(title = "5. Distribution of Returns", subtitle = "Compare spread/outliers") + theme_minimal()
ggsave("project_plots/q1_5_boxplot.png", plot = p5, width = 8, height = 6)

# Diagnostic plots (for checking normality assumptions)

# Linearity
png("project_plots/diag_linearity.png", width = 800, height = 600)
plot(clean_data$spy, clean_data$btc,
     main = "BTC vs SPY Log Returns",
     xlab = "SPY Log Returns", ylab = "BTC Log Returns",
     pch = 20, col = rgb(0,0,0,0.2))
abline(lm(btc ~ spy, clean_data), col = "red")
dev.off()

# Q-Q Plot SPY
png("project_plots/diag_qq_spy.png", width = 800, height = 600)
qqnorm(clean_data$spy, main = "Q-Q SPY")
qqline(clean_data$spy, col = "red")
dev.off()

# Q-Q Plot BTC
png("project_plots/diag_qq_btc.png", width = 800, height = 600)
qqnorm(clean_data$btc, main = "Q-Q BTC")
qqline(clean_data$btc, col = "red")
dev.off()

# ACF for BTC
png("project_plots/diag_acf_btc.png", width = 800, height = 600)
acf(clean_data$btc, main = "ACF BTC")
dev.off()


print("Plots saved.")

# We create a function so we can easily repeat the test 3 times
analyze_correlation <- function(df, label) {
  # Handling Method: Complete Case Analysis (na.omit)
  # We drop the rows where data is missing.
  df_clean <- na.omit(df)
  test <- cor.test(df_clean$spy, df_clean$btc, method = "pearson")
  return(data.frame(Scenario = label, N = nrow(df_clean), Correlation = round(test$estimate, 4), P_Value = format.pval(test$p.value, digits=3)))
}

# Baseline residuals
res_base <- analyze_correlation(clean_data, "1. Baseline")

# MCAR 
# Delete 15% of SPY data, I guess this is like a packet loss in downloading the data from the internet
# This implies no bias, just loss of sample size


mcar_data <- clean_data; mcar_data$spy[sample(nrow(mcar_data), 0.15*nrow(mcar_data))] <- NA
res_mcar <- analyze_correlation(mcar_data, "2. MCAR")

# MNAR
# For MNAR, we delete 80% of the bottom 10% of SPY returns

# Identify the bottom 10% of SPY returns (crashes)
mnar_data <- clean_data; thresh <- quantile(mnar_data$spy, 0.10, na.rm=TRUE)
mnar_data$spy[which(mnar_data$spy < thresh)[sample(sum(mnar_data$spy < thresh, na.rm=TRUE), 0.8*sum(mnar_data$spy < thresh, na.rm=TRUE))]] <- NA
res_mnar <- analyze_correlation(mnar_data, "3. MNAR")

print(rbind(res_base, res_mcar, res_mnar))
print("Below is the correlation coefficient between BTC and SPY:")
print(cor(clean_data$spy, clean_data$btc))