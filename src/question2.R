setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # Set working directory to current file location

library(quantmod)
library(tidyverse)
library(broom)
library(lmtest)   
library(sandwich) 
library(car)
library(reshape2) 

set.seed(42)
options(scipen = 999)
dir.create("project_plots", showWarnings = FALSE)

start_date <- "2023-11-01" 

native_coins <- c("ETH-USD", "SOL-USD", "XRP-USD", "ADA-USD", "AVAX-USD", 
                  "DOT-USD", "LINK-USD", "POL-USD", "LTC-USD", "ATOM-USD")
meme_coins <- c("DOGE-USD", "SHIB-USD", "PEPE-USD", "WIF-USD", "BONK-USD", 
                "FLOKI-USD", "MEME-USD", "ELON-USD", "BABYDOGE-USD", "SAFEMOON-USD")
btc_ticker <- "BTC-USD"

# Load local data
print("Loading stored CSV instead of downloading.")
all_df <- read.csv("stored_q2_data.csv")
all_df$date <- as.Date(all_df$date)

# Filter for Continuity
limit <- 0.05 * nrow(all_df)
keep_cols <- colSums(is.na(all_df)) < limit
all_df <- all_df[, keep_cols]
print(paste("Coins kept:", ncol(all_df)-1))

btc_col_name <- grep("BTC", colnames(all_df), value = TRUE)
if (length(btc_col_name) == 0) stop("BTC missing")

btc_data_clean <- all_df %>% select(date, all_of(btc_col_name)) %>% rename(log_return_btc = all_of(btc_col_name))
stacked_data <- all_df %>% select(-all_of(btc_col_name)) %>% pivot_longer(cols = -date, names_to = "coin", values_to = "log_return_coin")

final_data <- stacked_data %>%
  left_join(btc_data_clean, by = "date") %>%
  mutate(
    clean_ticker = gsub("\\.USD|-USD", "", coin),
    is_memecoin = ifelse(coin %in% meme_coins | clean_ticker %in% gsub("-USD","",meme_coins), 1, 0),
    is_memecoin = factor(is_memecoin, levels = c(0, 1), labels = c("Native", "Memecoin"))
  ) %>%
  na.omit() %>% 
  filter(is.finite(log_return_coin) & is.finite(log_return_btc))

print(paste("Final N:", nrow(final_data)))

# Plots
print("Generating 8 Plots...")
graphics.off()

growth_data <- all_df %>% pivot_longer(cols = -date, names_to = "coin_ticker", values_to = "log_ret") %>%
  na.omit() %>% group_by(coin_ticker) %>% arrange(date) %>% mutate(growth_of_1 = exp(cumsum(log_ret))) %>% ungroup() %>%
  mutate(clean_ticker = gsub("\\.USD|-USD", "", coin_ticker),
         Category = case_when(coin_ticker == btc_ticker ~ "Bitcoin", coin_ticker %in% meme_coins | clean_ticker %in% gsub("-USD","",meme_coins) ~ "Memecoin", TRUE ~ "Native"))

# Native coins growth
nat_sub <- growth_data %>% filter(Category %in% c("Bitcoin", "Native"))
p1 <- ggplot(nat_sub, aes(x=date, y=growth_of_1, color=coin_ticker)) +
  geom_line(size=0.8) + geom_line(data=filter(nat_sub, coin_ticker==btc_ticker), color="black", size=1.5) +
  labs(title="1. Native vs BTC Growth", y="Growth of $1") + theme_minimal()
ggsave("project_plots/q2_1_native_growth.png", plot=p1, width=12, height=7)

# Meme coins growth
mem_sub <- growth_data %>% filter(Category %in% c("Bitcoin", "Memecoin"))
p2 <- ggplot(mem_sub, aes(x=date, y=growth_of_1, color=coin_ticker)) +
  geom_line(size=0.8) + geom_line(data=filter(mem_sub, coin_ticker==btc_ticker), color="black", size=1.5) +
  labs(title="2. Meme vs BTC Growth", y="Growth of $1") + theme_minimal()
ggsave("project_plots/q2_2_meme_growth.png", plot=p2, width=12, height=7)

# Native scatter grid
p3 <- final_data %>% filter(is_memecoin=="Native") %>%
  ggplot(aes(x=log_return_btc, y=log_return_coin, color=coin)) + geom_point(alpha=0.3) + geom_smooth(method="lm", color="black") +
  facet_wrap(~coin, scales="free_y") + labs(title="3. Native Scatter Grid") + theme_bw() + theme(legend.position="none")
ggsave("project_plots/q2_3_native_grid.png", plot=p3, width=12, height=8)

# Meme scatter grid
p4 <- final_data %>% filter(is_memecoin=="Memecoin") %>%
  ggplot(aes(x=log_return_btc, y=log_return_coin, color=coin)) + geom_point(alpha=0.3) + geom_smooth(method="lm", color="black") +
  facet_wrap(~coin, scales="free_y") + labs(title="4. Meme Scatter Grid") + theme_bw() + theme(legend.position="none")
ggsave("project_plots/q2_4_meme_grid.png", plot=p4, width=12, height=8)

# Interaction overlay
p5 <- ggplot(final_data, aes(x=log_return_btc, y=log_return_coin, color=is_memecoin)) +
  geom_smooth(method="lm", se=TRUE, size=1.5) +
  labs(title="5. Interaction Effect: Comparing Slopes", subtitle="Comparing Leverage (Beta)", x="BTC Return", y="Coin Return") +
  theme_minimal()
ggsave("project_plots/q2_5_interaction_overlay.png", plot=p5, width=8, height=6)

# Correlation heatmap
cormat <- cor(all_df %>% select(-date)); melted <- melt(cormat)
p6 <- ggplot(melted, aes(x=Var1, y=Var2, fill=value)) + geom_tile() +
  scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0.5, limit=c(0,1)) +
  theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1)) + labs(title="6. Heatmap")
ggsave("project_plots/q2_6_heatmap.png", plot=p6, width=10, height=8)

# Distribution violin
p7 <- ggplot(final_data, aes(x=is_memecoin, y=log_return_coin, fill=is_memecoin)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") +
  labs(title="7. Return Distribution", subtitle="Wider violin = Fatter Tails") + theme_minimal()
ggsave("project_plots/q2_7_violin.png", plot=p7, width=8, height=6)

# BASELINE MODEL (unchanged)
baseline_model <- lm(log_return_coin ~ log_return_btc * is_memecoin, data = final_data)

# Diagnostics
png("project_plots/q2_diag_residuals_vs_fitted.png", width=900, height=700)
plot(baseline_model, which = 1)
dev.off()

png("project_plots/q2_diag_qqplot.png", width=900, height=700)
qqPlot(baseline_model, main="Normal Q-Q Plot")
dev.off()

png("project_plots/q2_diag_scale_location.png", width=900, height=700)
plot(baseline_model, which = 3)
dev.off()

png("project_plots/q2_diag_acf.png", width=900, height=700)
acf(baseline_model$residuals, main="ACF of Residuals")
dev.off()

# Interaction term extraction
run_mlr <- function(df, label) {
  df <- na.omit(df); if (length(levels(factor(df$is_memecoin))) < 2) return(NULL)
  model <- lm(log_return_coin ~ log_return_btc * is_memecoin, data = df)
  res <- tidy(model) %>% filter(grepl(":", term))
  return(data.frame(Scenario=label, Beta=round(res$estimate, 4), P=format.pval(res$p.value, digits=3)))
}

# Baseline / MCAR / MNAR
res_base <- run_mlr(final_data, "1. Baseline")
mcar <- final_data; meme_i <- which(mcar$is_memecoin=="Memecoin"); mcar$log_return_coin[sample(meme_i, 0.2*length(meme_i))] <- NA
res_mcar <- run_mlr(mcar, "2. MCAR")
mnar <- final_data; thr <- quantile(mnar$log_return_coin, 0.05); crash <- which(mnar$is_memecoin=="Memecoin" & mnar$log_return_coin < thr)
mnar$log_return_coin[sample(crash, 0.9*length(crash))] <- NA
res_mnar <- run_mlr(mnar, "3. MNAR")

print("--- Simulation Results ---")
print(rbind(res_base, res_mcar, res_mnar))

# HAC Standard Errors
robust_test <- coeftest(baseline_model, vcov = vcovHAC(baseline_model))
int_row <- robust_test[grep(":", rownames(robust_test)), ]

print("--- Comparison of Interaction Term ---")
print(paste("Original P-Value:", format.pval(summary(baseline_model)$coefficients[4,4], digits=3)))
print(paste("Robust P-Value:  ", format.pval(int_row[4], digits=3)))
