library(yfR)
library(TDA)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(igraph)

# The following preliminary script replicates the results in "Topology Data Analysis of Critical Transitions in Financial Networks" by Marian Gidea

# Tickers (DJIA as of Feb 19, 2008, excluding GM)
tickers <- c("MMM", "AA", "AXP", "AIG", "T", "BAC", "BA", "CAT", "CVX", "C", "KO", "DD", 
             "DIS", "XOM", "GE", "HD", "HPQ", "IBM", "INTC", "JNJ", "JPM", "MCD", 
             "MRK", "MSFT", "PFE", "PG", "UTX", "VZ", "WMT")

# Closing prices
data <- yf_get(
  tickers = tickers,
  first_date = "2004-01-01",
  last_date = "2008-09-30",
  freq_data = "daily"
)


prices <- data %>%
  select(ticker, ref_date, price_adjusted) %>%
  pivot_wider(names_from = ticker, values_from = price_adjusted) %>%
  arrange(ref_date)

# Arithmetic returns
returns <- prices %>%
  mutate(across(-ref_date, ~ (.x - lag(.x)) / lag(.x)))

# Drop first row with NA returns
returns <- returns[-1, ]

# Rolling distance matrices (window T=15, sampled every 10 days)
T_window <- 15
sample_step <- 10
start_idx <- T_window + 1  # First valid window
indices <- seq(start_idx, nrow(returns), by = sample_step)

dates <- returns$ref_date[indices] 
dist_sub <- list()
dist_super <- list()

for (idx in indices) {
  window <- returns[(idx - T_window):idx, -1]  # T+1 returns
  corr <- cor(window, use = "pairwise.complete.obs")
  dist <- sqrt(2 * (1 - corr))  # For sub-level
  dist_sub[[length(dist_sub) + 1]] <- dist
  
  dist_super[[length(dist_super) + 1]] <- 2 - dist  # For super-level (dual)
}

# Now we generate persistence diagrams (dim 0 and 1 combined)
max_dim <- 1
max_scale <- 2.0  # Range [0,2]

diags_sub <- list()
diags_super <- list()

for (i in 1:length(dist_sub)) {
  # Sub-level P.D.
  diags_sub[[i]] <- ripsDiag(
    X = as.matrix(dist_sub[[i]]),
    maxdimension = max_dim,
    maxscale = max_scale,
    dist = "arbitrary",
    printProgress = FALSE
  )$diagram
  
  # Super-level P.D.
  diags_super[[i]] <- ripsDiag(
    X = as.matrix(dist_super[[i]]),
    maxdimension = max_dim,
    maxscale = max_scale,
    dist = "arbitrary",
    printProgress = FALSE
  )$diagram
}

# Plotting function for persistence diagram (dim 0 black, dim 1 red, infinity ♦)
plot_pers_diag <- function(diag, title) {
  plot.diagram(diag, diagLim = c(0, 2), main = title)
}

# Generate and display persistence diagrams (early vs. late, like Figures 3 and 5)
num_windows <- length(diags_sub)
early_indices <- 1:round(0.3 * num_windows)
late_indices <- round(0.7 * num_windows):num_windows

# Generate Sub-level diagrams (like Figure 3)
par(mfrow = c(2, 3))  # 2 rows, 3 columns
for (i in early_indices[1:3]) {
  plot_pers_diag(diags_sub[[i]], paste("Sub-level Early:", dates[i]))
}
for (i in late_indices[1:3]) {
  plot_pers_diag(diags_sub[[i]], paste("Sub-level Late:", dates[i]))
}

# Generate Super-level diagrams (like Figure 5)
par(mfrow = c(2, 3))  # Reset for new plot
for (i in early_indices[1:3]) {
  plot_pers_diag(diags_super[[i]], paste("Super-level Early:", dates[i]))
}
for (i in late_indices[1:3]) {
  plot_pers_diag(diags_super[[i]], paste("Super-level Late:", dates[i]))
}

# Step 6: Compute Wasserstein distances (p=2) to reference (first diagram)
ref_sub <- diags_sub[[1]]
ref_super <- diags_super[[1]]

dist_sub_dim0 <- sapply(2:length(diags_sub), function(i) {
  wasserstein(ref_sub, diags_sub[[i]], p = 2, dimension = 0)
})
dist_sub_dim1 <- sapply(2:length(diags_sub), function(i) {
  wasserstein(ref_sub, diags_sub[[i]], p = 2, dimension = 1)
})

dist_super_dim0 <- sapply(2:length(diags_super), function(i) {
  wasserstein(ref_super, diags_super[[i]], p = 2, dimension = 0)
})
dist_super_dim1 <- sapply(2:length(diags_super), function(i) {
  wasserstein(ref_super, diags_super[[i]], p = 2, dimension = 1)
})

# Plot distances like Figures 4 and 6
time_seq <- 1:length(dist_sub_dim0)  # Time index (starting from second window)
plot_data <- data.frame(
  time = rep(time_seq, 4),
  distance = c(dist_sub_dim0, dist_sub_dim1, dist_super_dim0, dist_super_dim1),
  type = rep(c("Sub Dim0 (Fig 4 Left)", "Sub Dim1 (Fig 4 Right)", "Super Dim0 (Fig 6 Left)", "Super Dim1 (Fig 6 Right)"), each = length(time_seq))
)

ggplot(plot_data, aes(x = time, y = distance)) +
  geom_line() +
  facet_wrap(~ type, scales = "free_y", nrow = 2) +
  labs(title = "Wasserstein Distances Between Persistent Diagrams",
       x = "Time Index (sampled every 10 trading days)",
       y = "Wasserstein Distance (p=2)") +
  theme_minimal()