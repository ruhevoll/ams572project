# AMS 572 Project Question 2

# ***Please run question1.R before running this file.***
#EV : I'm assuming this means the expected value of the model's error
#So, we want to test if EV = 0, i.e., if the model is unbiased and fairly priced. I think
#the below code answers that question.

# When I run the code, I get a 95% CI of [-0.195088, 0.362646], which contains zero. So, we fail to reject
# H_0: EV = 0. 

# Add model residuals to dataframe
df_analysis <- df_merged %>%
  mutate(residuals = model$residuals)


# Bootstrap parameters


block_size <- 12  # 12-month periods
n_obs <- nrow(df_analysis)
n_replications <- 10000 # Sample size 


# List of all possible 12-month blocks, added to the dataframe
all_blocks <- list()
for (i in 1:(n_obs - block_size + 1)) {
  all_blocks[[i]] <- df_analysis[i:(i + block_size - 1), ]
}


# Simulation


n_total_blocks <- length(all_blocks)
bootstrap_evs <- numeric(n_replications) # Stores EV for each simulation

for (i in 1:n_replications) {
  # Randomly sample ONE 12-month block (with replacement)
  sampled_block <- all_blocks[[sample(1:n_total_blocks, 1, replace = TRUE)]]
  
  # Calculate the EV (mean residual) for that 12-month period
  bootstrap_evs[i] <- mean(sampled_block$residuals)
}


# Results
# We're testing H_0: EV = 0, H_1: EV =/= 0.


cat("Bootstrap parameters\n")
cat("Replications:", n_replications, "\n")
cat("Block Size:", block_size, "months\n\n")

mean_ev <- mean(bootstrap_evs) # Mean EV
ci_lower <- quantile(bootstrap_evs, 0.025) # Lower CI for EV
ci_upper <- quantile(bootstrap_evs, 0.975) # Upper CI for EV

cat("Average 12-Month EV:", round(mean_ev, 6), "\n")
cat("95% EV Confidence Interval: [", round(ci_lower, 6), ", ", round(ci_upper, 6), "]\n\n")

# Histogram of the bootstrap results
hist(bootstrap_evs, breaks = 50,
     main = "Bootstrap Distribution of 12-Month EV",
     xlab = "Mean Residual ($)",
     col = "lightblue")
# Add lines for the CI
abline(v = ci_lower, col = "red", lty = 2, lwd = 2)
abline(v = ci_upper, col = "red", lty = 2, lwd = 2)
abline(v = 0, col = "blue", lty = 1, lwd = 2)
legend("topright", legend = c("95% CI", "EV = 0"), col = c("red", "blue"), lty = c(2, 1))