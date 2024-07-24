library(glmnet)
library(ggplot2)
library(MASS)
library(tilting)

n_reps <- 1000
n <- 100
sigma = 1
alpha <- 0.05

sparsity <- 0.5

get_width_ratio <- function(var, p, n=100, rho=0.5, sigma=1, alpha=0.05, signal=0.2, sparsity=0.5) {
  c <- sqrt(var)
  
  # correlation matrix for sampling X
  Sigma_X <- matrix(0, p, p) + rho
  diag(Sigma_X) <- 1
  
  # Check coverages of methods
  dt_coverage <- 0
  dt_widths <- c()
  
  # For Oracle I&C
  biases <- c()
  
  for (rep in seq(n_reps)) {
    set.seed(rep)
    
    # Simulate X, Y
    X <- mvrnorm(n, rep(0, p), Sigma_X)
    X <- scale(X, center=FALSE, scale=col.norm(X))

    # Simulate and thin data
    if (signal == 0) {
      betas <- rep(0, p)
      Y <- rnorm(n, mean = 0, sd = sigma)
    } else {
      betas <- c( rexp(floor(p * sparsity), rate = 1 / signal), rep(0, ceiling(p * ( 1 - sparsity )) ))
      Y <- X %*% betas + rnorm(n, mean = 0, sd = sigma)
    }

    noise <- rnorm(n, mean = 0, sd = sigma)
    
    Y_train <- Y + c * noise
    Y_test <- Y - 1 / c * noise
    
    # Apply lasso to Y_train
    lasso_fit <- cv.glmnet(X, Y_train, nfolds=3)
    coeffs <- coef(lasso_fit, s = "lambda.min")[2:(p+1)]
    i_all <- which(coeffs != 0)
    if (length(i_all) == 0) {
      # If doesn't select anything, then we "cover" the target correctly
      if (runif(1) >= alpha) {
        dt_coverage <- dt_coverage + 1 / n_reps
      }
      next
    }
    i_winner <- which.max(abs(coeffs[i_all]))
    
    # Compute "true_beta"
    Z <- X[,i_all]
    lm_subset <- lm((X %*% betas) ~ Z)
    beta_winner <- coef(lm_subset)[i_winner + 1]
    
    # Fit naive to Y
    lm_naive <- lm(Y ~ X[,i_all])
    beta_winner_naive <- coef(lm_naive)[i_winner + 1]

    # Get bias of the infer_widen point estimate
    bias <- abs(beta_winner_naive - beta_winner)
    biases <- c(biases, bias)
    
    # Fit DT to Y
    lm_dt <- lm(Y_test ~ X[,i_all])

    # Since \sigma^2 is known
    beta_winner_dt <- coef(lm_dt)[i_winner + 1]
    scale_sd <- diag(solve(t(cbind(rep(1, n), X[,i_all])) %*% cbind(rep(1, n), X[,i_all])))[i_winner + 1]
    ci_dt <- c(
      beta_winner_dt - qnorm(1 - alpha / 2) * sigma * sqrt(scale_sd) * sqrt(1 + 1 / var),
      beta_winner_dt + qnorm(1 - alpha / 2) * sigma * sqrt(scale_sd) * sqrt(1 + 1 / var)
        )
    dt_widths <- c(dt_widths, ci_dt[2] - ci_dt[1])
    if ( beta_winner < ci_dt[2] & beta_winner > ci_dt[1]  ) {
      dt_coverage <- dt_coverage + 1 / n_reps
    }
  }

  runs <- length(dt_widths)
  coverages <- seq(1, runs) / runs
  dt_width <- mean(dt_widths)
  oracle_width <- 2*quantile(biases, 1-alpha)  # dt_width is two-sided
  
  return(oracle_width / dt_width)
}

p_range <- c(25, 50, 100, 200)
var_range <- c(1/4, 1/2, 1, 2, 4, 8) # 10^(seq(-2, 2, length.out = 5)) 
grid_values <- expand.grid(var = var_range, p = p_range)
signals <- c(0.33)
rho <- 0.5

for (signal in signals) {
  print(paste("Signal: ", signal))
  grid_values$ratio <- mapply(get_width_ratio, var=grid_values$var, p=grid_values$p, rho=rho, signal=signal, n=n)
  
  grid_values %>% mutate(
    temp_ratio = ratio,
    ratio = ifelse(is.na(ratio), Inf, ratio)
  )
  
  plot_df <- grid_values %>% mutate(
    temp_ratio = ratio, 
    ratio = ifelse(is.na(ratio), Inf, ratio),
    ratio = ifelse(!is.na(temp_ratio) & temp_ratio == Inf, NaN, ratio),
    ratio = log2(ratio)
  )
  g <- ggplot(
    # Forces AS width of Inf due to computational limits to be plotted as black,
    # vs any NAs plotted grey
    plot_df,
    aes(x = p, y = var, fill = ratio)) +
    geom_tile() +
    scale_fill_distiller(
      palette = "RdBu",
      na.value="black",
      limits=c(-1,1)*max(abs(plot_df$ratio), na.rm=TRUE),
      # limits=c(-4, 4),
      # trans="log2",
      direction = -1
    ) +
    labs(
      x = "p",
      y = "Added noise",
      fill = "Width Ratio\n(I&W / S&C) "
    ) +
    scale_y_continuous(trans="log2") +
    scale_x_continuous(trans='log2') +
    # scale_x_log10() +
    theme_bw() +
    theme(
      text = element_text(size = 11),
    )
  
  ggsave(paste0("figures/vignette_3/vignette-3_width_ratio_heatmap_rho=", rho, "_signal=", signal, ".png"), width = 3.5, height = 2.2, unit = "in")
}

