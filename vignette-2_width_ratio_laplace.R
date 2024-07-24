#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Construct 2d grid of n, epsilon values. Compare widths
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

library(ggplot2)
library(dplyr)

NAN_KEY <- 99999

delta_range <- seq(0, 1, length.out = 100)
hyperparam_grid <- expand.grid(delta = delta_range)

n <- 100
n_reps <- 250
sparsity <- 0.5

get_as_width <- function(scale, delta, p, sigma=1, alpha=0.05) {
  if (delta > 1 | delta < 0) {
    return(Inf)
  }
  # c <- sqrt(eps / ( 1 - eps ))
  eta <- 2 * qnorm(1 - alpha * delta / 2 / p) / scale
  
  # Compute base AS width, two sided!!!
  as_width <- sigma * qnorm(1 - alpha * (1 - delta) * exp(-eta) / 2)
  
  return(as_width)
}

get_width_ratio <- function(var, p, rho=0.5, sigma=1, alpha=0.05, signal=0.2) {
  # Compute c
  scale <- sqrt(var / 2)
  # Compute base AS width
  as_width <- 2 * min( mapply(get_as_width, scale, hyperparam_grid$delta, p) )
  if (as_width == Inf) {
    return(Inf)
  }
  if (as_width == NAN_KEY) {
    return(NA)
  }
  
  # Compute DF width
  Sigma_X <- matrix(0, p, p) + rho
  diag(Sigma_X) <- 1
  
  df_widths <- c()
  
  for (rep in seq(n_reps)) {
    set.seed(rep)
    
    # Simulate X, Y
    X <- mvrnorm(n, rep(0, p), Sigma_X)
    X <- scale(X, center=FALSE, scale=col.norm(X))
    
    # Simulate and fission data
    if (signal == 0) {
      betas <- rep(0, p)
      Y <- rnorm(n, mean = 0, sd = sigma)
    } else {
      betas <- c( rexp(floor(p * sparsity), rate = 1 / signal), rep(0, ceiling(p * ( 1 - sparsity )) ))
      Y <- X %*% betas + rnorm(n, mean = 0, sd = sigma)
    }
    noise <- rlaplace(p, scale=scale)
    
    Y_train <- t(X) %*% Y + noise
    Y_test <- t(X) %*% Y
    
    Sigma_Y <- sigma^2 * t(X) %*% X

    # Select winner
    i_hat <- which.max(Y_train)
    
    # Approach #5: Lee & Taylor correct...
    lb = -20
    ub = 20
    
    c <- Sigma_Y[,i_hat] # eqtn 5.3
    z <- Y_test - c * Y_test[i_hat] # eqtn 5.2
    
    # polyhedron matrices
    A <- do.call(rbind, list(
      rep(0, p),
      diag(p),
      -diag(p)
    ))
    A[, i_hat] <- -1
    b <- -A %*% noise
    
    A <- A * sign(Y_train[i_hat])
    b <- b * sign(Y_train[i_hat])
    
    vmin <- (b - A %*% z) / (A %*% c) # eqtn 5.4
    if (all((A %*% c) > 0)) {
      vmin <- -Inf
    } else {
      vmin <- max(vmin[(A %*% c) < 0])
    }
    
    vmax <- (b - A %*% z) / (A %*% c) # eqtn 5.5
    if (all((A %*% c) < 0)) {
      vmax <- Inf
    } else {
      vmax <- min(vmax[(A %*% c) > 0])
    }
    
    upper <- optimize({
      function(mu) (
        ptruncnorm( q = Y_test[i_hat], a = vmin, b = vmax, mean = mu, sd = sigma) - (alpha / 2))^2
    }, lower=lb, upper=ub)$minimum
    lower <- optimize({
      function(mu) (
        ptruncnorm( q = Y_test[i_hat], a = vmin, b = vmax, mean = mu, sd = sigma) - ( 1 - alpha / 2))^2
    }, lower=lb, upper=ub)$minimum
    
    df_width <- upper - lower
    df_widths <- c(df_widths, df_width)
  }
  
  
  return(as_width / mean(df_widths))
}

p_range <- c(25, 50, 100, 200, 400)
var_range <- 10^(seq(-2, 2, length.out = 5))

rho_range <- c(0.5)
signal_range <- c(0)

for (rho in rho_range) {
  for (signal in signal_range) {

    grid_values <- expand.grid(var = var_range, p = p_range)
    grid_values$ratio <- mapply(get_width_ratio, grid_values$var, grid_values$p, rho=rho, signal=signal)
    
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
      scale_y_log10() +
      scale_x_continuous(trans='log2') +
      # scale_x_log10() +
      theme_bw() +
      theme(
        text = element_text(size = 11),
      )
    
    ggsave(paste0("figures/vignette_2/vignette-2_laplace_width_ratio_rho=", rho, "_signal=", signal, ".png"), width = 4, height = 2, unit = "in")
  }
}
