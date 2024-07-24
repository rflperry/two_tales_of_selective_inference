library(ggplot2)
library(MASS)
library(tilting)
library(VGAM)
library(truncnorm)
library(tmvtnorm)

get_as_width <- function(scale, delta, p, sigma=1, alpha=0.05) {
  if (delta > 1 | delta < 0) {
    return(Inf)
  }
  eta <- 2 * qnorm(1 - alpha * delta / 2 / p) / scale

  # Compute base AS width
  as_width <- sigma * qnorm(1 - alpha * (1 - delta) * exp(-eta) / 2)
  
  return(as_width)
}

n_reps <- 1000
n <- 100
p <- 100
sigma = 1
sparsity <- 0.5

alphas <- c(0.0125, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)

delta_range <- seq(0, 1, length.out = 20)
# Generate grid of parameter combinations
grid_values <- expand.grid(delta = delta_range)

# rho_range <- c(0.1, 0.5, 0.9)
eps_range <- c(0.25, 0.75)
signal_range <- c(0, 0.2, 0.14)

rho <- 0
# eps <- 0.5
# signal <- 0.14

for (signal in signal_range) {
  for (eps in eps_range) {
    # For AS and naive
    alpha <- 0.05

    # Laplacian var = 2 (scale)^2, so preserve wrt epsilon
    scale <- sqrt( (1 - eps) / 2 / eps )
    
    # correlation matrix for sampling X
    Sigma_X <- matrix(0, p, p) + rho
    diag(Sigma_X) <- 1
    
    # AS width
    as_width <- min( mapply(get_as_width, scale, grid_values$delta, p, sigma, alpha) )
    
    # Naive CI width
    naive_width <- sigma * qnorm(1 - alpha/2)
    
    # Check coverages of methods
    df_covered <- c()
    as_coverage <- 0
    naive_coverage <- 0
    df_widths_alpha <- vector("list", length(alphas))
    df_coverages_alpha <- vector("list", length(alphas))
    
    # For optimal infer-and-widen
    biases <- c()

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
      i_hat <- which.max( abs( Y_train ) )
      mu_i <- (t(X) %*% X %*% betas)[i_hat]
      
      # Check naive coverage
      if ( mu_i < Y_test[i_hat] + naive_width & mu_i > Y_test[i_hat] - naive_width  ) {
        naive_coverage <- naive_coverage + 1 / n_reps
      }
      
      # Get bias of the infer_widen point estimate
      bias <- abs(mu_i - Y_test[i_hat])
      biases <- c(biases, bias)
      
      # Check if AS covered
      if ( mu_i < Y_test[i_hat] + as_width & mu_i > Y_test[i_hat] - as_width  ) {
        as_coverage <- as_coverage + 1 / n_reps
      }
      
      # Approach: Lee & Taylor corrected
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
      
      for (alpha_i in seq(length(alphas))) {
        alpha <- alphas[alpha_i]
        try({
          upper <- optimize({
            function(mu) (
              ptruncnorm( q = Y_test[i_hat], a = vmin, b = vmax, mean = mu, sd = sigma) - (alpha / 2))^2
          }, lower=lb, upper=ub)$minimum
          lower <- optimize({
            function(mu) (
              ptruncnorm( q = Y_test[i_hat], a = vmin, b = vmax, mean = mu, sd = sigma) - ( 1 - alpha / 2))^2
          }, lower=lb, upper=ub)$minimum
          
          
          
          df_widths_alpha[[alpha_i]] <- c(df_widths_alpha[[alpha_i]], upper - lower)
          
          # Check if DF covered
          if ( mu_i > lower & mu_i < upper  ) {
            df_coverages_alpha[[alpha_i]] <- c(df_coverages_alpha[[alpha_i]], 1)
          } else {
            df_coverages_alpha[[alpha_i]] <- c(df_coverages_alpha[[alpha_i]], 0)
          }
        })
      }
      
    }
    
    df_widths <- c()
    for (alpha_i in seq(length(alphas))) {
      df_widths <- c(df_widths, mean(df_widths_alpha[[alpha_i]]))
    }
    
    df_coverages <- c()
    for (alpha_i in seq(length(alphas))) {
      df_coverages <- c(df_coverages, mean(df_coverages_alpha[[alpha_i]]))
    }
    
    coverages <- seq(1, n_reps) / (n_reps)

    g <- ggplot() +
      geom_line(
        data = data.frame(
          sorted_points = 2 * sort(abs(biases)),
          percentages = coverages
        ),
        aes(x = sorted_points, y = percentages, color = "Oracle I&W"),
        linetype = "solid") +
      geom_hline(yintercept = 0.95, linetype = "dotted", color = "grey") +
      geom_point(data = data.frame(x = 2 * as_width, y = as_coverage), aes(x, y, color = "Infer-and-widen"), shape = 15, size=2) +
      geom_line(data = data.frame(x = df_widths, y = df_coverages), aes(x, y, color = "Split-and-condition"), linetype = 'dashed') +
      geom_point(data = data.frame(x = 2 * naive_width, y = naive_coverage), aes(x, y, color = "Classic"), shape = 17, size=2) +
      scale_color_manual(name = "",
                         values = c("Infer-and-widen" = "blue", "Split-and-condition" = "forestgreen", "Classic" = "red", "Oracle I&W" = "black"),
                         labels = c("Infer-and-widen" = "Infer-and-widen", "Split-and-condition" = "Split-and-condition", "Classic" = "Classic", "Oracle I&W" = "Oracle I&W")
      ) +
      labs(x = "CI Width", y = "Coverage", colour="") +
      theme_bw() +
      theme(text = element_text(size = 11)) +
      theme(legend.position = "right") +
      # scale_x_continuous(limits = c(0, 10)) +
      scale_y_continuous(limits = c(0, 1.01), n.breaks=3) +
      guides(color = guide_legend(override.aes = list(
        shape = c("Classic" = 17, "Infer-and-widen" = 15, "Oracle I&W" = NA, "Split-and-condition" = NA),
        linetype = c("Infer-and-widen" = 0, "Classic" = 0, "Oracle I&W" = 1, "Split-and-condition" = 2)
      ))) +
      coord_flip()
    
    ggsave(paste0("figures/vignette_2/vignette-2_oracle-curves_laplace_eps=", eps, "_rho=", rho, "_signal=", signal, '.png'), width = 3.25, height = 1.75, unit = "in")
  }
}
