library(ggplot2)
library(VGAM)
library(truncnorm)

fission_data <- function(Y, rep, scale) {
  noise <- rlaplace(length(Y), scale=scale)
  fission <- list()
  fission$Y_train <- Y + noise
  fission$Y_test <- Y
  return( fission )
}

get_as_width <- function(nu, scale, alpha=0.05) {
  h <- 2 * qnorm(1 - nu / (2*n) )
  eta <- h / scale
  
  # Compute base AS width, 1-sided
  as_width <- sigma * qnorm(1 - (alpha - nu) * exp(-eta) / 2)

  return(as_width)
}

n_reps <- 10000
n <- 100
sigma = 1
alphas <- c(0.0125, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)

nu_range <- seq(0, 0.05, length.out = 20)

eps_range <- c(0.25, 0.75)

for (eps in eps_range) {
  
  alpha <- 0.05
  
  # Laplacian var = 2 (scale)^2, so preserve wrt epsilon
  scale <- sqrt( (1 - eps) / 2 / eps )
  
  # AS width
  as_width <- min( mapply(get_as_width, nu_range, scale, alpha) )
  
  # Naive CI width
  naive_width <- sigma * qnorm(1 - alpha/2)
  
  # Check coverages of methods
  df_coverage <- 0
  as_coverage <- 0
  naive_coverage <- 0
  df_widths_alpha <- vector("list", length(alphas))
  df_coverages_alpha <- vector("list", length(alphas))

  # For optimal infer-and-widen
  biases <- c()

  for (rep in seq(n_reps)) {
    set.seed(rep)
    
    # Simulate and fission data
    Y <- rnorm(n, mean = 0, sd = sigma)
    fission <- fission_data(Y, scale=scale)
    
    # Select winner
    i_hat <- which.max(fission$Y_train)
    
    # Check naive coverage
    if ( 0 < Y[i_hat] + naive_width & 0 > Y[i_hat] - naive_width  ) {
      naive_coverage <- naive_coverage + 1 / n_reps
    }
    
    # Get bias of the infer_widen point estimate
    bias <- fission$Y_test[i_hat]
    biases <- c(biases, bias)
    
    # Check if AS covered
    if ( 0 < Y[i_hat] + as_width & 0 > Y[i_hat] - as_width  ) {
      as_coverage <- as_coverage + 1 / n_reps
    }

    for (alpha_i in seq(length(alphas))) {
      alpha <- alphas[alpha_i]

      # Get DF width and check if covered
      try({
        if (fission$Y_train[i_hat] < Y[i_hat]) {
          lb = -10
          ub = 10
          b <- optimize( {
            function(mu) (
              ptruncnorm( q = Y[i_hat], a = fission$Y_train[i_hat], mean = mu, sd = sigma) - (alpha / 2))^2
            }, lower=lb, upper=ub)$minimum
          a <- optimize( {
            function(mu) (
              ptruncnorm( q = Y[i_hat], a = fission$Y_train[i_hat], mean = mu, sd = sigma) - ( 1 - alpha / 2))^2
            }, lower=lb, upper=ub)$minimum
          
          a <- a + sigma^2 / scale
          b <- b + sigma^2 / scale
        } else {
          lb = -10
          ub = 10
          b <- optimize({
            function(mu) (
              ptruncnorm( q = Y[i_hat], b = fission$Y_train[i_hat], mean = mu, sd = sigma) - (alpha / 2))^2
            }, lower=lb, upper=ub)$minimum
          a <- optimize({
            function(mu) (
              ptruncnorm( q = Y[i_hat], b = fission$Y_train[i_hat], mean = mu, sd = sigma) - ( 1 - alpha / 2))^2
            }, lower=lb, upper=ub)$minimum
          
          a <- a - sigma^2 / scale
          b <- b - sigma^2 / scale
        }
  
        df_widths_alpha[[alpha_i]] <- c(df_widths_alpha[[alpha_i]], b - a)

        # Check if DF covered
        if ( 0 > a & 0 < b  ) {
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
  df_width <- mean(df_widths)

  # Create the plot
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
                       values = c("Infer-and-widen" = "blue", "Split-and-condition" = "forestgreen", "Classic" = "red", "Oracle I&W" = "black")
                       # labels = c("Infer-and-widen" = "Infer-and-widen", "Split-and-condition" = "Split-and-condition", "Classic" = "Classic", "Oracle I&W" = "Oracle I&W")
    ) +
    labs(x = "CI Width", y = "Coverage", colour="") +
    theme_bw() +
    theme(text = element_text(size = 11)) +
    theme(legend.position = "right") +
    # scale_x_continuous(limits = c(0, 10)) +
    scale_y_continuous(limits = c(0, 1), n.breaks=3) +
    guides(color = guide_legend(override.aes = list(
      shape = c("Classic" = 17, "Infer-and-widen" = 15, "Oracle I&W" = NA, "Split-and-condition" = NA),
      linetype = c("Infer-and-widen" = 0, "Classic" = 0, "Oracle I&W" = 1, "Split-and-condition" = 2)
    ))) +
    coord_flip()
  ggsave(paste0("figures/vignette_1/vignette-1_oracle-curves_laplace_eps=", eps, '.png'), width = 3.25, height = 1.75, unit = "in")
}