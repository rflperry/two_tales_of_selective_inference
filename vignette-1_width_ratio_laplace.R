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

get_as_width <- function(nu, scale, n, alpha=0.05) {
  h <- 2 * qnorm(1 - nu / (2*n) )
  eta <- h / scale
  as_width <- 2 * sigma * qnorm(1 - (alpha - nu) * exp(-eta) / 2)
  return(as_width)
}

get_best_width_ratio <- function(var, n, alpha = 0.05) {
  
  # AS width
  nu_range <- seq(0, 0.05, length.out = 20)
  
  # scale <- sqrt( (1 - eps) / 2 / eps )
  scale <- sqrt( var / 2 )
  as_width <- min( mapply(get_as_width, nu_range, scale, n, alpha) )
  if (as_width == Inf) {
    return(NaN)
  }
  
  
  # (mean) DF width
  df_widths <- c()
  for (rep in seq(n_reps)) {
    set.seed(rep)
    
    # Simulate and fission data
    Y <- rnorm(n, mean = 0, sd = sigma)
    fission <- fission_data(Y, scale=scale)
    
    # Select winner
    i_hat <- which.max(fission$Y_train)

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
      df_widths <- c(df_widths, b - a)
    })
  }
  
  df_width <- mean(df_widths)
  return( as_width / df_width )
}


#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Construct 2d grid of n, epsilon values. Compare widths
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

library(ggplot2)
library(dplyr)


nu_range <- seq(0, 0.05, length.out = 20)
n_reps <- 100
sigma = 1
n_range <- 10^seq(1, 6)
var_range <- 2^(-seq(-4, 4, length.out=9))

grid_values <- expand.grid(var = var_range, n = n_range)
grid_values$ratio <- mapply(get_best_width_ratio, grid_values$var, grid_values$n)
grid_values <- grid_values %>% mutate(
  new_var = var + 1
)

g <- ggplot(
  # Forces AS width of Inf due to computational limits to be plotted as black,
  # vs any NAs plotted grey
  grid_values %>% mutate(
    temp_ratio = ratio, 
    ratio = ifelse(is.na(ratio), Inf, ratio),
    ratio = ifelse(!is.na(temp_ratio) & temp_ratio == Inf, NaN, ratio)
  ),
  aes(x = n, y = var, fill = ratio)) +
  geom_tile() +
  scale_fill_distiller(palette = "Reds", direction = 1, na.value="black") +
  labs(
    x = "n",
    y = "Added noise",
    fill = "Width Ratio\n(I&W / S&C)"
  ) +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() +
  theme(text = element_text(size = 11))
# theme_minimal()

ggsave("figures/vignette_1/vignette-1_laplace_width_ratio.png", width = 4, height = 2, unit = "in")
