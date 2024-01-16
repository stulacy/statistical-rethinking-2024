library(rethinking)

# Question 1 --------------------------------------------------------------
# Obtain the posterior of p given 3W and 11L
grid_approx_globe <- function(grid_size, prior, W, L) {
    # R code 2.3
    p_grid <- seq(from=0, to=1, length.out=grid_size)
    likelihood <- dbinom(W, size=W+L, prob=p_grid)
    unstd.posterior <- likelihood * prior(p_grid)
    posterior <- unstd.posterior / sum(unstd.posterior)
    plot(p_grid, posterior, type="b", xlab="p", ylab="P(p|D)")
    mtext(sprintf("%d points", grid_size))
    return(posterior)
}

# I'll use a uniform prior here, but here's some others shown in chapter 2
uniform_prior <- function(p) rep(1, length(p))
step_prior <- function(p) ifelse(p < 0.5, 0, 1)
sharp_prior <- function(p) exp(-5*abs(p-0.5))

grid_size <- 1e3
post <- grid_approx_globe(grid_size, uniform_prior, 3, 11)

# Question 2 --------------------------------------------------------------
# Compute the posterior predictive for the next 5 tosses
p_grid <- seq(0, 1, length.out=grid_size)
samples <- sample(p_grid, prob=post, size=1e5, replace=TRUE)
posterior_predictive <- rbinom(1e5, size=5, prob=samples)
simplehist(posterior_predictive)

# Question 3 --------------------------------------------------------------
# Observe W=7 but L (and hence N) are unknown
# If p=0.7, which is posterior for N?
# Same principle as calculating posterior for p, N is just another parameter 
# of the binomial
# Use grid approximation but this time the grid is for N
N_grid <- seq(from=0, to=50)  # Sensible-ish size grid
W <- 7
likelihood <- dbinom(W, size=N_grid, prob=0.7)
unstd.posterior <- likelihood * uniform_prior(N_grid)
posterior <- unstd.posterior / sum(unstd.posterior)
plot(N_grid, posterior, type="b", xlab="p", ylab="P(N|D)")
mtext(sprintf("%d points", grid_size))

# MAP estimate of N is 9, but from above plot 10 looks to have similar probability
N_grid[which.max(posterior)]
# Yes 10 has equal probability 
posterior[N_grid %in% c(8:11)]

# Can also obtain summary statistics about the posterior by sampling
# By playing with the seed can get MAP = 10 or 11
# But median is generally always 10
set.seed(17)
samples <- sample(N_grid, prob=posterior, size=1e6, replace=TRUE)
chainmode(samples, adj=0.01)  # This gives MAP as 10, why is it different to posterior?
mean(samples)
median(samples)
