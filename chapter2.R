# Contains notes, workings, and end-of-chapter question answers

# 2.3 grid approximation --------------------------------------------------

grid_approx_globe <- function(grid_size, prior, W, L) {
    # R code 2.3
    p_grid <- seq(from=0, to=1, length.out=grid_size)
    likelihood <- dbinom(W, size=W+L, prob=p_grid)
    unstd.posterior <- likelihood * prior(p_grid)
    posterior <- unstd.posterior / sum(unstd.posterior)
    plot(p_grid, posterior, type="b", xlab="p", ylab="P(p|D)")
    mtext(sprintf("%d points", grid_size))
}

uniform_prior <- function(p) rep(1, length(p))

grid_approx_globe(20, uniform_prior, 6, 3)
grid_approx_globe(5, uniform_prior, 6, 3)

step_prior <- function(p) ifelse(p < 0.5, 0, 1)
sharp_prior <- function(p) exp(-5*abs(p-0.5))
grid_approx_globe(50, step_prior, 6, 3)
grid_approx_globe(50, sharp_prior, 6, 3)

# 2.6 quadratic approximation ---------------------------------------------

library(rethinking)
globe.qa <- quap(
    alist(
        W ~ dbinom(W+L, p), 
        p ~ dunif(0, 1)
    ),
    data=list(W=6, L=3)
)
precis(globe.qa)

# Comparison with analytical
W <- 6
L <- 3
curve(dbeta(x, W+1, L+1), from=0, to=1)       # analytical, from conjugate beta-binomial
                                              # beta(x+a, n-x+Beta) with prior Beta(1,1),
                                              # where Beta(1,1) is flat prior
curve(dnorm(x, 0.67, 0.16), lty=2, add=TRUE)  # From quadratic approx

# 2.8 MCMC ----------------------------------------------------------------

mcmc_globe <- function(n_samples, p_0) {
    p <- rep(NA, n_samples)
    p[1] <- p_0
    W <- 6
    L <- 3
    for (i in 2:n_samples) {
        p_new <- rnorm(1, p[i-1], 0.1)  # How to choose variance?
        if (p_new < 0) p_new <- abs(p_new)  # Truncated normal within 0-1, I imagine because
                                            # binomial p has to be [0, 1]
        if (p_new > 1) p_new <- 2 - p_new
        q0 <- dbinom(W, W+L, p[i-1])        # likelihood of previous and new samples
        q1 <- dbinom(W, W+L, p_new)
        p[i] <- ifelse(runif(1) < q1/q0, p_new, p[i-1])  # Seem to have bias towards new sample?
                                                         # if q1 > q0 then take new value (ratio > 1)
                                                         # if q1 < q0 then take new value at rate of q1/q0 
                                                         # i.e. if q1 = 0.5 and q0 = 0.3, then take q1
                                                         # if q1 = 0.3 and q0 = 0.5 then have 60% chance to take q1
    }
    dens(p, xlim=c(0, 1))
    curve(dbeta(x, W+1, L+1), lty=2, add=TRUE)  # analytical
}

mcmc_globe(10, 0.5)
mcmc_globe(100, 0.5)
mcmc_globe(1000, 0.5)
mcmc_globe(10000, 0.5)

# Practice ----------------------------------------------------------------

# 2M1.
qa <- function(Win, Lin) {
    precis(quap(
        alist(
            W ~ dbinom(W+L, p), 
            p ~ dunif(0, 1)
        ),
        data=list(W=Win, L=Lin)
    ))
}
grid_approx_globe(100, uniform_prior, 3, 0)
grid_approx_globe(100, uniform_prior, 3, 1)
grid_approx_globe(100, uniform_prior, 5, 2)

# 2M2.
new_prior <- function(p) ifelse(p < 0.5, 0, 0.8)
grid_approx_globe(100, new_prior, 3, 0)
grid_approx_globe(100, new_prior, 3, 1)
grid_approx_globe(100, new_prior, 5, 2)

# 2M3.
# have just 2 parameters, 0.7 and 1
# So can use the logic from grid approximation
p_lands <- c(0.3, 1)
likelihoods <- dbinom(1, size=1, prob=p_lands)
# The uniform prior doesn't change anything here
unstd.posterior <- likelihoods * uniform_prior(p_lands)
# Standardize to sum to 1 and we get 0.23
unstd.posterior / sum(unstd.posterior)

# 2M4.
# card (conjecture) | ways to produce data (B) | normalised 
# ~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~ |~~~~~~~~~~~~
#     BB            |       2                  |   2/3
#     BW            |       1                  |   1/3
#     WW            |       0                  |    0

# 2M5
# card (conjecture) | ways to produce data (B) | normalised, summing first 2 gives 4/5
# ~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~ |~~~~~~~~~~~~
#     BB            |       2                  |   2/5
#     BB            |       2                  |   2/5
#     BW            |       1                  |   1/5
#     WW            |       0                  |    0

# 2M6
# But now with uneven priors
# card (conjecture) | prior | ways to produce data (B) | product  | normalised 
# ~~~~~~~~~~~~~~~~~~| ~~~~~~| ~~~~~~~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~|~~~~~~~~~~~
#     BB            |   1    |      2                  |     2    |    1/2
#     BW            |   2   |       1                  |     2    |    1/2
#     WW            |   3   |       0                  |     0    |     0

# 2M7
# card (conjecture) | ways to produce data (BW) | product  | normalised 
# ~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~ |~~~~~~~~~~~| ~~~~~~~~~~
#     BB            |       2 x 3              |  6        |  3/4
#     BW            |       1 x 2              |  2        |  1/4
#     WW            |       0 x 1              |  0        |   0

# 2H1
# binomial with 1 trial and candidate ps are 0.1 and 0.2
# assume even prior ('equally common in the wild')
p_species <- c(0.1, 0.2)
likelihoods <- dbinom(1, size=1, prob=p_species)
# The uniform prior doesn't change anything here, just showing for completeness
unstd.posterior <- likelihoods * uniform_prior(p_species)
# Posterior is 1/3 species A and 1/3 species B
posterior <- unstd.posterior / sum(unstd.posterior)
posterior
# Probability of twins is probability twin given A OR (+) probability twin given B
# 1/6
(1/3) * 0.1 + (2/3) * 0.2

# 2H2
# Just posterior from above, i.e. 1/3 
# Have I missed something? Need to do this step in order to answer previous question

# 2H3
# Use 1/3 and 2/3 as new priors and update posterior
# Get 0.36 as prob that have A
likelihood_new <- dbinom(0, size=1, prob=p_species)
posterior_unstd <- likelihood_new * posterior
posterior_new <- posterior_unstd / sum(posterior_unstd)
posterior_new

# 2H4
# P(A|+ve test) = P(+ve test | A) * P(A) / P(+ve test)
# P(+ve test) = P(+ve test |A) * P(A) + P(+ve test | B) * P(B)
p_pve_test <- 0.8 * 0.5 + 0.65 * 0.5
p_a_pve_test <- 0.8 * 0.5 / p_pve_test
# 0.551
p_a_pve_test

# part 2, rather than use 0.5 as priors, use 0.36 and 0.64
p_pve_test <- 0.8 * 0.36 + 0.65 * 0.64
p_a_pve_test <- 0.8 * 0.36 / p_pve_test
# 0.409, so the probability has lowered
p_a_pve_test
