
# 3.1 sampling from grid-approximate posterior ----------------------------
grid_approx_globe <- function(grid_size, prior, W, L) {
    # R code 2.3
    p_grid <- seq(from=0, to=1, length.out=grid_size)
    likelihood <- dbinom(W, size=W+L, prob=p_grid)
    unstd.posterior <- likelihood * prior(p_grid)
    unstd.posterior / sum(unstd.posterior)
}

N <- 1e3
p_grid <- seq(0, 1, length.out=N)
posterior <- grid_approx_globe(N, function(x) rep(1, length(x)), 6, 3)
# Sample from posterior of the parameters using the posterior as weights (never used the prob argument before)
# So the posterior is the PDF of the posterior, whereas the sampling here gives the actual values, which is what MCMC
# gives (MCMC just *samples* the posterior, doesn't return it directly).
samples <- sample(p_grid, prob=posterior, size=1e6, replace=TRUE)

library(rethinking)
# Then plotting the density of the samples gives an estimate of the actual posterior
dens(samples)


# 3.2 sampling to summarize -----------------------------------------------

# Obtain posterior P(p<0.5), can do this by integrating the area under the posterior from 0-0.5
# integrating = adding up the probabilities at each grid point
sum(posterior[p_grid < 0.5])
# Can also achieve same answer by looking at proportion of samples below 0.5
mean(samples < 0.5)

# Can get intervals, i.e. P(p > 0.5, p < 0.75)
mean(samples > 0.5 & samples < 0.75)

# Or can get the value for specified quantiles
quantile(samples, 0.8) 
quantile(samples, c(0.1, 0.9))
quantile(samples, c(0.025, 0.975))  # Classic 95% CI.
# Mentions that when have asymmetric CIs (compatibility interval rather than credible), the 
# PIs (Percentile Intervals) aren't that useful
# i.e.

p_grid <- seq(0, 1, length.out=N)
posterior <- grid_approx_globe(N, function(x) rep(1, length(x)), 3, 0)
samples <- sample(p_grid, prob=posterior, size=1e6, replace=TRUE)
dens(samples)
# 50% CI, or PI
PI(samples, prob=0.5)  # function from rethinking package
quantile(samples, c(0.25, 0.75)) # manual

# Highest Posterior Density Interval (HPDI) is narrowest interval that has x% of the data
HPDI(samples, prob=0.5)

# For symmetric distributions there's less of a difference
posterior <- grid_approx_globe(N, function(x) rep(1, length(x)), 6, 3)
samples <- sample(p_grid, prob=posterior, size=1e6, replace=TRUE)
dens(samples)
PI(samples, prob=0.8)
HPDI(samples, prob=0.8)

PI(samples, prob=0.95)
HPDI(samples, prob=0.95)

# MAP (aka mode) can be gotten from the posterior directly
# Going back to the 3n3 binomial
posterior <- grid_approx_globe(N, function(x) rep(1, length(x)), 3, 0)
samples <- sample(p_grid, prob=posterior, size=1e6, replace=TRUE)

p_grid[which.max(posterior)]
# Or can be approximated from the samples
chainmode(samples, adj=0.01)  # function from rethinking

# Can also get posterior mean and median which isn't the same (depending on distribution's symmetry)
mean(samples)
median(samples)

# Different loss functions imply different point estimates
# median minimises loss abs(x - x_hat)
# So if predict p=0.5, then loss here is:
#   - abs(0.5 - p_grid) gets the absolute distance from my chosen point for every value in the grid
#   - multiplying by posterior weights them, and then sum is integrating. But I can't grok how this is the loss at 0.5
# We want to minimise this, find the value which is the 'closest' to everywhere else on average
sum(posterior * abs(0.5 - p_grid))
# So firstly get loss for every point in grid
loss <- sapply(p_grid, function(x) sum(posterior * abs(x - p_grid)))
# And this is minimized at the median
p_grid[which.min(loss)]
median(samples)

# quadratic loss is minimized at the mean
quad_loss <- sapply(p_grid, function(x) sum(posterior * (x - p_grid)**2))
# And this is minimized at the median
p_grid[which.min(quad_loss)]
mean(samples)

# 3.3 sampling to simulate prediction -------------------------------------
# NB sampling from the LIKELIHOOD here, not the posterior!
simplehist(rbinom(1e5, size=9, prob=0.7))
simplehist(rbinom(1e3, size=9, prob=0.7))
simplehist(rbinom(1e2, size=9, prob=0.7))

# Approaching normal, or skewed?
simplehist(rbinom(1e5, size=100, prob=0.7))
simplehist(rbinom(1e5, size=100, prob=0.5))

simplehist(rbinom(1e5, size=100, prob=0.9))

# posterior predictive is data generated from the likelihood weighted by posterior
# (in Stan terms this is done by sampling data from the likelihood using all parameters, i.e.
#  the posterior alphas, betas, sigmas etc...)
# For globe tossing can generate it by drawing from likelihood but weighted by the posterior, rather than using a specific value
# Here using 6W, 3L
posterior <- grid_approx_globe(N, function(x) rep(1, length(x)), 6, 3)
samples <- sample(p_grid, prob=posterior, size=1e6, replace=TRUE)
pp <- rbinom(1e4, size=9, prob=samples)
# Samples is samples from posterior of p
hist(samples)
# Posterior predictive is distribution of actual data, rather than parameters
hist(pp)

# How does passing vector of probs into rbinom work?
# Passing one value means take 1e4 draws from a Binomial distribution of (9, 0.3)
likelihood <- rbinom(1e4, size=9, prob=0.3)
hist(likelihood)
# Passing 2 probs...
# What is this doing?
# It sounds like the probabilities are recycled to be the same length as the number of draws
likelihood <- rbinom(1e4, size=9, prob=c(0.3, 0.9))
hist(likelihood)

# So if we ask for 4 draws but only provide 2 probs, then it will use the first prob, then second, then first, then second
# Which appears to be the case
rbinom(4, size=9, prob=c(0.1, 0.9))
# If have more probs than results it just seems to use the first N
rbinom(1, size=9, prob=c(0.1, 0.9))
# So in the posterior predictive sampling it is doing 1 draw from the likelihood for each sample
# from the posterior, which seems intituitive
# The main idea is that once you've got your posterior, sample the data using all of it, rather than 
# a point estimate

# Practice ----------------------------------------------------------------
p_grid <- seq(0, 1, length.out=1e3)
posterior <- grid_approx_globe(1e3, function(x) rep(1, length(x)), 6, 3)
set.seed(100)
samples <- sample(p_grid, prob=posterior, size=1e4, replace=TRUE)

# 3E1
# 0.04%
# Why don't these give the same result?!
# RNG, rerunning with a finer grid and more samples gives answers closer to 0.08%
# i.e. the posterior method is more reliable
mean(samples < 0.2) * 100 
sum(posterior[p_grid < 0.2]) * 100

# 3E2
# 11.16% or 12.03%
mean(samples > 0.8) * 100 
sum(posterior[p_grid > 0.8]) * 100

# 3E3
# 88.8%
mean(samples > 0.2 & samples < 0.8) * 100 
sum(posterior[p_grid > 0.2 & p_grid < 0.8])*100

# 3E4
# 0.519
quantile(samples, 0.2)

# 3E5
# 0.756
quantile(samples, 0.8)

# 3E6
# 0.509 - 0.774
HPDI(samples, 0.66)

# 3E7
# 0.503 - 0.770
PI(samples, 0.66)

# 3M1
posterior <- grid_approx_globe(1e3, function(x) rep(1, length(x)), 8, 7)

# 3M2
# 0.329 - 0.717
samples_3m2 <- sample(p_grid, prob=posterior, size=1e4, replace=TRUE)
dens(samples_3m2)
HPDI(samples_3m2, 0.90)

# 3M3
# 14.44%
pp_3m3 <- rbinom(1e4, size=15, prob=samples_3m2)
dens(pp_3m3)
mean(pp_3m3 == 8) * 100

# 3M4
# 17.51%
# what is the interpretation of this number?
# given a posterior about the earth obtained from a sample
# what would be the probability of observing a new sample?
pp_3m4 <- rbinom(1e4, size=9, prob=samples_3m2)
mean(pp_3m4 == 6) * 100

# 3M5
step_prior <- function(x) ifelse(x < 0.5, 0, 1)
# 3M2b
posterior_step <- grid_approx_globe(1e3, step_prior, 8, 7)
samples_step <- sample(p_grid, prob=posterior_step, size=1e4, replace=TRUE)
# Far steeper distribution, nothing below 0.5!
plot(posterior_step)
dens(samples_step, col='red', xlim=c(0, 1))
dens(samples_3m2, add=TRUE)
# Far narrower HPDI, lower boundary is at 0.5 now, but upper boundary unchanged
HPDI(samples_3m2, 0.90)
HPDI(samples_step, 0.90)

# 3M3b
# Have higher number of earth touching draws now, although perhaps less of
# a difference than I'd expect from not having any p < 0.5
# There is a bigger difference at 9, 10, 11 / 15, 6-8 are around the middle
# And a higher percentage of seeing 8 draws at 15.25% vs 14.51%
pp_step <- rbinom(1e4, size=15, prob=samples_step)
dens(pp_step, xlim=c(0, 16), col='red')
dens(pp_3m3, add=TRUE)
mean(pp_3m3 == 8) * 100
mean(pp_step == 8) * 100

# 3M4b
# More of a difference in the predictive posterior to get 6/9
# 23.22% vs 17.98%
pp_step <- rbinom(1e4, size=9, prob=samples_step)
mean(pp_step == 6) * 100
mean(pp_3m4 == 6) * 100

# 3M6
# We continually sample data from a system with known p 
# and observe its 99% PI get smaller over time
n_tosses <- function(p_grid, prior, p, width) {
    i <- 0
    while(TRUE) {
        i <- i + 1
        # sample 1 toss
        obs <- rbinom(1, 1, prob=p)
        # calculate likelihood of this observation
        like <- dbinom(obs, 1, p_grid)
        # calculate posterior
        posterior <- like * prior
        posterior <- posterior / sum(posterior)
        
        # Sample from posterior to get PI width
        # Potentially possible to do from the posterior directly
        samples <- sample(p_grid, size=1e4, replace=TRUE, prob=posterior)
        width <- diff(PI(samples, 0.99))
        
        # Update prior to be old posterior
        prior <- posterior
        
        if (width < 0.05) break
    }
    i
}
# Run this simulation a number of times and get final answer
# Around 2,450 times
p_grid <- seq(0, 1, length.out=1e3)
res <- mclapply(1:100, function(i) n_tosses(p_grid, rep(1, length(p_grid)), 0.6, 0.05), mc.cores = 8)
res <- unlist(res)
hist(res)

# 3H1
library(rethinking)
data(homeworkch3)
p_grid <- seq(0, 1, length.out=1e3)
# Prob of any birth being a boy
# so have 111 boys / 200 trials
n_boys <- sum(birth1) + sum(birth2)
n_boys
posterior <- grid_approx_globe(1e3, function(x) rep(1, length(x)), n_boys, 200-n_boys)
plot(p_grid, posterior)
# Want MAP, which gives 0.55
p_grid[which.max(posterior)]
# chainmode gives roughly same answer
samples <- sample(p_grid, prob=posterior, size=1e6, replace=TRUE)
chainmode(samples, adj=0.01)

# 3H2
# 0.526 - 0.572
HPDI(samples, 0.5)
# 0.496 - 0.608
HPDI(samples, 0.89)
# 0.475 - 0.627 
HPDI(samples, 0.97)

# 3H3
# Posterior predictive, see how well simulating future data from the posterior
# matches the observed
# Yep seems reasonable
res <- rbinom(1e4, 200, prob=samples)
dens(res)
abline(v=n_boys, col='red')

# 3H4
# Posterior predictive again, but this time on a new sample of 
# the first cohort
# Looks slightly less likely than the total cohort but still not unreasonable
# Model predicts more boys - is there a difference between the 2 groups?
res <- rbinom(1e4, 100, prob=samples)
dens(res)
abline(v=sum(birth1), col='red')

# Yep birth2 has slightly more
sum(birth1)
sum(birth2)

# So looking at the likelihood of the second cohort (blue) and we can
# see this time the model is _under_estimating
abline(v=sum(birth2), col='blue')

# 3H5
# Posterior predictive again, this time looking at the subgroup of boys following girls
# So my model would predict far fewer boys to families with girls first, then we actually had
n_firstborn_girls <- sum(!birth1)
n_boys_following_girls <- sum(birth2[!birth1])
res <- rbinom(1e4, n_firstborn_girls, prob=samples)
dens(res)
abline(v=n_boys_following_girls, col='red')

