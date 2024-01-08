library(rethinking)
data("Howell1")
d <- Howell1
d2 <- d[ d$age >= 18, ]

# 4.3.3 Grix approximation of the posterior distribution ------------------

mu.list <- seq(from=150, to=160, length.out = 100)
sigma.list <- seq(from=7, to=9, length.out = 100)
post <- expand.grid(mu=mu.list, sigma=sigma.list)
post$LL <- sapply(1:nrow(post), function(i) {
    # Log-Likelihood over all data points (sum log = product raw)
    sum(dnorm(d2$height, post$mu[i], post$sigma[i], log=TRUE))
})
# Add on log priors (= multiplying by raw probs)
post$prod <- post$LL + dnorm(post$mu, 178, 20, log=TRUE) + dunif(post$sigma, 0, 50, log=TRUE)
# Normalise to 0-1, max log prob will be exp(0) = 1. 
# ! Why isn't this instead squashing so the sum = 1?!
post$prob <- exp(post$prod - max(post$prod))

contour_xyz(post$mu, post$sigma, post$prob)
image_xyz(post$mu, post$sigma, post$prob)

# Sampling from posterior, sample mu and sigma proportionally to their posterior prob
sample.rows <- sample(1:nrow(post), size=1e4, replace=TRUE, prob=post$prob)
sample.mu <- post$mu[sample.rows]
sample.sigma <- post$sigma[sample.rows]
plot(sample.mu, sample.sigma, cex=0.5, pch=16, col=col.alpha(rangi2, 0.1))

# why isn't mu's distribution continuous?!
# Is it because grid size is 'only' 100 points?
dens(sample.mu)
dens(sample.sigma)

# Show what happens wen posterior isn't gaussian
d3 <- sample(d2$height, size=20)
mu.list <- seq(from=150, to=160, length.out = 200)
sigma.list <- seq(from=7, to=9, length.out = 200)
post2 <- expand.grid(mu=mu.list, sigma=sigma.list)
post2$LL <- sapply(1:nrow(post2), function(i) {
    sum(dnorm(d3, post2$mu[i], post2$sigma[i], log=TRUE))
})
post2$prod <- post2$LL + dnorm(post2$mu, 178, 20, log=TRUE) + dunif(post2$sigma, 0, 50, log=TRUE)
post2$prob <- exp(post2$prod - max(post2$prod))
sample.rows <- sample(1:nrow(post2), size=1e4, replace=TRUE, prob=post2$prob)
sample.mu2 <- post2$mu[sample.rows]
sample.sigma2 <- post$sigma[sample.rows]
plot(sample.mu2, sample.sigma2, cex=0.5, pch=16, col=col.alpha(rangi2, 0.1))
dens(sample.sigma2)  # Poor normal fit!

# 4.3.5 Finding the posterior distribution with quap ----------------------

flist <- alist(
    height ~ dnorm(mu, sigma),
    mu ~ dnorm(170, 20),
    sigma ~ dunif(0, 50)
)
m4.1 <- quap(flist, data=d2)
m4.1
precis(m4.1)

# Can change the optimisation initial conditions
# Although doesn't really change much here
start <- list(
    mu = mean(d2$height),
    sigma=sd(d2$height)
)
m4.1 <- quap(flist, data=d2, start=start)
m4.1

# More informative prior on mu, much narrower
# Look how the SD of mu is = to the prior sd
# And Sigma's mean has become much larger to accomodate
# the reduced variance of mu
m4.2 <- quap(alist(
    height ~ dnorm(mu, sigma),
    mu ~ dnorm(170, 0.1),
    sigma ~ dunif(0, 50)
), data=d2)
precis(m4.2)

# Under the hood quap uses a multivariate gaussian, as it assumes
# all posteriors are gaussian. can see cov matrix
vcov(m4.1)
# And can extract the variances and correlation matrix
diag(vcov(m4.1))
cov2cor(vcov(m4.1))

# Sampling from multi-dimensional posterior
# Get n-column data frame
post <- extract.samples(m4.1, n=1e4)
dim(post)
head(post)
# And can see samples distribution, which is very close to posterior MAP
precis(post)

# 4.4 Linear prediction ---------------------------------------------------
# Adding in predictor covariate
# Clear relationship between height and weight
plot(d2$height, d2$weight)

# Prior predictive
# I.e. sample from prior, use samples in model form to generate samples of the outcome
# Useful sanity check of assumptions! 
# Here, the heights can go negative with low weight,
# and with the highest weights can go far above the highest ever recorded
# As well as having potential negative relationships - counter to our domain knowledge
set.seed(2971)
N <- 100
a <- rnorm(N, 178, 20)
b <- rnorm(N, 0, 10)
plot(NULL, xlim=range(d2$weight), ylim=c(-100, 400), 
     xlab="weight", ylab="height")
abline(h=0, lty=2)
abline(h=272, lty=1, lwd=0.5)
mtext("b ~ dnorm(0, 10)")
xbar <- mean(d2$weight)
for (i in 1:N) {
    curve(a[i] + b[i] * (x-xbar),
          from=min(d2$weight), to=max(d2$weight),
          add=TRUE,
          col=col.alpha("black", 0.2)
    )
}

# Try beta prior LogNormal(0, 1)
# Enforces positive relationship
# NB: why does Stan recommend half-normal
# or cauchy or exponential rather than log-normal
# for scale terms?
# Exponential is similar to log-normal, but has
# thinner tails, which could explain it
# This is a really useful reminder of the importance
# of checking priors rather than just trying new ones
# in Stan until the MCMC converges
b <- rlnorm(1e4, 0, 1)
dens(b, xlim=c(0, 5), adj=0.1)
d <- rexp(1e4, 1)
dens(d, xlim=c(0, 5), adj=0.1, add=TRUE, col="red")

# Much better!
set.seed(2971)
N <- 100
a <- rnorm(N, 178, 20)
b <- rlnorm(N, 0, 1)
plot(NULL, xlim=range(d2$weight), ylim=c(-100, 400), 
     xlab="weight", ylab="height")
abline(h=0, lty=2)
abline(h=272, lty=1, lwd=0.5)
mtext("b ~ dnorm(0, 10)")
xbar <- mean(d2$weight)
for (i in 1:N) {
    curve(a[i] + b[i] * (x-xbar),
          from=min(d2$weight), to=max(d2$weight),
          add=TRUE,
          col=col.alpha("black", 0.2)
    )
}

# 4.4.2 Finding the posterior distribution --------------------------------

xbar <- mean(d2$weight)
m4.3 <- quap(
    alist(
        height ~ dnorm(mu, sigma),
        mu <- a + b*(weight-xbar),
        a ~ dnorm(178, 20),
        b ~ dlnorm(0, 1),
        sigma ~ dunif(0, 50)
    ),
    data=d2
)
m4.3
precis(m4.3)

# Or can parameterize with logb instead
# Both of these ways get the benefit of forcing positive relationship,
# while still getting an interpretable coefficient
# I guess the trade-off is the first method gives you samples directly
# on the original scale, while with this one you have more flexibility
# to use whatever distribution you want, but you have to exponentiate
# the samples to understand them
m4.3b <- quap(
    alist(
        height ~ dnorm(mu, sigma),
        mu <- a + exp(log_b)*(weight-xbar),
        a ~ dnorm(178, 20),
        log_b ~ dnorm(0, 1),
        sigma ~ dunif(0, 50)
    ),
    data=d2
)
m4.3b
# Get the exact same results for alpha & sigma,
# but log_b is different as expected
precis(m4.3b)
# But can convert back to regular scale
exp(-0.1)

# 4.4.3 Interpreting the posterior distribution ---------------------------
round(vcov(m4.3), 3)
pairs(m4.3)

# Plotting posterior (NB: NOT posterior predictive, as not accounting for
# uncertainty in observations, sigma)
# Also not full posterior as just using MAP/mean parameter estimates
plot(height ~ weight, data=d2, col=rangi2)
post <- extract.samples(m4.3)
a_map <- mean(post$a)
b_map <- mean(post$b)
curve(a_map + b_map*(x-xbar), add=TRUE)

# Showing uncertainty from parameter estimates rather than just using
# MAP estimates. Starting off by subsetting to 10 data points
plot_regression_subset <- function(N) {
    dN <- d2[1:N, ]
    mN <- quap(
        alist(
            height ~ dnorm(mu, sigma),
            mu <- a + b * (weight - mean(weight)),
            a ~ dnorm(178, 20),
            b ~ dlnorm(0, 1),
            sigma ~ dunif(0, 50)
        ),
        data=dN
    )
    # Sample 20 sets of parameters from posterior and plot
    post <- extract.samples(mN, n=20)
    plot(dN$weight, dN$height,
         xlim=range(d2$weight), ylim=range(d2$height),
         col=rangi2, xlab="weight", ylab="height")
    mtext(concat("N = ", N))
    # Plot regression lines
    for (i in 1:nrow(post)) {
        curve(
            post$a[i] + post$b[i] * (x-mean(dN$weight)),
            col=col.alpha("black", 0.3), add=TRUE
        )
    }
}
# As the amount of data points increases, the uncertainty in the parameters decreases
# Could also see this in the marginal densities of a and b, but those don't take the covariance
# into account, hence plotting regression lines
plot_regression_subset(10)
plot_regression_subset(50)
plot_regression_subset(100)
plot_regression_subset(352)

# Plotting uncertainty through intervals and contours
post <- extract.samples(m4.3)
# This posterior contains the correlation structure from a and b
# TODO: When to look at posterior of linear predictor, and posterior predictive 
# (the same thing but including sigma?)
# I guess look at linear predictor first as a sanity check for model structure
# and THEN look at posterior predictive
mu_at_50 <- post$a + post$b * (50 - xbar)
dens(mu_at_50, col=rangi2, lwd=2, xlab="mu | weight at 50")

# out of my interest compare to posterior predictive
# have to resort to ggplot as my base R plotting skills are lacking
# But can see clearly that the posterior mu has far less variance, even
# when we're conditioning on weight at 50kg
# The difference is that the posterior is the distribution over
# MEAN HEIGHTS conditioned on 50
# Whereas posterior predictive is distribution over ALL HEIGHTS
pp <- rnorm(1e4, mean=mu_at_50, sd=post$sigma)
library(tidyverse)
tibble(post=mu_at_50, post_pred=pp) |>
    mutate(i=row_number()) |>
    pivot_longer(-i) |>
    ggplot(aes(x=value, colour=name)) +
        geom_density()

# Back to the book...
PI(mu_at_50, prob=0.89)

# Want to repeat this for every datapoint, not just for weight = 50kg
# Could do this in Stan in generated quantities
mu <- link(m4.3)
dim(mu)

# Although rather than generating this for _every point in the original dataset_, 
# we want to calculate this for a known set of data points (again can be done in Stan
# through generated quantities)
# This shows the uncertainty in mu gets larger at the extremes (makes sense, fewer data points)
# similar to the effect shown from plotting multiple regression lines
weights.seq <- seq(from=25, to=70, by=1)
mu <- link(m4.3, data=data.frame(weight=weights.seq))
dim(mu)
plot(height ~ weight, d2, type="n")  # Not actually plotting these! Just to get the canvas setup
for (i in 1:100) {
    points(weights.seq, mu[i, ], pch=16, col=col.alpha(rangi2, 0.1))
}

# Summarize the distribution of mu, calculating statistics for each weight value
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob=0.89)
str(mu.mean)
str(mu.PI)

# Tie all together in one plot
plot(height ~ weight, d2, col=col.alpha(rangi2, 0.5))  
# Plot MAP line, which is mean mu for each weight (in Gaussian approx MAP=mean=median)
lines(weights.seq, mu.mean)
# Shaded 89% PI region. Function from rethinking
# I should always show regions when I do my own plots!
shade(mu.PI, weights.seq)

# rethinking provides sim function to generate posterior predictive
sim.height <- sim(m4.3, data=list(weight=weights.seq))
# Has simulated heights, NOT distributions of plausible average heights mu
str(sim.height)
height.PI <- apply(sim.height, 2, PI, prob=0.89)

# Plot posterior of mu alongside posterior predictive
plot(height ~ weight, d2, col=col.alpha(rangi2, 0.5))  
# MAP line of mu (is basically the same for posterior and posterior predictive)
lines(weights.seq, mu.mean)
lines(weights.seq, apply(sim.height, 2, mean), col='red')

# 95% PIs, note how much extra variance is in the posterior predictive!
# Which makes sense, should contain 95% of all values?
shade(mu.PI, weights.seq)
shade(height.PI, weights.seq)

# 4.5 Curves from lines ---------------------------------------------------
data(Howell1)
d <- Howell1
d$weight_centered <- d$weight - mean(d$weight)
d$weight_standardized <- d$weight_centered / sd(d$weight)
plot(height ~ weight, d)
# Same relationship when centering weight
plot(height ~ weight_centered, d)
# and when standardized
plot(height ~ weight_standardized, d)

# Also create variable for squared term
d$weight_standardized_2 <- d$weight_standardized**2
# And this inverts the relationship, as weight with a Z of -1 and 1 are now equivalent
plot(height ~ weight_standardized_2, d)
# NB: STANDARDIZING PREDICTOR, NOT OUTCOME!
# If fit model in such a way that only allows for positive outcomes, then don't always need to
# transform outcome, depends on the problem.

# Fit polynomial model
m4.5 <- quap(
    alist(
        height ~ dnorm(mu, sigma),
        mu <- a + b1*weight_standardized + b2*weight_standardized_2,
        a ~ dnorm(178, 20),
        b1 ~ dlnorm(0, 1),
        b2 ~ dnorm(0, 1),
        sigma ~ dunif(0, 50)
    ),
    data=d
)
# b1 and b2 lose interpretability!
# intercept is no longer the mean height in dataset
precis(m4.5)

# Will plot posterior of mu to understand
weights.seq <- seq(from=-2.2, to=2, length.out=30)
pred_dat <- list(weight_standardized=weights.seq, weight_standardized_2=weights.seq**2)
mu <- link(m4.5, data=pred_dat)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob=0.89)

# Sample from posterior predictive too
sim.height <- sim(m4.5, data=pred_dat)
height.PI <- apply(sim.height, 2, PI, prob=0.89)

# That fits quite well!
plot(height ~ weight_standardized, d, col=col.alpha(rangi2, 0.5))
lines(weights.seq, mu.mean)
shade(mu.PI, weights.seq)
shade(height.PI, weights.seq)

# What about a log transform?
m4.5_log <- quap(
    alist(
        height ~ dnorm(mu, sigma),
        mu <- a + b1*log(weight),
        a ~ dnorm(178, 20),
        b1 ~ dlnorm(0, 1),
        b2 ~ dnorm(0, 1),
        sigma ~ dunif(0, 50)
    ),
    data=d
)
weights <- seq(min(d$weight), max(d$weight))
plot(height ~ log(weight), d)
mu <- link(m4.5_log, data=data.frame(weight=weights))
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob=0.89)
sim.height <- sim(m4.5_log, data=data.frame(weight=weights))
height.PI <- apply(sim.height, 2, PI, prob=0.89)

# Potentially slightly worse fit than polynomial due to the lower tail behaviour
# But above that extremity it fits quite well
plot(height ~ log(weight), d, col=col.alpha(rangi2, 0.5))
lines(log(weights), mu.mean)
shade(mu.PI, log(weights))
shade(height.PI, log(weights))
