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

# Show what happens when posterior isn't gaussian
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

# 4.5.1 Polynomial regression ---------------------------------------------------
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

# 4.5.2 Splines -----------------------------------------------------------
data(cherry_blossoms)
d <- cherry_blossoms
precis(d)

# Fit model with 15 knots
d2 <- d[complete.cases(d$doy), ]
num_knots <- 15
knot_list <- quantile(d2$year, probs=seq(0, 1, length.out=num_knots))
knot_list

# Use
library(splines)
B <- bs(d2$year, knots=knot_list[-c(1, num_knots)],
        degree=3, intercept=TRUE)
# How do we get 17 cols? 15 knots, - start & end = 13 knots, add intercept=??
head(B)

# The first column is intercept, so get 16 cols without
# Why use intercept if fit it in regression model?
B2 <- bs(d2$year, knots=knot_list[-c(1, num_knots)],
        degree=3, intercept=FALSE)
head(B2)

plot(NULL, xlim=range(d2$year), ylim=c(0, 1), xlab="year", ylab="basis")
for (i in 1:ncol(B)) {
    lines(d2$year, B[, i])
}

# Without intercept lose that starting basis
plot(NULL, xlim=range(d2$year), ylim=c(0, 1), xlab="year", ylab="basis")
for (i in 1:ncol(B2)) {
    lines(d2$year, B2[, i])
}

# What if kept the first knot, is that the same as an intercept?
B3 <- bs(d2$year, knots=knot_list[-c(num_knots)],
        degree=3, intercept=FALSE)
# The basis functions look the same!
# And plot looks the same!
head(B3)
head(B)
plot(NULL, xlim=range(d2$year), ylim=c(0, 1), xlab="year", ylab="basis")
for (i in 1:ncol(B3)) {
    lines(d2$year, B3[, i])
}


# Finally, what happens if ask for intercept and don't remove anything from knot list?
B4 <- bs(d2$year, knots=knot_list[-c(num_knots)],
        degree=3, intercept=TRUE)
# Get an extra column
dim(B4)
# Which is all zeros... an intercept should be 1?
head(B4)
table(B4[, 1])

# Now fit the model
m4.7 <- quap(
    alist(
        D ~ dnorm(mu, sigma),
        mu <- a + B %*% w,
        a ~ dnorm(100, 10),
        w ~ dnorm(0, 10),
        sigma ~ dexp(1)
    ),
    data=list(D=d2$doy, B=B),
    start=list(w=rep(0, ncol(B)))
)

# Plot MAP posterior basis effects (modelled weights * basis functions)
post <- extract.samples(m4.7)
w <- apply(post$w, 2, mean)
plot(NULL, xlim=range(d2$year), ylim=c(-6, 16), xlab="year", ylab="basis*weight")
for (i in 1:ncol(B)) {
    lines(d2$year, w[i] * B[, i])
}

# Add 97% PI for mu (i.e. posterior, not posterior predictive)
mu <- link(m4.7)
mu_PI <- apply(mu, 2, PI, 0.97)
plot(d2$year, d2$doy, col=col.alpha(rangi2, 0.3), pch=16)
shade(mu_PI, d2$year, col=col.alpha("black", 0.5))

# Problems ----------------------------------------------------------------
# 4E1: yi ~ N(mu, sigma) is likelihood
# 4E2: There are 2 parameters (mu and sigma)
# 4E3: p(u, sigma | y) = p(sigma, u | y) * p(sigma, u) / p(y)
# 4E4: ui = alpha + beta*xi is the linear model
# 4E5: There are 3 parameters: alpha, beta, sigma

# 4M1
# Simulated observed y from the prior
# yi ~ N(mu, sigma)
# u ~ N(0, 10)
# sigma ~ Exp(1)
mu_sim <- rnorm(1e4, 0, 10)
sigma_sim <- rexp(1e4, 1)
y_sim <- rnorm(1e4, mu_sim, sigma_sim)
dens(y_sim)

# 4M2
# quap implementation
m_4m2 <- quap(
    alist(
        y ~ dnorm(mu, sigma),
        mu ~ dnorm(0, 10),
        sigma ~ dexp(1)
    ),
    data = list(y=rnorm(1e3, 2, 5))
)
# Enough data so that it found the true parameters!
precis(m_4m2)
# If use fewer data points get less accurate estimates with greater variance
quap(
    alist(
        y ~ dnorm(mu, sigma),
        mu ~ dnorm(0, 10),
        sigma ~ dexp(1)
    ),
    data = list(y=rnorm(1e2, 2, 5))
) |> precis()

# 4M3
# yi ~ N(mu, sigma)
# mui = a + b*xi
# a ~ N(0, 10)
# b ~ Unif(0, 1)
# sigma ~ Exp(1)

# 4M4
# height ~ N(mu, sigma)
# mui = alpha + beta * (year - year_0)
# alpha ~ N(130, 40)  # should cover most heights, not sure if students means primary/secondary/FE
# beta ~ N(0, 5)  # No idea how much people grow in 1 year, but 95% +/- 10cm feels about right
# sigma ~ Exp(1)  # Standard scale prior
# Let's simulate from prior
a <- rnorm(1e4, 130, 40)
b <- rnorm(1e4, 0, 5)
years <- 1:3
year_0 <- years[1]

# High prior has resulted in zero negative heights
# But there are some people approaching the tallest known height ever!
plot(NULL, xlim=c(0, 3), ylim=c(-100, 400), 
     xlab="Years ", ylab="height")
abline(h=0, lty=2)
abline(h=272, lty=1, lwd=0.5)
mtext("")
for (i in 1:N) {
    curve(a[i] + b[i] * (x-year_0),
          from=1, to=3,
          add=TRUE,
          col=col.alpha("black", 0.2)
    )
}

# Will tighten intercept a bit
a <- rnorm(1e4, 130, 30)
b <- rnorm(1e4, 0, 5)
years <- 1:3
year_0 <- years[1]

# High prior has resulted in zero negative heights
# But there are some people approaching the tallest known height ever!
plot(NULL, xlim=c(0, 3), ylim=c(-100, 400), 
     xlab="Years ", ylab="height")
abline(h=0, lty=2)
abline(h=272, lty=1, lwd=0.5)
mtext("")
for (i in 1:N) {
    curve(a[i] + b[i] * (x-year_0),
          from=1, to=3,
          add=TRUE,
          col=col.alpha("black", 0.2)
    )
}

# 4M5
# Change to log normal to enforce positive relationship between year and growth
# Should also put a positive constraint on the likelihood or alpha?
# Have quite a wide alpha prior that could get negative heights
# beta ~ lognormal(0, 1) 
# Now only have positive growth!
b <- rlnorm(1e4, 0, 1)
plot(NULL, xlim=c(0, 3), ylim=c(-100, 400), 
     xlab="Years ", ylab="height")
abline(h=0, lty=2)
abline(h=272, lty=1, lwd=0.5)
mtext("")
for (i in 1:N) {
    curve(a[i] + b[i] * (x-year_0),
          from=1, to=3,
          add=TRUE,
          col=col.alpha("black", 0.2)
    )
}

# One thing I find hard to get my head around is interpreting log-normal priors
# So this doesn't have much weight above ten
# So what does '1' mean exactly? The variance of log-beta? Hardly interpretable!
# Often I just want a coefficient that must be positive and I know roughly what values
# it should take and don't want to have to think about log-scales
# But I guess if you standardize your data then having a lot of weight around 0-1 is fine
dens(rlnorm(1e4, 0, 1))

# 4M6
# Could use sigma ~ Unif(0, 64), but that is a very vague prior

# 4M7
# Original model
d2 <- Howell1[ Howell1$age >= 18, ]
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

round(vcov(m4.3), 3)
pairs(m4.3)

# Without centering first
# First time it threw a convergence warning, first time I've seen that!
m4.3_noncentered <- quap(
    alist(
        height ~ dnorm(mu, sigma),
        mu <- a + b*weight,
        a ~ dnorm(178, 20),
        b ~ dlnorm(0, 1),
        sigma ~ dunif(0, 50)
    ),
    data=d2
)
# Alpha is very different now, because it has a different interpretation
# Now alpha is the height when someone has 0 weight, i.e. doesn't exist
# Previously it was the height of an average weighted person
# The other parameters are essentially the same
# beta is the same as it's still the effect of a 1kg weight increase on height
# and sigma still is the model uncertainty
precis(m4.3)
precis(m4.3_noncentered)

# But look how much higher the correlations are now!
round(vcov(m4.3), 3)
round(vcov(m4.3_noncentered), 3)

# Can see this graphically especially between a and b
# As a gets larger, b gets smaller
# As the weight of a non-existent person increases, 
# the height that each additional kg of weight adds is less
# But why isn't this true for the centered version?
pairs(m4.3)
pairs(m4.3_noncentered)

pp_orig <- sim(m4.3)
pp_noncentered <- sim(m4.3_noncentered)

# But doesn't effect posterior predictive
dens(pp_orig)
dens(pp_noncentered, col='red', add=TRUE)

# Nor the posterior mu
dens(link(m4.3))
dens(link(m4.3_noncentered), col='red', add=TRUE)

# So does it matter if we have correlated parameters if it doesn't affect
# the predictions? Or is it more that having correlated parameters makes sampling harder?

# M4.8
# Splines with increasing knots
# Width of the prior on weights
d2 <- cherry_blossoms[complete.cases(cherry_blossoms$doy), ]
fit_basis_funcs <- function(x, n_knots) {
    knot_list <- quantile(x, probs=seq(0, 1, length.out=n_knots))
    B <- bs(x, knots=knot_list[-c(1, n_knots)],
            degree=3, intercept=TRUE)
}

knots_15 <- fit_basis_funcs(d2$year, 15)
knots_30 <- fit_basis_funcs(d2$year, 30)

plot(NULL, xlim=range(d2$year), ylim=c(0, 1), xlab="year", ylab="basis")
for (i in 1:ncol(knots_15)) {
    lines(d2$year, knots_15[, i])
}

# And for 30 knots, as expected more basis functions
plot(NULL, xlim=range(d2$year), ylim=c(0, 1), xlab="year", ylab="basis")
for (i in 1:ncol(knots_30)) {
    lines(d2$year, knots_30[, i])
}

# Imagine it will make predictions more flexible
# Now fit the model
mod_15 <- quap(
    alist(
        D ~ dnorm(mu, sigma),
        mu <- a + B %*% w,
        a ~ dnorm(100, 10),
        w ~ dnorm(0, 10),
        sigma ~ dexp(1)
    ),
    data=list(D=d2$doy, B=knots_15),
    start=list(w=rep(0, ncol(knots_15)))
)
mod_30 <- quap(
    alist(
        D ~ dnorm(mu, sigma),
        mu <- a + B %*% w,
        a ~ dnorm(100, 10),
        w ~ dnorm(0, 10),
        sigma ~ dexp(1)
    ),
    data=list(D=d2$doy, B=knots_30),
    start=list(w=rep(0, ncol(knots_30)))
)

# Plot MAP posterior basis effects (modelled weights * basis functions)
post <- extract.samples(mod_15)
w <- apply(post$w, 2, mean)
plot(NULL, xlim=range(d2$year), ylim=c(-6, 16), xlab="year", ylab="basis*weight")
for (i in 1:ncol(knots_15)) {
    lines(d2$year, w[i] * knots_15[, i])
}

# And for 30 knots can really see each basis function, but are only 3 still effecting
# the spline at any point? Definitely looks like 4 at some points below, but maybe that's just
# when intercept is in play?
post <- extract.samples(mod_30)
w <- apply(post$w, 2, mean)
plot(NULL, xlim=range(d2$year), ylim=c(-6, 16), xlab="year", ylab="basis*weight")
for (i in 1:ncol(knots_30)) {
    lines(d2$year, w[i] * knots_30[, i])
}

# And showing PI of mu
mu <- link(mod_15)
mu_PI <- apply(mu, 2, PI, 0.97)
plot(d2$year, d2$doy, col=col.alpha(rangi2, 0.3), pch=16)
shade(mu_PI, d2$year, col=col.alpha("black", 0.5))

# Yep much wiggler
mu <- link(mod_30)
mu_PI <- apply(mu, 2, PI, 0.97)
plot(d2$year, d2$doy, col=col.alpha(rangi2, 0.3), pch=16)
shade(mu_PI, d2$year, col=col.alpha("black", 0.5))

# Widening variance on basis weights 
# From N(0, 10) to N(0, 20)
# Will stick with 15 knots
# My guess is this will make the spline less smooth
mod_15_wide <- quap(
    alist(
        D ~ dnorm(mu, sigma),
        mu <- a + B %*% w,
        a ~ dnorm(100, 10),
        w ~ dnorm(0, 20),
        sigma ~ dexp(1)
    ),
    data=list(D=d2$doy, B=knots_15),
    start=list(w=rep(0, ncol(knots_15)))
)

# Plot MAP posterior basis effects (modelled weights * basis functions)
# Means are the same as this posterior basis effects plot shows
post <- extract.samples(mod_15_wide)
w <- apply(post$w, 2, mean)
plot(NULL, xlim=range(d2$year), ylim=c(-6, 16), xlab="year", ylab="basis*weight")
for (i in 1:ncol(knots_15)) {
    lines(d2$year, w[i] * knots_15[, i])
}

# And showing PI of mu
# Barely made any difference!
mu <- link(mod_15_wide)
mu_PI <- apply(mu, 2, PI, 0.97)
plot(d2$year, d2$doy, col=col.alpha(rangi2, 0.3), pch=16)
shade(mu_PI, d2$year, col=col.alpha("black", 0.5))

# What about narrowing it to N(0, 2)?
# This made it smoother, contrary to my expectation
# Is it that having a narrow prior means the data has less chance
# to actually find a sensible coefficient so the effect stays closer to zero
mod_15_narrow <- quap(
    alist(
        D ~ dnorm(mu, sigma),
        mu <- a + B %*% w,
        a ~ dnorm(100, 10),
        w ~ dnorm(0, 2),
        sigma ~ dexp(1)
    ),
    data=list(D=d2$doy, B=knots_15),
    start=list(w=rep(0, ncol(knots_15)))
)

# And showing PI of mu
# Barely made any difference!
mu <- link(mod_15_narrow)
mu_PI <- apply(mu, 2, PI, 0.97)
plot(d2$year, d2$doy, col=col.alpha(rangi2, 0.3), pch=16)
shade(mu_PI, d2$year, col=col.alpha("black", 0.5))

# 4H1
# Use m4.3, model of height on weight with narrowish priors
# Should this use posterior or posterior predictive?
mu <- link(m4.3, data=data.frame(weight=c(46.95, 43.72, 64.78, 32.59, 54.63)))
apply(mu, 2, mean)
apply(mu, 2, PI, 0.89)

# 4H2
d2 <- Howell1[ Howell1$age < 18, ]
nrow(d2)

# Fit same model as m4.3 but with different priors
# Height is around 60-160! so just put a uniform prior in this range
hist(d2$height)
xbar <- mean(d2$weight)
h2 <- quap(
    alist(
        height ~ dnorm(mu, sigma),
        mu <- a + b*(weight-xbar),
        a ~ dunif(60, 160),
        b ~ dlnorm(0, 1),
        sigma ~ dunif(0, 50)
    ),
    data=d2
)

# a) The expected height for an average weight person is 108cm
# For every 1 increase in kg, there is an associated height increase of 2.72cm
# So 10kg is 27.2cm
precis(h2)

# b) Plot
plot(height ~ weight, data=d2, col=rangi2)
post <- extract.samples(h2)
a_map <- mean(post$a)
b_map <- mean(post$b)
curve(a_map + b_map*(x-xbar), add=TRUE)

new_weights <- seq(min(d2$weight), max(d2$weight))
mu <- link(h2, data=data.frame(weight=new_weights))
mu_PI <- apply(mu, 2, PI, 0.89)
pp <- sim(h2, data=data.frame(weight=new_weights))
pp_PI <- apply(pp, 2, PI, 0.89)

plot(d2$weight, d2$height, col=col.alpha(rangi2, 0.3), pch=16)
shade(mu_PI, new_weights, col=col.alpha("black", 0.3))
shade(pp_PI, new_weights, col=col.alpha("black", 0.3))

# c) 
# The most concerning fact is the linear fit for a demonstrably non-linear relationship
# Could try a log-transform of weight, or piecewise linear with a break at ~30kg
# Or even adjust for age, as there's a clear confounder here

# 4H3
# Using log transforms of the full dataset
d <- Howell1
# Very large distribution!
# Will use wide prior
hist(d$height)
xbar_log <- mean(log(d$weight))

h3 <- quap(
    alist(
        height ~ dnorm(mu, sigma),
        mu <- a + b*(log(weight) - xbar_log),
        a ~ dunif(40, 200),
        b ~ dnorm(0, 2.5),
        sigma ~ dunif(0, 50)
    ),
    data=d
)

# So now a 1 increase in log(weight), which means an increase from 
# xbar_log 3.44233 to 4.44233 (31.2597 to 84.9727kg) is associated with a 
# 46cm height increase
# Intercept still has same interpretation
precis(h3)

# b) Plot it, and that looks much better!
mu <- link(h3, data=data.frame(weight=new_weights))
mu_PI <- apply(mu, 2, PI, 0.89)
pp <- sim(h3, data=data.frame(weight=new_weights))
pp_PI <- apply(pp, 2, PI, 0.89)

plot(d2$weight, d2$height, col=col.alpha(rangi2, 0.3), pch=16)
shade(mu_PI, new_weights, col=col.alpha("black", 0.3))
shade(pp_PI, new_weights, col=col.alpha("black", 0.3))

# 4H4
# Prior predictive for parabolic polynomial
# Parabolic polynomial regression
# Modify priors on alpha, beta1, beta2 so that model looks sensible
data(Howell1)
d <- Howell1
d$weight_centered <- d$weight - mean(d$weight)
d$weight_standardized <- d$weight_centered / sd(d$weight)
plot(height ~ weight, d)
# Also create variable for squared term
d$weight_standardized_2 <- d$weight_standardized**2

plot_priors <- function(prior_a, prior_b1, prior_b2, N=100) {
    sim_a <- prior_a(N)
    sim_b1 <- prior_b1(N)
    sim_b2 <- prior_b2(N)
    plot(NULL, xlim=c(0, 70), ylim=c(50, 180), 
         xlab="weight", ylab="height")
    for (i in 1:N) {
        curve({
                x_std <- (x - mean(d$weight)) / sd(d$weight)
                sim_a[i] + sim_b1[i] * x_std + sim_b2[i] * x_std**2
            },
              from=min(d$weight), to=max(d$weight),
              add=TRUE,
              col=col.alpha("black", 0.2)
        )
    }
}

# The default priors that we had
# Don't look very good!
# Want to see some kind of initial increase that then peters off
plot_priors(
    function(N) rnorm(N, 178, 20), 
    function(N) rlnorm(N, 0, 1), 
    function(N) rnorm(N, 0, 1)
)

# Will want to lower alpha now we're not just looking at adults
# How to enforce this?
# Think 1 SD should be ~10cm
# This has some decent shapes but it doesn't tailor off and has some negative parabolas
# According to https://rcompanion.org/handbook/I_11.html:
# The plateau value is calculated as a + b * clx –0.5 * b * clx
# Where a = intercept, b = slope, clx = quadratic parameter
# The quadratic coefficient is calculated as –0.5 * b / clx.
# Want to plateau around 170, so 170 = 140 + (10 * clx) - (0.5 * 10 * clx)
# So 30 = 10clx - 5clx
# clx = 6
# So the quadratic coef = -0.5 * 10 / 6 = -30
# That looks reasonable ish
plot_priors(
    function(N) rnorm(N, 140, 20), 
    function(N) rnorm(N, 10, 5), 
    function(N) rnorm(N, -0.83, 2)
)

# 4H5
# Predict blossom date from march temperature using either linear, polynomial, or spline
data(cherry_blossoms)
d <- cherry_blossoms

d2 <- d[complete.cases(d$doy, d$temp), ]
# Potentially non-linear, but tbh none look good
plot(d2$temp, d2$doy)

# Fit linear model
# Use wide priors as since I don't know anything about data would have to look 
# at it first to form prior - bad practice!
h5_linear <- quap(
    alist(
        doy ~ dnorm(mu, sigma),
        mu <- a + b * temp,
        a ~ dnorm(0, 10),
        b ~ dnorm(0, 5),
        sigma ~ dexp(1)
    ),
    data=d2
)
# At 0C would expect to have doy at 118th day of year
# Every +1C corresponds to 2.3 days earlier
precis(h5_linear)

# Plot linear fit
# Not a great fit, a lot of variance!
temp_seq <- seq(min(d2$temp), max(d2$temp), length.out=50)
plot_h5 <- function(mod) {
    mu <- link(mod, data=data.frame(temp=temp_seq))
    mu_PI <- apply(mu, 2, PI, 0.97)
    pp <- sim(mod, data=data.frame(temp=temp_seq))
    pp_PI <- apply(pp, 2, PI, 0.97)
    plot(d2$temp, d2$doy)
    shade(mu_PI, temp_seq, col=col.alpha("black", 0.3))
    shade(pp_PI, temp_seq, col=col.alpha("black", 0.3))
}
plot_h5(h5_linear)

# Try quadratic model
h5_quad <- quap(
    alist(
        doy ~ dnorm(mu, sigma),
        mu <- a + b * temp + d * temp**2,
        a ~ dnorm(0, 10),
        b ~ dnorm(0, 5),
        d ~ dnorm(0, 5),
        sigma ~ dexp(1)
    ),
    data=d2
)
# Less easy to interpret model now
precis(h5_quad)

# Plot quadratic fit
# Potentially better, but potentially too steep a drop-off
plot_h5(h5_quad)

# Fit spline model - not very wiggly dataset so will go for 3 knotes
num_knots <- 3
knot_list <- quantile(d2$temp, probs=seq(0, 1, length.out=num_knots))
B <- bs(d2$temp, knots=knot_list[-c(1, num_knots)],
        degree=3, intercept=TRUE)

# Now fit the model
h5_poly <- quap(
    alist(
        doy ~ dnorm(mu, sigma),
        mu <- a + B %*% w,
        a ~ dnorm(0, 10),
        w ~ dnorm(0, 10),
        sigma ~ dexp(1)
    ),
    data=data.frame(doy=d2$doy, B=B),
    start=list(w=rep(0, ncol(B)))
)

# Plot MAP posterior basis effects (modelled weights * basis functions)
# Doesn't really help provide insight...
post <- extract.samples(h5_poly)
w <- apply(post$w, 2, mean)
plot(NULL, xlim=range(d2$temp), ylim=c(-6, 16), xlab="temp", ylab="basis*weight")
for (i in 1:ncol(B)) {
    lines(d2$temp, w[i] * B[, i])
}

B_seq <- bs(temp_seq, knots=knot_list[-c(1, num_knots)],
            degree=3, intercept=TRUE)
mu <- link(h5_poly, data=list(B=B_seq))
mu_PI <- apply(mu, 2, PI, 0.97)
pp <- sim(h5_poly, data=list(B=B_seq))
pp_PI <- apply(pp, 2, PI, 0.97)
plot(d2$temp, d2$doy)
shade(mu_PI, temp_seq, col=col.alpha("black", 0.3))
# Can't plot posterior predictive as Sim returns 250 points vs
# the 50 that actually provide
# Link returns 50 though and basically fits a linear fit
shade(pp_PI, temp_seq, col=col.alpha("black", 0.3))

# In which case my conclusion is that a linear fit is sufficient for this data
# But low correlation so still not a very good fit anyway
cor(d2$temp, d2$doy)

# 4H6
# That's this model
num_knots <- 15
knot_list <- quantile(d2$year, probs=seq(0, 1, length.out=num_knots))
B <- bs(d2$year, knots=knot_list[-c(1, num_knots)],
        degree=3, intercept=TRUE)
m4.7 <- quap(
    alist(
        D ~ dnorm(mu, sigma),
        mu <- a + B %*% w,
        a ~ dnorm(100, 10),
        w ~ dnorm(0, 10),
        sigma ~ dexp(1)
    ),
    data=list(D=d2$doy, B=B),
    start=list(w=rep(0, ncol(B)))
)

mu <- link(m4.7)
mu_PI <- apply(mu, 2, PI, 0.97)
plot(d2$year, d2$doy, col=col.alpha(rangi2, 0.3), pch=16)
shade(mu_PI, d2$year, col=col.alpha("black", 0.5))

plot_h6 <- function(prior_a, prior_w, N=100) {
    sim_a <- prior_a(N)
    sim_w <- replicate(num_knots+2, n_prior_w(N))
    plot(NULL, xlim=range(d2$year), ylim=c(50, 150),
         xlab="Year", ylab="Date of first blossom")
    for (i in 1:N) {
        y <- sim_a[i] + B %*% sim_w[i, ]
        lines(d2$year, y[, 1], col=col.alpha("black", 0.2))
    }
    
}

# Default priors
plot_h6(
    function(N) rnorm(N, 100, 10),
    function(N) rnorm(N, 0, 10),
)

# What if tighten variance on weights
# Possibly smoother? Hard to tell
plot_h6(
    function(N) rnorm(N, 100, 10),
    function(N) rnorm(N, 0, 1),
)

# What if I make them positive?
# Hard to tell again...
plot_h6(
    function(N) rnorm(N, 100, 10),
    function(N) rnorm(N, 10, 2),
)

# Increasing variance on weights seem to mae them subtley more positive
plot_h6(
    function(N) rnorm(N, 100, 10),
    function(N) rnorm(N, 10, 10),
)

# Really can't tell much difference at all...
plot_h6(
    function(N) rnorm(N, 100, 10),
    function(N) rnorm(N, 500, 10),
)

# 4H7 fit the spline without an intercept, instead using the first basis function
# If I use the Basis functions that had:
# INTERCEPT but kept the first basis function we get a junk intercept (because the 
# basis function is constant 0)
h7_1 <- quap(
    alist(
        D ~ dnorm(mu, sigma),
        mu <- B %*% w,
        w ~ dnorm(0, 10),
        sigma ~ dexp(1)
    ),
    data=list(D=d2$doy, B=B4),
    start=list(w=rep(0, ncol(B4)))
)
precis(h7_1, depth=2)

# If we use INTERCEPT but removing first knot (this was model we used)
h7_2 <- quap(
    alist(
        D ~ dnorm(mu, sigma),
        mu <- B %*% w,
        w ~ dnorm(0, 10),
        sigma ~ dexp(1)
    ),
    data=list(D=d2$doy, B=B),
    start=list(w=rep(0, ncol(B)))
)
precis(h7_2, depth=2)
# No point repeating this for when INTERCEPT argument is FALSE

# So neither of these seem to be the intercept
# Surely the answer isn't just add a column of 1s to the basis functions?
B5 <- B
B5 <- cbind(B5, 1)
# Nope that hasn't worked, probably due to violating the constraint
h7_3 <- quap(
    alist(
        D ~ dnorm(mu, sigma),
        mu <- B %*% w,
        w ~ dnorm(0, 10),
        sigma ~ dexp(1)
    ),
    data=list(D=d2$doy, B=B5),
    start=list(w=rep(0, ncol(B5)))
)
precis(h7_3, depth=2)
