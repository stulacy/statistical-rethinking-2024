library(rethinking)
library(tidyverse)

# 7.1.1 More parameters (almost) always improves fit ----------------------

sppnames <- c("afarensis", "africanus", "habilis", "boisei", "rudolfensis", "ergaster", "sapiens")
brainvolcc <-c(438, 452, 612, 521, 752, 871, 1350)
masskg <- c(37.0, 35.5, 34.5, 41.5, 55.5, 61.0, 53.5)
d <- tibble(species=sppnames, brain=brainvolcc, mass=masskg)
d

# Standardize, NB: keep brain 0-1 so always positive (could use log transform)
d$mass_std <- scale(d$mass)
d$brain_std <- d$brain / max(d$brain)
d$brain_log <- log(d$brain)

# Fit linear model first, wanting to see if different species have different brain volumes
# after adjusting for size
# Deliberately flat priors
# Q: Why using log-normal for sigma rather than exp?!
m7.1 <- quap(
    alist(
        brain_std ~ dnorm(mu, exp(log_sigma)),
        mu <- a + b * mass_std,
        a ~ dnorm(0.5, 1),
        b ~ dnorm(0, 10),
        log_sigma ~ dnorm(0, 1)  # log(x) ~ norm() is the same as x ~ log-norm()
    ),
    data=d
)

# Calculating R2 manually
R2_calc <- function(mod, outcome) {
    # Obtain posterior predictive
    s <- sim(mod, refresh=0)
    # Get residuals of MEAN estimates. In Bayesian world should really take full posterior
    # into account but R2 is only about MEAN estimates
    r <- colMeans(s) - outcome
    resid_var <- var2(r)
    # NB: 'var' in R is the Frequentist estimator so it has the adjusted denominator
    resid_var
    var(r)
    outcome_var <- var2(outcome)
    1 - resid_var / outcome_var
}

set.seed(12)
# Pretty poor!
R2_calc(m7.1, d$brain_std)
    
# Now will fit 5 more models in increasing complexity (i.e. higher degree poly)
m7.2 <- quap(
    alist(
        brain_std ~ dnorm(mu, exp(log_sigma)),
        mu <- a + b[1] * mass_std + b[2] * mass_std^2,
        a ~ dnorm(0.5, 1),
        b ~ dnorm(0, 10),  # Can pass vector
        log_sigma ~ dnorm(0, 1)  # log(x) ~ norm() is the same as x ~ log-norm()
    ),
    data=d,
    start=list(b=rep(0, 2))  # Here is where B is defined as having length 2
)

m7.3 <- quap(
    alist(
        brain_std ~ dnorm(mu, exp(log_sigma)),
        mu <- a + b[1] * mass_std + b[2] * mass_std^2 + b[3] * mass_std^3,
        a ~ dnorm(0.5, 1),
        b ~ dnorm(0, 10),  # Can pass vector
        log_sigma ~ dnorm(0, 1)  # log(x) ~ norm() is the same as x ~ log-norm()
    ),
    data=d,
    start=list(b=rep(0, 3))  # Here is where B is defined as having length 3
)

m7.4 <- quap(
    alist(
        brain_std ~ dnorm(mu, exp(log_sigma)),
        mu <- a + b[1] * mass_std + b[2] * mass_std^2 + b[3] * mass_std^3 + 
                  b[4] * mass_std^4,
        a ~ dnorm(0.5, 1),
        b ~ dnorm(0, 10),  # Can pass vector
        log_sigma ~ dnorm(0, 1)  # log(x) ~ norm() is the same as x ~ log-norm()
    ),
    data=d,
    start=list(b=rep(0, 4))  # Here is where B is defined as having length 4
)

m7.5 <- quap(
    alist(
        brain_std ~ dnorm(mu, exp(log_sigma)),
        mu <- a + b[1] * mass_std + b[2] * mass_std^2 + b[3] * mass_std^3 + 
                  b[4] * mass_std^4 + b[5] * mass_std^5,
        a ~ dnorm(0.5, 1),
        b ~ dnorm(0, 10),  # Can pass vector
        log_sigma ~ dnorm(0, 1)  # log(x) ~ norm() is the same as x ~ log-norm()
    ),
    data=d,
    start=list(b=rep(0, 5))  # Here is where B is defined as having length 5
)

# NB: For the final model we hardcode the variance, will see why shortly
m7.6 <- quap(
    alist(
        brain_std ~ dnorm(mu, 0.001),
        mu <- a + b[1] * mass_std + b[2] * mass_std^2 + b[3] * mass_std^3 + 
                  b[4] * mass_std^4 + b[5] * mass_std^5 + b[6] * mass_std^6,
        a ~ dnorm(0.5, 1),
        b ~ dnorm(0, 10)  # Can pass vector
    ),
    data=d,
    start=list(b=rep(0, 6))  # Here is where B is defined as having length 5
)

# Want to plot posterior predictive throughout the x-space
plot_mod <- function(mod) {
    r2 <- R2_calc(mod, d$brain_std)
    mass_seq <- seq(from=min(d$mass_std), to=max(d$mass_std), length.out=100)
    l <- link(mod, data=list(mass_std=mass_seq))
    mu <- colMeans(l)
    ci <- apply(l, 2, PI)
    d |>
        ggplot(aes(x=mass_std)) +
            geom_ribbon(aes(x=mass_seq, ymin=min, ymax=max), data=tibble(mass_seq, 
                                                                         min=ci[1, ],
                                                                         max=ci[2, ]),
                        alpha=0.5) +
            geom_line(aes(x=mass_seq, y=mu), data=tibble(mass_seq, mu)) +
            geom_point(aes(y=brain_std)) +
            geom_text(aes(x=x, y=y, label=lab), vjust=1,
                      tibble(x=-0.5, y=Inf, lab=sprintf("R2 = %.2f", r2)))
}

library(patchwork)
mods <- list(m7.1, m7.2, m7.3, m7.4, m7.5, m7.6)
plts <- lapply(mods, plot_mod)

# Can see that R2 always improves when add parameter, is this always true or just for this specific
# dataset?
# Think always true!
# I.e. as number of parameters approaches number of points can get a perfect fit, as seen in the last
# plot where have 6 parameters (well 6 + intercept=7) for 7 data points
# This for the 6 degree poly is why had to fix variance, otherwise would tend to zero
# Also get negative brain size! (NB: thought we'd made this impossible? actually no as had
# scaled 0-1 but still possible to get LP < 0)
wrap_plots(plts, tag_level='new')

# Minimum Descriptive Length (MDL) sounds like AIC, way of compressing data into as few parameters
# that adequately describe the data

# 7.2 Entropy and accuracy ----------------------------------------------------

# R code 7.14
# Manually calculating log-pointwise-predictive-density (lppd)
set.seed(1)
logprob <- sim(m7.1, ll=TRUE, n=1e4)
foo <- sim(m7.1, ll=FALSE, n=1e4)
# So sim gives samples for each datapoint
summary(foo)
# ll=TRUE gives log probabilities instead
summary(logprob)
# If exponentiate get probabilities > 1, which I think is valid because they are densities
summary(exp(logprob))
n <- ncol(logprob)
ns <- nrow(logprob)
# Calculate sum probs per data point (column) across all samples
f <- function(i) log_sum_exp(logprob[, i]) - log(ns)
lppd_man <- sapply(1:n, f)
lppd_man
# The same as the lppd function (from rethinking)
# BUT! What exactly is this doing? I.e. how can it be calculated for new data points?
lppd(m7.1, n=1e4)
# Log sum exp is used to calculate the sum of log values more stabley
# Say have lots of probabilities
x <- runif(5)
# That access on a log scale for ease of computation
logx <- log(x)
# If we want to calculate the log of their sum it can be unreliable to directly use
# on the original scale
log(sum(x))
# So instead use LSE on log scale
log_sum_exp(log(x))

# 7.2.5 Scoring the right data --------------------------------------------
set.seed(1)
# Calculate total lppd overall all (7) data points per model
# Notice how it gets higher with the more parameters, that's not ideal! Just
# the same flaw as R2
sapply(mods, function(x) sum(lppd(x)))

# 7.4 Predicting Predictive Accuracy --------------------------------------
# Overthinking: WAIC calculations
# R Code 7.19
# Going to code WAIC from scratch
data(cars)
m <- quap(
    alist(
        dist ~ dnorm(mu, sigma),
        mu <- a + b*speed,
        a ~ dnorm(0, 100),
        b ~ dnorm(0, 10),
        sigma ~ dexp(1)
    ), data=cars
)
set.seed(94)
post <- extract.samples(m, n=1e3)
# Calculate log-likelihood of each observation i and each posterior sample s
n_samples <- 1e3
logprob <- sapply(1:n_samples, function(s) {
    mu <- post$a[s] + post$b[s] * cars$speed
    dnorm(cars$dist, mu, post$sigma[s], log=TRUE)
})
# Get one log-likelihood per point per posterior sample
# Could easily run this over a new dataset! Doesn't have to be the training set
dim(logprob)
# Calculate lppd (the Bayesian Deviance) by getting the log of the sum of all probabilties
# per posterior sample
# Can do this using the LSE
log_sum_exp(logprob[1, ])

# Just to make it clear why this is used, if we instead work on a raw probability scale
prob <- sapply(1:n_samples, function(s) {
    mu <- post$a[s] + post$b[s] * cars$speed
    dnorm(cars$dist, mu, post$sigma[s], log=FALSE)
})

prob[1:5, 1:5]
exp(logprob[1:5, 1:5])

# And calculate the log of the sums (as per the eqn on page 210 in the Overthinking box)
# We get the same answer
# The LSE is only necessary because we want to work with logprobs to have more computational
# stability
log(sum(prob[1, ]))

# So getting the sum over samples, and then subtracting the log of number samples,
# which is the same as dividing by the number of samples shown in the eqn on P 210
n_cases <- nrow(cars)
lppd <- sapply(1:n_cases, function(i) log_sum_exp(logprob[i, ]) - log(n_samples))
# This gives us the lppd per observation, can get the overall lppd using sum(lppd)

# Can calculate the penalty term per obs, which is variance across samples for each observation
pWAIC <- sapply(1:n_cases, function(i) var(logprob[i, ]))
pWAIC

# Then combine these two components to get WAIC:
# Pretty much the same, with the only difference being simulation variance from the 
# extract.samples function call
-2 * (sum(lppd) - sum(pWAIC))
WAIC(m)

# Can compute standard error of WAIC using the observation level lppd and pWAIC
waic_vec <- -2 * (lppd - pWAIC)
# Again, this agrees with the std_error in WAIC function
sqrt(n_cases * var(waic_vec))

# 7.5.1 Model mis-selection -----------------------------------------------
# Refit models m6 
N <- 100
set.seed(909)
height <- rnorm(N, 10, 2)
leg_prop <- runif(N, 0.4, 0.5) # leg as proportion of height
leg_left <- leg_prop*height + rnorm(N, 0, 0.02) # add error so legs aren't same size exactly
leg_right <- leg_prop*height + rnorm(N, 0, 0.02)
d <- tibble(height, leg_left, leg_right)
m6.1 <- quap(
    alist(
        height ~ dnorm(mu, sigma),
        mu <- a + bl*leg_left + br*leg_right,
        a ~ dnorm(10, 100),
        bl ~ dnorm(2, 10),
        br ~ dnorm(2, 10),
        sigma ~ dexp(1)
    ),
    data=d
)
set.seed(71)
N <- 100  # Number of plants
h0 <- rnorm(N, 10, 2)  # Initial heights
treatment <- rep(0:1, each=N/2)
table(treatment)  # evenly split
fungus <- rbinom(N, size=1, prob=0.5 - treatment*0.4)
table(fungus) # fungus is more likely in untreated plants
h1 <- h0 + rnorm(N, 5-3*fungus) # final height is a function of fungus, NOT treatment
d <- tibble(h0=h0, h1=h1, treatment=treatment, fungus=fungus)

m6.6 <- quap(
    alist(
        h1 ~ dnorm(mu, sigma),
        mu <- h0*p,
        p ~ dlnorm(0, 0.25),
        sigma ~ dexp(1)
    ), data=d
)

m6.7 <- quap(
    alist(
        h1 ~ dnorm(mu, sigma),
        mu <- h0*p,
        p <- a + bT*treatment + bF*fungus,
        a ~ dlnorm(0, 0.25),
        bT ~ dnorm(0, 0.5),
        bF ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), data=d
)

m6.8 <- quap(
    alist(
        h1 ~ dnorm(mu, sigma),
        mu <- h0*p,
        p <- a + bT*treatment,
        a ~ dlnorm(0, 0.25),
        bT ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), data=d
)
# 
set.seed(11)
WAIC(m6.1)
set.seed(77)
compare(m6.6, m6.7, m6.8, func=WAIC)

# Or can use PSIS, but get basically the same results
compare(m6.6, m6.7, m6.8, func=PSIS)

# Standard error is uncertainty around WAIC/PSIS
# If want uncertainty around how different models are, need to use dSE
# This is calculated from the observation-wise estimates
set.seed(91)
waic_m6.7 <- WAIC(m6.7, pointwise=TRUE)$WAIC
# Use pointwise=TRUE to get the estimate per observation
length(waic_m6.7)
# If it's FALSE just get the sum
length(WAIC(m6.7, pointwise=FALSE)$WAIC)

waic_m6.8 <- WAIC(m6.8, pointwise=TRUE)$WAIC
n <- length(waic_m6.7)
diff_m6.8_m6.7 <- waic_m6.7 - waic_m6.8
# This is the standard error of the difference
sqrt(n*var(diff_m6.8_m6.7))
# So the difference between M6.8 and M6.7 is 41 +/- 10, so the difference is
# greater than zero at least!

# There's also a plot method, very useful!
# Solid points are in-sample deviance
# Hollow points are WAIC/PSIS
# The triangle line shows the standard error of the difference between 2 models (not really sure how this works)
plot(compare(m6.6, m6.7, m6.8, func=PSIS))

# The compare table only shows the difference between the best model and any others
# But it does calculate all pairwise differences
set.seed(93)
# So here the intercept only model (6.6) has only has a std error in difference of 4.8 
# vs the treatment only model (6.8), and their WAIC is only 3 units apart, so this does
# cross zero
cmp <- compare(m6.6, m6.7, m6.8, func=PSIS)
cmp@dSE

# 7.5.2 Outliers and other illusions --------------------------------------
# Refit Waffle data as remember had the Mormon state outliers
data("WaffleDivorce")
d <- WaffleDivorce

# standardize
d$A <- scale(d$MedianAgeMarriage)
d$D <- scale(d$Divorce)
d$M <- scale(d$Marriage)

# Regress divorce rate on age
m5.1 <- quap(
    alist(
        D ~ dnorm(mu, sigma),
        mu <- a + bA * A,
        a ~ dnorm(0, 0.2),
        bA ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), 
    data=d
)

m5.2 <- quap(
    alist(
        D ~ dnorm(mu, sigma),
        mu <- a + bM * M,
        a ~ dnorm(0, 0.2),
        bM ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ),
    data=d
)

m5.3 <- quap(
    alist(
        D ~ dnorm(mu, sigma),
        mu <- a + bM*M + bA*A,
        a ~ dnorm(0, 0.2),
        bM ~ dnorm(0, 0.5),
        bA ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ),
    data=d
)

# m5.1: D ~ A
# m5.2: D ~ M
# m5.3: D ~ A + M (M now has no effect)
set.seed(24071847)
compare(m5.1, m5.2, m5.3, func=PSIS)  # Raises warning
# Here the model with just age is preferred, which makes sense as Marriage offers nothing
# else

set.seed(24071847)
PSIS_m5.3 <- PSIS(m5.3, pointwise=TRUE)  # Raises warning
set.seed(24071847)
WAIC_m5.3 <- WAIC(m5.3, pointwise=TRUE)
# PSIS and WAIC have values for identifying outliers, or strong influences
# This is k in PSIS (when > 0.5 the warning is raised), and the WAIC penalty
# Can think of WAIC penalty as # effective parameters, so this one obs is so
# influential it's like adding 2 parameters 
# Here can see that Idaho (top right) is a massive outlier in both of these
# ME is the other outlier
plot(PSIS_m5.3$k, WAIC_m5.3$penalty, xlab="PSIS Pareto k", 
     ylab="WAIC penalty", col=rangi2, lwd=2)

# Refit with Student with 2DOF for thicker tails to reduce outliers effects
m5.3t <- quap(
    alist(
        D ~ dstudent(2, mu, sigma),
        mu <- a + bM*M + bA*A,
        a ~ dnorm(0, 0.2),
        bM ~ dnorm(0, 0.5),
        bA ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ),
    data=d
)

# No outliers, far lower WAIC, no Pareto close to 0.5!
set.seed(24071847)
PSIS_m5.3t <- PSIS(m5.3t, pointwise=TRUE)  # Raises warning
set.seed(24071847)
WAIC_m5.3t <- WAIC(m5.3t, pointwise=TRUE)
plot(PSIS_m5.3t$k, WAIC_m5.3t$penalty, xlab="PSIS Pareto k", 
     ylab="WAIC penalty", col=rangi2, lwd=2)

# The t model has given age a stronger effect!
# This is because Idaho had a very low divorce rate for its age and was
# shrinking this a bit. With more outliers allowed, the coefficient can be stronger
plot(coeftab(m5.3, m5.3t))

# Problems ----------------------------------------------------------------

# 7E1: 3 motivations for information entropy (quantification of uncertainty in data):
# 1. continuous
# 2. should increase with number of events
# 3. additive

# 7E2: coin lands up 70% of the time
# definition is -sum(p_i * log(p_i)) across all possible events i
# Here 0.250
-sum(0.7 * log(0.7))

# 7E3: 1.38
probs <- c(0.2, 0.25, 0.25, 0.3)
-sum(probs * log(probs))

# 7E4: 1.084
# Mathematically not possible! But according to the Overthinking box on p206
# just drop these events
probs <- c(0.3, 0.3, 0.3)
-sum(probs * log(probs))

# 7M1: Definitions of AIC and WAIC including their assumptions
# AIC: -2lppd + 2p (p = # parameters)
#   Assumes:
#      - priors are flat (or overwhelmed by likelihood) (not realistic, espec for multilevel)
#      - posterior is approx multivariate Gaussian
#      - sample size N is much greater than # parameters p
# WAIC: -2(lppd - penalty) where penalty is based on observation level variance of logprob
#   Assumes:
#      - Can't see anything
#   So is AIC when have flat priors, posterior is approx multi Gaussian, N >> p

# 7M2
# Model selection vs model comparison
# Comparison is just that, comparing models. Selection is a decision making process
# to use a single one. Which isn't really best practice. All models can tell us 
# useful information that should be considered in context, rather than relying solely
# upon a single 'best' model

# 7M3
# When using information criteria, why must all models be fit to the same data?
# Because the lppd is the SUM over all observations, see Overthinking on p210
# not the AVG, so with more observations it will be on a different scale

# 7M4
# What happens to effective number of parameters as a prior becomes more concentrated?
# Effective # parameters

# Let's use the m5.3 model that has the Idaho outlier
plot(PSIS_m5.3$k, WAIC_m5.3$penalty, xlab="PSIS Pareto k", 
     ylab="WAIC penalty", col=rangi2, lwd=2)

# And then relax the priors
m5.3_flat <- quap(
    alist(
        D ~ dnorm(mu, sigma),
        mu <- a + bM*M + bA*A,
        a ~ dnorm(0, 2),
        bM ~ dnorm(0, 5),
        bA ~ dnorm(0, 5),
        sigma ~ dexp(1)
    ),
    data=d
)

m5.3_narrow <- quap(
    alist(
        D ~ dnorm(mu, sigma),
        mu <- a + bM*M + bA*A,
        a ~ dnorm(0, 0.1),
        bM ~ dnorm(0, 0.25),
        bA ~ dnorm(0, 0.25),
        sigma ~ dexp(1)
    ),
    data=d
)

PSIS_narrow <- PSIS(m5.3_narrow, pointwise=TRUE) 
PSIS_flat <- PSIS(m5.3_flat, pointwise=TRUE) 
WAIC_narrow <- WAIC(m5.3_narrow, pointwise=TRUE) 
WAIC_flat <- WAIC(m5.3_flat, pointwise=TRUE) 

PSISs <- list("default"=PSIS_m5.3, "narrow"=PSIS_narrow, "flat"=PSIS_flat)
WAICs <- list("default"=WAIC_m5.3, "narrow"=WAIC_narrow, "flat"=WAIC_flat)
PSISs_eff <- map_dfr(PSISs, function(x) x$k, .id="model")
WAICs_eff <- map_dfr(WAICs, function(x) x$penalty, .id="model")
colnames(PSISs_eff) <- paste0("PSIS_", colnames(PSISs_eff))
colnames(WAICs_eff) <- paste0("WAIC_", colnames(WAICs_eff))

# Plot!
cbind(
    d,
    PSISs_eff,
    WAICs_eff
) |>
    select(Loc, contains("PSIS"), contains("WAIC")) |>
    pivot_longer(-Loc, names_pattern="(.+)_(.+)", names_to=c("metric", "model")) |>
    pivot_wider(names_from=metric, values_from=value) |>
    mutate(model = factor(model, levels=c("narrow", "default", "flat"))) |>
    ggplot(aes(x=PSIS, y=WAIC)) +
        geom_point() +
        facet_wrap(~model) +
        theme_bw()

# Making the priors flatter makes the effective parameters more
# I.e. increases the penalties

# 7M5
# Informative priors reduce the effect of the likelihood on the posterior
# So that provided the priors don't align with the data, it pulls
# the posterior away from being solely affected by the likelihood

# 7M6
# Overly informative prior means the prior have a dominating affect over the
# likelihood and no matter the data aren't able to have much effect.

# 7H1