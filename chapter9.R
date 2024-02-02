library(rethinking)
library(tidyverse)

# 9.1 Good King Markov and his island kingdom -----------------------------
# Metropolis algorithm
metropolis <- function(n_samples, starting) {
    positions <- rep(0, n_samples)
    current <- starting
    for (i in 1:n_samples) {
        # record current position
        positions[i] <- current
        # flip coin to move either side
        proposal <- current + sample(c(-1, 1), size=1)
        # ensure proposal is within bounds
        proposal <- min(proposal, 10)
        proposal <- max(proposal, 1)
        # prob of moving is ratio. If proposal > current then move
        # with 100% prob
        prob_move <- proposal/current
        current <- ifelse(runif(1) < prob_move, proposal, current)
    }
    positions
}

results <- metropolis(1e5, 3)
# Seemingly random sampling
plot(1:100, results[1:100])
# But are clearly sampling in proportion to island size as desired
plot(table(results))

# 9.2 Metropolis algorithms -----------------------------------------------
# Metropolis is grandfather algorithm
# Metropolis-Hastings allows for asymmetric proposals
# Gibbs is a variant of MH that is more efficient in some way by using 'adaptive proposals' 
#   based on conjugacy
# Metropolis and Gibbs can get stuck in high-dimensional problems and/or with correlated parameters

# Simulation to show Concentration of Measure, that as dimensions increase, the bulk of
# the probability mass lies further and further from the mode
sim_conc_measure <- function(n_d, n_samples=1e3) {
    x <- rmvnorm(n_samples, rep(0, n_d), diag(n_d))
    rad_dist <- function(x) sqrt(sum(x**2))
    # Calculate radial distance for each point
    sapply(1:n_samples, function(i) rad_dist(x[i, ]))
}
dens(sim_conc_measure(1000), xlim=c(0, 34))
dens(sim_conc_measure(100),  add=TRUE)
dens(sim_conc_measure(10), add=TRUE)
dens(sim_conc_measure(1), add=TRUE)

# Hamiltonian Monte Carlo -------------------------------------------------
# AKA Hybrid Monte Carlo (HMC)
# Much more efficient than Metropolis or Gibbs but much computationally expensive
# Can handle very complex models

# Example HMC for 2D Gaussin dataset:
# xi ~ N(ux, 1)
# yi ~ N(uy, 1)
# ux ~ N(0, 0.5)
# uy ~ N(0, 0.5)

# Firstly need function to calculate neg-log-prob at data
# q are the parameters, the unknown values getting sampled
U <- function(q, a=1, b=1, k=0, d=1) {
    # NB: Seems to use y and x from global scope rather than function args
    muy <- q[1]
    mux <- q[2]
    # Need the sum on y and x as have multiple data points to evaluate likelihood over
    # Where have just 1 muy and mux
    U <- sum(dnorm(y, muy, 1, log=TRUE)) + sum(dnorm(x, mux, 1, log=TRUE)) +
        dnorm(muy, a, b, log=TRUE) + dnorm(mux, k, d, log=TRUE)
    return(-U) # Work with negative log probs!
}

# Need a corresponding gradient function, which is relatively straight forward for Gaussian
# likelihoods: dlogN(y | a, b) / da = (y - a) / (b^2)
# Not sure how this would change if sigma was a parameter and not just mean
U_gradient <- function(q, a=0, b=1, k=0, d=1) {
    muy <- q[1]
    mux <- q[2]
    G1 <- sum(y - muy) + (a - muy)/b**2 # dU/dmuy
    G2 <- sum(x - mux) + (k - mux)/d**2 # dU/dmux
    return(c(-G1, -G2))  # Again negative log prob!
}

# Create some dummy data
set.seed(7)
y <- rnorm(50)
x <- rnorm(50)
y <- as.numeric(scale(y))
x <- as.numeric(scale(x))

# Now run! Using the rethinking::HMC2 function

library(shape)  # For plotting
# red means dH > 0.1? What does that mean?
Q <- list()
Q$q <- c(-0.1, 0.2)  # Starting values of muy and muyx. Wrapped in a list because it's used between HMC iterations
pr <- 0.3
plot(NULL, ylab="muy", xlab="mux", xlim=c(-pr, pr), ylim=c(-pr, pr))
step <- 0.03
L <- 28  # Working example, set to 28 to see U-turns
n_samples <- 9
path_col <- col.alpha("black", 0.5)
points(Q$q[1], Q$q[2], pch=4, col="black")
for (i in 1:n_samples) {
    Q <- HMC2(U, U_gradient, step, L, Q$q)
    if (n_samples < 10) {
        for (j in 1:L) {
            K0 <- sum(Q$ptraj[j,]**2) / 2 # Kinetic energy
            lines(Q$traj[j:(j+1), 1], Q$traj[j:(j+1), 2], col=path_col, lwd=1+2*K0)
        }
        points(Q$traj[1:(L+1), ], pch=16, col="white", cex=0.35)
        Arrows(Q$traj[L, 1], Q$traj[L, 2], Q$traj[L+1, 1], Q$traj[L+1, 2],
               arr.length=0.35, arr.adj=0.7)
        text(Q$traj[L=1, 1], Q$traj[L+1, 2], i, cex=0.8, pos=4, offset=0.4)
    }
    points(Q$traj[L+1, 1], Q$traj[L+1, 2], pch=ifelse(Q$accept==1, 16, 1),
           col=ifelse(abs(Q$dH)>0.1, 'red', 'black'))
}

# Easy HMC: ulam ----------------------------------------------------------

# Reload rugged dataset
data(rugged)
d <- rugged
d$log_gdp <- log(d$rgdppc_2000)
dd <- d[complete.cases(d$rgdppc_2000), ]
dd$log_gdp_std <- dd$log_gdp / mean(dd$log_gdp)
dd$rugged_std <- dd$rugged / max(dd$rugged)
dd$cid <- ifelse(dd$cont_africa==1, 1, 2)

m8.3 <- quap(
    alist(
        log_gdp_std ~ dnorm(mu, sigma),
        mu <- a[cid] + b[cid]*(rugged_std - 0.215),
        a[cid] ~ dnorm(1, 0.1),
        b[cid] ~ dnorm(0, 0.3),
        sigma ~ dexp(1)
    ),
    data=dd
)
precis(m8.3, depth=2)

# For ulam need to have transferred everything up front, but have already done that here
# Best practice to make clean dataset (as a list too)
dat_slim <- list(
    log_gdp_std=dd$log_gdp_std,
    rugged_std = dd$rugged_std,
    cid = as.integer(dd$cid)
)


m9.1 <- ulam(
    alist(
        log_gdp_std ~ dnorm(mu, sigma),
        mu <- a[cid] + b[cid]*(rugged_std - 0.215),
        a[cid] ~ dnorm(1, 0.1),
        b[cid] ~ dnorm(0, 0.3),
        sigma ~ dexp(1)
    ),
    data=dat_slim, chains=1  # Only difference between quap and ulam is adding in chains
)
# This then translates into Stan and compiles
# Do get a warning about scale parameter being 0, but just once
# Very similar estimates to quap!
precis(m9.1, depth=2)

# Can sample over multiple chains and in parallel
# Seems to always re-compile model rather than identify that executable is up-to-date
# like Stan does
m9.1 <- ulam(
    alist(
        log_gdp_std ~ dnorm(mu, sigma),
        mu <- a[cid] + b[cid]*(rugged_std - 0.215),
        a[cid] ~ dnorm(1, 0.1),
        b[cid] ~ dnorm(0, 0.3),
        sigma ~ dexp(1)
    ),
    data=dat_slim, chains=4, cores=4
)
# same results, look at the ess_bulk! Higher than the 2,000 (500 x 4 chains) samples we have because of very 
# efficient sampling (better than uncorrelated...)
precis(m9.1, depth=2)

# 'show' gives general overview of the MCMC
show(m9.1)

# Graphical summary
pairs(m9.1)

# Trace plot shows very bad initial behaviour, but since this is during warm-up
# this is to be expected and thus isn't that useful
traceplot(m9.1)

# Trim argument doesn't seem to be working...
# Should remove the first x points
traceplot(m9.1, trim=500)

# Can also see rank histogram - want to be even ranking over sample space
trankplot(m9.1)

# 9.5.3 Taming a wild chain -----------------------------------------------
# Example of very flat prior
y <- c(-1, 1)
set.seed(11)
m9.2 <- ulam(
    alist(
        y ~ dnorm(mu, sigma),
        mu <- alpha,
        alpha ~ dnorm(0, 1000),   # Flaaaaaaaaaaat
        sigma ~ dexp(0.0001)      # Flaaaaaaaaaaat
    ),
    data=list(y=y), chains=3
)

# Get divergence warning
# Estimates are far from expected (mean 0 would be reasonable, sd of 0.5)
# And very poor ess and rhat!
precis(m9.2)

# Book says can access Stan model through @stanfit, but there's no such slot!
# But can easily see poor traces
traceplot(m9.2)

set.seed(11)
m9.3 <- ulam(
    alist(
        y ~ dnorm(mu, sigma),
        mu <- alpha,
        alpha ~ dnorm(0, 10),   # Weakly informative
        sigma ~ dexp(1)      # Weakly informative
    ),
    data=list(y=y), chains=3
)

# Much more reasonable!
# And much higher ess and lower rhat!
precis(m9.3)
# Much healthier exploration, just from removing the space away from unrealistic values
traceplot(m9.3)

# 9.5.4 Non-identifiable parameters ---------------------------------------
# Remember when had 2 leg height values in a model and it struggled because it
# was their sum that had the impact, not them individually, but there's an infinite
# way to get them to sum up to the effect size
# Now see how this plays out in MCMC
set.seed(41)
# Have data with N(0, 1)
y <- rnorm(100, 0, 1)

# But assume model
# y ~ N(a1 + a2, 0)
# So there's an infinite number of ways that a1 + a2 = 0
# Also exacerbate this by using flat priors
set.seed(384)
m9.4 <- ulam(
    alist(
        y ~ dnorm(mu, sigma),
        mu <- a1 + a2,
        a1 ~ dnorm(0, 1000),
        a2 ~ dnorm(0, 1000),
        sigma ~ dexp(1)
    ),
    data=list(y=y), chains=3
)
# Strange a1 & a2 values with massive uncertainties!
# Also very poor n_eff and Rhat!
# And a lot of transitions high maximum treedepth limit
precis(m9.4)

# Try weakly informative priors
m9.5 <- ulam(
    alist(
        y ~ dnorm(mu, sigma),
        mu <- a1 + a2,
        a1 ~ dnorm(0, 10),
        a2 ~ dnorm(0, 10),
        sigma ~ dexp(1)
    ),
    data=list(y=y), chains=3
)
# Much better ess, rhat, far fewer max treedepth!
precis(m9.5)

# And much improved traceplot!
traceplot(m9.4)
traceplot(m9.5)

# Problems ----------------------------------------------------------------

# 9E1: Metropolis assumes symmetric proposals

# 9E2: Gibbs can make more efficient proposals, but these require conjugacy
# Doesn't do well with high-dimensional problems or correlated parameters

# 9E3: HMC can't handle discrete parameters (remember issues with mixture models in Stan)

# 9E4: n_eff is how effective sampling was, for some posteriors the proposals are so efficient
# that can end up with n_eff > actual samples

# 9E5: Want Rhat = 1

# 9E6: Want effectively random noise with constant variance and mean

# 9E7: Want flat histogram with all chains having similar shapes

# 9M1: 
# Changing sigma ~ dexp(0, 1) to dunif(0, 1) and seeing effect
# Using rugged terrain dataset, refit model M9.1
d <- rugged
d$log_gdp <- log(d$rgdppc_2000)
dd <- d[complete.cases(d$rgdppc_2000), ]
dd$log_gdp_std <- dd$log_gdp / mean(dd$log_gdp)
dd$rugged_std <- dd$rugged / max(dd$rugged)
dd$cid <- ifelse(dd$cont_africa==1, 1, 2)

dat_slim <- list(
    log_gdp_std=dd$log_gdp_std,
    rugged_std = dd$rugged_std,
    cid = as.integer(dd$cid)
)

m9.1 <- ulam(
    alist(
        log_gdp_std ~ dnorm(mu, sigma),
        mu <- a[cid] + b[cid]*(rugged_std - 0.215),
        a[cid] ~ dnorm(1, 0.1),
        b[cid] ~ dnorm(0, 0.3),
        sigma ~ dexp(1)
    ),
    data=dat_slim, chains=4, cores=4
)

# And now fit with the new prior
m_9m1 <- ulam(
    alist(
        log_gdp_std ~ dnorm(mu, sigma),
        mu <- a[cid] + b[cid]*(rugged_std - 0.215),
        a[cid] ~ dnorm(1, 0.1),
        b[cid] ~ dnorm(0, 0.3),
        sigma ~ dunif(0, 1)
    ),
    data=dat_slim, chains=4, cores=4, log_lik = TRUE
)

# No real impact on sigma at all!
plot(coeftab(m9.1, m_9m1))

# Quite a difference in prior shapes!
dens(rexp(1e3, 1))
dens(runif(1e3, 0, 1))

# So why no impact on sigma?
# Could be because there's easily enough data that prior doesn't impact as much
nrow(d)

# 9M2
# Changing continent slope from N(0, 0.3) to Exp(0.3)
# Will have massive difference! Exp(0.3) for a start constrains positive,
# so that the African continent slope will now be forced to be positive
m_9m2 <- ulam(
    alist(
        log_gdp_std ~ dnorm(mu, sigma),
        mu <- a[cid] + b[cid]*(rugged_std - 0.215),
        a[cid] ~ dnorm(1, 0.1),
        b[cid] ~ dexp(0.3),
        sigma ~ dunif(0, 1)
    ),
    data=dat_slim, chains=4, cores=4, log_lik = TRUE
)
# Sure enough, look at that effect that's now positive on b[2]
plot(coeftab(m_9m1, m_9m2))

# And sure enough the first model has a much better WAIC
compare(m_9m1, m_9m2)

# 9M3
# Investigating the impact of warmup on n_eff
# Will use m9.1, with the default warmup it gets ~2.5k n_eff
precis(m_9m1, depth=2)

# Default settings are 2k samples, and 1k warmup
# Let's try 100, 200, 500, 1500
warmups <- list("10"=10, "100"=100, "200"=200, "500"=500)
mods <- lapply(warmups, function(w) {
  ulam(
    alist(
        log_gdp_std ~ dnorm(mu, sigma),
        mu <- a[cid] + b[cid]*(rugged_std - 0.215),
        a[cid] ~ dnorm(1, 0.1),
        b[cid] ~ dnorm(0, 0.3),
        sigma ~ dunif(0, 1)
    ),
    data=dat_slim, chains=4, cores=4, warmup = w
    )  
})
mods$`1000` <- m9.1

# Interestingly it's not necessarily linear...
# So we actually get the best ESS_bulk from using 200 warmup samples
# Even 100 seems to often outperform 1,000!
# What does this imply for using Stan?!
# Also ofc 10 samples is stupidly bad
map_dfr(mods, function(x) {
    foo <- precis(x, depth=2)
    as_tibble(foo) |> 
        mutate(param=row.names(foo))
}, .id="warmup") |>
    mutate(warmup = as.factor(as.numeric(warmup))) |>
    filter(warmup != 10) |>
    ggplot(aes(x=param, y=ess_bulk, colour=warmup)) +
        geom_point() +
        theme_minimal()

# Plotting the relationship directly
# Shows rougly that for b[1] and sigma n_eff increases with warmup
# But for a[1], a[2], b[2], it gets worse
# I wouldn't interpret this causally, but it's far less straightforward than I'd expect
map_dfr(mods, function(x) {
    foo <- precis(x, depth=2)
    as_tibble(foo) |> 
        mutate(param=row.names(foo))
}, .id="warmup") |>
    mutate(warmup = as.numeric(warmup)) |>
    filter(warmup != 10) |>
    ggplot(aes(x=warmup, y=ess_bulk, colour=param)) +
        stat_smooth(method="lm", se = FALSE) +
        geom_point() +
        theme_minimal()

# 9H1
# Cauchy is thicker tailed than normal
# It's the resulting distribution from the ratio of 2 normally distributed rvs with mean=0
# Looks from the code like it's just sampling from the prior with no likelihood to update
mp <- ulam(
    alist(
        a ~ dnorm(0, 1),
        b ~ dcauchy(0, 1)
    ),
    data=list(y=1), chains=1
)

# Now have actually fit it can inspect posterior
post <- extract.samples(mp)
# Yep see cauchy's much thicker tails and lower peak
tibble(i=rep(1:500, 2), value=c(post$a, post$b), 
       dist=rep(c('normal', 'cauchy'), each=500)) |>
    ggplot(aes(x=value, colour=dist)) +
        geom_density() +
        theme_minimal()

# The traceplot for the cauchy looks quite different to the normal
# But mostly due to the thicker tails!
# I don't follow the hint about the Cauchy having no expected value, as it looks like 0 here
traceplot(mp)

# 9H2
data("WaffleDivorce")
d <- WaffleDivorce

# standardize
d$A <- scale(d$MedianAgeMarriage)
d$D <- scale(d$Divorce)
d$M <- scale(d$Marriage)

# D ~ A
m5.1 <- ulam(
    alist(
        D ~ dnorm(mu, sigma),
        mu <- a + bA * A,
        a ~ dnorm(0, 0.2),
        bA ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), 
    data=list(D=d$D, A=d$A), chains=4, cores=4, log_lik = TRUE
)

# D ~ M
m5.2 <- ulam(
    alist(
        D ~ dnorm(mu, sigma),
        mu <- a + bM * M,
        a ~ dnorm(0, 0.2),
        bM ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ),
    data=list(D=d$D, M=d$M), chains=4, cores=4, log_lik = TRUE
)

# D ~ M + A
m5.3 <- ulam(
    alist(
        D ~ dnorm(mu, sigma),
        mu <- a + bM*M + bA*A,
        a ~ dnorm(0, 0.2),
        bM ~ dnorm(0, 0.5),
        bA ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ),
    data=list(D=d$D, M=d$M, A=d$A), chains=4, cores=4, log_lik = TRUE
)
# Here the most predictive model is the one with just age 
# It's not surprising that marriage alone doesn't add much (m5.3)
# But even though M is a confounder would still think it'd help predictions
# But evidently not enough to overcome having the additional parameter
compare(m5.1, m5.2, m5.3)

# 9H3
# Leg length example
N <- 100
set.seed(909)
height <- rnorm(N, 10, 2)
leg_prop <- runif(N, 0.4, 0.5) # leg as proportion of height
leg_left <- leg_prop*height + rnorm(N, 0, 0.02) # add error so legs aren't same size exactly
leg_right <- leg_prop*height + rnorm(N, 0, 0.02)
d <- tibble(height, leg_left, leg_right)

# Before had highly correlated parameters for the 2 leg length slopes, their
# sum was the simulated slope (2.2), but were very anti-correlated
m_9h3_1 <- ulam(
    alist(
        height ~ dnorm(mu, sigma),
        mu <- a + bl*leg_left + br*leg_right,
        a ~ dnorm(10, 100),
        bl ~ dnorm(2, 10),
        br ~ dnorm(2, 10),
        sigma ~ dexp(1)
    ),
    data=d, chains=4, cores=4,
    start=list(a=10, bl=0, br=0.1, sigma=1),
    log_lik=TRUE
)
# Here because we start bl < br, it helps make br higher
# I think the starting points are just to help the model fit
# Although even then we got max treedepth warnings (4%) and n_eff is low
precis(m_9h3_1)

# If we don't provide the starting point we get 56% hitting
# max treedepth!
m_9h3_nostart <- ulam(
    alist(
        height ~ dnorm(mu, sigma),
        mu <- a + bl*leg_left + br*leg_right,
        a ~ dnorm(10, 100),
        bl ~ dnorm(2, 10),
        br ~ dnorm(2, 10),
        sigma ~ dexp(1)
    ),
    data=d, chains=4, cores=4
)
# And similar n_eff
precis(m_9h3_nostart)

# Anyway back to the question. Now are going to do the same
# but with a constraint on one leg to be positive
# Which shouldn't make a huge difference normally since the leg effect should always be positive
# But because it's anti-correlated with the other leg we might find strange behaviour
# NB: internally this positive constraint is handled with a log transform
# I don't think it enforces a linear-log relationship, but I'd like to look into that
# for my own understanding
m_9h3_2 <- ulam(
    alist(
        height ~ dnorm(mu, sigma),
        mu <- a + bl*leg_left + br*leg_right,
        a ~ dnorm(10, 100),
        bl ~ dnorm(2, 10),
        br ~ dnorm(2, 10),
        sigma ~ dexp(1)
    ),
    data=d, chains=4, cores=4,
    constraints=list(br="lower=0"),
    start=list(a=10, bl=0, br=0.1, sigma=1),
    log_lik=TRUE
)
precis(m_9h3_2)

# Annoyingly Stan won't let me extract samples from the posterior for either model
# Saying that '$ operator not defined for this S4 class'
# Which makes me think there's been a cmdstanr update that isn't reflected in rethinking
#extract.samples(m_9h3_1)  # Errors!
#extract.samples(m_9h3_2)  # Errors!

# The precis function still works so will need to use that for inference
# With no constraints, both leg effects are in the range [-4, 6], and assumedly
# are very anti-correlated (i.e. if one is 6, the other is -4)
precis(m_9h3_1)

# The positive constraint on br has given it a minimum of 0.37,
# which now means the maximum of bl is 1.60 (as 1.60 + 0.37 ~= 2.2)
# So yes positively constraining br is implicitly constraining bl
precis(m_9h3_2)

# 9H4
# Doesn't work because of the same S4 error...
compare(m_9h3_1, m_9h3_2)

# 9H5
metropolis_2 <- function(n_samples, starting, populations) {
    positions <- rep(0, n_samples)
    current <- starting
    for (i in 1:n_samples) {
        positions[i] <- current
        proposal <- current + sample(c(-1, 1), size=1)
        proposal <- min(proposal, length(populations))
        proposal <- max(proposal, 1)
        # THIS LINE IS THE ONLY CHANGE!
        # Rather than take ratio of proposed island number / current island number,
        # explicitly use their populations
        prob_move <- populations[proposal]/populations[current]
        current <- ifelse(runif(1) < prob_move, proposal, current)
    }
    positions
}

# Firstly test with the default assumption that populations are correlated with island number
results <- metropolis_2(1e5, 3, seq(1, 10)*100)
# Seemingly random sampling
plot(1:100, results[1:100])
# But are clearly sampling in proportion to island size as desired
plot(table(results))

# Now try again with random populations
pops <- rnorm(10, 50000, 1000)
results2 <- metropolis_2(1e5, 3, pops)
# Again seemingly random sampling
plot(1:100, results2[1:100])
cts <- table(results2)
# And sure enough the number of visits looks proportional to their populations
tibble(island = 1:10, n_visits=as.numeric(cts), population=pops) |>
    ggplot(aes(x=population, y=n_visits)) +
        geom_point() +
        stat_smooth(method="lm") +
        theme_minimal()

# 9H6
# Globe tossing data and model using Metropolis
# The islands are possible parameters of p
# So at each step get the likelihood of the current point
metropolis_globe <- function(n_samples, W, L, grid_size, prior, starting) {
    p_grid <- seq(0, 1, length.out=grid_size)
    positions <- rep(0, n_samples)
    current <- starting
    for (i in 1:n_samples) {
        positions[i] <- current
        proposal <- current + sample(c(-1, 1), size=1)
        proposal <- min(proposal, grid_size)
        proposal <- max(proposal, 1)
        # TODO ideally would use prior
        likelihood_proposal <- dbinom(W, size=W+L, prob=p_grid[proposal])
        likelihood_current <- dbinom(W, size=W+L, prob=p_grid[current])
        prob_move <- likelihood_proposal / likelihood_current
        current <- ifelse(runif(1) < prob_move, proposal, current)
    }
    positions
}

post_globe <- metropolis_globe(1e5, 5, 2, 100, function(x) 1, 1)
# That actually looks quite reasonable!
tibble(p=seq(0, 1, length.out=100), i=1:100) |>
    inner_join(tibble(i=post_globe), by="i") |>
    count(p) |>
    mutate(prop = n / sum(n)) |>
    ggplot(aes(x=p, y=prop)) +
        geom_col(colour="black", fill="white") +
        theme_minimal() +
        labs(x="P", y="N samples", title="Sampling p from W5 2L using Metropolis")

# Now compare it to what I get using grid search
grid_approx_globe <- function(grid_size, prior, W, L) {
    # R code 2.3
    p_grid <- seq(from=0, to=1, length.out=grid_size)
    likelihood <- dbinom(W, size=W+L, prob=p_grid)
    unstd.posterior <- likelihood * prior(p_grid)
    unstd.posterior / sum(unstd.posterior)
}
res_grid <- grid_approx_globe(100, function(x) rep(1, length(x)), 5, 2)
dens_grid <- tibble(
    p = seq(0, 1, length.out=100),
    prob = res_grid
) |>
  mutate(algorithm="Grid Approximation") 
    
dens_metropolis <- tibble(p=seq(0, 1, length.out=100), i=1:100) |>
    inner_join(tibble(i=post_globe), by="i") |>
    count(p) |>
    mutate(prob = n / sum(n)) |>
    select(-n) |>
    mutate(algorithm="Metropolis")

# Not bad approximation from Metropolis!
dens_grid |> 
    rbind(dens_metropolis) |>
    ggplot(aes(x=p, y=prob)) +
        geom_col(colour="black", fill="white") +
        facet_wrap(~algorithm) +
        theme_minimal() +
        labs(x="P", y="Posterior probability", title="Sampling p from W5 2L")

# 9H7
# HMC for globe tossing
# Firstly need function to calculate neg-log-prob at data
Water <- 5
Land <- 2
U <- function(q) {
    # NB: Seems to use y and x from global scope rather than function args
    p <- q[1]
    # constrain!
    p <- min(p, 1)
    p <- max(p, 0)
    # Return negative log-prob
    -1 * dbinom(Water, Water+Land, prob=p, log=TRUE)
}

# Need the gradient of binomial likelihood
# https://math.stackexchange.com/q/651976
U_gradient <- function(q) {
    p <- q[1]
    log_prob <- (Water/p) - Land/(1-p)
    -1 * log_prob
}

# And here's the code to run it
Q <- list()
Q$q <- c(0.5)  # Starting values of p
pr <- 0.3
step <- 0.03
L <- 11  # Working example, set to 28 to see U-turns
n_samples <- 1e4
# NB: Had to give up with 1e5 samples as was taking too long
results <- tibble(
    iteration=numeric(n_samples),
    p=numeric(n_samples),
    accepted=integer(n_samples)
)
for (i in 1:n_samples) {
    Q <- HMC2(U, U_gradient, step, L, Q$q)
    results$iteration[i] <- i
    results$p[i] <- Q$traj[L+1, 1]
    results$accepted[i] <- Q$accept
}

# Accepted the vast majority of samples!
results |>
    count(accepted)

# Now have a continuous parameter rather than discrete so can't easily
# join in on the other algorithms, would need to calculate kernel density first
# Although annoyingly it looks although this hasn't worked - it's still sampling
# around the starting value
raw_density <- density(results |> filter(accepted==1) |> pull(p))
dens_hmc <- tibble(
    p = raw_density$x,
    prob = raw_density$y
) |>
    mutate(algorithm="HMC") |>
    mutate(prob = prob / sum(prob))

# Shape isn't quite right, far less certain around the estimate, even if mode is in the right place
dens_grid |> 
    rbind(dens_metropolis) |>
    rbind(dens_hmc) |>
    ggplot(aes(x=p, y=prob)) +
        geom_col(colour="black", fill="white") +
        facet_wrap(~algorithm) +
        theme_minimal() +
        labs(x="P", y="Posterior probability", title="Sampling p from W5 2L")

# Just try once more with larger L
Q <- list()
Q$q <- c(0.5)  # Starting values of p
pr <- 0.3
step <- 0.03
L <- 25  # Working example, set to 28 to see U-turns
n_samples <- 1e4
results <- tibble(
    iteration=numeric(n_samples),
    p=numeric(n_samples),
    accepted=integer(n_samples)
)
for (i in 1:n_samples) {
    Q <- HMC2(U, U_gradient, step, L, Q$q)
    results$iteration[i] <- i
    results$p[i] <- Q$traj[L+1, 1]
    results$accepted[i] <- Q$accept
}

raw_density <- density(results |> filter(accepted==1) |> pull(p))
dens_hmc <- tibble(
    p = raw_density$x,
    prob = raw_density$y
) |>
    mutate(algorithm="HMC") |>
    mutate(prob = prob / sum(prob))

# Not a real difference!
# I wonder if the problem is because I only ran it for 1e4 samples rather than 1e5?
dens_grid |> 
    rbind(dens_metropolis) |>
    rbind(dens_hmc) |>
    ggplot(aes(x=p, y=prob)) +
        geom_col(colour="black", fill="white") +
        facet_wrap(~algorithm) +
        theme_minimal() +
        labs(x="P", y="Posterior probability", title="Sampling p from W5 2L")


# Things to ask -----------------------------------------------------------
#   - General rule of thumb for scaling vars? i.e. when to Z-transform, when to divide / max, 
#     when to divide / mean, when to just center?
#   - When using WAIC to compare, should focus on weight or dWAIC and dSE? The former always
#     seems to have much stronger differences than the latter
