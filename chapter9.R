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
