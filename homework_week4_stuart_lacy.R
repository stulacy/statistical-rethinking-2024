library(rethinking)
library(dagitty)
library(tidyverse)

# Question 1 --------------------------------------------------------------
# 1. Fit models m6.9 and m6.10 and compare using PSIS and WAIC
# Which is expected to make better predictions, and which better causal inference?
d <- sim_happiness(seed=1977, N_years=1000)
d2 <- d[d$age > 17, ]
d2$A <- (d2$age - 18) / (65 - 18)
d2$mid <- d2$married + 1

m6.9 <- ulam(
    alist(
        happiness ~ dnorm(mu, sigma),
        mu <- a[mid] + bA*A,
        a[mid] ~ dnorm(0, 1),
        bA ~ dnorm(0, 2),
        sigma ~ dexp(1)
    ), data=d2, chains=4, cores=4, log_lik=TRUE
)

m6.10 <- ulam(
    alist(
        happiness ~ dnorm(mu, sigma),
        mu <- a + bA*A,
        a ~ dnorm(0, 1),
        bA ~ dnorm(0, 2),
        sigma ~ dexp(1)
    ), data=d2, chains=4, cores=4, log_lik=TRUE
)

# m6.9 (conditioning on marriage) is far more predictive of happiness than
# just using age! Huge dWAIC (388) with small dSE (35) and has weight=1
compare(m6.9, m6.10)

# HOWEVER! Just because it is more predictive doesn't mean that is provides
# us with causal inference. And here marriage is a collider so it opens an
# association between age and happiness which doesn't exist
# So if we want to know the direct causal effect of age on happiness
# we must only look at m6.10, where it is correctly 0.
# If we used the more predictive model where age has a slope of -0.76
# we would infer that age is negatively associated with happiness which 
# isn't true.
coeftab(m6.9, m6.10)

# Question 2 --------------------------------------------------------------
# Urban fox dataset
# Which combination of variables best predicts body weight?
# What causal interpretation can we assign each coefficient from the best 
# scoring model?
data(foxes)
d <- foxes
d$F <- scale(d$avgfood)
d$A <- scale(d$area)
d$G <- scale(d$groupsize)
d$W <- scale(d$weight)

# This is the model remember
# So will fit a regression on W with every permutation of A, F, G 
# Will use quap for speed
dag_1 <- dagitty("dag {
    A -> F -> G -> W <- F
                   }")
coordinates(dag_1) <- list(x=c(F=0, A=1, W=1, G=2),
                               y=c(F=1, A=0, W=2, G=1))
drawdag(dag_1)

m2_1 <- quap(
    alist(
        W ~ dnorm(mu, sigma),
        mu <- a + bF * F + bG * G + bA * A,
        a ~ dnorm(0, 0.2),
        c(bF, bG, bA) ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), 
    data=d
)

m2_2 <- quap(
    alist(
        W ~ dnorm(mu, sigma),
        mu <- a + bF * F + bG * G,
        a ~ dnorm(0, 0.2),
        c(bF, bG) ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), 
    data=d
)

m2_3 <- quap(
    alist(
        W ~ dnorm(mu, sigma),
        mu <- a + bG * G + bA * A,
        a ~ dnorm(0, 0.2),
        c(bG, bA) ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), 
    data=d
)

m2_4 <- quap(
    alist(
        W ~ dnorm(mu, sigma),
        mu <- a + bF * F + bA * A,
        a ~ dnorm(0, 0.2),
        c(bF, bA) ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), 
    data=d
)

m2_5 <- quap(
    alist(
        W ~ dnorm(mu, sigma),
        mu <- a + bF * F,
        a ~ dnorm(0, 0.2),
        c(bF) ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), 
    data=d
)

m2_6 <- quap(
    alist(
        W ~ dnorm(mu, sigma),
        mu <- a + bA * A,
        a ~ dnorm(0, 0.2),
        c(bA) ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), 
    data=d
)

m2_7 <- quap(
    alist(
        W ~ dnorm(mu, sigma),
        mu <- a + bG * G,
        a ~ dnorm(0, 0.2),
        c(bG) ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), 
    data=d
)

#   1 F + G + A
#   2 F + G
#   3 G + A
#   4 F + A
#   5 F
#   6 A
#   7 G
# All the models with G in (1, 2, 3, 7) are the top 4,
# although the model with just G isn't quite as predictive, 
# it fares much better with either A or F (and remember A causes F).
# The the most predictive model is the one with F, G, and A.
# This leads to the interpretation that G has the biggest causal effect on W
# although it's mediated through F (and indirectly through A), so knowing these
# adds additional predictive ability once G is known.
# F and A on their own is actually the least predictive model.
compare(m2_1, m2_2, m2_3, m2_4, m2_5, m2_6, m2_7)

# Question 3 --------------------------------------------------------------
data("Dinosaurs")
d <- Dinosaurs
# Have 32 dinosaurs with 4 columns:
#  - age
#  - mass
#  - species (categorical)
#  - species_id (looks like just a numeric index for species)
nrow(d)
summary(d)

# Yes sp_id is an index variable for species
table(d$species, d$sp_id)

# Model growth as function of age for at least 1 species
# Makes sense to use varying-slopes & intercepts (different species 
# both have different average weight and different growth rates)
# Will try a linear model and at least one other non-linear version

# Firstly data scaling
# Want to have 0 as reference here, so will use 0-1 transform
dens(d$age) 
dens(d$age / max(d$age))
# But mass is massively skewed
dens(d$mass) 
# Log transform makes it better
dens(log(d$mass))
# Log transform with / max gives us outputs in 0-1.... nearly!
dens(log(d$mass) / max(log(d$mass)))
# There is one point < 0 so will use / max to get output in 0-1 and hope the skewed nature is ok
dens(d$mass / max(d$mass))
d$A <- d$age / max(d$age)
d$M <- d$mass / max(d$mass)

# Now fit first linear model with varying intercepts and slopes
# Prior for intercept is that for average age, what is the expected mass?
# Since this is z-score of log, should be in -2, 2, centered around 0
# So N(0, 1) should be ok, although perhaps too wide given the lack of data
# For slopes, again centered on 0 with max rise/run being 4/4 = 1, so N(0, 0.5) should be ok
# NB: If get massive posterior PIs on the coefficients might want to use multi-level model
m_3_linear <- ulam(
    alist(
        M ~ dnorm(mu, sigma),
        mu <- a[sp_id] + b[sp_id] * A,
        a[sp_id] ~ dnorm(0, 1),
        b[sp_id] ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), data=d, chains=4, cores=4, log_lik=TRUE
)

# Species 1 looks like an outlier! Has the lowest average mass but highest rate of growth
plot(precis(m_3_linear, depth=2))

# Now plot posterior predictions to see functional form
ndata <- expand_grid(
    A = seq(0, 1, length.out=100),
    sp_id=d$sp_id
)

post <- link(m_3_linear, data=ndata)
mu <- colMeans(post)
mu_PI <- apply(post, 2, PI)
sim_PI <- apply(sim(m_3_linear, data=ndata), 2, PI)

# The main issue here is that the masses are so different per species, we ideally should standardize on a per-species basis
# As otherwise it's hard to tell how well the model fit is
ndata |>
    mutate(mu=mu,
           mu_lower = mu_PI[1, ],
           mu_upper = mu_PI[2, ],
           sim_lower = sim_PI[1, ],
           sim_upper = sim_PI[2, ],
           ) |>
    inner_join(d |> select(sp_id, species), by="sp_id", relationship="many-to-many") |>
    ggplot(aes(x=A, y=mu, colour=species, fill=species)) +
        #geom_ribbon(aes(ymin=sim_lower, ymax=sim_upper), alpha=0.3) +   # Had to remove as the simulated PIs don't look sensible, it looks like they're in the wrong order
        geom_ribbon(aes(ymin=mu_lower, ymax=mu_upper), alpha=0.5) +
        geom_line() +
        geom_point(aes(y=M),data=d) +
        facet_wrap(~species) +
        guides(colour="none", fill="none") +
        theme_bw()

# So now normalise by species for mass. Age should be ok...
d <- d |> 
    group_by(species) |>
    mutate(M = mass / max(mass)) |>
    ungroup()
    
# Now refit
m_3_linear_2 <- ulam(
    alist(
        M ~ dnorm(mu, sigma),
        mu <- a[sp_id] + b[sp_id] * A,
        a[sp_id] ~ dnorm(0, 1),
        b[sp_id] ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), data=d, chains=4, cores=4, log_lik=TRUE
)

# That looks more realistic, different intercepts and different growth rates
plot(precis(m_3_linear_2, depth=2))

# Now plot posterior predictions to see functional form
ndata <- expand_grid(
    A = seq(0, 1, length.out=100),
    sp_id=d$sp_id
)

post <- link(m_3_linear_2, data=ndata)
mu <- colMeans(post)
mu_PI <- apply(post, 2, PI)
#sim_PI <- apply(sim(m_3_linear, data=ndata), 2, PI)

# The scaling is better now, but some issues, namely it predicts negative mass quite often, noticeably for the Massoponylus dinosaur
ndata |>
    mutate(mu=mu,
           mu_lower = mu_PI[1, ],
           mu_upper = mu_PI[2, ],
           sim_lower = sim_PI[1, ],
           sim_upper = sim_PI[2, ],
           ) |>
    inner_join(d |> select(sp_id, species), by="sp_id", relationship="many-to-many") |>
    ggplot(aes(x=A, y=mu, colour=species, fill=species)) +
        #geom_ribbon(aes(ymin=sim_lower, ymax=sim_upper), alpha=0.3) +   # Had to remove as the simulated PIs don't look sensible, it looks like they're in the wrong order
        geom_ribbon(aes(ymin=mu_lower, ymax=mu_upper), alpha=0.5) +
        geom_line() +
        geom_hline(yintercept=0, linetype="dashed") +
        geom_point(aes(y=M),data=d) +
        facet_wrap(~species) +
        guides(colour="none", fill="none") +
        theme_bw()

# So now will try a proper logistic growth model
# Logistic growth model is dN / dt = r * N(K -N) / K
# N = K /(1 + A * exp(-Kt))
# Where K is carrying capacity and r is rate of growth, N is current value
# As shown here
# https://web.ma.utexas.edu/users/m408s/CurrentWeb/LM9-4-3.php

# Use uninformative priors as I have no idea what values to expect here!
# Ah this gives some bad diagnostics, max tree depth and divergences...
m_3_logistic <- ulam(
    alist(
        M ~ dnorm(mu, sigma),
        mu <- K / (1 + B * exp(-r*A)),
        K ~ dnorm(1, 0.1), # K is the maximum mass so would expect around 1
        B ~ dnorm(0, 0.5),
        r ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), data=d |> filter(sp_id == 3) |> select(M, A), chains=4, cores=4, log_lik=TRUE
)

post <- link(m_3_logistic, data=ndata |> filter(sp_id == 3))
mu <- colMeans(post)
mu_PI <- apply(post, 2, PI)

# Something seems to have gone wrong here... this should produce a smooth line, not something that weirdly varies at high values
ndata |>
    filter(sp_id == 3) |>
    mutate(mu=mu,
           mu_lower = mu_PI[1, ],
           mu_upper = mu_PI[2, ],
           ) |>
    inner_join(d |> select(sp_id, species), by="sp_id", relationship="many-to-many") |>
    ggplot(aes(x=A, y=mu, colour=species, fill=species)) +
        geom_ribbon(aes(ymin=mu_lower, ymax=mu_upper), alpha=0.5) +
        geom_line() +
        geom_point(aes(y=M),data=d |> filter(sp_id == 3)) +
        facet_wrap(~species) +
        guides(colour="none", fill="none") +
        theme_bw()

# Just plot for one set of coefs
# Yep these are lines, not curves!
coefs <- extract.samples(m_3_logistic)
A_seq <- seq(0, 1, length.out=100)
foo <- map_dfr(1:20, function(i) {
    y <- coefs$K[i] / (1 + coefs$B[i] * exp(-coefs$r[i]*A_seq))
    tibble(A=A_seq, M=y, i=i)
})
foo |>
    ggplot(aes(x=A, y=M, group=i)) +
        geom_line()
# Going to give up now as have run out of time!
# Will try and pick this up later
