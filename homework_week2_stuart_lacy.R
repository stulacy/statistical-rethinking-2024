library(rethinking)
library(dagitty)
data("Howell1")
d <- Howell1[ Howell1$age < 13, ]
summary(d)



# Question 1 --------------------------------------------------------------
# Estimate causal association between age & weight given:
# A -> H, H->W, A->W
dag_q1 <- dagitty("dag{ W <- A -> H -> W }")
coordinates(dag_q1) <- list(x=c(A=0, H=1, W=1), y=c(A=0, H=0, W=1))
drawdag(dag_q1)

# Can find the direct causal effect without needing to condition on anything as there's no backdoor paths
# into A
adjustmentSets(dag_q1, exposure="A", outcome="W")

# Will standardize variables before modelling
d$A <- scale(d$age)
d$H <- scale(d$height)
d$W <- scale(d$weight)
# Age has a strong effect, with 1stdev of change associated with ~0.9 stdevs difference in weight
quap(
    alist(
        W ~ dnorm(mu, sigma),
        mu <- a + bA * A,
        a ~ dnorm(0, 0.2),
        c(bA) ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), 
    data=d
) |> precis() |> plot()

# To estimate the total causal effect of this dag need 2 models: 1 for each path
m2 <- quap(
    alist(
        # Both A->W and H->W
        W ~ dnorm(mu_2, sigma_1),
        mu_2 <- a_2 + bA_2 * A + bH * H,
        a_2 ~ dnorm(0, 0.2),
        bA_2 ~ dnorm(0, 0.5),
        bH ~ dnorm(0, 0.5),
        sigma_1 ~ dexp(1),
        # Direct A->H
        H ~ dnorm(mu_1, sigma_2),
        mu_1 <- a_1 + bA_1 * A,
        a_1 ~ dnorm(0, 0.2),
        bA_1 ~ dnorm(0, 0.5),
        sigma_2 ~ dexp(1)
    ), 
    data=d
)
precis(m2)

# Now to simulate from this we need a range of ages
# Age ranges from -1.4 to 1.7
# Will use these limits as -1.4 corresponds to an age of 0 years and 1.7 corresponds to 12 years
summary(d$A)
A_seq <- seq(min(d$A), max(d$A), length.out=50)
sim_dat <- data.frame(A = A_seq)

# Ensure we simulate H first
s <- sim(m2, data=sim_dat, vars=c("H", "W"))

# Marking -----------------------------------------------------------------
# I misunderstood the question. Rather than forming a full generative model
# of the DAG, it was just looking for a simulation that represents these
# relationships.
# i.e.
N <- 100
age <- runif(N, 1, 12)
# Coefficients don't matter, the order of execution is the focus here
height <- rnorm(N, age * 0.7, 3) 
weight <- rnorm(N, age * 0.3 + height * 0.6, 2)
# Age vs height
plot(age, height)
# Age has direct effect on weight
plot(age, weight)
# Height has stronger direct effect on weight
plot(height, weight)


# Question 2 --------------------------------------------------------------
# Estimate the total causal effect of each year of growth on weight
# Using the above simulation

# Plot counter-factual
# Quite a strong line!
plot(sim_dat$A, colMeans(s$W), ylim=c(-2, 2), type='l',
     xlab="manipulated A", ylab="counterfactual W")
shade(apply(s$W, 2, PI), sim_dat$A)
mtext("Total counterfactual effect of A on W")
# Average gradient of 0.89 stdevs in weight per stdev of age
lm(colMeans(s$W) ~ sim_dat$A)
# Back into normal units this is ~9 years is estimated to correspond to ~20kg of weight increase
1 * sd(d$age) + mean(d$age)
0.89 * sd(d$weight) + mean(d$weight)

# NB: this is NOT the same as adding the 2 age coefficients from the model
precis(m2)
# 0.07 + 0.91

# Marking -----------------------------------------------------------------
# I misunderstood the question for the same reason as in Chapter 11's 
# Problems, that I misunderstood how to estimate a total causal effect,
# thinking that a full joint model and looking at counterfactuals
# was the correct way. Whereas in fact a single model of W ~ A is needed
adjustmentSets(dag_q1, "A", "W", effect="total")
# I'll use the same priors as I had for my original answer
# Will be interesting to see if the same result is reached
m2_corrected <- quap(
    alist(
        W ~ dnorm(mu, sigma),
        mu <- a + bA * A,
        a ~ dnorm(0, 0.2),
        bA ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), 
    data=d
)
# This gives an effect of 0.89 per SD of age on SD of weight
# Which is actually what my counter-factual method identified, but
# obviously doing it directly is far more straight forward.
# The counter-factual method is more flexible however and allows us to ask,
# well, counterfactuals.
precis(m2_corrected)

# Converting it back into original measurement units to compare with the
# given answer gives 9.22 years is associated with 19.8kg
# Or 1 year
age_sd <- 1 * sd(d$age) + mean(d$age)
weight_change_per_sd <- 0.89 * sd(d$weight) + mean(d$weight)
age_sd
weight_change_per_sd
# Or ~ 2.14kg per year, which is actually higher than the 1.38kg/year
# that the solutions has, although they didn't use any standardization
# which I don't think should change the effect size this much
weight_change_per_sd / age_sd
# I tried again with a wider prior and it didn't change the age coefficient

# Question 3 --------------------------------------------------------------
data("Oxboys")
d <- Oxboys
nrow(d)  # 234 observations
length(unique(d$Subject)) # 26 boys
length(unique(d$Occasion)) # 9 occasions
# Age looks standardized-ish (mean != 0?!)
# Height is in cm
summary(d)

# Want to model the rate of change of height with each occasion
# Will want hierarchical intercept on each boy (their initial height)
# And a common growth factor
# I think the question wants a separate growth factor for each occasion,
# i.e. assuming non-linear growth
# But I'll try with a single factor first in 2 ways:

# Attempt 1: Hierarchical intercept:
# Give each boy their own intercept which is their initial height
d$Occasion_0 <- d$Occasion - 1 # Zero-index Occasion
set.seed(17)  # Failed to converge a few times, I guess because have so many intercepts
m3_1 <- quap(
    alist(
        height ~ dnorm(mu, sigma),
        mu <- a[Subject] + bO * Occasion_0,
        a[Subject] ~ dnorm(140, 20),
        bO ~ dlnorm(0, 1),  # Because need positive constraint
        sigma ~ dexp(1)
    ), data=d
)
# Get an estimate of 1.63 cm per occasion!
precis(m3_1, depth=2)

# Attempt 2: centered data
library(tidyverse)
# Center height so that it is relative to the first measurement
d <- d |>
    group_by(Subject) |> 
    mutate(height_centered = height - min(height)) |>
    ungroup() 
d

set.seed(17)
m3_2 <- quap(
    alist(
        height_centered ~ dnorm(mu, sigma),
        mu <- bO * Occasion_0,
        a ~ dnorm(0, 0.2),  # Expect a = 0, but putting intercept in just in case (it will converge without it)
        bO ~ dlnorm(0, 1),  # Because need positive constraint
        sigma ~ dexp(1)
    ), data=d |> filter(Occasion_0 > 0) # Don't want the first row since this is just height = 0
)
# Get the same result: 1.61 cm per occasion
# Slightly higher sigma
precis(m3_2)

# Attempt 3: Separate growth parameters
# I think the question is asking for 
# height_delta_i = beta_i
# Where i is the occasion number
# So although beta is a rate factor because the objective is a rate, it's modelled as an intercept here
# Need a vector of height differences, easy!
d <- d |>
    arrange(Subject, Occasion) |>
    group_by(Subject) |>
    mutate(height_diff = height - lag(height)) |>
    ungroup()

set.seed(3) # Had to do more tuning seed until found a fit - not a great sign for the model
m3_3 <- quap(
    alist(
        height_diff ~ dnorm(mu, sigma),
        mu <- a[Occasion_0],
        a[Occasion_0] ~ dlnorm(0, 1),  # Because need positive constraint
        sigma ~ dexp(1)
    ), data=d |> filter(Occasion_0 > 0)  # Likewise don't want first row as don't have a delta
)
# The estimates do seem to vary, i.e. the 4th (5th in 1-index) Occasion seems to be associated with a 
# much slower growth than the average ~1.6, while the 6th (7th in 1-index) Occasion seems to be
# much higher, with 89% PI over 2
precis(m3_3, depth=2) |> plot()

# Compute the posterior of the total growth over all 9 occasions
post <- extract.samples(m3_3)
# Have 1 row per sample and 8 columns representing 8 growth points
dim(post$a)
# So take row sum to get 1000 samples of total growth
total_growth <- rowSums(post$a)
dens(total_growth)
# Mean growth of 13cm
mean(total_growth)
# Median is the same - I guess because quap assumes Gaussian posterior?
median(total_growth)
# Likewise the MAP is 13.05
chainmode(total_growth)
# Range of 11.3 - 14.4
range(total_growth)

# Marking -----------------------------------------------------------------
# I went into way more depth than the question asked.
# They didn't want specific growth-factors per occasions,
# or to model each subject with a different intercept
# It just wanted the 'pooled' average inter-occasion growth
# i.e.
m4 <- quap(
    alist(
        height_diff ~ dlnorm(mu, sigma),
        mu <- alpha,
        alpha ~ dnorm(0, 0.1),
        sigma ~ dexp(3)
    ), data=d |> filter(!is.na(height_diff))
)
# So the average height growth is 0.29 in log-normal
precis(m4)

# Looking at the exponent of this to get back to the measurement scale
dens(sim(m4))
# the median is 1.33cm per occasion
median(sim(m4))
# The mean is 1.62, which matches up with my models 1 and 2 (which were indirectly asking each occasions contribution to
# total height)
mean(sim(m4))

# Can simulate the total growth over 8 occasions
# By simulating 8 sums
post <- extract.samples(m4)
total <- sapply(1:nrow(post), function(i) {
    sum(rlnorm(8, post$alpha[i], post$sigma[i]))
})
dens(total)

# Median growth of 12.6
median(total)
# Mean growth of 13cm, which is the same as my 
# model that allowed each occasion to have a differing rate
mean(total)

# NB: My models had the same median and mean as my growth factor
# despite also having the log-normal prior to enforce positivity
# How come I got a Gaussian posterior but the solutions model
# has log-normal? I believe because quap assumes every posterior
# is Gaussian, including that of my growth rate (with log-normal)
# prior. HOWEVER, the solutions model has the growth-rate as an
# OUTCOME, so it strictly doesn't have a posterior, just a
# POSTERIOR PREDICTIVE which we run outside of quap and so is
# log-normal
