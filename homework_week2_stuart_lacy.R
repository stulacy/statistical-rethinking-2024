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