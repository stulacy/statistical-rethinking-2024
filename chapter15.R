library(rethinking)
library(tidyverse)
library(dagitty)
library(patchwork)

# Monty Hall Problem ------------------------------------------------------
# But with toast!
sim_pancake <- function() {
    # These are the 3 pancakes with their 2 sides (BB, BU, UU)
    sides <- matrix(c(1,1,1, 0,0,0), 2, 3)
    # Randomly select a pancake
    pancake <- sides[, sample(1:3, 1)]
    # Randomly select the order in which the pancake appears (i.e. face up is first
    # face down is second)
    sample(pancake)
}

# This means our pancake is Unburnt faceup, and Burnt face down
set.seed(17)
sim_pancake()

# Simulated 1e4 pancakes
pancakes <- replicate(1e4, sim_pancake())
# Get the face up and downs
up <- pancakes[1, ]
down <- pancakes[2, ]

# Get proportion of all flips THAT ARE BURNT FACE UP (i.e. the condition) which are burnt face down
# Near 2/3s
sum(up == 1 & down == 1) / sum(up == 1)

# 15.1 Measurement Error --------------------------------------------------
data(WaffleDivorce)
d <- as_tibble(WaffleDivorce)
d

# Marriage rate and divorce rate are actually measured with measurement error!
plot(d$Divorce ~ d$MedianAgeMarriage, ylim=c(4, 15),
     xlab="Median age marriage", ylab="Divorce rate")

# Standard errors
# We already have standard errors in the dataset, so show +/- 1 SE
for (i in 1:nrow(d)) {
    ci <- d$Divorce[i] + c(-1, 1) * d$Divorce.SE[i]
    x <- d$MedianAgeMarriage[i]
    lines(c(x,x), ci)
}

# There is a negative exponential relationship between SE and Population size
# It's actually likely a 1/Sqrt(N) relationship, as this is the normalising factor in the SE
d |>
    ggplot(aes(x=Population, y=Divorce.SE)) +
        geom_point()

# Yep, 1/sqrt(N) 
d |>
    mutate(y = 1 / sqrt(Population)) |>
    ggplot(aes(x=Population, y=y)) +
        geom_point()


dlist <- list(
    D_obs = standardize(d$Divorce),
    D_sd = d$Divorce.SE / sd(d$Divorce),  # Interesting standardization... I guess to keep it on the
                                          # same scale as D_obs, but there's no centering needed
    M = standardize(d$Marriage),
    A = standardize(d$MedianAgeMarriage),
    N = nrow(d)
)

m15.1 <- ulam(
    alist(
        D_obs ~ dnorm(D_true, D_sd),
        vector[N]:D_true ~ dnorm(mu, sigma),
        mu <- a + bA*A + bM*M,
        a ~ dnorm(0, 0.2),
        bA ~ dnorm(0, 0.5),
        bM ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), data=dlist, chains=4, cores=4
)

# The bA slope has shrunk from its original estimate (before added measurement errors)
# Measurement errors were increasing the association between divorce and age at marriage
precis(m15.1, depth=2)

# Measurement error on both outcome and predictor
# same basic steps
dlist <- list(
    D_obs = standardize(d$Divorce),
    D_sd = d$Divorce.SE / sd(d$Divorce),  # Interesting standardization... I guess to keep it on the
    M_obs = standardize(d$Marriage),
    M_sd = d$Marriage.SE / sd(d$Marriage),  # Interesting standardization... I guess to keep it on the
                                          # same scale as D_obs, but there's no centering needed
    A = standardize(d$MedianAgeMarriage),
    N = nrow(d)
)

m15.2 <- ulam(
    alist(
        D_obs ~ dnorm(D_true, D_sd),
        M_obs ~ dnorm(M_true, M_sd),
        vector[N]:D_true ~ dnorm(mu, sigma),
        # What is this [i] doing here?! Won't compile without it
        mu <- a + bA*A + bM*M_true[i],
        # Need to provide a prior on true covariate! Can just put weakly informative 
        # This might be what I'm missing from my York regression analysis
        vector[N]:M_true ~ dnorm(0, 1),
        a ~ dnorm(0, 0.2),
        bA ~ dnorm(0, 0.5),
        bM ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), data=dlist, chains=4, cores=4
)

# No real difference on coefficients
precis(m15.1, pars=c("bA", "bM"))
precis(m15.2, pars=c("bA", "bM"))

post <- extract.samples(m15.2)
D_true <- apply(post$D_true, 2, mean)
M_true <- apply(post$M_true, 2, mean)

# Really nice visualisation of shrinkage!
# I should again try this for my York regression analysis
plot(dlist$M_obs, dlist$D_obs, pch=16, col=rangi2,
     xlab="Marriage rate (std)", ylab="Divorce rate (std)")
points(M_true, D_true)
for (i in 1:nrow(d)) {
    lines(c(dlist$M_obs[i], M_true[i]), c(dlist$D_obs[i], D_true[i]))
}

# 15.2 Missing Data -------------------------------------------------------
# Simulate missing homework
N <- 100
S <- rnorm(N)
H <- rbinom(N, size=10, inv_logit(S))

# If the (missing) outcome H is d-separated from the missing-causing exposure D,
# then can proceed with analysis
# Can condition on other variables as usual to achieve this

# Missing at random - no action taken
D <- rbern(N)
Hm <- H
Hm[D==1] <- NA

# Missing influenced by the predictor
D <- ifelse(S > 0, 1, 0)
Hm <- H
Hm[D==1] <- NA

# Missing influenced by an unobserved variable, which also influences the outcome
set.seed(501)
N <- 1e3
X <- rnorm(N) # Latent variable
S <- rnorm(N) # Exposure
H <- rbinom(N, size=10, inv_logit(2+S-2*X)) # Outcome is a function of both
D <- ifelse(X > 1, 1, 0) # Missingness is a function of the latent parameter
Hm <- H
Hm[D==1] <- NA

# What happens if we don't have missingness (but still have latent variable X influencing H)
dat_list <- list(
    H=H, S=S
)

m15.3 <- ulam(
    alist(
        H ~ binomial(10, p),
        logit(p) <- a + bS*S,
        a ~ normal(0, 1),
        bS ~ normal(0, 0.5)
    ), data=dat_list, chains=4, cores=4
)
# True coefficient of S should be 1, but we have X involved hence it's reduced
precis(m15.3)

# Now fit the model using the non-missing data only
dat_list0 <- list(H=H[D==0], S=S[D==0])
m15.4 <- ulam(
    alist(
        H ~ binomial(10, p),
        logit(p) <- a + bS*S,
        a ~ normal(0, 1),
        bS ~ normal(0, 0.5)
    ), data=dat_list0, chains=4, cores=4
)
# Actually get less biased result! Intercept is closer to 2 as well
precis(m15.4)

# Final situation, the outcome influences missingness
N <- 100
S <- rnorm(N)
H <- rbinom(N, size=10, inv_logit(S))
D <- ifelse(H < 5, 1, 0)
Hm <- H
Hm[D==1] <- NA

# No way we can block the causal path here

# 15.2.2 Imputing primates ------------------------------------------------
data(milk)
d <- milk |> as_tibble()
d <- d |>
    mutate(
        neocortex.prop = neocortex.perc / 100, # % to decimal
        logmass=log(mass),
        K=standardize(kcal.per.g),
        B=standardize(neocortex.prop),
        M=standardize(logmass)
    )

m15.5 <- ulam(
    alist(
        K ~ dnorm(mu, sigma),
        mu <- a + bB * B + bM*M,
        B ~ dnorm(nu, sigma_B),
        c(a, nu) ~ dnorm(0, 0.5),
        c(bB, bM) ~ dnorm(0, 0.5),
        sigma_B ~ dexp(1),
        sigma ~ dexp(1)
    ), data=d |> select(K, B, M), chains=4, cores=4
)
precis(m15.5, depth=2)

m15.6 <- ulam(
    alist(
        K ~ dnorm(mu, sigma),
        mu <- a + bB * B + bM*M,
        B ~ dnorm(nu, sigma_B),
        c(a, nu) ~ dnorm(0, 0.5),
        c(bB, bM) ~ dnorm(0, 0.5),
        sigma_B ~ dexp(1),
        sigma ~ dexp(1)
    ), data=d |> select(K, B, M) |> filter(!is.na(K), !is.na(B), !is.na(M)), chains=4, cores=4
)

# Why were we modelling B there in the complete case analysis? Let's remove that likelihood
m15.7 <- ulam(
    alist(
        K ~ dnorm(mu, sigma),
        mu <- a + bB * B + bM*M,
        c(a) ~ dnorm(0, 0.5),
        c(bB, bM) ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), data=d |> select(K, B, M) |> filter(!is.na(K), !is.na(B), !is.na(M)), chains=4, cores=4
)
# Compare intercept & slopes to complete case model
# m15.6 & m15.7 are identical, which is reassuring as we didn't need to model the likelihood of B
# But since these parameters weren't used in the K likelihood it didn't affect the slopes
# Between the imputation and complete case model, there is greater uncertainty in the complete case model!
# That is surprising! Maybe just from having more data available?
plot(coeftab(m15.5, m15.6, m15.7), pars=c("a", "bB", "bM"))

post <- extract.samples(m15.5)
B_impute_mu <- apply(post$B_impute, 2, mean)
B_impute_ci <- apply(post$B_impute, 2, PI)

# Plot B vs K
# Blue are data, black are imputed values
plot(d$B, d$K, pch=16, col=rangi2,
     xlab="Bain (std)", ylab="K (std)")
miss_idx <- which(is.na(d$B))
Ki <- d$K[miss_idx]
points(B_impute_mu, Ki)
for (i in 1:12) lines(B_impute_ci[, i], rep(Ki[i], 2))

# Plot M vs B
# blue are data, black are imputed
# Not bad, but the imputation is missing that B is correlated with M through
# U. So a better imputation model would take this into account
plot(d$M, d$B, pch=16, col=rangi2,
     xlab="M (std)", ylab="B (std)")
Mi <- d$M[miss_idx]
points(Mi, B_impute_mu)
for (i in 1:12) lines(rep(Mi[i], 2), B_impute_ci[, i])

m15.7 <- ulam(
    alist(
        # K as function of B and M
        K ~ normal(mu, sigma),
        mu <- a + bB * B_merge + bM*M,
        
        # M and B correlation
        MB ~ multi_normal(c(muM, muB), Rho_BM, Sigma_BM),
        matrix[29, 2]:MB <<- append_col(M, B_merge),
        
        # define B_merge as mix of obs and imputed
        vector[29]:B_merge <- merge_missing(B, B_impute),
        
        # priors
        c(a, muB, muM) ~ normal(0, 0.5),
        c(bB, bM) ~ normal(0, 0.5),
        sigma ~ dexp(1),
        Rho_BM ~ lkj_corr(2),
        Sigma_BM ~ dexp(1)
    ), data=d |> select(B, M, K), chains=4, cores=4
)
precis(m15.7, depth=3, pars=c("bM", "bB", "Rho_BM"))

# Plot again and can see that the imputed values are more in line with the data
post <- extract.samples(m15.7)
B_impute_mu <- apply(post$B_impute, 2, mean)
B_impute_ci <- apply(post$B_impute, 2, PI)

# Plot B vs K
# Blue are data, black are imputed values
plot(d$B, d$K, pch=16, col=rangi2,
     xlab="Bain (std)", ylab="K (std)")
miss_idx <- which(is.na(d$B))
Ki <- d$K[miss_idx]
points(B_impute_mu, Ki)
for (i in 1:12) lines(B_impute_ci[, i], rep(Ki[i], 2))

# Plot M vs B
# blue are data, black are imputed
# Not bad, but the imputation is missing that B is correlated with M through
# U. So a better imputation model would take this into account
plot(d$M, d$B, pch=16, col=rangi2,
     xlab="M (std)", ylab="B (std)")
Mi <- d$M[miss_idx]
points(Mi, B_impute_mu)
for (i in 1:12) lines(rep(Mi[i], 2), B_impute_ci[, i])

# merge_missing is a custom function that Richard wrote to combine observed data
# and parameters to impute
# It also identifies missing values for you, but it's not that hard to do this yourself
# and pass it into the model
stancode(m15.7)

# 15.2.3 Where is your god now? -------------------------------------------
data("Moralizing_gods")
Moralizing_gods <- as_tibble(Moralizing_gods)
# 5 cols:
#  - policty (location)
#  - year 
#  - population (log scale)
#  - moralizing_gods (binary, 0/1 with some NAs, often meaning lack of written evidence)
#  - writing (did the society have literary?)
Moralizing_gods

# Have 30 locations with varying numbers of measurements
Moralizing_gods |>
    count(polity)

# Very imbalanced classes, although a huge amount of missing data
Moralizing_gods |>
    count(moralizing_gods)

# Question: Does belief in moralizing gods increase the rate of population growth?
# Let's plot as a time-series
# Can see that moralizing gods seem to be associated with both increased time, as well as higher population
Moralizing_gods |>
    mutate(moralizing_gods = factor(ifelse(is.na(moralizing_gods), 'Missing', 
                                    ifelse(moralizing_gods == 0, 'No', 'Yes')),
                                    levels=c("Missing", "No", "Yes"))) |>
    ggplot(aes(x=year, y=population, colour=moralizing_gods)) +
        geom_point(shape=1) +
        theme_classic() +
        theme(
            legend.position = "bottom"
        )
# But it's not missing at random, it's missing when we don't have writing available
# Most of the missing data is when writing wasn't available
Moralizing_gods |>
    count(moralizing_gods, writing) 

# 15.3 Categorical errors and discrete absences ---------------------------
set.seed(9)
N_houses <- 100L
alpha <- 5
beta <- -3
k <- 0.5
r <- 0.2

cat <- rbern(N_houses, k)
# 50/50 house has a cat
table(cat)
# Reduce expected notes by 3 if have cat present
notes <- rpois(N_houses, alpha + beta*cat)
table(notes)
# Generate missing data in 20% of cases
R_C <- rbern(N_houses, r)
cat_obs <- cat
# Missing data code?
cat_obs[R_C == 1] <- -9L
table(cat_obs)
dat <- list(notes=notes, cat=cat_obs, RC=R_C, N=as.integer(N_houses))

# For non-missing data we have just the standard Poisson likelihood
# For the missing data we have: 
# Pr(cat) * Pr(notes | cat)   +   Pr(no cat) * Pr(notes | no cat)
# NB: this is averaging over the possible outcomes!
# Can give custom likelihoods to stan
# Think about how can apply this to missing Covid PCRs? It's not the same though
# as we know they didn't have the test
m15.8 <- ulam(
    alist(
        # cat known to be present or absent
        notes|RC==0 ~ poisson(lambda),
        log(lambda) <- a + b*cat,
        
        # unknown whether there was a cat or not
        notes|RC==1 ~ custom(
            log_sum_exp(  # LSE to sum log probs
                log(k) + poisson_lpmf(notes | exp(a + b)), # Pr(cat) * Pr(notes | cat), using log prob
                log(1-k) + poisson_lpmf(notes | exp(a))
            )
        ),
        
        # Priors
        a ~ normal(0, 1),
        b ~ normal(0, 0.5),
        
        # Only estimate the cat presence probability from non-missing data!
        cat|RC==0 ~ bernoulli(k),
        k ~ beta(2, 2)
    ),
    data=dat, chains=4, cores=4
)
# converges well!
precis(m15.8)

# Now to generate the probabilities of cats being present for the missing values
# We don't need to change the DAG for this or to explicitly model this
# as a function of anything else in order to calculate it
# Instead can just use Bayes rule using the notes sung!
# Again this is a nice trick! Can remember this for the Covid PCR test model where want to get
# probability of a test without making it an outcome
# Pr(cat | notes) =  Pr(notes | cat) * Pr(cat)  (likelihood which is our posterior predictive and prior) 
#                    -------------------------
#         Pr(notes | cat)*Pr(cat) + Pr(notes| no cat) * Pr(no cat)  (normalizing constant, averaging over outcomes)
m15.9 <- ulam(
    alist(
        # cat known to be present or absent
        notes|RC==0 ~ poisson(lambda),
        log(lambda) <- a + b*cat,
        
        # unknown whether there was a cat or not
        notes|RC==1 ~ custom(
            log_sum_exp(  # LSE to sum log probs
                log(k) + poisson_lpmf(notes | exp(a + b)), # Pr(cat) * Pr(notes | cat), using log prob
                log(1-k) + poisson_lpmf(notes | exp(a))
            )
        ),
        
        # Priors
        a ~ normal(0, 1),
        b ~ normal(0, 0.5),
        
        # Only estimate the cat presence probability from non-missing data!
        cat|RC==0 ~ bernoulli(k),
        k ~ beta(2, 2),
        
        # Imputed values
        gq> vector[N]:PrC1 <- exp(lpC1) / (exp(lpC1) + exp(lpC0)),
        gq> vector[N]:lpC1 <- log(k) + poisson_lpmf(notes[i] | exp(a+b)),
        gq> vector[N]:lpC0 <- log(1-k) + poisson_lpmf(notes[i] | exp(a))
    ),
    data=dat, chains=4, cores=4
)
# same estimates as before, only difference is the generated quantities
precis(m15.9)

# Problems ----------------------------------------------------------------

# 15E1
# Add measurement error on P
# Ti ~ Poisson(ui)
# log(ui) = a + b*logPreali
# a ~ normal(0, 1.5)
# b ~ normal(0, 1)
# Pi ~ normal(Preali, sigmaPi)

# 15E2
# Have missing values on P, allow for imputation
# Similar idea, put a prior on Pi (here the prior is directly on Pi and Pi is
# used in the model, whereas before the prior was on Pi but Preali was used in the
# the linear predictor
# Ti ~ Poisson(ui)
# log(ui) = a + b*logPi
# a ~ normal(0, 1.5)
# b ~ normal(0, 1)
# Pi ~ normal(0, 1)

# 15M1
# Assuming talking about the model on page 506, It assumes the missing data were generated through 
# a varying-intercept model (i.e. there's a hierarchical prior on the missing data prior)

# 15M2
# Refit the imputation model with Beta
# As B is [0, 1]
d <- milk |> 
    as_tibble() |>
    mutate(
        neocortex.prop = neocortex.perc / 100, # % to decimal
        logmass=log(mass),
        K=standardize(kcal.per.g),
        B=standardize(neocortex.prop),
        M=standardize(logmass)
    )
m15.5 <- ulam(
    alist(
        K ~ dnorm(mu, sigma),
        mu <- a + bB * B + bM*M,
        B ~ dnorm(nu, sigma_B),
        c(a, nu) ~ dnorm(0, 0.5),
        c(bB, bM) ~ dnorm(0, 0.5),
        sigma_B ~ dexp(1),
        sigma ~ dexp(1)
    ), data=d |> select(K, B, M), chains=4, cores=4
)
# What prior to use for Beta?
# Want to not work with standardised B otherwise lose the constraint
# Quite high!
dens(d$neocortex.prop)
# What's a genericy flat beta?
# Let's go 2,2
dens(rbeta(1e4, 1, 1))
dens(rbeta(1e4, 2, 2), add=TRUE, col='blue')
dens(rbeta(1e4, 10, 10), add=TRUE, col="red")

# But of course need a prior on the hyper-parameters
# What's an approriate distribution here?
# Can you have a negative shape or scale?
# Nope, they must be positive
# So could use an exponential
#dens(rbeta(1e4, -2, 2))

# But if use exp(1) on the 2 parameters end up getting too close to 0!
b1 <- rexp(1e4, 1)
b2 <- rexp(1e4, 1)
dens(rbeta(1e4, b1, b2))

# Can try half-normal instead
# That's more like it!
# How to use half-normal in ulam?
# NB: Used tighter variance on these parameters else got 0 and 1
b1 <- abs(rnorm(1e4, 2, 0.25))
b2 <- abs(rnorm(1e4, 2, 0.25))
dens(rbeta(1e4, b1, b2))
summary(rbeta(1e4, b1, b2))

# Can't get this to fit using half-normal - it sounds like it's spitting out negatives so maybe the constraint isn't working
#m15m2 <- ulam(
#    alist(
#        K ~ dnorm(mu, sigma),
#        mu <- a + bB * B + bM*M,
#        B ~ dbeta(b1, b2),
#        c(a) ~ dnorm(0, 0.5),
#        c(bB, bM) ~ dnorm(0, 0.5),
#        sigma ~ dexp(1),
#        c(b1, b2) ~ normal(2, 0.25)
#    ), data=d |> select(K, B=neocortex.prop, M), chains=4, cores=4,
#    constraints=list(b1="lower=0", b2="lower=0")
#)

# Will try log-normal
# Has a fair bit of skew but should be usable
dens(rlnorm(1e4, 0, 0.5))
# This is what I wanted to try and replicate
dens(abs(rnorm(1e4, 2, 0.5)))

# That looks acceptable, can try it in the model
b1 <- rlnorm(1e4, 2, 0.5)
b2 <- rlnorm(1e4, 2, 0.5)
dens(rbeta(1e4, b1, b2))
summary(rbeta(1e4, b1, b2))

# No that still fails, getting errors about initialisation between -2, 2 which seems to cause negative beta values
# But no idea why they are being initialised in this way?
# Maybe I need to put a constraint on B so that it can't be initialised to negative?
# No that also fails...
#m15m2 <- ulam(
#    alist(
#        K ~ dnorm(mu, sigma),
#        mu <- a + bB * B + bM*M,
#        B ~ dbeta(b1, b2),
#        c(a) ~ dnorm(0, 0.5),
#        c(bB, bM) ~ dnorm(0, 0.5),
#        sigma ~ dexp(1),
#        c(b1, b2) ~ dlnorm(0, 0.5)
#    ), data=d |> select(K, B=neocortex.prop, M), chains=4, cores=4,
#    constraint=list(B="lower=0")
#)

# Is it because B is too big?
# I.e. what if we center it?
# Nope also still fails!
# Not really sure why it's failing
mean(d$neocortex.prop, na.rm=TRUE)
m15m2 <- ulam(
    alist(
        K ~ dnorm(mu, sigma),
        mu <- a + bB * (B-0.676) + bM*M,
        B ~ dbeta(b1, b2),
        c(a) ~ dnorm(0, 0.5),
        c(bB, bM) ~ dnorm(0, 0.5),
        sigma ~ dexp(1),
        c(b1, b2) ~ dlnorm(0, 0.5)
    ), data=d |> select(K, B=neocortex.prop, M), chains=4, cores=4,
    constraint=list(B="lower=0")
)

# Put constraint on B_impute
# Looking at the Stan code, there's a B_impute parameter which will have a wide prior by default
# Maybe this is the problem?
# Maybe I should just roll it myself?
stancode(m15m2)

code_m15m2 <- "
data{
     int N;
     int N_miss;
     int N_obs;
     vector[N] K;
     vector[N] M;
     array[N_obs] int B_obs;
     array[N_obs] int B_obs_idx;
     array[N_miss] int B_miss_idx;
}
parameters{
     real a;
     real bM;
     real bB;
     real<lower=0> sigma;
     real<lower=0> b2;
     real<lower=0> b1;
     vector<lower=0>[N_miss] B_impute;
}
model{
    vector[N] mu;
    B_impute ~ beta( b1 , b2 );
    b1 ~ lognormal( 0 , 0.5 );
    b2 ~ lognormal( 0 , 0.5 );
    sigma ~ exponential( 1 );
    bB ~ normal( 0 , 0.5 );
    bM ~ normal( 0 , 0.5 );
    a ~ normal( 0 , 0.5 );
    B_obs ~ beta( b1 , b2 );
    for ( i in 1:N_miss ) {
        mu[B_miss_idx[i]] = a + bB * (B_impute[i] - 0.676) + bM * M[B_miss_idx[i]];
    }
    for ( i in 1:N_obs ) {
        mu[B_obs_idx[i]] = a + bB * (B_obs[B_obs_idx[i]] - 0.676) + bM * M[B_obs_idx[i]];
    }
    K ~ normal( mu , sigma );
}"

dlist_miss <- list(
    N=nrow(d),
    N_miss=sum(is.na(d$neocortex.prop)),
    N_obs=sum(!is.na(d$neocortex.prop)),
    K=d$K,
    M=d$M,
    B_obs=d$neocortex.prop[!is.na(d$neocortex.prop)],
    B_obs_idx = which(!is.na(d$neocortex.prop)),
    B_miss_idx = which(is.na(d$neocortex.prop))
)

# Still not sure why this isn't working!
# Still saying that B_impute is outside of [0, 1], but how can that be the case
# with this beta prior?
# NB: I tried commenting out the B_obs ~ dbeta... likelihood in case this was the problem
# and it still failed
# Give up for now. Sadly that's the first medium question I've had to give up on
#m12m2 <- stan(model_code=code_m15m2, data=dlist_miss, chains=4, cores=4)

# 15M3
# Repeat divorce measurement errors but double standard errors
d <- as_tibble(WaffleDivorce)

# Measurement error on both outcome and predictor
# This time DOUBLING standard error
dlist <- list(
    D_obs = standardize(d$Divorce),
    D_sd = d$Divorce.SE*2 / sd(d$Divorce),  # Interesting standardization... I guess to keep it on the
    M_obs = standardize(d$Marriage),
    M_sd = d$Marriage.SE*2 / sd(d$Marriage),  # Interesting standardization... I guess to keep it on the
                                          # same scale as D_obs, but there's no centering needed
    A = standardize(d$MedianAgeMarriage),
    N = nrow(d)
)

m15m3 <- ulam(
    alist(
        D_obs ~ dnorm(D_true, D_sd),
        M_obs ~ dnorm(M_true, M_sd),
        vector[N]:D_true ~ dnorm(mu, sigma),
        # What is this [i] doing here?! Won't compile without it
        mu <- a + bA*A + bM*M_true[i],
        # Need to provide a prior on true covariate! Can just put weakly informative 
        # This might be what I'm missing from my York regression analysis
        vector[N]:M_true ~ dnorm(0, 1),
        a ~ dnorm(0, 0.2),
        bA ~ dnorm(0, 0.5),
        bM ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), data=dlist, chains=4, cores=4
)
# Get divergence warnings for a start!
# And 4/4 had low E-BFMI
# Extremely poor convergence on sigma
# This is like what I had with my model
precis(m15m3, depth=2)

# Try non-centering
# Yep that works!
# NB: I hadn't assigned a group prior to M_true so that is ok
m15m3 <- ulam(
    alist(
        D_obs ~ dnorm(D_Z * sigma + mu, D_sd),
        M_obs ~ dnorm(M_true, M_sd),
        # It's the D_obs prior using D_true which is itself assigned a prior
        # that's the problem. D_sd is ok as this is data
        # Want D_true non-centered
        #vector[N]:D_true ~ dnorm(mu, sigma),
        vector[N]: D_Z ~ normal(0, 1),
        # What is this [i] doing here?! Won't compile without it
        mu <- a + bA*A + bM*M_true[i],
        # Need to provide a prior on true covariate! Can just put weakly informative 
        # This might be what I'm missing from my York regression analysis
        vector[N]:M_true ~ dnorm(0, 1),
        a ~ dnorm(0, 0.2),
        bA ~ dnorm(0, 0.5),
        bM ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), data=dlist, chains=4, cores=4
)
# Much better!
precis(m15m3, depth=2)

# Quite a difference on bM!
# So now the effect of marriage is much stronger, despite having less certainty 
# in our measurements, this is because now have different effects on the regression
# from different size states that were previously much more confident
plot(coeftab(m15.2, m15m3), pars=c("bA", "bM"))

# 15M4
# Fit Pipe X->Y->Z
# Just use the same effect size for now rather than complicating it
N <- 1e3
x <- rnorm(N)
y <- rnorm(N, x)
z <- rnorm(N, y)
# So X on Y should be 0.5, and here it's 0.43
# I don't know how you'd predict how Z should influence Y, Y->Z should be 0.5
# But here we get 0.5 from the other direction. I guess the association should
# be the same strength as linear regression alone doesn't provide CI, that's what the DAG
# + theory is for. So maybe it's not unexpected to have this
# We are conditioning on a child of the outcome
# Y ~ X is the 'correct' causal model, so we are conditioning on Z too as the descendant
# of Y. Not as bad as conditioning on the outcome but close
summary(lm(y ~ z + x))

# 15M5
# Compare posterior of missing cat presence (PrC1) to simulated
# How good is model at inferring missing data?
# How can change model to improve precision of these estimates?
set.seed(9)
N_houses <- 100L
alpha <- 5
beta <- -3
k <- 0.5
r <- 0.2

cat <- rbern(N_houses, k)
# 50/50 house has a cat
table(cat)
# Reduce expected notes by 3 if have cat present
notes <- rpois(N_houses, alpha + beta*cat)
table(notes)
# Generate missing data in 20% of cases
R_C <- rbern(N_houses, r)
cat_obs <- cat
# Missing data code?
cat_obs[R_C == 1] <- -9L
table(cat_obs)
dat <- list(notes=notes, cat=cat_obs, RC=R_C, N=as.integer(N_houses))

# Here will impute the missing values
# Pr(cat | notes) =  Pr(notes | cat) * Pr(cat)  (likelihood which is our posterior predictive and prior) 
#                    -------------------------
#         Pr(notes | cat)*Pr(cat) + Pr(notes| no cat) * Pr(no cat)  (normalizing constant, averaging over outcomes)
m15.9 <- ulam(
    alist(
        # cat known to be present or absent
        notes|RC==0 ~ poisson(lambda),
        log(lambda) <- a + b*cat,
        
        # unknown whether there was a cat or not
        notes|RC==1 ~ custom(
            log_sum_exp(  # LSE to sum log probs
                log(k) + poisson_lpmf(notes | exp(a + b)), # Pr(cat) * Pr(notes | cat), using log prob
                log(1-k) + poisson_lpmf(notes | exp(a))
            )
        ),
        
        # Priors
        a ~ normal(0, 1),
        b ~ normal(0, 0.5),
        
        # Only estimate the cat presence probability from non-missing data!
        cat|RC==0 ~ bernoulli(k),
        k ~ beta(2, 2),
        
        # Imputed values
        gq> vector[N]:PrC1 <- exp(lpC1) / (exp(lpC1) + exp(lpC0)),
        gq> vector[N]:lpC1 <- log(k) + poisson_lpmf(notes[i] | exp(a+b)),
        gq> vector[N]:lpC0 <- log(1-k) + poisson_lpmf(notes[i] | exp(a))
    ),
    data=dat, chains=4, cores=4
)
precis(m15.9)

post <- extract.samples(m15.9, pars="PrC1")
post$PrC1 |> dim()
# Yes the imputations are quite good!
# The houses that do have cats have higher predicted probability
# than those that don't
tibble(
    pred_mu=colMeans(post$PrC1),
    pred_lower=apply(post$PrC1, 2, PI)[1, ],
    pred_upper=apply(post$PrC1, 2, PI)[2, ],
    actual=cat,
    house=1:N_houses
) |>
    ggplot(aes(x=actual, y=pred_mu, group=actual)) +
        geom_boxplot()

# How can we improve precision of these estimates? Not really sure, unless
# there's any predictors we're missing?

# 15M6
# dog eat homework
# simulate data for each
# can the estimate S->H be recovered?

# 1: completely at random as the missing-causing exposure D
# is d-separated from the outcome H
dag_1 <- dagitty("dag { 
 H[u];
 S -> H -> H_obs <- D;
}")
drawdag(dag_1)
# Dagitty can do this for us
dseparated(dag_1, "D", "H")
# No need to adjust
adjustmentSets(dag_1, "S", "H_obs", effect="total")
adjustmentSets(dag_1, "S", "H_obs", effect = "direct")

N <- 1e4
S <- rnorm(N)
H <- rnorm(N, -3 + 0.8*S)
D <- rbern(N)
Hm <- H
Hm[D==1] <- NA
# Yep can retrieve the relationship
lm(Hm ~ S)

# 2: exposure causes the missingness
dag_2 <- dagitty("dag { 
 H[u];
 D <- S -> H -> H_obs <- D;
}")
drawdag(dag_2)
# They aren't d-separated
dseparated(dag_2, "H", "D")
# Don't need to condition on anything to get the total effect
adjustmentSets(dag_2, "S", "H_obs", effect="total")
# But _can't_ get the direct effect
adjustmentSets(dag_2, "S", "H_obs", effect = "direct")

# But still are able to recover effects
N <- 1e4
S <- rnorm(N)
H <- rnorm(N, -3 + 0.8*S)
D <- rbern(N, inv_logit(0.2 + 0.8*S))
Hm <- H
Hm[D==1] <- NA
# Yep can retrieve the relationship
lm(Hm ~ S)

# 3: Missing influenced by an unobserved variable, 
# which also influences the outcome
dag_3 <- dagitty("dag { 
 H[u];
 X[u];
 S -> H -> H_obs <- D;
 H <- X -> D;
}")
drawdag(dag_3)
# They aren't d-separated
dseparated(dag_3, "H", "D")
# Likewise can't get the direct effect but nothing to condition on
# to get the total effect
adjustmentSets(dag_3, "S", "H_obs", effect="total")
adjustmentSets(dag_3, "S", "H_obs", effect = "direct")

# Can now get an estimate of S reasonably well
N <- 1e4
X <- rnorm(N) # Latent variable
S <- rnorm(N) # Exposure
H <- rnorm(N, -3 + 0.8*S - 2 * X)
D <- rbern(N, inv_logit(0.2 + 1.2*X))
Hm <- H
Hm[D==1] <- NA
lm(Hm ~ S)

# 4: outcome influences missingness
dag_4 <- dagitty("dag { 
 H[u];
 S -> H -> H_obs <- D;
 H -> D;
}")
drawdag(dag_4)
# They aren't d-separated
dseparated(dag_4, "H", "D")
# And again can't get the direct effect
adjustmentSets(dag_4, "S", "H_obs", effect="total")
adjustmentSets(dag_4, "S", "H_obs", effect = "direct")

N <- 100
S <- rnorm(N)
H <- rnorm(N, -3 + 0.8*S)
D <- rbern(N, inv_logit(0.2 + 1.2*H))
Hm <- H
Hm[D==1] <- NA
# Ah no the estimate is way off here
lm(Hm ~ S)

# 15H1
data("elephants")
d <- elephants |> as_tibble()
d
# Predict number of matings as function of age
m15h1 <- ulam(
    alist(
        # cat known to be present or absent
        MATINGS ~ poisson(lambda),
        log(lambda) <- a + b*AGE,
        
        # Priors
        a ~ normal(0, 1),
        b ~ normal(0, 1)
    ),
    data=elephants |> mutate(AGE=as.numeric(scale(AGE))), chains=4, cores=4
)
# Converged ok
precis(m15h1)

# Now with measurement error
# Need standardised sds, 5 years is 0.76
5 / sd(d$AGE)
m15h1_err <- ulam(
    alist(
        # cat known to be present or absent
        MATINGS ~ poisson(lambda),
        log(lambda) <- a + b*AGE_REAL,
        AGE_OBS ~ normal(AGE_REAL, 0.76),
        
        # Need prior on age real
        AGE_REAL ~ normal(0, 1),
        
        # Priors
        a ~ normal(0, 1),
        b ~ normal(0, 1)
    ),
    data=elephants |> mutate(AGE=as.numeric(scale(AGE))), chains=4, cores=4
)
# Also compiles OK, could have non-centered if needed
precis(m15h1_err)

# The difference is the massive uncertainty in the slope and intercept estimates for the
# model with error, I guess because the solution space has now increased massively
plot(coeftab(m15h1, m15h1_err))

# 15H2
# Now increase standard error until posterior for coefficient on age reaches 0
# Double standard error to 10 years
# Which is a lot given that ages vary between 25-55 years
dens(d$AGE)
10 / sd(d$AGE)  # 1.52 standardized
m15h2_10 <- ulam(
    alist(
        # cat known to be present or absent
        MATINGS ~ poisson(lambda),
        log(lambda) <- a + b*AGE_REAL,
        AGE_OBS ~ normal(AGE_REAL, 1.52),
        
        # Need prior on age real
        AGE_REAL ~ normal(0, 1),
        
        # Priors
        a ~ normal(0, 1),
        b ~ normal(0, 1)
    ),
    data=elephants |> mutate(AGE=as.numeric(scale(AGE))), chains=4, cores=4
)
# Also compiles OK, could have non-centered if needed
precis(m15h2_10)

# The difference is the massive uncertainty in the slope and intercept estimates for the
# model with error, I guess because the solution space has now increased massively
# That basically means 0 age effect
plot(coeftab(m15h1, m15h1_err, m15h2_10))

# 15H3
set.seed(100)
x <- c(rnorm(10), NA)
y <- c(rnorm(10, x), 100)
d <- tibble(x, y)

# Firstly fit with just complete cases
m15h3_complete <- ulam(
    alist(
        y ~ dnorm(mu, sigma),
        mu <- a + b * x,
        a ~ dnorm(0, 100),
        b ~ dnorm(0, 100),
        sigma ~ dexp(1)
    ), data=d |> filter(!is.na(x)), chains=4, cores=4
)
m15h3_impute <- ulam(
    alist(
        y ~ dnorm(mu, sigma),
        mu <- a + b * x,
        x ~ dnorm(0, 1),
        a ~ dnorm(0, 100),
        b ~ dnorm(0, 100),
        sigma ~ dexp(1)
    ), data=d, chains=4, cores=4
)
# Complete cases converged ok
precis(m15h3_complete)
# But awful Rhat on the imputed b!
precis(m15h3_impute, depth=2)
# Beta's effect changes from effectively 0 (complete cases) 
# to 1.45 (with high CI) in the imputed dataset
# This is because the final value is 10x higher and therefore has huge leverage
# It wasn't included in the first model and hence the slope being so low
# Not least because the imputed b is constrained to be small rather than around 100
# itself by its prior

# Looking in more detail
post <- extract.samples(m15h3_impute)
# Observe the bi-modal values for a/b/sigma/x_impute
# Where are these from?
# Looks like it stems from x_impute having 2 modes, which then forces 2 different
# intercept/slope/variances
pairs(post)

# x impute is either negative or positive
dens(post$x_impute)
# Likewise for slope or intercept
dens(post$b)

# Why does the imputed value have 2 modes then?
# With a y of 100, x should be on a similar order of magnitude but the prior constrains it
# So end up in the current scenario
# Let's look at both scenarios
impute_pos <- post$x_impute > 0
impute_neg <- post$x_impute < 0

d_pos <- d |>
    mutate(x = ifelse(is.na(x),
                      mean(post$x_impute[impute_pos]),
                      x
                      ))
x_pos <- seq(min(d_pos$x)-0.2, max(d_pos$x)+0.2, length.out=100)
y_pos <- sapply(x_pos, function(x) post$a[impute_pos] + post$b[impute_pos]*x)
pos_smooth <- tibble(
    x=x_pos,
    y=colMeans(y_pos),
    y_lower=apply(y_pos, 2, PI)[1, ],
    y_upper=apply(y_pos, 2, PI)[2, ]
)

p_pos <- d_pos |>
    ggplot(aes(x=x, y=y)) +
        geom_ribbon(aes(ymin=y_lower, ymax=y_upper), data=pos_smooth, alpha=0.3) +
        geom_line(data=pos_smooth, colour="orange", linewidth=1) +
        geom_point() +
        theme_classic()

d_neg <- d |>
    mutate(x = ifelse(is.na(x),
                      mean(post$x_impute[impute_neg]),
                      x
                      ))
x_neg <- seq(min(d_neg$x)-0.2, max(d_neg$x)+0.2, length.out=100)
y_neg <- sapply(x_neg, function(x) post$a[impute_neg] + post$b[impute_neg]*x)
neg_smooth <- tibble(
    x=x_neg,
    y=colMeans(y_neg),
    y_lower=apply(y_neg, 2, PI)[1, ],
    y_upper=apply(y_neg, 2, PI)[2, ]
)

p_neg <- d_neg |>
    ggplot(aes(x=x, y=y)) +
        geom_ribbon(aes(ymin=y_lower, ymax=y_upper), data=neg_smooth, alpha=0.3) +
        geom_line(data=neg_smooth, colour="orange", linewidth=1) +
        geom_point() +
        theme_classic()

# The only way the model can handle having such an extreme y-value with a limited x value
# is through having an extreme slope, either positive or negative are equally likely
p_pos | p_neg

# 15H4
data("Primates301")
d <- Primates301 |> as_tibble()
cc <- complete.cases(d$brain, d$body)
B <- d$brain[cc]
M <- d$body[cc]
B <- B / max(B)
M <- M / max(M)
# Assume std errors are 10% of measurement
Bse <- B*0.1
Mse <- M*0.1
dat_list <- list(B=B, M=M, Bse=Bse, Mse=Mse, N=length(B))
m15h4 <- ulam(
    alist(
        B ~ dnorm(B_true, Bse),
        M ~ dnorm(M_true, Mse),
        vector[N]:B_true ~ dlnorm(mu, sigma),
        mu <- a + b*log(M_true[i]),
        vector[N]:M_true ~ dnorm(0, 1),
        a ~ normal(0, 1),
        b ~ normal(0, 1),
        sigma ~ exponential(1)
    ), data=dat_list, chains=4, cores=4,
    start=list(M_true=dat_list$M, B_true=dat_list$B)
)
# HMC ran ok
precis(m15h4)


# How do coefficients compare with no measurement error?
m15h4_noerr <- ulam(
    alist(
        B ~ dlnorm(mu, sigma),
        mu <- a + b*log(M),
        a ~ normal(0, 1),
        b ~ normal(0, 1),
        sigma ~ exponential(1)
    ), data=dat_list, chains=4, cores=4
)
# Very little difference
precis(m15h4)
precis(m15h4_noerr)

# Plot predictions
# Can't extract posterior!
# post <- extract.samples(m15h4)       # DOESN'T WORK
# post <- extract.samples(m15h4_noerr) # Works
# Both have same class, but m15h4 says "$ operator not defined for S4 class"
class(m15h4)
class(m15h4_noerr)

# postcheck also fails!
postcheck(m15h4)

# 15H5
# Lots of missing data (301 rows in total)
colSums(is.na(d)) 

d <- Primates301 |> as_tibble()
# Now only consider complete M, so will have missing B
cc <- complete.cases(d$body)
B <- d$brain[cc]
M <- d$body[cc]
B <- B / max(B, na.rm=T)
M <- M / max(M)
length(B)
sum(is.na(B))


# What's the DAG?
# MCAR
dag_1 <- dagitty("dag {
  B[u];
  M -> B;
  B -> B_obs;
  Miss -> B_obs;
}")
drawdag(dag_1)
# This implies the following:
#  - M is independent of Miss
impliedConditionalIndependencies(dag_1)
# M looks like it's more commonly found at smaller masses
# Although not a huge difference
# But potentially enough to invalidate these assumptions
tibble(B,M) |>
    mutate(B_miss = as.factor(is.na(B))) |>
    ggplot(aes(x=M, fill=B_miss)) +
        geom_density(alpha=0.3)
tibble(B,M) |>
    mutate(B_miss = as.factor(is.na(B))) |>
    ggplot(aes(x=B_miss, y=M)) +
        geom_boxplot(alpha=0.3)

# MAR
dag_2 <- dagitty("dag {
  B[u];
  M -> B;
  B -> B_obs;
  M -> Miss;
  Miss -> B_obs;
}")
drawdag(dag_2)
# This has no implied conditional independencies so can't test it empirically
# Which means that this could be a valid missing model for our data as it doesn't require
# that M and Missingness are independent
impliedConditionalIndependencies(dag_2)

dag_3 <- dagitty("dag {
  X[u];
  B[u];
  M -> B;
  B -> B_obs;
  Miss -> B_obs;
  X -> B;
  X -> Miss;
}")
drawdag(dag_3)
# This also implies M is independent of Miss which seems unlikely
impliedConditionalIndependencies(dag_3)

# MNAR
dag_4 <- dagitty("dag {
  B[u];
  M -> B;
  B -> Miss;
  B -> B_obs;
  Miss -> B_obs;
}")
drawdag(dag_4)
# Also no implied conditional independencies
# So this is also a possibility
impliedConditionalIndependencies(dag_4)

# The difference comes down to whether I think that BRAIN SIZE or MASS
# causes the missingness. In reality it's probably an unobserved confounder that
# causes both, i.e. something that jointly determines species mass AND brain size,
# is responsible for whether these primates are less likely to have data
# Although out of the 2, I'd say mass directly makes more likely to be missing as physically
# harder to locate smaller mammals

# Impute values for Brain Size
# Don't think need to change model as already have a likelihood on the outcome B
# which will become a prior when missing

# Getting the same problem as 15M2 where a positive constrained variable (B here)
# is having negative values, since the stan code is creating a B_impute variable
# that by default has a wide prior
# Here can get it working by giving starting values for the missing B values (not same as prior!)
m15h5_impute <- ulam(
    alist(
        B ~ dlnorm(mu, sigma),
        mu <- a + b*log(M),
        a ~ normal(0, 1),
        b ~ normal(0, 1),
        sigma ~ exponential(1)
    ), data=list(B=B, M=M), chains=4, cores=4,
    start = list(B_impute=rep(0.5, 56)) 
)

# Compare to model with complete values
m15h5_complete <- ulam(
    alist(
        B ~ dlnorm(mu, sigma),
        mu <- a + b*log(M),
        a ~ normal(0, 1),
        b ~ normal(0, 1),
        sigma ~ exponential(1)
    ), data=list(B=B[!is.na(B)], M=M[!is.na(B)]), chains=4, cores=4
)

# Very little difference in slope!
precis(m15h5_impute)
precis(m15h5_complete)

# Why is that? Plot imputed and actual values
# Can't do that again because of stanfit error

# 15H6
data(WaffleDivorce)
d <- as_tibble(WaffleDivorce)

dlist <- list(
    D_obs = standardize(d$Divorce),
    D_sd = d$Divorce.SE / sd(d$Divorce),  # Interesting standardization... I guess to keep it on the
    M_obs = standardize(d$Marriage),
    M_sd = d$Marriage.SE / sd(d$Marriage),  # Interesting standardization... I guess to keep it on the
                                          # same scale as D_obs, but there's no centering needed
    A = standardize(d$MedianAgeMarriage),
    N = nrow(d)
)

# Want to model DAG
# A -> M -> D
# A -> D
# Already have D ~ M + A
# Now need A -> M (M ~ A)
m15h6 <- ulam(
    alist(
        D_obs ~ dnorm(D_true, D_sd),
        M_obs ~ dnorm(M_true, M_sd),
        vector[N]:D_true ~ dnorm(mu, sigma),
        # What is this [i] doing here?! Won't compile without it
        mu <- a + bA*A + bM*M_true[i],
        # Need to provide a prior on true covariate! Can just put weakly informative 
        # This might be what I'm missing from my York regression analysis
        vector[N]:M_true ~ dnorm(nu, tau),
        nu <- a2 + bA2*A,
        c(a, a2) ~ dnorm(0, 0.2),
        c(bA, bA2) ~ dnorm(0, 0.5),
        bM ~ dnorm(0, 0.5),
        sigma ~ dexp(1),
        tau ~ dexp(1)
    ), data=dlist, chains=4, cores=4
)
# Compare to m15.2
# Has subtle change on bA and bM
# If can't extract posterior than can't do much further
# But ideally would look at the TRUE values themselves and see if they have changed
# between models, as this is the one of the primary outputs of the model
# And as seen in an earlier question, the imputed values DRIVE the rest of the model
# (slopes, intercepts, variances), so start all diagnoses with them
precis(m15.2)
precis(m15h6)

# 15H7
val <- seq(8)
freqs <- c(18, 19, 22, NA, NA, 19, 20, 22)
# Need prior on total number of spins and marginalize over the unknown total
# Use dirichlet for "not every value is equally likely, but none of the values are twice
# as likely as any other"
# Can't work out how to start this, trying to think of a way to phrase it in the same
# manner as the cat example: Pr(Ni) = Pr(cati) * Pr(Ni|cati) + Pr(nocati) * Pr(Ni | nocati)
# But all I have is Pr(dice=4) + Pr(dice=5)...
# I also think this requires multinomial which we haven't seen before
# Likely possible to turn into categorical but makes it more complicated

# Ask Richard --------------------------------------------------------------
# Measurement error leading to perfect model


