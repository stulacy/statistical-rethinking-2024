library(rethinking)
library(tidyverse)

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


# Ask Richard --------------------------------------------------------------
# Measurement error leading to perfect model


