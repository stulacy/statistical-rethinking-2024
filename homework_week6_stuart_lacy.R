library(rethinking)
library(tidyverse)
library(dagitty)
library(ggridges)


# Question 1 --------------------------------------------------------------
predictive_prior <- function(mu, sigma, rate, N=100) {
    a_bar <- rnorm(N, mu, sigma)
    sigma_a <- rexp(N, rate)
    a_int <- rnorm(a_bar, sigma_a)
    tibble(sim=1:N, prob=inv_logit(a_int))
}

models <- list(
    "N(0, 1)  Exp(1)"=predictive_prior(0, 1, 1),
    "N(0, 1)  Exp(10)"=predictive_prior(0, 1, 10),
    "N(0, 1)  Exp(0.1)"=predictive_prior(0, 1, 0.1),
    "N(0, 10)  Exp(1)"=predictive_prior(0, 10, 1),
    "N(1, 1)  Exp(1)"=predictive_prior(1, 1, 1)
)
# The biggest difference seems to be from decreasing sigma from Exp(1) to Exp(0.1), which
# shifts the bulk of the mass to around 1
# Whereas increasing Exp(1) to Exp(10) doesn't have quite so much of a difference, slightly shifting
# the mass to the left. I don't understand why changing sigma changes the central point of the mass,
# I thought it would just increase its spread
# Changing the a_bar normal distribution has expected effects, increasing the mean from 0 to 1
# pushes the mass to the right, while increasing the variance to N(0, 10) spreads the bulk out more
# in the way that I thought changing the Exp prior would.
bind_rows(models, .id="prior") |>
    ggplot(aes(x=prob, y=prior)) +
        geom_density_ridges(alpha=0.5) +
        theme_classic() +
        labs(x="Probability of survival", y="Density") +
        scale_fill_discrete("") +
        theme(legend.position = "bottom")


# Question 2 --------------------------------------------------------------
# Standard varying intercepts model
data(reedfrogs)
d <- reedfrogs |> as_tibble()
d$tank <- 1:nrow(d)
m13.2 <- ulam(
    alist(
        S ~ dbinom(N, p),
        logit(p) <- a[tank],
        a[tank] ~ dnorm(a_bar, sigma),
        a_bar ~ dnorm(0, 1.5),
        sigma ~ dexp(1)
    ), data=d |> select(S=surv, N=density, tank),
    chains=4, log_lik=TRUE, cores=4
)

# Add predation and size to varying intercepts
# Predation is binary so can add as index or indicator var
table(d$pred)
# Likewise for size
table(d$size)
d$pid <- as.integer(d$pred)
d$sid <- as.integer(d$size)

# Want to fit 4 models:
#   - Just size
#   - Just predation
#   - Size & predation
#   - Size & predation & Interaction
m13m1_a <- ulam(
    alist(
        S ~ dbinom(N, p),
        logit(p) <- a[tank] + s[sid],
        a[tank] ~ dnorm(a_bar, sigma),
        s[sid] ~ dnorm(0, 1),
        a_bar ~ dnorm(0, 1.5),
        sigma ~ dexp(1)
    ), data=d |> select(S=surv, N=density, tank, sid),
    chains=4, log_lik=TRUE, cores=4
)

m13m1_b <- ulam(
    alist(
        S ~ dbinom(N, p),
        logit(p) <- a[tank] + pred[pid],
        a[tank] ~ dnorm(a_bar, sigma),
        pred[pid] ~ dnorm(0, 1),
        a_bar ~ dnorm(0, 1.5),
        sigma ~ dexp(1)
    ), data=d |> select(S=surv, N=density, tank, pid),
    chains=4, log_lik=TRUE, cores=4
)

m13m1_c <- ulam(
    alist(
        S ~ dbinom(N, p),
        logit(p) <- a[tank] + pred[pid] + s[sid],
        a[tank] ~ dnorm(a_bar, sigma),
        pred[pid] ~ dnorm(0, 1),
        s[sid] ~ dnorm(0, 1),
        a_bar ~ dnorm(0, 1.5),
        sigma ~ dexp(1)
    ), data=d |> select(S=surv, N=density, tank, pid, sid),
    chains=4, log_lik=TRUE, cores=4
)

m13m1_d <- ulam(
    alist(
        S ~ dbinom(N, p),
        logit(p) <- a[tank] + pred[pid] + s[sid] + pred[pid] * s[sid],
        a[tank] ~ dnorm(a_bar, sigma),
        pred[pid] ~ dnorm(0, 1),
        s[sid] ~ dnorm(0, 1),
        a_bar ~ dnorm(0, 1.5),
        sigma ~ dexp(1)
    ), data=d |> select(S=surv, N=density, tank, pid, sid),
    chains=4, log_lik=TRUE, cores=4
)

# Want to compare on the inferred variation across tanks, this is sigma
# So it starts off very high at 1.62
# Doesn't decrease when adding size, but when adding predation it shrinks a lot
# It then doesn't change much when adding predation & size, or them + their interaction
# So adding predation removes some of the inter-tank variability
# So predation is likely a key predictor of survival, and also varies between tank
coeftab(m13.2, m13m1_a, m13m1_b, m13m1_c, m13m1_d)

# There's minimal predictive difference, although the posteriors for the coefficients changed
# a lot. This isn't surprising, the additional variation between tanks that is caused by the predation
# was accounted for by the wider between-group variance (sigma) when predation wasn't included in
# the model. This doesn't effect the model predictions.
compare(m13.2, m13m1_a, m13m1_b, m13m1_c, m13m1_d)

# This suggests the following DAG:
# I don't think size affects anything
dag2 <- dagitty("dag {
  Tank;
  Predation;
  Survival;
  Tank -> Survival;
  Predation -> Tank;
  Predation -> Survival;
}")
coordinates(dag2) <- list(x=c(Predation=0, Survival=2, Tank=1),
                          y=c(Predation=1, Survival=1, Tank=0))
drawdag(dag2)


# Question 3 --------------------------------------------------------------
data(Trolley)
m13h2_noid <- ulam(
    alist(
        R ~ dordlogit(phi, cutpoints),
        phi <- a + bA*A + bC*C + BI*I,  # I personally find the fully explicit model clearer than separating it like this
        BI <- bI + bIA*A + bIC*C,
        a ~ dnorm(0, 1),
        c(bA, bI, bC, bIA, bIC) ~ dnorm(0, 0.5),
        cutpoints ~ dnorm(0, 1.5)
    ),
    data=Trolley |> select(R=response, A=action, I=intention, C=contact), 
    chains=4,
    cores=4,
    log_lik=TRUE
)

# NB: Using non-centered parameterisation as got low neff without it
m13h2_hierarchical <- ulam(
    alist(
        R ~ dordlogit(phi, cutpoints),
        phi <- a_bar + sigma_a * a_z[id] + bA*A + bC*C + BI*I,  # I personally find the fully explicit model clearer than separating it like this
        a_z[id] ~ dnorm(0, 1),
        BI <- bI + bIA*A + bIC*C,
        c(bA, bI, bC, bIA, bIC) ~ dnorm(0, 0.5),
        a_bar ~ dnorm(0, 1),
        sigma_a ~ dexp(1),
        cutpoints ~ dnorm(0, 1.5),
        gq> vector[id]:a <<- a_bar + a_z*sigma_a
    ),
    data=Trolley |> mutate(id = as.integer(id)) |> select(id, R=response, A=action, I=intention, C=contact), 
    chains=4, cores=4, log_lik = TRUE
)

# Adding in varying intercept for user id results in a massive predictive change
# suggesting that there is a significant amount of between-user variance in terms
# of their baseline moral views
compare(m13h2_noid, m13h2_hierarchical, func = PSIS)

# There are 331 unique individuals in this dataset
# The model has a_bar = 0.8, but sigma_a is relatively high at ~=1.9, suggesting a large
# amount of inter-individual variance
# NB: There are relatively few neff for sigma, although rhat is ~okay at 1.02...
# In my first attempt with a centered parameterisation a_bar had Rhat of 1.6 and
# < 100 n_eff
precis(m13h2_hierarchical, pars=c("a_bar", "sigma_a"))

# Let's run some counter factual analysis by holding all inputs constant except
# the individual, we'll use the baseline model with no action, intent, or contact
# Get get the distribution across individuals holding everything else constant
# For both models
preds_hier <- link(m13h2_hierarchical, 
     data = tibble(
         A=0, C=0, I=0, id=unique(as.integer(Trolley$id))
     )
)
preds_default <- link(m13h2_noid,
     data = tibble(
         A=0, C=0, I=0
     )
)

preds_hier_summary <- tibble(
    id=1:331, 
    mu = colMeans(preds_hier$phi),
    lower=apply(preds_hier$phi, 2, PI)[1, ],
    upper=apply(preds_hier$phi, 2, PI)[2, ]
)

# Order
id_ord <- preds_hier_summary |>
    arrange(mu) |>
    pull(id)
preds_hier_summary$id <- factor(preds_hier_summary$id, levels=id_ord)
pooled_mean <- colMeans(preds_default$phi)
pooled_lower <- apply(preds_default$phi, 2, PI)[1, ]
pooled_upper <- apply(preds_default$phi, 2, PI)[2, ]

# Plot!
# Huge amount of inter-person variance
# 64% of people have an average mu > 0, i.e. more morally permissible, while 36% are < 0
# Looks relatively normally distrubted, although some largish tails on the positive end
preds_hier_summary |>
    ggplot(aes(x=id)) +
        geom_point(aes(y=mu)) +
        geom_errorbar(aes(ymin=lower, ymax=upper)) +
        theme_classic()

preds_hier_summary |>
    summarise(mean(mu > 0))

