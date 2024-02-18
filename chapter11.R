library(rethinking)
library(tidyverse)
library(dagitty)

# 11.1 Binomial regression ------------------------------------------------

data("chimpanzees")
d <- chimpanzees
# 4 conditions:
# 2 binary vars: prosoc_left and condition
#   - two food items on right and no partner
#   - two food items on left and no partner
#   - two food items on right and partner present
#   - two food items on left and partner present
# Create 1 indicator var to store these cases
d$treatment <- 1 + d$prosoc_left + 2*d$condition
# Yep that's right
# But rather than have binary outcome for pulled left/right and dependent variable for whether
# the food is on left or right, why not just have the outcome as "pulled lever the same side as person"?
# Ah I guess because that only works when there was a person!
# 7 chimps in total but with different trials
xtabs(~treatment + prosoc_left + condition, d)

# Demonstrating flat priors in logistic regression
# y ~ Binom(1, p)
# logit(p) = a
# a ~ N(0, w)
# Here trying with w=10, so a wide prior for OLS sure but nothing crazy
m11.1 <- quap(
    alist(
        pulled_left ~ dbinom(1, p),
        logit(p) <- a,
        a ~ dnorm(0, 10)
    ), data=d
)
# Let's sample from prior
set.seed(1999)
prior <- extract.prior(m11.1, n=1e4)
# These are on logit scale!
head(prior$a)
p <- inv_logit(prior$a)
# Wow, very bi-modal. "A flat prior in the logit space is not a flat prior in the outcome
# probability space"
dens(p, adj=0.1)

p <- quap(
    alist(
        pulled_left ~ dbinom(1, p),
        logit(p) <- a,
        a ~ dnorm(0, 1.5)
    ), data=d
) |>
    extract.prior(n=1e4) 
# Much more 'flat' with N(0, 1.5)
p$a |>
    inv_logit() |>
    dens(adj=0.1)

# Now for priors for treatment effects, which are also intercepts (one per treatment type)
# Example: N(0, 10)
m11.2 <- quap(
    alist(
        pulled_left ~ dbinom(1, p),
        logit(p) <- a + b[treatment],
        a ~ dnorm(0, 1.5),
        b[treatment] ~ dnorm(0, 10)
    ), data=d
)

set.seed(1999)
prior <- extract.prior(m11.2, n=1e4)
p <- sapply(1:4, function(k) inv_logit(prior$a + prior$b[, k]))

# These priors are the priors for each treatment
# We're more interested in the prior DIFFERENCES between treatments
# But just for now can see that these have the same bathtub shape as with
# a ~ N(0, 10)
colnames(p) <- paste0("treatment_", 1:4)
as_tibble(p) |>
    mutate(i=row_number()) |>
    pivot_longer(-i) |>
    ggplot(aes(x=value, colour=name)) +
        geom_density()

# So difference between first 2 treatments:
# Is again bi-modal!
dens(abs(p[, 1] - p[, 2]), adj=0.1)

# Try again with N(0, 0.5)
m11.3 <- quap(
    alist(
        pulled_left ~ dbinom(1, p),
        logit(p) <- a + b[treatment],
        a ~ dnorm(0, 1.5),
        b[treatment] ~ dnorm(0, 0.5)
    ), data=d
)
set.seed(1999)
prior <- extract.prior(m11.3, n=1e4)
p <- sapply(1:4, function(k) inv_logit(prior$a + prior$b[, k]))
# Now mean difference is 10%
mean(abs(p[, 1] - p[, 2]))
# And visually looks a lot more in line with expectations
dens(abs(p[, 1] - p[, 2]), adj=0.1)

# Now fit with HMC
data_list <- list(
    pulled_left = d$pulled_left,
    actor=d$actor,
    treatment = as.integer(d$treatment)
)

# NB: adding log_lik allows for calculation of PSIS/WAIC straight from the fitted object
m11.4 <- ulam(
    alist(
        pulled_left ~ dbinom(1, p),
        logit(p) <- a[actor] + b[treatment],  # NB: varying intercept for each actor
        a[actor] ~ dnorm(0, 1.5),
        b[treatment] ~ dnorm(0, 0.5)
    ), data=data_list, chains=4, log_lik=TRUE
)

# Different signs on alphas!
# Also quite different treatment effects...
# Weak SEs on treatment effects, and not the highest ESS, although rhat is low
# Remember b[2] and b[4] are when the food was on the left, but the partner was only in b[4],
# which I guess makes sense as b[4] >> b[2]
# Also b[3] and b[4] have opposite signs, which could make sense because in b[3] the food was on
# the right, so this could mean that the chimp pulled the right lever
precis(m11.4, depth=2)

# Let's look graphically, firstly at the intercepts
post <- extract.samples(m11.4)
p_left <- inv_logit(post$a)
# Nb: chimp 2 nearly always pulled the left lever regardless!
plot(precis(as.data.frame(p_left)), xlim=c(0, 1))

# Now treatment effects
labs <- c("Food Right/No partner", "Food Left/No partner", "Food Right/Partner", "Food Left/Partner")
# Yep see those massive SEs
plot(precis(m11.4, depth=2, pars="b"), labels=labs)

# We want evidence that chimps chose the prosocial (i.e. food on left) option more when a partner is present
# So want to compare first row with third, and second with fourth
diffs <- list(
    db13 = post$b[, 1] - post$b[, 3],
    db24 = post$b[, 2] - post$b[, 4]
)
# The contrast with the prosocial (db24) is actually quite low, espec compared to db13.
# The interpretation of this is that chimps were more likely to use the anti-social move when there was a 
# partner present than when there wasn't, then they were more likely to use the pro-social move
# when a partner was present vs absent (NB: I might have missed some detail in the methodology as
# this isn't Richard's interpretation. He uses db13 as slight evidence chimps were social!)
# Actually I see why, looking at the 4 treatment log-odds, the log-odds for Food Right/Partner
# are more negative than Food Right / No partner (i.e. more likely to pull right lever when partner
# was present)
# So in the contrast plot, the chimps were more likely to pull left when there was no partner, i.e.
# more likely to pull right when there was a partner
# And for b24, the chimps were actually more likely to pull left when there was no partner!
# For a 'successful' treatment in b24 woudl expect to see _negative_ result
# The scale is log-odds of pulling the left-lever
plot(precis(diffs))

# Posterior predictive check
# Summarise proportions of left pulls for each actor in each treatment and plot vs predictions
# My life for the tidyverse...
pl <- by(d$pulled_left, list(d$actor, d$treatment), mean)
dim(pl)

ndata <- expand_grid(actor=unique(d$actor), treatment=unique(d$treatment))
post <- link(m11.4, data=ndata)
# Is there any difference between sim and link for binomial?
# For Gaussian sim added in the observation noise (sigma), but here we
# have a 1-param distribution with no noise
post2 <- sim(m11.4, data=ndata)
# Ah yes of course, sim returns an output simulated from the actual distribution
# (Binomial here), whereas link just returns the linear predictor.
# So even though there's only 1 parameter to Binomial (ignoring N...), the draw
# from the Binomial adds some uncertainty that sim takes into account
# Returns p
dens(post[, 1])
# Returns discrete 0/1
dens(post2[, 1])
table(post2[, 1]) 

# So rather than comparing proportions from data to linear predictor p,
# I'll compare proportions from data to proportions of simulated outcomes
ndata$propSim <- colMeans(post2)
d_plt <- d |>
    group_by(actor, treatment) |>
    summarise(propActual = mean(pulled_left)) |>
    ungroup() |>
    inner_join(ndata, by=c("actor", "treatment"))

# The actual data looks right, but the simulated values are all far lower than the book shows
# Or at least are around 50% on average!
d_plt |>
    pivot_longer(c(propActual, propSim)) |>
    mutate(hand = ifelse(treatment %% 2 == 1, 'Right', 'Left'),
           partner = ifelse(treatment > 2, "Present", "Absent")) |>
    ggplot(aes(x=partner, y=value, colour=hand)) +
        geom_point() +
        geom_line(aes(group=hand)) +
        facet_grid(name ~ actor) +
        geom_hline(yintercept=0.5, linetype="dashed") +
        theme_bw() +
        theme(legend.position="bottom") +
        labs(x="Partner Present?", y="Proportion pulled left lever")

# Try again using the simulated proportions from the link!
# Yep that's consistent with the book
# So basically in the real data there's no strong effect between actors
# And neither is there in the simulated data really
ndata$propSim <- colMeans(post)
d |>
    group_by(actor, treatment) |>
    summarise(propActual = mean(pulled_left)) |>
    ungroup() |>
    inner_join(ndata, by=c("actor", "treatment")) |>
    pivot_longer(c(propActual, propSim)) |>
    mutate(hand = ifelse(treatment %% 2 == 1, 'Right', 'Left'),
           partner = ifelse(treatment > 2, "Present", "Absent")) |>
    ggplot(aes(x=partner, y=value, colour=hand)) +
        geom_point() +
        geom_line(aes(group=hand)) +
        facet_grid(name ~ actor) +
        geom_hline(yintercept=0.5, linetype="dashed") +
        theme_bw() +
        theme(legend.position="bottom") +
        labs(x="Partner Present?", y="Proportion pulled left lever")

# Now fit a model with index variables for side and whether a partner was present
# (I think I commented before saying this is the natural way I'd think about it)
# Richard says this isn't the approach we took because we don't just care about
# the effect of a partner being present, it's the *interaction* of a partner present
# _and_ the food on the pro-social side that makes the difference
# But could encode this in a model explicitly?
d$side <- d$prosoc_left + 1 # Right 1, left 2
d$cond <- d$condition + 1   # No partner 1, partner 2
dat_list2 <- list(
    pulled_left = d$pulled_left,
    actor=d$actor,
    side=d$side,
    cond=d$cond
)

m11.5 <- ulam(
    alist(
        pulled_left ~ dbinom(1, p),
        logit(p) <- a[actor] + bs[side] + bc[cond],
        a[actor] ~ dnorm(0, 1.5),
        bs[side] ~ dnorm(0, 0.5),
        bc[cond] ~ dnorm(0, 0.5)
    ),
    data=dat_list2, chains=4, log_lik = TRUE
)

# Model without interaction is slightly preferred here in terms of predictive accuracy
# But since we're interested in treatment effecs, we need m11.4
# NB: Both models have the same number of parameters! 11! (7 for actors, 
# then 4 for treatment/side combinations either separately or in 1)
compare(m11.5, m11.4, func=PSIS)

# Can get loglik
foo <- extract.samples(m11.4, clean=TRUE)
bar <- extract.samples(m11.4, clean=FALSE)
# By default just get the posterior samples
names(foo)
# But if have asked for log_lik=TRUE in ulam, then get lp__ and log_lik fields
# NB: log_lik is defined as a generated quantity using binomial_lpmf (or equivalent)
names(bar)
# log_lik is the log prob of each observation at each sample
dim(bar$log_lik)
# lp__ is the log-prob of each sample (is this summed over obs?)
dim(bar$lp__)
# Not directly the sums...
head(rowSums(bar$log_lik))
head(bar$lp__)

# 11.1.2 Relative shark and absolute deer ---------------------------------
# Rather than looking at effect on the absolute scale (i.e. how likely are we to pull a level),
# can also look at relative effects. This is Proportional Odds for Logistic Regression
# Simply exponentiate parameter
post <- extract.samples(m11.4)
# This is relative odds of pulling lever on left when person present vs absent
# I.e. 8% reduction in pulling left lever
# Always useful to look at both relative and absolute effects!
mean(exp(post$b[, 4]-post$b[, 2]))

# 11.1.3 Aggregated binomial: Chimpanzees again, condensed ----------------
# Organise data so that each row is a collection of binaries, only useful when
# have multiple binary outcomes per subject
d_aggregated <- aggregate(
    d$pulled_left,
    list(
        treatment=d$treatment, actor=d$actor,
        side=d$side, cond=d$cond
    ),
    sum
)
colnames(d_aggregated)[5] <- "left_pulls"
# So now have each unique combination of actor, side, treatment, condition, and number of left pulls
d_aggregated

# Can get same inference as before
dat <- with(d_aggregated, 
            list(
                left_pulls=left_pulls,
                treatment=treatment,
                actor=actor,
                side=side,
                cond=cond
            ))

m11.6 <- ulam(
    alist(
        left_pulls ~ dbinom(18, p),  # 18 combinations of side/treatment/condition
        logit(p) <- a[actor] + b[treatment],  # NB: varying intercept for each actor
        a[actor] ~ dnorm(0, 1.5),
        b[treatment] ~ dnorm(0, 0.5)
    ), data=dat, chains=4, log_lik=TRUE
)
# Very similar - if not identical - coefficents!
plot(coeftab(m11.4, m11.6))
# But very diferent PSIS! I guess because different number of observations
# General advice: If going to use WAIC or PSIS, use logistic regression rather than aggregated!
compare(m11.6, m11.4, func=PSIS)

# Aggregated binomial: Graduate school admisisons -------------------------
# Above had 18 trials per subject - not always the case!
data(UCBadmit)
d <- UCBadmit
# Data from 6 depts for male/female applications
# Research question is to investigate any gender bias
# 1 row is dept/sex, and want to find the p for male and female across
# all depts. As ever will use indicator var for sex, not index
# Admitted_dept_sex ~ Binom(N_applied_dept_sex, p_dept_sex)
# logit(p_dept_sex) ~ alpha_sex
# alpha_sex ~ N(0, 1.5)
dat_list <- list(
    admit=d$admit,
    applications=d$applications,
    gid=ifelse(d$applicant.gender=='male', 1, 2)
)
m11.7 <- ulam(
    alist(
        admit ~ dbinom(applications, p),
        logit(p) <- a[gid],
        a[gid] ~ dnorm(0, 1.5)
    ), 
    data=dat_list, chains=4
)
# Good ess and rhat
# So female has lower log odds
# With index variable would exponentiate to get the relative log-odds
# but can't do that directly with indicator variable. Instead want to calculate
# contrasts!
precis(m11.7, depth=2)
exp(-0.22)
exp(-0.83)

post <- extract.samples(m11.7)
diff_a <- post$a[, 1] - post$a[, 2]
diff_p <- inv_logit(post$a[, 1]) - inv_logit(post$a[, 2])
# So the absolute difference is 0.61 log odds
precis(list(diff_a = diff_a, diff_p =diff_p))
# And relative difference is 0.14 on probability scale

# Posterior predictive checks
# Can use the `postcheck` function
# Admittedly, this isn't the most useful plot in the world...
postcheck(m11.7)

# Model doesn't fit very well!
# The dept with cases 11 & 12 is very below the modelled predictions, while depts 1 & 2 
# (cases 1-2 and 3-4) are much higher
# Also note that males are the first points in each pair and women the second.
# There are only 2 depts where women have lower admittence rates than men: cases 5-6, and 8-9
# Yet the average has women lower. Is this Simpson's paradox?
# The explanation isn't really shown on this plot, but basically although women had a higher _rate_
# of admittance to dept A, B etc..., they didn't apply to these in big numbers!
# Instead they applied to the depts with lower admittance rates (like F) and so while on a dept-by-dept
# basis it looks like there's no bias (or indeed a slight anti-male one), looking at the whole dataset
# shows the full story.

# Instead want to ask a better question, like "what is the average difference in probability of admission
# between men and women within depts?" I.e. taking dept-baseline into account
# So each probability is a function of sex & dept!
# Admitted_dept_sex ~ Binom(N_applied_dept_sex, p_dept_sex)
# logit(p_dept_sex) ~ alpha_sex + beta_dept
# alpha_sex ~ N(0, 1.5)
# beta_sex ~ N(0, 1.5)
dat_list$dept_id <- rep(1:6, each=2)

m11.8 <- ulam(
    alist(
        admit ~ dbinom(applications, p),
        logit(p) <- a[gid] + delta[dept_id],
        a[gid] ~ dnorm(0, 1.5),
        delta[dept_id] ~ dnorm(0, 1.5)
    ), 
    data=dat_list, chains=4, iter=4000  # Bumped up to 4,000 samples!
)
# Not great ess, rhat is ok.
# I imagine this is why we bulked up iters
precis(m11.8, depth=2)

# What happens if we use the default 2,000 samples?
# Would expect even fewer ess
# And yep, ~190
ulam(
    alist(
        admit ~ dbinom(applications, p),
        logit(p) <- a[gid] + delta[dept_id],
        a[gid] ~ dnorm(0, 1.5),
        delta[dept_id] ~ dnorm(0, 1.5)
    ), 
    data=dat_list, chains=4
) |>
    precis(depth=2)

# Back to the model results, and now we can see that the male acceptance rate is a little lower
# than that for womens
# Can also see the large difference in acceptance rate per dept
precis(m11.8, depth=2)

# So we compare on contrasts
post <- extract.samples(m11.8)
diff_a <- post$a[, 1] - post$a[, 2]
diff_p <- inv_logit(post$a[, 1]) - inv_logit(post$a[, 2])
# So the absolute difference is 0.10 log odds, or 2% probability
precis(list(diff_a = diff_a, diff_p =diff_p))

# Now to complete the story and show the application numbers per dept
# And yep women applied to the tougher depts
d |>
    group_by(dept) |>
    mutate(n_dept = sum(applications)) |>
    ungroup() |>
    mutate(prop_applied = applications / n_dept * 100) |>
    select(applicant.gender, prop_applied, dept) |>
    pivot_wider(names_from=dept, values_from=prop_applied)

# Much better preds now!
postcheck(m11.8)

# This is why, as dept is a confounder
# gender influences dept, but dept also affects acceptance rate
# so to estimate the direct effect G -> A need to condition on D
# Which m11.8 does
dag_dept <- dagitty("dag { G -> A <- D <- G }")
drawdag(dag_dept)

# As shown here
adjustmentSets(dag_dept, exposure="G", outcome="A", effect = "direct")
# But remember that 'total' causal effect is different!
# Here I would fit a model of both D ~ G, and A ~ G + D and then hold G constant
# and see the counterfactual on A
# NB: I think I've really misunderstand total causal effect here, stemming from
# Chapter 5 where it looked at the 'total COUNTERFACTUAL effect', which did this, fit
# a joint model of both Confounder ~ Exposure, and Outcome ~ Exposure + Confounder
# And simulate Confounder using Exposure and then Outcome from both.
# I mistook this to mean that this is how you calculate the TOTAL CAUSAL EFFECT
# when in fact all you have to do is fit the model as described by adjustmentSets
adjustmentSets(dag_dept, exposure="G", outcome="A", effect = "total")
adjustmentSets(dag_dept, exposure="G", outcome="A", effect = "direct")

# 11.2 Poisson regression -------------------------------------------------
data(Kline)
d <- Kline
# 10 island states with population, number toolsets, and contact rates between pops
# Where contact rate is ordinal
d
d$P <- scale(log(d$population))
d$contact_id <- ifelse(d$contact == 'high', 2, 1)
# Want to model total tools as a function of:
#   - log population size (already believe it's the magnitude that matters)
#   - contact rate (more contact rate = more tools)
#   - should be an interaction between population and contact, i.e. if have high 
#     population then high contact should have greater effect
# So model:
#  T_i ~ Poisson(lambda_i)
#  log(lambda_i) = a_cid[i] + b_cid[i] * log(P_i)
# Then need priors on a_cid and b_cid
# a_cid is contact ordinal variable as intercept
# b_cid is contact ordinal variable as multiplier by population
# No other intercept needed, if were using index variables would need one
# for the 'default/zeroth' contact state

# For priors want to experiment to see what the log transform does
# I.e. if had y_i ~ Poisson(l_i)
#              log(l_i) = a
#                     a ~ N(0, 10)
# Stupidly big values, with most of weight at a tiny fraction 
foo <- exp(rnorm(1e3, 0, 10))
dens(foo)
# I.e. 52% of values are < 1! 
# This is because we're centering on 0, so ~50% of values will be < 0 in 
# measurement scale, or < 1 in log scale.
# Then have massive tail from the high variance
# Ideally want something with a flatter peak and far less tail
mean(foo < 0.5)

# If reduce the variance we get something that looks a lot better
# But the mode is still around 1, ideally would be a bit higher
foo <- exp(rnorm(1e3, 0, 0.5))
dens(foo)

# Richard suggests using mean 3, which gives expected value around 20
# Expected value is exp(u + sigma^2 / 2)
foo <- exp(rnorm(1e3, 3, 0.5))
dens(foo)

# For the coefficient of log population will again start with N(0, 10) in combination
# with the selected intercept
# Plotting on STANDARDIZED X-AXIS which is the unit of modelling
plot_prior <- function(a_mu, a_sigma, b_mu, b_sigma, N=100) {
    a <- rnorm(N, a_mu, a_sigma)
    b <- rnorm(N, b_mu, b_sigma)
    plot(NULL, xlim=c(-2, 2), ylim=c(0, 100))
    for (i in 1:N) {
        curve(exp(a[i] + b[i] * x), add=TRUE, col=grau())
    }
}

# Firstly with N(0, 10) for beta
# Want to see that at least get sensible values out on the y-axis
# And that there's a general linear increase from x=-2 (-2SDs) to x=+2 (+2SDs)
# But instead they blow up around mean logPop!
plot_prior(3, 0.5, 0, 10)

# Massively reducing variance helps contain the priors to the measurement scale (i.e. most
# of these outcomes are plausible values)
plot_prior(3, 0.5, 0, 0.2)

# What happens if we change the mean here?
# Get explosive trends that aren't realisic! Don't expect the relationship between
# log pop and tools to be that strong!
plot_prior(3, 0.5, 2, 0.2)

# Maybe N(1, 0.2) isn't so bad!
# But actually since we've already logged pop we should probably try and stick to linear
# If we were working with raw population this might be more appropriate
plot_prior(3, 0.5, 1, 0.2)

# Have plotted priors on standardized x-axis
# But population size has natural zero and so it helps to see that on plot
# So now plot on raw log scale (NB: doesn't mean won't still use z-transform for model,
# just choosing here the plotting scale)
# That doesn't look too bad, a bit less linear that had hoped but most of lines are
# linearish. Quite a few negative trends but don't want to force that relationship (would rather
# model find it)
# Still a few too many datasets with > 100 tools (roughly the max we'd expect to see)
a <- rnorm(100, 3, 0.5)
b <- rnorm(100, 0, 0.2)
x_seq <- seq(from=log(100), to=log(200000), length.out=100)
lambda <- sapply(x_seq, function(x) exp(a + b*x))
plot(NULL, xlim=range(x_seq), ylim=c(0, 500))
for (i in 1:100) {
    lines(x_seq, lambda[i, ], col=grau(), lwd=1.5)
}

# Now plotting on natural population scale
# Can see the log relationship clearly here
# The Poisson link functions means we have log-linear relationship
# And the logging of population turns this into log-log
plot(NULL, xlim=range(exp(x_seq)), ylim=c(0, 500))
for (i in 1:100) {
    lines(exp(x_seq), lambda[i, ], col=grau(), lwd=1.5)
}

# Going to fit the interaction model & intercept only (1 intercept, not varying by contact id)
dat <- list(
    T = d$total_tools,
    P = d$P,
    cid = d$contact_id
)

m11.9 <- ulam(
    alist(
        T ~ dpois(lambda),
        log(lambda) <- a,
        a ~ dnorm(3, 0.5)
    ), data=dat, chains=4, log_lik=TRUE
)

# Getting the same Stan error I've had before!
# Because trying to do element-wise multiplication, need to use 2 vectors
# whereas here have vector (b[cid]) & array(P)
# Will just have to follow along with book without replicating unless can
# tweak the Stan code easily
m11.10 <- ulam(
    alist(
        T ~ dpois(lambda),
        log(lambda) <- a[cid] + b[cid] * P,
        a[cid] ~ dnorm(3, 0.5),
        b[cid] ~ dnorm(0, 0.2)
    ), data=dat, chains=4, log_lik=TRUE
)

# Much better PSIS for the interaction model, but look at the warning about k!
# Means we have some outliers
# Also note the effective parameters (pPSIS) is higher for the intercept only model!
compare(m11.9, m11.10, func=PSIS)

k <- PSIS(m11.10, pointwise=TRUE)$k
plot(dat$P, dat$T, xlab="Log pop (std)", ylab="totla tools",
     col=rangi2, pch=ifelse(dat$cid == 1, 1, 16), lwd=2,
     ylim=c(0, 75), cex=1+normalize(k))
# Set up horizontal axis values to compare predictions at
ns <- 100
P_seq <- seq(from=min(dat$P)-0.15, to=max(dat$P)+0.15, length.out=ns)

# Predictions for low contact
lambda <- link(m11.10, data=tibble(P=P_seq, cid=1))
lmu <- colMeans(lambda)
lci <- apply(lambda, 2, PI)
lines(P_seq, lmu, lty=2, lwd=1.5)
shade(lci, P_seq, xpd=TRUE)

# Predictions for high contact
lambda <- link(m11.10, data=tibble(P=P_seq, cid=2))
lmu <- colMeans(lambda)
lci <- apply(lambda, 2, PI)
lines(P_seq, lmu, lty=1, lwd=1.5)
shade(lci, P_seq, xpd=TRUE)

# Note how the model thinks that high contact regions are likely to have lower tools at higher pop
# This isn't logical at all! It's just that there's no data to say otherwise and no constraint on 
# the model
# Hawaii is also a massive outlier in Population and is a low-contact island

# Scientific model!
dat2 <- list(T=d$total_tools, P=d$population, cid=d$contact_id)

# Here have no log link, and because not using linear model there is no intercept, i.e. 0 
# population = 0 tools
# Can see this as if P is 0 then lambda is zero as there's no additive intercept
m11.11 <- ulam(
    alist(
        T ~ dpois(lambda),
        lambda <- exp(a[cid]) * P^b[cid]/g,
        a[cid] ~ dnorm(1, 1),
        b[cid] ~ dexp(1),
        g ~ dexp(1)
    ), data=dat2, chains=4, log_lik=TRUE
)

# Slightly better WAIC for the scientific model but not by much
# But the predictions make much more sense!
# Both in terms of individually, but having zero intercept is more realistic
compare(m11.10, m11.11)

postcheck(m11.10)
postcheck(m11.11)

# 11.2.3 Example: Exposure and the offset ---------------------------------
# Sample counts over 30 days at a rate of 1.5 per day
num_days <- 30
y <- rpois(num_days, 1.5)
# Second monastery gives weekly totals, so 4 weeks at 3.5 manuscripts per week
# How can we adjust to compare on same scale?
num_weeks <- 4
y_new <- rpois(num_weeks, 3.5)
# To compare in same model need to add the log of exposure (time delta)
y_all <- c(y, y_new)
exposure <- c(rep(1, num_days), rep(7, num_weeks))  # Units are days
monastery <- c(rep(0, num_days), rep(1, num_weeks))
d <- tibble(y=y_all, days=exposure, monastery=monastery)
d$log_days <- log(d$days)

# NB: Using index variable here so need 1 intercept + index variable
# Otherwise would use 2 indicator variable intercepts for each monastery
m11.12 <- quap(
    alist(
        y ~ dpois(lambda),
        log(lambda) <- log_days + a + b*monastery,
        a ~ dnorm(0, 1),
        b ~ dnorm(0, 1)
    ), data=d
)

# To compute posterior of lambda in each monastery (which is what we want), sample from posterior
# and use the linear model but without the offset
# Don't use offset again as predictions are on daily scale
# See, the offset isn't returned from the model as it's observed
precis(m11.12)
post <- extract.samples(m11.12)
lambda_old <- exp(post$a)
lambda_new <- exp(post$a + post$b)

# So the rate of writing is far quicker in the weekly monastery
# Don't really need a model to tell this, as 3.5 books a week is 0.5 books a day,
# which is < 1.5 books a day from the first monastery
precis(data.frame(lambda_old, lambda_new))


# Multinomial -------------------------------------------------------------
# 11.3.1 Predictors matched to outcomes -----------------------------------
# career choice example with simulated data
N <- 500              # number of individuals
income <- c(1, 2, 5)  # expected income of each of the 3 career options
score <- 0.5 * income # scores for each career, based on income
p <- softmax(score[1], score[2], score[3]) # Convert scores to prob

# Simulate choice
career <- rep(NA, N)
set.seed(34302) 
for (i in 1:N) {
    career[i] <- sample(1:3, size=1, prob=p)
}
# Obviously have most career 3s
table(career)

# I don't really get this example. Surely it would make more sense to simulate an 
# evenish distribution of the 3 career choices and have individual level incomes (which
# are correlated with career).
# Here, the distribution of careers is correlated with income, and it's hard to understand
# why the 'income' (that maps income of each career) is known? What's the point of passing
# these values in? Instead shouldn't we have the known incomes of all N individuals and
# we have to estimate these 3 values that map career choice and income?
code_m11.13 <- "
data {
    int N;  // Number of individuals
    int K;  // Number of possible careers
    int career[N];  // outcome
    vector[K] career_income;  // Number of individuals
}

parameters {
    vector[K-1] a;  // Intercepts
    real<lower=0> b; // Association of income with choice
}

model {
    vector[K] p;  // Linear predictor post link
    vector[K] s;  // Linear predictor pre link
    a ~ normal(0, 1);
    b ~ normal(0, 0.5);
    s[1] = a[1] + b*career_income[1];
    s[2] = a[2] + b*career_income[2];
    s[3] = 0; // pivot NB: Don't understand what this is doing?! Surely should be a free parameter?
    p = softmax(s);
    career ~ categorical(p);
}
"

dat_list <- list(N=N, K=3, career=career, career_income=income)
# This is using the stan function from rstan
m11.13 <- stan(model_code=code_m11.13, data=dat_list, chains=4)
precis(m11.13, depth=2)

# Hard to interpret coefs on original scale so need to look at predictions
post <- extract.samples(m11.13)

# Setup logit scores
# Here showing counterfactual example where income of second career doubles
s1 <- with(post, a[, 1] + b*income[1])
s2_orig <- with(post, a[, 2] + b*income[2])
s2_new <- with(post, a[, 2] + b*income[2]*2)

# calculate probabilities
p_orig <- sapply(1:length(post$b), function(i) {
    softmax(c(s1[i], s2_orig[i], 0))   # Remember models are fitted with third value as 0
})
p_new <- sapply(1:length(post$b), function(i) {
    softmax(c(s1[i], s2_new[i], 0))
})
# Calculate contrast
p_dif <- p_new[2, ] - p_orig[2, ]
# 5% increase in probability of choosing the career when income is doubled
precis(p_dif)

# 11.3.2 Predictors matched to observations -------------------------------
# Now consider where each outcome has a different lp form with different predictors
# Still modelling career choice but now using each person's family income
N <- 500              # number of individuals
family_income <- runif(N)
b <- c(-2, 0, 2)      # Different coef for each event type.
career <- rep(NA, N)
for (i in 1:N) {
    score <- 0.5 * (1:3) + b*family_income[i]
    p <- softmax(score[1], score[2], score[3])
    career[i] <- sample(1:3, size=1, prob=p)
}

code_m11.14 <- "
data {
    int N;  // Number of individuals
    int K;  // Number of possible careers
    int career[N];  // outcome
    real family_income[N];  // Number of individuals
}

parameters {
    vector[K-1] a;  // Intercepts
    vector[K-1] b;  // coefficients on income
}

model {
    vector[K] p;  // Linear predictor post link
    vector[K] s;  // Linear predictor pre link
    a ~ normal(0, 1.5);
    b ~ normal(0, 1);
    for (i in 1:N) {
        for (j in 1:K-1) {
            s[j] = a[j] + b[j] * family_income[i];
        }
        s[K] = 0; // pivot
        p = softmax(s);
        career[i] ~ categorical(p);
    }
}
"

dat_list <- list(N=N, K=3, career=career, family_income=family_income)
m11.14 <- stan(model_code=code_m11.14, data=dat_list, chains=4)
# Again need to convert to probs to understand
precis(m11.14, depth=2)

# Can't run sim on stanfit!
#sim(m11.14, data=list(N=5, K=3, family_income=runif(5)))
# Or link! Guess that's why we manually built predictions from posterior in previous example
#link(m11.14, data=list(N=5, K=3, family_income=runif(5)))

# 11.3.3 Multinomial in disguise as Poisson -------------------------------
# Final way to fit multinomial is as a series of Poisson likelihoods
d <- UCBadmit

# Binomial model of admissions with single intercept
m_binom <- quap(
    alist(
        admit ~ dbinom(applications, p),
        logit(p) <- a,
        a ~ dnorm(0, 1.5)
    ),
    data=d
)

# Poisson model of admission and rejection rate
dat <- list(admit=d$admit, rej=d$reject)
m_pois <- ulam(
    alist(
        admit ~ dpois(lambda1),
        rej ~ dpois(lambda2),
        log(lambda1) <- a1,
        log(lambda2) <- a2,
        c(a1, a2) ~ dnorm(0, 1.5)
    ), data=dat, chains=3, cores=3
)

# What is the overall acceptance rate?
# From binomial model: (using posterior mean for simplicity)
inv_logit(coef(m_binom))

# From Poisson we get the same answer!
k <- coef(m_pois)
exp(k['a1']) / (exp(k['a1']) + exp(k['a2']))

# This example was for binary, but same principle holds for multinomial
# Can be computationally quicker to run


# Problems ----------------------------------------------------------------
# 11E1: 
# prob to log-odds is logit
logit(0.35)
# NB: logit(0.5) = 0
logit(0.5)

# 11E2: Inverse logit is log odds to prob
inv_logit(3.2)

# 11E3: Exponentiating gives proportional odds!
# So a five-fold increase, which is ginormous
exp(1.7)

# 11E4: Offsets are used if the rows contain counts over different time/spatial-scales
# effectively to normalise them

# 11M1: Good question... something to do with how the Binomial likelihood behaves with
# multiple trials? Ideally would write it out

# 11M2: 1.7 implies an exp(1.7) = 5.47 multiplicative change in outcome

# 11M3: Because we want the p parameter of a Binomial to be [0, 1]

# 11M4: Because we want the rate parameter of a Poisson to be > 0

# 11M5: We have a counting process where the rate can't be > 1.
# However, this doesn't mean we never observe counts > 1
# So where would the rate be constrained < 1?

# 11M6: Binomial have max entropy when each outcome must have 1 of 2
# outcomes, and probability is constant (with respect to predictor vars)
# Poisson is when output is discrete positive value with constant rate

# 11M7: 
m11.4_quap <- quap(
    alist(
        pulled_left ~ dbinom(1, p),
        logit(p) <- a[actor] + b[treatment],  # NB: varying intercept for each actor
        a[actor] ~ dnorm(0, 1.5),
        b[treatment] ~ dnorm(0, 0.5)
    ), data=data_list
)
# There aren't really any differences... I'm not sure what the question is expecting
plot(coeftab(m11.4, m11.4_quap))

m11.4_wide_ulam <- ulam(
    alist(
        pulled_left ~ dbinom(1, p),
        logit(p) <- a[actor] + b[treatment],  # NB: varying intercept for each actor
        a[actor] ~ dnorm(0, 1.5),
        b[treatment] ~ dnorm(0, 0.5)
    ), data=data_list, chains=4, cores=4
)

m11.4_wide_quap <- quap(
    alist(
        pulled_left ~ dbinom(1, p),
        logit(p) <- a[actor] + b[treatment],  # NB: varying intercept for each actor
        a[actor] ~ dnorm(0, 1.5),
        b[treatment] ~ dnorm(0, 0.5)
    ), data=data_list
)

# Still no real differences...
plot(coeftab(m11.4_wide_quap, m11.4_wide_ulam))

# 11M8
# Dropping Hawaii and refitting
d <- Kline |> filter(culture != 'Hawaii')
d$P <- scale(log(d$population))
d$contact_id <- ifelse(d$contact == 'high', 2, 1)

# Going to fit the interaction model & intercept only (1 intercept, not varying by contact id)
dat <- list(
    T = d$total_tools,
    P = d$P,
    cid = d$contact_id
)
m11.10_sans_hawaii <- ulam(
    alist(
        T ~ dpois(lambda),
        log(lambda) <- a[cid] + b[cid] * P,
        a[cid] ~ dnorm(3, 0.5),
        b[cid] ~ dnorm(0, 0.2)
    ), data=dat, chains=4, log_lik=TRUE
)

# Can't compare directly to the original m11.10 as have different outcomes
k <- PSIS(m11.10_sans_hawaii, pointwise=TRUE)$k

plot_kline_posterior <- function(mod) {
    plot(dat$P, dat$T, xlab="Log pop (std)", ylab="totla tools",
         col=rangi2, pch=ifelse(dat$cid == 1, 1, 16), lwd=2,
         ylim=c(0, 75), cex=1+normalize(k))
    # Set up horizontal axis values to compare predictions at
    ns <- 100
    P_seq <- seq(from=min(dat$P)-0.15, to=max(dat$P)+0.15, length.out=ns)
    
    # Predictions for low contact
    lambda <- link(mod, data=tibble(P=P_seq, cid=1))
    lmu <- colMeans(lambda)
    lci <- apply(lambda, 2, PI)
    lines(P_seq, lmu, lty=2, lwd=1.5)
    shade(lci, P_seq, xpd=TRUE)
    
    # Predictions for high contact
    lambda <- link(mod, data=tibble(P=P_seq, cid=2))
    lmu <- colMeans(lambda)
    lci <- apply(lambda, 2, PI)
    lines(P_seq, lmu, lty=1, lwd=1.5)
    shade(lci, P_seq, xpd=TRUE)
}

plot_kline_posterior(m11.10_sans_hawaii)
plot_kline_posterior(m11.10)

# Low contact looks better, but high still has outliers!
# At least high contact is always predicted to be higher tool usage than 
# low tool countries now and the low-contact line is a lot flatter

# Still have Tonga and Yap with high k
d |>
    mutate(k=k) |>
    arrange(desc(k))

# It looks like these were always present to some extent (Tonga was highest with 0.92)
PSIS(m11.10, pointwise=TRUE)$k

# Actually looks ok? Hawaii was within the PIs before
# Case 6 (Trobriand) still poorly predicted
postcheck(m11.10)
postcheck(m11.10_sans_hawaii)

# Now see the scientific model sans Hawaii
dat2 <- list(T=d$total_tools, P=d$population, cid=d$contact_id)

m11.11_sans_hawaii <- ulam(
    alist(
        T ~ dpois(lambda),
        lambda <- exp(a[cid]) * P^b[cid]/g,
        a[cid] ~ dnorm(1, 1),
        b[cid] ~ dexp(1),
        g ~ dexp(1)
    ), data=dat2, chains=4, log_lik=TRUE
)

# The scientific model is again better by WAIC
compare(m11.10_sans_hawaii, m11.11_sans_hawaii)

# 11H1: Use PSIS to compare chimpanzee model (m11.4) with simpler models
# m11.1: intercept only
# m11.2: intercept + treatment (wide prior)
# m11.3: intercept + treatment (tighter prior)
# m11.4: varying intercept + treatment (tighter prior)
# The model with actor specific intercept fares far better
# Implies that each actor has a very different baseline probability of pulling
# left lever
compare(m11.1, m11.2, m11.3, m11.4_quap)

# 11H2: 
library(MASS)
data("eagles")
eagles <- as_tibble(eagles)
# Each row is a record of different victim-pirate pairs
# 5 columns:
#   - y: number of successful attempts
#   - n: total number of attempts
#   - P: whether pirate was large or not
#   - V: whether victim was large
#   - A: whether pirate was adult
eagles
# Convert factors into dummy vars
eagles <- eagles |>
    mutate(
        P = as.integer(factor(P, levels=c('S', 'L')))-1,
        V = as.integer(factor(V, levels=c('S', 'L')))-1,
        A = as.integer(factor(A, levels=c('I', 'A')))-1
    )

# Fit model with indicator variables first (so baseline intercept)
m11h2_quap <- quap(
    alist(
        y ~ dbinom(n, p),
        logit(p) <- a + bP*P + bV*V + bA*A,
        a ~ dnorm(0, 1.5),
        c(bP, bV, bA) ~ dnorm(0, 0.5)
    ), data=eagles
)

m11h2_ulam <- ulam(
    alist(
        y ~ dbinom(n, p),
        logit(p) <- a + bP*P + bV*V + bA*A,
        a ~ dnorm(0, 1.5),
        c(bP, bV, bA) ~ dnorm(0, 0.5)
    ), data=eagles, chains=4, cores=4
)

# a) Is the quadratic approximation OK? Yep agrees very well with ulam!
plot(coeftab(m11h2_quap, m11h2_ulam))

# b) Interpret estimates: If pirate is big, victim is small, and pirate is adult
# are all associated with higher chance of success
# plot posterior predictions. want to see the predicted probability of
# success with its 89% interval for each row, as well as the predicted success
# count and its 89% interval
# The first part is asking for the link value, and the second for sim
mu <- link(m11h2_quap)
mu_mean <- colMeans(mu)
mu_PI <- apply(mu, 2, PI)

preds <- sim(m11h2_quap)
preds_mean <- colMeans(preds)
preds_PI <- apply(preds, 2, PI)

# Now how to plot?
# Rather than plot for each row I'd probably sweep across all permutations of inputs
# But plotting as requested and we don't seem to have enough variability in the model
# I.e. cases 2&4 had a perfect success rate, with p > the PI for these cases, and
# the actual number of values just on the top of the interval (this should be rarer!)
# Likewise case 7 had a 0% success rate which was far below the p PI, and even the 
# outcome y was well below the PI
# What's the difference in information between p and y?
# p is the average prediction, while y allows for the variance in the Binomial
# How you would use their information content differently I'm not sure
# I guess you'd want to get the mean as close as possible (low bias) rather than
# try to optimise for y (low variance)
eagles2 <- eagles |>
    mutate(
        p_raw = y / n,
        p_mod = mu_mean,
        p_lower = mu_PI[1, ],
        p_upper = mu_PI[2, ],
        y_mod = preds_mean,
        y_lower = preds_PI[1, ],
        y_upper = preds_PI[2, ],
        i = row_number()
    ) |>
    rename(y_raw = y) |>
    pivot_longer(c(starts_with("y_"), starts_with("p_")),
                 names_pattern = "([yp])_([a-z]+)",
                 names_to=c("var", "type")) |>
    pivot_wider(names_from=type, values_from=value) 

eagles2 |>
    ggplot(aes(x=i)) +
        geom_point(aes(y=raw), colour="orange", size=2) +
        geom_point(aes(y=mod), colour="steelblue", size=2) +
        geom_errorbar(aes(ymin=lower, y=mod, ymax=upper), colour="steelblue") +
        facet_wrap(~var, scales="free", ncol=2)

# Let's see if this is to do with any of the inputs
# So just plotting the same data but rearranged to see if we can see a pattern 
# I will just use p
# NB: The dataset has every 8 combination of predictors so can just all plot together
# This suggests that an interaction is in order:
# It looks like Pirate has a smaller effect when the victim is Big, as would expect
# (thus interaction?)
# It also seems like Victim size and Adult have an interaction too, as when
# have a small adult but infant predator have 25%-100% success rate depending
# on whether predator is big
# But when have adult predator success rates jump to near 100% regardless of predator size
eagles2 |>
    filter(var == 'p') |>
    mutate(
        A = factor(A, levels=c(0, 1), labels=c("Young", "Adult")),
        P = factor(P, levels=c(0, 1), labels=c("Small Pirate", "Big Pirate")),
        V = factor(V, levels=c(1, 0), labels=c("Big Victim", "Small Victim"))
    ) |>
    ggplot(aes(x=P)) +
        geom_point(aes(y=raw), colour="orange", size=2) +
        geom_point(aes(y=mod), colour="steelblue", size=2) +
        geom_errorbar(aes(ymin=lower, y=mod, ymax=upper), colour="steelblue") +
        facet_grid(A~V)

# c) Been asked to model interaction between pirate size and age, which is not
# the interaction I was going for
# I will actually try all 2-way interactions and see which is most predictive
m11h2_AV <- quap(
    alist(
        y ~ dbinom(n, p),
        logit(p) <- a + bP*P + bV*V + bA*A + bAV*A*V,
        a ~ dnorm(0, 1.5),
        c(bP, bV, bA, bAV) ~ dnorm(0, 0.5)
    ), data=eagles
)
m11h2_AP <- quap(
    alist(
        y ~ dbinom(n, p),
        logit(p) <- a + bP*P + bV*V + bA*A + bAP*A*P,
        a ~ dnorm(0, 1.5),
        c(bP, bV, bA, bAP) ~ dnorm(0, 0.5)
    ), data=eagles
)
m11h2_VP <- quap(
    alist(
        y ~ dbinom(n, p),
        logit(p) <- a + bP*P + bV*V + bA*A + bVP*V*P,
        a ~ dnorm(0, 1.5),
        c(bP, bV, bA, bVP) ~ dnorm(0, 0.5)
    ), data=eagles
)

# This gives highest weight to the Size-Victim model, which is what I predicted
# from my interpretation. Then the Victim-Predator size interaction is next, and the
# suggested interaction is last.
compare(m11h2_quap, m11h2_AP, m11h2_AV, m11h2_VP)

# But just thinking about Pirate size and age, this means that if a pirate is big,
# it has the same likelihood of success whether it's an adult or immature

# 11H3
data("salamanders")
salamanders <- as_tibble(salamanders)
# Have 47 rows where each row is a count from a different area
# Each row has 4 columns:
#  - SITE (id)
#  - SALAMAN (count)
#  - PCTCOVER (continuous)
#  - FORESTAGE (continuous)
# Model using quap
# What priors?
# Firstly, what transforms?
# Keep outcome the same
# But what about pctcover?
# Bi-modal (no log) between 0-100
# Will keep with natural 0, so 0-1
dens(salamanders$PCTCOVER)
salamanders$P <- salamanders$PCTCOVER / max(salamanders$PCTCOVER)
# Max number of salamanders is 13, min is 0, with median of 1
# So need low rate!
summary(salamanders)

# So for priors then, hard to really intuit out like for a normal
# So instead just plot and get a feel for that way
# For the intercept, don't want more than 20 counts, or 50 at push 
# Start with N(0, 1)
foo <- exp(rnorm(1e3, 0, 1))
# Very tight and very large tail!
dens(foo)
# Half the values are < 1, and 94% are < 5
# Which actually could be ok in this dataset
mean(foo < 1)
mean(foo < 5)
# But perhaps could do with more variance as only 1% of values are > 10
mean(foo > 10)
# And max is only 20
max(foo)

# Setting sd=2 has added too much weight to the tails!
# Has given us 273 as max counts, no way we'd ever get that many from our study!
foo <- exp(rnorm(1e3, 0, 2))
dens(foo)
mean(foo < 1)
mean(foo < 5)
mean(foo > 10)
max(foo)

# With 1.5 it's still giving far too high outliers
foo <- exp(rnorm(1e3, 0, 1.5))
dens(foo)
mean(foo < 1)
mean(foo < 5)
mean(foo > 10)
max(foo)

# I've gone for 1.2, it has slightly too high tails (would allow for some crazy outliers)
# but allows for several predictions over 10
# Ideally I think this would be a zero-inflated model as there's a lot of zeros
foo <- exp(rnorm(1e3, 0, 1.2))
dens(foo)
mean(foo < 1)
mean(foo < 5)
mean(foo > 10)
max(foo)

# Now to identify prior for slope
plot_prior_sal <- function(a_mu, a_sigma, b_mu, b_sigma, N=100) {
    a <- rnorm(N, a_mu, a_sigma)
    b <- rnorm(N, b_mu, b_sigma)
    plot(NULL, xlim=c(0, 2), ylim=c(0, 100))
    for (i in 1:N) {
        curve(exp(a[i] + b[i] * x), add=TRUE, col=grau())
    }
}
# Start off with N(0, 1)
# Allows for some crazy responses
plot_prior_sal(0, 1.2, 0, 1)

# N(0, 0.5)
# That's more sedate, that will do!
plot_prior_sal(0, 1.2, 0, 0.5)

# Now to fit the model
m11h3_quap <- quap(
    alist(
        SALAMAN ~ dpois(lambda),
        log(lambda) <- a + bP * P,
        a ~ dnorm(0, 1.2),
        bP ~ dnorm(0, 0.5)
    ), data=salamanders |> dplyr::select(SALAMAN, P)
)

m11h3_ulam <- ulam(
    alist(
        SALAMAN ~ dpois(lambda),
        log(lambda) <- a + b * P,
        a ~ dnorm(0, 1.2),
        b ~ dnorm(0, 0.5)
    ), data=salamanders |> dplyr::select(SALAMAN, P),
    chains=4, cores=4
)
# Yep quap is OK
# So as forest cover increases, get higher chance of finding salamanders
plot(coeftab(m11h3_quap, m11h3_ulam))

# a) Expected counts and 89% vs P
ndata <- tibble(P=seq(0, 1, length.out=100))
mu <- link(m11h3_quap, data=ndata)
mu_mean <- colMeans(mu)
mu_PI <- apply(mu, 2, PI)
preds <- sim(m11h3_quap, data=ndata)
preds_mean <- colMeans(preds)
preds_PI <- apply(preds, 2, PI)
# Massive uncertainty!
# Also why isn't the line continuous? I guess just due to sampling from Poisson?
# If plotted vs link would be smooth function as that's deterministic
# Just a poor fit in general!
# Looks very bimodal as saw in the original counts
# If have < 75% cover then it's basically one model with its own intercept
# And above 75% cover I'd use a different rate, although there's a fair amount of 
# variance there too
ndata |>
    mutate(
        mu=mu_mean,
        mu_lower = mu_PI[1, ],
        mu_upper = mu_PI[2, ],
        pred=preds_mean,
        pred_lower = preds_PI[1, ],
        pred_upper = preds_PI[2, ]
    ) |>
    ggplot(aes(x=P)) +
        geom_ribbon(aes(ymin=pred_lower, ymax=pred_upper), alpha=0.3) +
        geom_ribbon(aes(ymin=mu_lower, ymax=mu_upper), alpha=0.3) +
        geom_line(aes(y=mu)) +
        geom_point(aes(y=SALAMAN), data=salamanders)

# b) Now add FORESTAGE
# This looks skewed so will do log then scale
dens(salamanders$FORESTAGE)
# Just need to add 1 as had a zero reading
dens(scale(log(salamanders$FORESTAGE+1)))
salamanders$F <- scale(log(salamanders$FORESTAGE+1))

# Will use the same priors
# Maybe widen them a bit
m11h3_quap_b <- quap(
    alist(
        SALAMAN ~ dpois(lambda),
        log(lambda) <- a + bP * P + bF * F,
        a ~ dnorm(0, 1.2),
        c(bP, bF) ~ dnorm(0, 0.5)
    ), data=salamanders |> dplyr::select(SALAMAN, P, F)
)
# Forest age has minimal difference
# Looks like having a forest makes a bigger difference than its age
plot(coeftab(m11h3_quap, m11h3_quap_b))

m11h3_quap_c <- quap(
    alist(
        SALAMAN ~ dpois(lambda),
        log(lambda) <- a + bP * P + bF * F + bPF*P*F,
        a ~ dnorm(0, 1.2),
        c(bP, bF, bPF) ~ dnorm(0, 0.5)
    ), data=salamanders |> dplyr::select(SALAMAN, P, F)
)

# Is there an interaction?
# I.e. if have older forests with a lot of coverage do we see a lot more salamanders?
# Not really, negative if anything!
plot(coeftab(m11h3_quap, m11h3_quap_b, m11h3_quap_c))
# Yep simplest model favoured here, age of forest doesn't seem to impact much
compare(m11h3_quap, m11h3_quap_b, m11h3_quap_c)

# Try a stratified model
salamanders <- salamanders |>
    mutate(coverage=as.integer(PCTCOVER > 0.75)+1)
m11h3_quap_d <- quap(
    alist(
        SALAMAN ~ dpois(lambda),
        log(lambda) <- a[coverage] + bP[coverage] * P,
        a[coverage] ~ dnorm(0, 1.5),
        bP[coverage] ~ dnorm(0, 0.5)
    ), data=salamanders |> dplyr::select(SALAMAN, coverage, P=PCTCOVER)
)
precis(m11h3_quap_d, depth=2)
# Rates of 1 below and 1.3 above, so barely any difference.
# I find that surprising, I was expecting a much higher slope for the second part
exp(0.03)

# Ah, this is the most predictive model by WAIC though
compare(m11h3_quap, m11h3_quap_b, m11h3_quap_c, m11h3_quap_d)

# Does it look a better fit at all?
# Nope, basically flat!
# So how does it have better WAIC?
ndata <- ndata |>
    mutate(coverage=as.integer(P > 0.75)+1)
mu <- link(m11h3_quap_d, data=ndata)
mu_mean <- colMeans(mu)
mu_PI <- apply(mu, 2, PI)
preds <- sim(m11h3_quap_d, data=ndata)
preds_mean <- colMeans(preds)
preds_PI <- apply(preds, 2, PI)
ndata |>
    mutate(
        mu=mu_mean,
        mu_lower = mu_PI[1, ],
        mu_upper = mu_PI[2, ],
        pred=preds_mean,
        pred_lower = preds_PI[1, ],
        pred_upper = preds_PI[2, ]
    ) |>
    ggplot(aes(x=P)) +
        geom_ribbon(aes(ymin=pred_lower, ymax=pred_upper), alpha=0.3) +
        geom_ribbon(aes(ymin=mu_lower, ymax=mu_upper), alpha=0.3) +
        geom_line(aes(y=mu)) +
        geom_point(aes(y=SALAMAN), data=salamanders)

# 11H4
data("NWOGrants")
NWOGrants <- NWOGrants |> as_tibble()
NWOGrants |>
    head()
# Have dept/discipline, gender, applications, and number of awards
# "What are the total and indirect causal effects of gender on grant awards?"
# Consider mediation path (aka pipe) through discipline
dag <- dagitty("dag { 
  awards <- gender -> discipline -> awards;
}")
drawdag(dag)
# Can estimate direct effects as A ~ G + D,
# Total effects as A ~ G
# Then indirect effect is the difference
adjustmentSets(dag, exposure="gender", outcome="awards", effect="direct")
adjustmentSets(dag, exposure="gender", outcome="awards", effect="total")

NWOGrants <- NWOGrants |>
    mutate(
        gid = as.numeric(factor(gender, levels=c('m', 'f'))),
        did = as.numeric(discipline)
    )

m11h4_direct <- quap(
    alist(
        awards ~ dbinom(applications, p),
        logit(p) <- a[did] + b[gid],
        a[did] ~ dnorm(0, 1.5),  # Using same priors as UBCadmit
        b[gid] ~ dnorm(0, 1.5)
    ),
    data=NWOGrants
)
# Quite high variance between disciplines
plot(precis(m11h4_direct, depth=2))
# Let's look at gender contrasts
post_direct <- extract.samples(m11h4_direct)
# Slightly higher acceptance rate for men, although 2.5% probability difference
# in absolute terms is minute
gender_diff_direct <- inv_logit(post_direct$b[, 1]) - inv_logit(post_direct$b[, 2])
plot(precis(tibble(direct=gender_diff_direct)))
# That's the direct effect anyway!

# For indirect we need the total as well as the direct, so need to model
# A ~ G.
m11h4_total <- quap(
    alist(
        awards ~ dbinom(applications, p),
        logit(p) <- b[gid],
        b[gid] ~ dnorm(0, 1.5)
    ),
    data=NWOGrants
)
post_total <- extract.samples(m11h4_total)
# So will calculate the gender contrast for each model separately, then take
# the difference to get the indirect effect
gender_diff_total <- inv_logit(post_total$b[, 1]) - inv_logit(post_total$b[, 2])
# Total difference is ~3% too... that seems wrong?!
plot(precis(tibble(total=gender_diff_total)))

# Indirect effect is thus the difference
# And is ~ 0.5%. So most of the way in which gender affects awards being granted
# is directly, rather than through the choice of discipline
# However, this difference is low in absolute terms at 2.5%
plot(precis(tibble(indirect=gender_diff_total - gender_diff_direct)))

# M11H5
dag2 <- dagitty("dag { 
  career[u];
  awards <- gender -> discipline -> awards;
  career -> discipline;
  career -> awards;
}")
coordinates(dag2) <- list(x=c(gender=0, awards=1, discipline=1, career=2),
                          y=c(gender=0, discipline=0, awards=1, career=0))

# When you condition in discipline it's a collider so it opens the backdoor path through
# career stage to awards
drawdag(dag2)

# 11H6
data("Primates301")
Primates301 <- as_tibble(Primates301)
# 301 rows where each row is species
# Want to investigate how brain size is associatd with social learning
Primates301

# a) Model number of observations of social learning as function of log brain size
# Will need to z-scale
dens(log(Primates301$brain))
Primates301$log_brain <- as.numeric(scale(log(Primates301$brain)))

# Very long tail!
dens(Primates301$social_learning)
summary(Primates301$social_learning)

# Just use weakly informative priors
# Remember that for the intercept, half the values are below 1 if have
# normal centered on 0, which seems reasonable for this data
# It shouldn't matter too much here since we have so much data

# Here's the possible priors using N(0, 1.5) for intercept and
# N(0, 1) for slope
# Can see will easily get some very high counts if needed
plot_prior(0, 1.5, 0, 1)

summary(Primates301)
m11h6 <- quap(
    alist(
        social_learning ~ dpois(lambda),
        log(lambda) <- a + bBrain * log_brain,
        a ~ dnorm(0, 1.5),
        bBrain ~ dnorm(0, 1)
    ), data=Primates301 |> filter(!is.na(log_brain), !is.na(social_learning))
)
# Fits ok
# Very high coefficient on brain!
# Intercept: When we're at an average log-brain size
# The expected number of social learnings is exp(-1.3), or 0.27
# But for every SD increase in log brain size we get 20 fold increase in # social learnings
# So the range of social learnings is around 0 for an average or below brain size,
# but then increases up to around 2*20*0.27 (10)
plot(precis(m11h6))

# b) Add log research effect
dens(log(Primates301$research_effort))
Primates301$log_research <- as.numeric(scale(log(Primates301$research_effort)))

m11h6_b <- quap(
    alist(
        social_learning ~ dpois(lambda),
        log(lambda) <- a + bBrain * log_brain + bResearch * log_research,
        a ~ dnorm(0, 1.5),
        c(bBrain, bResearch) ~ dnorm(0, 1)
    ), data=Primates301 |> filter(!is.na(log_brain), !is.na(social_learning),
                                  !is.na(log_research))
)
# Ah that makes a huge difference in the brain parameter
# So now brain has a far lower influence on the number of social learnings
# But of course the amount of research interest should play a role in the number of 
# social learnings observed. What we really want to model is the PROPORTION of 
# social learnings per research effort, rather than number
plot(coeftab(m11h6, m11h6_b))

# c)
# I would say B -> R, (as researchers are more interested in primates with big brains)
#             B -> L 
#             R -> L
# For the latter 2 there is evidence after conditioning on the other, that each variable
# does have a direct causal effect on L (although R is stronger by far) 
# Now to quickly test B->R
# For priors, max rise/run is 4/4 = 1 (since using both z-scores)
# And a reasonable prior on slope is could be either way so center on 0,
# and this would allow for 1 SD change having 2 SD change in output
# Again, the data is big enough it won't make a huge difference
m11h6_c <- quap(
    alist(
        log_research ~ dnorm(mu, sigma),
        mu <- a + b * log_brain,
        a ~ dnorm(0, 1),
        b ~ dnorm(0, 1),
        sigma ~ dexp(1)
    ), data=Primates301 |> filter(!is.na(log_research), !is.na(log_brain))
)
# Yep, there is a moderate (0.4) association between brain size and research
# So I'll stick with the mediation model
plot(precis(m11h6_c))

     