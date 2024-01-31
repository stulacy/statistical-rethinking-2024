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

# TODO reproduce this plot myself rather than using base R for the data wrangling + plotting
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
mean(foo < 1)

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
