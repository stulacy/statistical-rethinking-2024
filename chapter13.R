library(tidyverse)
library(rethinking)
library(patchwork)

# 13.1 Multilevel tadpoles ------------------------------------------------

data(reedfrogs)
d <- reedfrogs |> as_tibble()
d

# Standard model - one intercept for each tank, but all having the same
# prior
d$tank <- 1:nrow(d)
m13.1 <- ulam(
    alist(
        S ~ dbinom(N, p),
        logit(p) <- a[tank],
        a[tank] ~ dnorm(0, 1.5)
    ), data=d |> select(S=surv, N=density, tank),
    chains=4, log_lik=TRUE, cores=4
)

# Hierarchical prior on intercept
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

# Despite having more parameters, the hierarchical model has a better WAIC
# Notice how it has a lower pWAIC (i.e. # effective parameters) than the
# unpooled model, despite having more parameters
compare(m13.1, m13.2)

# The hierarchical prior is (mean values) N(1.35, 1.62)
# This is a regularizing prior as it's quite strong, but NB:
# this prior has been LEARNED FROM THE DATA rather than hardcoded
# This results in more regularization than the unpooled model and thus
# more effective parameters
precis(m13.2, pars = c("a_bar", "sigma"))

# Plot posterior predictions
post <- extract.samples(m13.2)
# Compute mean intercept for each tank, and convert to probability scale
d$propsurv.est <- logistic(apply(post$a, 2, mean))

# Display raw survival proportions
plot(d$propsurv, ylim=c(0, 1), pch=16, xaxt="n", xlab="tank", ylab="proportion survival",
     col=rangi2)
axis(1, at=c(1, 16, 32, 48), labels=c(1, 16, 32, 48))

# Overlay posterior means
points(d$propsurv.est)

# Mark posterior mean probability across tanks
abline(h=mean(inv_logit(post$a_bar)), lty=2)

# Vertical dividers
abline(v=16.5, lwd=0.5)
abline(v=32.5, lwd=0.5)
text(8, 0, "small tanks")
text(16+8, 0, "medium tanks")
text(32+8, 0, "large tanks")

# NB: 
#   - 1. More shrinkage in smaller tanks (i.e. where have less data)
#   - 2. More shrinkage the further an observed value is from a_bar

# Plot posterior
# One way of doing this is plotting the density of the hierarchical prior
# And NOT just plotting the density of this prior, but the uncertainty in
# an individual intercept drawn from it
# I.e. drawing X new tanks intercept, but rather than just doing a single
# draw to get the intercept, we'll use the full distribution and see the density
x_new <- seq(-3, 4, length.out=1e3)
map_dfr(1:100, function(i) {
    a_dens <- dnorm(x_new, post$a_bar[i], post$sigma[i])
    tibble(i=i, x=x_new, dens=a_dens)
}) |>
    ggplot(aes(x=x, y=dens, group=i)) +
        geom_line() +
        theme_classic()

# Reminder: since we have the prior distribution of a_i ~ N(a_bar, sigma)
# To see the distribution of a_i we could either analytically look at the density
# across a range of x-values, or we could look at the density of a random sample
# drawn from the distribution. Note the former is far quicker!
plot(seq(0, 20, length.out=1e3), dnorm(seq(0, 20, length.out=1e3), 10, 3))
dens(rnorm(1e3, 10, 3))

# Now take the second approach and draw straight from the hierarchical prior distribution and
# take the density
sim_tanks <- rnorm(8000, post$a_bar, post$sigma)
# This gives us a similar shape to the first plot
dens(sim_tanks, lwd=2, adj=0.1)
# Can also plot this on a probability scale, which somehow turns
# into a cumulative density curve
# How is there that much density at prob=0.99?
# Must be the non-linear nature of this transform
dens(inv_logit(sim_tanks), lwd=2, adj=0.1)

# Yep if we tweak the original plot to show results on the probability scale
# we get the same
map_dfr(1:100, function(i) {
    a_dens <- dnorm(x_new, post$a_bar[i], post$sigma[i])
    tibble(i=i, x=x_new, dens=a_dens)
}) |>
    mutate(x_prob=inv_logit(x)) |>
    ggplot(aes(x=x_prob, y=dens, group=i)) +
        geom_line() +
        theme_classic()

# So of these approaches to visualising the posterior, which is most useful?
# I like the first approach as it's simple to generate PIs through quantiles
# While the second doesn't show uncertainty since we're using the density
# to aggregate over our random draws, whereas the first approach uses the analytical
# pdf for each random draw, so can get uncertainty over density

# Varying effects and the underfitting/overfitting trade-off --------------

# Simulating data to see if model can recover hardocded parameters
# Really nice idea for debugging models!
# NB: DON'T INCLUDE PRIORS WHEN SIMULATING DATA (except for their hierarchies, definitely
# want to include a_i ~ N(a_bar, sigma) in our simulated dataset)
# Also need to fix values for these params, rather than distributions
a_bar <- 1.5
sigma <- 1.5
nponds <- 60
# Initial number of tadpoles in each pond
# Have 15 ponds of each 'size'
Ni <- as.integer(rep(c(5, 10, 25, 35), each=15))
set.seed(5005)
# Draw an intercept for each pond
# NB: this is on the LOG ODDS scale, as the link function converts this to prob
# Could do it in prob scale but then the model would retrieve parameters on the prob scale
a_pond <- rnorm(nponds, mean=a_bar, sd=sigma)
dsim <- tibble(pond=1:nponds, Ni=Ni, true_a=a_pond)
# Simulate number of survivors
dsim$Si <- rbinom(nponds, prob=logistic(dsim$true_a), size=dsim$Ni)

# 13.2.4 No pooling -------------------------------------------------------
# No need to fit model here, just proportion within each tank
dsim$p_nopool <- dsim$Si / dsim$Ni

# 13.2.5 Partial pooling --------------------------------------------------
m13.3 <- ulam(
    alist(
        Si ~ dbinom(Ni, p),
        logit(p) <- a_pond[pond],
        a_pond[pond] ~ dnorm(a_bar, sigma),
        a_bar ~ dnorm(0, 1.5),
        sigma ~ dexp(1)
    ), data=dsim |> select(Si, Ni, pond),
    chains=4, cores=4
)
precis(m13.3, depth=2)

# Add mean predicted survival probs to dataset
post <- extract.samples(m13.3)
dsim$p_partiallypooled <- colMeans(inv_logit(post$a_pond))

# Compare to the actual probs (we simulated log-odds)
dsim$prob_a <- inv_logit(dsim$true_a)

# Reproducing Fig 13.3 in ggplot
dsim_plt <- dsim |>
    pivot_longer(starts_with("p_"), names_pattern = "p_(.+)",
                 names_to="model", values_to="predicted") |>
    mutate(
        error = abs(predicted - prob_a)
    ) 

# Voila, the same as Fig 13.3.
# Incidentally I accidentally used the wrong seed and got very different resuls (unpooling fared better
# than partially pooled on more tank sizes)
# Main takeaway isn't that partially pooled is better, but that it does better on smaller tank sizes
dsim_plt |>
    ggplot(aes(x=pond, y=error)) +
        geom_point(aes(colour=model)) +
        geom_hline(aes(yintercept=y, colour=model),
                   data=dsim_plt |> group_by(model, Ni) |> summarise(y=mean(error))) +
        facet_wrap(~Ni, scales="free_x", nrow=1) +
        theme_classic() +
        labs(x="Pond", y="Absolute error") +
        theme(legend.position = "bottom")

# Can pass new data to compiled model and resample without having to recompile
# NB: if using vanilla Cmdstan you get this behaviour for free as it only
# compiles if the binary isn't up-to-date
a <- 1.5
sigma <- 1.5
nponds <- 60
Ni <- as.integer(rep(c(5, 10, 25, 35), each=15))
a_pond <- rnorm(nponds, mean=a, sd=sigma)
dsim <- tibble(pond = 1:nponds, Ni, true_a=a_pond)
dsim$Si <- rbinom(nponds, prob=inv_logit(dsim$true_a), size=dsim$Ni)
dsim$p_npool <- dsim$Si / dsim$Ni
newdat <- dsim |> select(pond, Si, Ni)

# Error! m13.3 doesn't have stanfit object!
#m13.3new <- stan(fit=m13.3@stanfit, data=newdat, chains=4)
# The example then proceeds by extracting the posterior using extract.samples
# and calculating the same predictions / errors as before
# But what is this Stan call doing? 

# 13.3 More than one type of cluster --------------------------------------
data("chimpanzees")
d <- as_tibble(chimpanzees)
d$treatment <- 1 + d$prosoc_left + 2*d$condition

set.seed(13)
m13.4 <- ulam(
    alist(
        pulled_left ~ dbinom(1, p),
        logit(p) <- a[actor] + g[block_id] + b[treatment],
        b[treatment] ~ dnorm(0, 0.5),
        a[actor] ~ dnorm(a_bar, sigma_a),
        # NB: Since adding 2 intercepts, can't add a hyper-parameter
        # for G's mean as would be non-identifiable
        g[block_id] ~ dnorm(0, sigma_g),
        # Hyper-parameters
        a_bar ~ dnorm(0, 1.5),
        sigma_a ~ dexp(1),
        sigma_g ~ dexp(1)
    ), data=d |> 
        mutate(block_id=block, treatment=as.integer(treatment)) |>
        select(pulled_left, actor, block_id, treatment),
    chains=4, cores=4, log_lik=TRUE
)

# Divergent transitions!
# A lot of parameters...
# Rhat is ok, but neff varies, notably sigma_g
precis(m13.4, depth=2)
# A is actor, G is Block, B is treatment
# Block seems to make little effect (both in terms of sigma_a vs sigma_g, and comparing
# the variability amongst a intercepts vs that amongst g intercepts (which are basically all 0))
# A is very varied - has partially pooling shrunk the intercepts a bit compared to
# the old unpooled model?
plot(precis(m13.4, depth=2))

# Now fit sans block and see if makes a difference
set.seed(14)
m13.5 <- ulam(
    alist(
        pulled_left ~ dbinom(1, p),
        logit(p) <- a[actor] + b[treatment],
        b[treatment] ~ dnorm(0, 0.5),
        a[actor] ~ dnorm(a_bar, sigma_a),
        # Hyper-parameters
        a_bar ~ dnorm(0, 1.5),
        sigma_a ~ dexp(1)
    ), data=d |> 
        mutate(treatment=as.integer(treatment)) |>
        select(pulled_left, actor, treatment),
    chains=4, cores=4, log_lik=TRUE
)

# M13.4 had 7 more parameters than m13.5 but only 2 more effective parameters (pWAIC)
# because the block parameters were effectively zero
# Minute WAIC difference
compare(m13.5, m13.4)

# Fit varying intercept on treatment
# NB: Again as this is another intercept can't define a hyper-parameter
# for mu, instead just sigma
set.seed(15)
m13.6 <- ulam(
    alist(
        pulled_left ~ dbinom(1, p),
        logit(p) <- a[actor] + g[block_id] + b[treatment],
        b[treatment] ~ dnorm(0, sigma_b),
        a[actor] ~ dnorm(a_bar, sigma_a),
        # NB: Since adding 2 intercepts, can't add a hyper-parameter
        # for G's mean as would be non-identifiable
        g[block_id] ~ dnorm(0, sigma_g),
        # Hyper-parameters
        a_bar ~ dnorm(0, 1.5),
        sigma_a ~ dexp(1),
        sigma_g ~ dexp(1),
        sigma_b ~ dexp(1)
    ), data=d |> 
        mutate(block_id=block, treatment=as.integer(treatment)) |>
        select(pulled_left, actor, block_id, treatment),
    chains=4, cores=4, log_lik=TRUE
)

# Very little difference in the treatment (b) estimates using either model
# This is because there is enough data in each model that partially pooling doesn't
# do much shrinkage (remember the tadpoles in tanks example)
coeftab(m13.4, m13.6)

# Can see sigma_b is also quite small - much smaller than actor sigma (a) although
# is larger than the block sigma (g)
plot(precis(m13.6))

# Very little infomration difference between them
# WAIC favours the simpler model, purely because there's little predictive difference
# but a fair difference in the number of parameters
compare(m13.4, m13.5, m13.6)

# 13.4 Divergent transitions and non-centered priors ----------------------
# Divergent samples are rejected, which isn't an issue in itself, but it means
# the posterior in that region is hard to sample, which could be emblematic of
# a wider problem
# Reparameterisation is best solution

# Having scale of one parameter depend on another is a fundamental part of multi-level
# modelling and can lead to these divergent transitions
m13.7 <- ulam(
    alist(
        v ~ normal(0, 3),
        x ~ normal(0, exp(v))
    ), data=list(N=1), chains=4, cores=4
)
# Notice the awful rhat and neff values!
precis(m13.7)

# Traceplot for x is awful too
traceplot(m13.7)

# Centered parameterisation: where distribution of one var is conditional on one or more other params
# e.g. v ~ N(0, 3), x ~ N(0, exp(v))
# Non-centered is the opposite, where a variable's distribution isn't conditional
# on any other var
# e.g. v ~ N(0, 3), z ~ N(0, 1), x = z*exp(v)
# Can think of this as z being the z-transform of x that we're going to sample
# and now want to convert it back to the original scale so apply the usual
# back-transform (z * sd) + mu
# NB: Use Generated Quantities to get back the parameter on the scale we're interested in
m13.7_nc <- ulam(
    alist(
        v ~ normal(0, 3),
        z ~ normal(0, 1),
        gq > real[1]: x <<- z*exp(v)
    ), data=list(N=1), chains=4, cores=4
)
# Much better sampling!
# And can see the expected values for x
precis(m13.7_nc)

# Can also try increasing the adapt_delta, which is a parameter that controls the 
# acceptance rate. Stan tunes the step size during the warm up to meet this, so by increasing it
# will get a smaller more cautious step size (at the cost of increased computational time)
set.seed(13)
m13.4b <- ulam(m13.4, chains=4, cores=4, control=list(adapt_delta=0.99))

# We've decreased the # of divergent transitions!
divergent(m13.4)
divergent(m13.4b)

# But still have relatively low neff and some poor Rhat
precis(m13.4b, depth=2)

# So fit a non-centered parameterisation
# Basically replace the individual level intercepts with a_bar + sigma_a * a_z
# Where a_z ~ N(0, 1)
# Then in generated quantities recreate the separate alpha intercepts
# NB: This step is optional, only do it if want to compare the intercepts
set.seed(13)
m13.4_nc <- ulam(
    alist(
        pulled_left ~ dbinom(1, p),
        logit(p) <- a_bar + a_z[actor]*sigma_a + g_z[block_id]*sigma_g + b[treatment],
        b[treatment] ~ dnorm(0, 0.5),
        a_z[actor] ~ dnorm(0, 1),
        g_z[block_id] ~ dnorm(0, 1),
        # Hyper-parameters are unchanged, since these aren't centered
        # i.e. functions of other parameters
        a_bar ~ dnorm(0, 1.5),
        sigma_a ~ dexp(1),
        sigma_g ~ dexp(1),
        gq> vector[actor]:a <<- a_bar + a_z*sigma_a,
        gq> vector[block_id]:g <<- g_z*sigma_g
    ), data=d |> 
        mutate(block_id=block, treatment=as.integer(treatment)) |>
        select(pulled_left, actor, block_id, treatment),
    chains=4, cores=4, log_lik=TRUE
)

# Much better n_eff now!
# NB: Generated Quantities also get an n_eff and rhat
precis(m13.4_nc, depth=2)

# Can plot n_eff directly
# Note that n_eff is better in the non-centered parameterisation for every parameter bar one!
models <- list("Centered"=m13.4, "Non-centered"=m13.4_nc)
map_dfr(models, function(mod) {
    precis <- precis(mod, depth=2)
    tibble(var = rownames(precis), 
           n_eff = precis$ess_bulk)
}, .id="model") |>
    pivot_wider(names_from=model, values_from=n_eff) |>
    ggplot(aes(x=Centered, y=`Non-centered`)) +
        geom_point() +
        geom_abline(intercept=0, slope=1) +
        theme_classic() +
        ylim(0, 2000) +
        xlim(0, 2000) 
        
# And that parameter which has a higher n_eff in the centered
# model is a_bar
map_dfr(models, function(mod) {
    precis <- precis(mod, depth=2)
    tibble(var = rownames(precis), 
           n_eff = precis$ess_bulk)
}, .id="model") |>
    group_by(var) |>
    dplyr::top_n(1, n_eff) |>
    filter(model == 'Centered')

# Also note that **NON-CENTERED ISN'T ALWAYS BETTER**
# It tends to help more where there is either:
#  - A cluster with low variation (blocks in m13.4)
#  - A cluster with a large number of values, but each value doesn't have many observations
# Otherwise there is no guarantee so it's worth trying both if a model is struggling

# Also also note can non-center non-Gaussians too
# x ~ Exp(lambda)  # Centered
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# z ~ Exp(1)       # Non-centered
# x = z*lambda

# 13.5 Multilevel posterior predictions -----------------------------------

# By design, multi-level models will never fully retrodict the training sample
# as the multi-level prior is 'adaptive regularization', so instead should fare better out-of-sample

# 2 Predictive situations:
#   - Predicting for units with the same clusters as in the training data
#   - Predicting for NEW clusters

# 13.5.1 Posterior prediction for same clusters ---------------------------
# Just use the unit level varying slope/intercept paramters, don't need to
# touch the group level parameter
chimp <- 2
d_pred <- list(
    actor=rep(chimp, 4),
    treatment=1:4,
    block_id=rep(1, 4)
)

# M13.4 is model with varying intercept for actor, block, but unpooled treatment
# Can either use link to get predictions
p <- link(m13.4, data=d_pred)
p_mu <- colMeans(p)
p_pi <- apply(p, 2, PI)

# Or manually from the posterior parameters
post <- extract.samples(m13.4)
# There are 7 actors
dim(post$a)

my_link <- function(treatment, actor=1, block_id=1) {
    # So just rebuild the linear predictor from the samples
    logodds <- post$a[, actor] + post$g[, block_id] + post$b[, treatment]
    # Return on the original scale
    inv_logit(logodds)
}

# Calculate predictions for every treatment, the same as using rethinking::link above
p_raw <- sapply(1:4, function(i) my_link(i, actor=2, block_id=1))
p_mu <- colMeans(p_raw)
p_pi <- apply(p_raw, 2, PI)

# 13.5.2 Posterior prediction for new clusters ----------------------------
# Generate predictions for a new unobserved _average_ actor
p_link_abar <- function(treatment) {
    # Deliberately ignoring block as we have a new unobserved actor
    # who doesn't have a block, and we'll omit it to keep the example
    # more concise. But it doesn't really matter since the block effect is ~=0
    logodds <- post$a_bar + post$b[, treatment]
    inv_logit(logodds)
}

p_raw <- sapply(1:4, function(i) p_link_abar(i))
p_mu <- colMeans(p_raw)
p_pi <- apply(p_raw, 2, PI)

p1 <- tibble(
    treatment=1:4,
    mu=p_mu,
    lower=p_pi[1, ],
    upper=p_pi[2, ]
) |>
    ggplot(aes(x=treatment, y=mu)) +
        geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
        scale_x_continuous(labels=c("R/N", "L/N", "R/P", "L/P")) +
        geom_point() +
        labs(x="Treatment", y="Proportion pulled left", title="Average actor") +
        geom_line() +
        ylim(0, 1) +
        theme_classic()
p1

# Now want to show the variation between actors, rather than just an 'average' actor
# Going to draw new actors from the hierarchical prior
a_sim <- rnorm(length(post$a_bar), post$a_bar, post$sigma_a)
# Now have 2,000 different intercepts
p_link_asim <- function(treatment) {
    logodds <- a_sim + post$b[, treatment]
    inv_logit(logodds)
}

p_raw_asim <- sapply(1:4, p_link_asim)
p_mu_asim <- colMeans(p_raw_asim)
p_pi_asim <- apply(p_raw_asim, 2, PI)

p2 <- tibble(
    treatment=1:4,
    mu=p_mu_asim,
    lower=p_pi_asim[1, ],
    upper=p_pi_asim[2, ]
) |>
    ggplot(aes(x=treatment, y=mu)) +
        geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
        scale_x_continuous(labels=c("R/N", "L/N", "R/P", "L/P")) +
        geom_point() +
        labs(x="Treatment", y="Proportion pulled left", title="Marginal of actor") +
        ylim(0, 1) +
        geom_line() +
        theme_classic()

# Much larger uncertainty when account for the variation between actors
# NB: p2 is termed the 'marginal of actor' as it averages over the uncertainty amongst actors
# They show different things and have different uses
# p1 shows the treatment effect for an 'average person' to understand how the treatment works
# Whereas p2 shows the variability amongst chimpanzees
p1 + p2

# Final option is to visualize both the treatment effect and variability amongst actors
# Although in my mind this is just a more cluttered version of p2...
# But rather than summarise the 2,000 simulated actors into a mean + PI, display a subset
df <- as_tibble(as.matrix(p_raw_asim))
colnames(df) <- paste0('T', 1:4)
df$sim <- 1:nrow(df)

p3 <- df |>
    pivot_longer(-sim, names_pattern="T([1-4])", names_to="treatment",
                 values_to="prop") |>
    mutate(treatment = as.integer(treatment)) |>
    filter(sim <= 100) |>  # HAVE TO LIMIT TO 100 ROWS ELSE TOO MUCH NOISE!
    ggplot(aes(x=treatment, y=prop, group=sim)) +
        scale_x_continuous(labels=c("R/N", "L/N", "R/P", "L/P")) +
        geom_line() +
        labs(x="Treatment", y="Proportion pulled left", title="Marginal of actor") +
        ylim(0, 1) +
        theme_classic()
# How is p3 useful? We lose data (can only show 100 regression lines) and it's cluttered
# Hard to tell where the bulk of the density is
# The book claims that this plot makes it easier to see the zig-zag treatment effect (already seen
# in ps 1 & 2 imo), and the variation amongst actors (already shown in p2 imo)
# Finally, one thing p3 does show is that for actors with very large intercepts, treatment has 
# very little effect, also for some reason those near the bottom
# I'd be interested to understand why this is the case...
p1 + p2 + p3
