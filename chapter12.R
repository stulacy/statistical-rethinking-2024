library(rethinking)
library(tidyverse)
library(dagitty)

# 12.1 Overdispersed Counts -----------------------------------------------


# 12.1.1 Beta-binomial ----------------------------------------------------

pbar <- 0.5
theta <- 5
# Density around 0.5 with theta=5
# NB: dbeta2 is from rethinking
curve(dbeta2(x, pbar, theta), from=0, to=1, xlab="Probability", ylab="Density")

# And it seems to be parameterised differently to the dbeta from base R!
curve(dbeta(x, pbar, theta), from=0, to=1, xlab="Probability", ylab="Density")

# When theta is 2 is flat
curve(dbeta2(x, pbar, 2), from=0, to=1, xlab="Probability", ylab="Density")

# Theta < 1 becomes bimodal at extremes
curve(dbeta2(x, pbar, 1), from=0, to=1, xlab="Probability", ylab="Density")

# Theta >> 2 becomes far steeper around pbar
curve(dbeta2(x, pbar, 50), from=0, to=1, xlab="Probability", ylab="Density")

# When used in a model the Beta-Binomial has the affect of assigning a different
# probability to each row. Could likely achieve similar by manually specifying a
# linear predictor for p as a function of row-specific variables
data(UCBadmit)
d <- UCBadmit
d$gid <- ifelse(d$applicant.gender == 'male', 1L, 2L)
dat <- list(A=d$admit, N=d$applications, gid=d$gid)

# The transformed parameters means minimum theta is 2
# So no bimodal prior!
m12.1 <- ulam(
    alist(
        A ~ dbetabinom(N, pbar, theta),
        logit(pbar) <- a[gid],
        a[gid] ~ dnorm(0, 1.5),
        transpars> theta <<- phi+2.0,
        phi ~ dexp(1)
    ), data=dat, chains=4, cores=4
)

# Compare gender contrast
# Very weak gender effect
# Now that allow each dept to have different p, we can see there is little gender
# effect after all, similar to what we got before when we explicitly included dept
# as a predictor
post <- extract.samples(m12.1)
post$da <- post$a[, 1] <- post$a[, 2]
precis(post, depth=2)

# We can see the distribution over p for given input variables
gid <- 2
curve(dbeta2(x, mean(logistic(post$a[, gid])),
             mean(post$theta)),
      from=0, to=1, ylab="Density", xlab="Probability admit", ylim=c(0, 3), lwd=2)

for (i in 1:50) {
    p <- logistic(post$a[i, gid])
    curve(dbeta2(x, p, theta), add=TRUE, col=col.alpha("black", 0.2))
}
mtext("Distribution of female admission rate")

# Much better fit now, all actual points are within CIs, although this isn't because
# the mean has shifted, but because it can handle additional variance
postcheck(m12.1)

# Compare this to the model that uses dept
dat$dept_id <- rep(1:6, each=2)

m11.8 <- ulam(
    alist(
        A ~ dbinom(N, p),
        logit(p) <- a[gid] + delta[dept_id],
        a[gid] ~ dnorm(0, 1.5),
        delta[dept_id] ~ dnorm(0, 1.5)
    ), 
    data=dat, chains=4, iter=4000  # Bumped up to 4,000 samples!
)

# See here the model that explicitly assigns a different mean prob by dept
# better handles the data!
# So only use over-dispersed model if we are seeing more variance than our model handles
# AND we don't have an explanatory variable
postcheck(m11.8)

# 11.1.2 Negative-binomial (gamma-Poisson) --------------------------------
data("Kline")
d <- Kline
d$P <- standardize(log(d$population))
d$contact_id <- ifelse(d$contact == 'high', 2L, 1L)
dat2 <- list(
    T=d$total_tools,
    P=d$population,
    cid=d$contact_id
)

m12.2 <- ulam(
    alist(
        T ~ dgampois(lambda, phi),
        lambda <- exp(a[cid])*P^b[cid]/g,
        a[cid] ~ dnorm(1, 1),
        b[cid] ~ dexp(1),
        g ~ dexp(1),
        phi ~ dexp(1)
    ),
    data=dat2, chains=4, log_lik=TRUE
)

# Compare to the regular Poisson model
m11.11 <- ulam(
    alist(
        T ~ dpois(lambda),
        lambda <- exp(a[cid]) * P^b[cid]/g,
        a[cid] ~ dnorm(1, 1),
        b[cid] ~ dexp(1),
        g ~ dexp(1)
    ), data=dat2, chains=4, log_lik=TRUE
)

# Less influenced by Hawaii outlier
postcheck(m12.2)
postcheck(m11.11)

# 12.2 Zero-inflated outcomes ---------------------------------------------
# Essentially have logistic regression model of probability of 0, then poisson
# model if non-zero (which can still give zeros ofc)
# Can use same predictors or separate in the 2 linear predictors (the first is for
# the probability of zero in binomial, and second is for lambda in Poisson)

# Simulate data
prob_drink <- 0.2 # 0n 20% of days get zero data as monks are drunk, separate from
                  # the Poisson model of # books written
rate_work <- 1    # average 1 manuscript per day when actually working
N <- 365          # time-scale is days so simulate 1 year

# simulate days monk drink
set.seed(365)
drink <- rbinom(N, 1, prob_drink)
table(drink)

# simulate manuscripts completed 
y <- (1-drink) * rpois(N, rate_work)
hist(y)

# Colour on the histogram the number of days with 0 books written due to drinking
simplehist(y)
zeros_drink <- sum(drink)
zeros_work <- sum(y==0 & drink==0)
zeros_total <- sum(y==0)
zeros_drink + zeros_work == zeros_total
lines(c(0, 0), c(zeros_work, zeros_total), lwd=4, col=rangi2)

m12.3 <- ulam(
    alist(
        y ~ dzipois(p, lambda),
        logit(p) <- ap,
        log(lambda) <- al,
        ap ~ dnorm(-1.5, 1),  # Prior on less < 0.5 chance that monks are drinking
        al ~ dnorm(1, 0.5)
    ),
    data=list(y=y), chains=4, cores=4
)

# Compare to a regular Poisson model
m12.3_pois <- ulam(
    alist(
        y ~ dpois(lambda),
        log(lambda) <- al,
        al ~ dnorm(1, 0.5)
    ),
    data=list(y=y), chains=4, cores=4
)

# Coefficients are a bit awkward due to the log scale
# But basically exp(0) = 1, so the zero-inflated model
# got the writing rate correct, and inv_logit(-1.31) ~= 21% so
# it also correctly estimated the monks' drinking habits (pun intended)
# The pure Poisson model estimates the writing rate at exp(-0.24) = 0.79,
# undershooting it because of the inflated zeros
coeftab(m12.3, m12.3_pois)

# 12.3 Ordered categorical outcomes ---------------------------------------
# Use *CUMULATIVE LINK FUNCTION*
# A long long long time ago I used a continuous Linear Predictor with latent
# thresholds separating categories
data("Trolley")
d <- Trolley
# 9,930 rows with 12 columns. data from 331 individuals, so have 30 conditions
# for each person
# Response is 1-7 ordinal with 1 being least morally permissible, and 7 being most
summary(d)

# Here are counts of outcomes (note the modal answer is the middle one!)
simplehist(d$response, xlim=c(1, 7), xlab="response")

# Want to get log-cumulative-odds
# Firstly get proportion of each response
pr_k <- table(d$response) / nrow(d)
# Now get cumulative probs
cum_pr_k <- cumsum(pr_k)

plot(1:7, cum_pr_k, type="b", xlab="response",
     ylab="cumulative proportion", ylim=c(0, 1))

# Finally conversion to log-cumulative-odds
# Want threshold intercepts for each outcome (similar to what I had in my 
# football modelling all those years ago, but I used dummy encoding rather than index)
logit <- function(x) log(x / (1-x))
lco <- logit(cum_pr_k)
lco
# Note have infinity for k=7, so therefore just need 6 intercepts, so this is the
# same as my football modelling
# Need to be able to convert between log-cumulative-odds (linear predictor scale
# to guarantee ordering) and probability scale (scale needed for likelihood)

# Richard's written a function to handle the likelihood for you so you don't need to
# write out the likelihood in terms of a categorical dist of K probabilities
# which are formed from K-1 thresholds and any external predictors
m12.4 <- ulam(
    alist(
        R ~ dordlogit(0, cutpoints),  # Later will replace 0 with predictors
        cutpoints ~ dnorm(0, 1.5)
    ), 
    data=list(R=d$response), chains=4, cores=4
)

# QUAP needs the starting values to show the order of the thresholds
m12.4_b <- quap(
    alist(
        R ~ dordlogit(0, cutpoints),  # Later will replace 0 with predictors
        cutpoints ~ dnorm(0, 1.5)
    ), 
    data=list(R=d$response), start=list(cutpoints=c(-2, -1, 0, 1, 2, 2.5))
)
# Both algorithms found the same values and note the ordering looks decent
# NB: Estimates are on log-cum-odds scale
precis(m12.4, depth=2)
precis(m12.4_b, depth=2)

# Can return back to cumulative probs using inverse logit
inv_logit(coef(m12.4))
# And these almost perfectly match the observed cumulative probs
cum_pr_k

# 12.3.3 Adding predictor variables ---------------------------------------
# Define logcum-odds of each response as a sum of its threshold + some linear factors
# In particular, define log_cum_odds(k) = alpha_k - theta_i
# theta_i = beta*x_i
# NB: the minus sign on theta_i is used because an increase in x_i with positive beta
# means the cumulative probability at this threshold is lower, so it places greater
# probability mass at higher thresholds.
# +ve Beta having this interpretation as "increases outcome probability" is more
# intuitive than the opposite!
pk <- dordlogit(1:7, 0, coef(m12.4))
# As an example, here are the posterior probabilities from the previous model
round(pk, 2)

# Mean value of 4.2
sum(pk * 1:7)

# If we were to subtract 0.5 from each coefficient
pk2 <- dordlogit(1:7, 0, coef(m12.4)-0.5)
# The probabilities below the average decrease and those at
# the top increase
rbind(
    round(pk, 2),
    round(pk2, 2)
)

# The mean value is now higher
sum(pk2 * 1:7)

dat <- list(
    R = d$response,
    A = d$action,
    I = d$intention,
    C = d$contact
)

m12.5 <- ulam(
    alist(
        R ~ dordlogit(phi, cutpoints),
        phi <- bA*A + bC*C + BI*I,  # I personally find the fully explicit model clearer than separating it like this
        BI <- bI + bIA*A + bIC*C,
        c(bA, bI, bC, bIA, bIC) ~ dnorm(0, 0.5),
        cutpoints ~ dnorm(0, 1.5)
    ),
    data=dat, chains=4, cores=4
)
# All the slopes are negative, so if a trolley problem features any of these it is considered to be
# less morally justifiable
precis(m12.5, depth=2)

# The interaction between Intent and Contact has the biggest effect by far!
# This is a bit surprising as Intent and Contact on their own don't make a huge difference
plot(precis(m12.5))

# Plotting posterior predictions requires a bit of thought
plot_vals <- function(A, C, I) {
    plot(NULL, type="n", xlab="Intention", ylab="Probability", 
         xlim=c(0, 1), ylim=c(0, 1), xaxp=c(0, 1, 1), yaxp=c(0, 1, 2))
    # Just want to compare Intention 0/1 (and thus still its interaction with Contact)
    # So hold the other values constant
    ndata <- expand_grid(A=A, C=C, I=I)
    phi <- link(m12.5, data=ndata)
    # Can get posterior values for BI (which is really just a transformed parameter rather than a link)
    # or phi
    names(phi)
    phi <- phi$phi
    # Phi has 2 columns for the 2 input values (intention = 0, intention = 1)
    dim(phi)
    # The first col is all zeros, so I'm guessing we just want the second
    summary(phi)
    
    # Plot first 50 predictions
    # Take care to compute cum prob for each outcome
    # Didn't necessarily need to use link above, could have derived phi from posterior
    # This doesn't look helpful at all...
    post <- extract.samples(m12.5)
    for (s in 1:50) {
        pk <- pordlogit(1:6, phi[s,], post$cutpoints[s, ])
        for (i in 1:6) lines(0:1, pk[, i], col=grau(0.1))
    }
    mtext(sprintf("A=%d C=%d", A, C))
}

# Note the steepest change is when switching from intention 0 to 1 when C = 1, hence the interaction's magnitude
plot_vals(0, 0, 0:1)
plot_vals(1, 0, 0:1)
plot_vals(0, 1, 0:1)

# Admittedly I don't find this plot that counter-intuitive, I'd rather have the 7 outcomes on the x-axis and plot
# using facets w/ colours
plot_vals_gg <- function(A, C, I) {
    ndata <- expand_grid(A=A, C=C, I=I)
    phi <- link(m12.5, data=ndata)
    phi <- phi$phi
    # Plot first 50 predictions
    # Take care to compute cum prob for each outcome
    # Didn't necessarily need to use link above, could have derived phi from posterior
    # This doesn't look helpful at all...
    post <- extract.samples(m12.5)
    res <- map_dfr(1:50, function(s) {
        pk <- pordlogit(1:6, phi[s,], post$cutpoints[s, ])
        colnames(pk) <- paste0("O", 1:6)
        pk <- as_tibble(pk)
        ndata |>
            cbind(pk)
    }, .id="Sample") |>
        pivot_longer(matches("O[1-6]"), names_pattern="O([1-6])", names_to="Outcome", values_to="Prob")
    res |>
        ggplot(aes(x=Outcome, y=Prob, colour=as.factor(I))) +
            geom_boxplot() +
            facet_grid(A~C, labeller = label_both)
}
# I prefer this style of plot personally
plot_vals_gg(A=0:1, C=0:1, I=0:1)

# An alternative that I also prefer to the line plots is to look at the simulated outcomes
# Here we can see that with I=0 have far fewer 1 responses (least morally permissible), but with I=1 (light blue)
# the proportion of 1 responses goes up drastically
ndata <- expand_grid(A=0, C=1, I=0:1)
s <- sim(m12.5, data=ndata)
simplehist(s, xlab="Response")

# 12.4 Ordered categorical predictors -------------------------------------
# When have ordered categorical inputs want the model to know the order without
# having to encode as numeric
# The 'edu' field contains (ordered!) information about educational status
# How would we adjust for this in a model of the Trolley problem?
summary(d$edu)

# I'd never recommend relevelling a factor implicitly like this, but if it works...
edu_levels <- c(6, 1, 8, 4, 7, 2, 5, 3)
d$edu_new <- edu_levels[d$edu]
summary(d$edu_new)

# For 8 levels will have 7 parameters (as the parameter represents the difference in the outcome
# from a change from the lower educational status to the current one)
# Therefore the lowest status will be absorbed into the intercept

# Basically will model 1 effect (the effect of the full predictor at level K), with K-1 fractions
# The fractions will be random parameters that sum to 1, so use a DIRICHLET distribution

# Will use a weak prior to get familiar with Dirichlet
# Setting alpha = 2, 2, 2, 2, 2, 2, 2
library(gtools)
set.seed(1805)
delta <- rdirichlet(10, alpha=rep(2, 7))
# Get 10 reps of 7 probs
dim(delta)
# Each row sums to 1
rowSums(delta)

# Can plot!
h <- 3
plot(NULL, xlim=c(1, 7), ylim=c(0, 0.4), xlab="index", ylab="probability")
for (i in 1:nrow(delta)) {
    lines(1:7, delta[i, ], type="b",
          pch=ifelse(i==h, 16, 1), lwd=ifelse(i==h, 4, 1.5),
          col=ifelse(i==h, "black", col.alpha("black", 0.7)))
}

# However, remember that we are setting delta_1 = 0, as the first category should be
# absorbed into the intercept
dat$E <- as.integer(d$edu_new)
dat$alpha <- rep(2, 7) # Delta prior

# Much slower to sample! Combination of complex model with 10k points
# TODO: Why is the simplex needed, if dirichlet always returns values that sum to 1?
# Also this is the first time we've just used normal rather than dnorm
m12.6 <- ulam(
    alist(
        R ~ ordered_logistic(phi, kappa),
        phi <- bE*sum(delta_j[1:E]) + bA*A + bI*I + bC*C,
        kappa ~ normal(0, 1.5),
        c(bA, bI, bC, bE) ~ normal(0, 1),
        vector[8]: delta_j <<- append_row(0, delta),
        simplex[7]: delta ~ dirichlet(alpha)
    ),
    data=dat, chains=4, cores=4
)
# Here bE has a NEGATIVE effect overall, so more education results in a tendency to
# give morally less justifiable scores
# But note that the effect of this is far lower than bC, bI, or bA (yet alone bIC which we didn't include)
precis(m12.6, depth=2, omit="kappa")

delta_labels <- c("Elem", "MidSch", "SHS", "HSG", "SCol", "Bach", "Mast", "Grad")
# Ignoring the pars argument and displaying all!
# But looking at values in book can see all negatively correlated
# Because they have the constraint of summing to 1
# Only SCol has very little additional effect
pairs(m12.6, pars="delta", labels=delta_labels)

# And can see the Stancode itself for fully understanding what this does
stancode(m12.6)

# Compare to model with education as continuous (rescaled to 0-1)
dat$edu_norm <- normalize(d$edu_new)
m12.7 <- ulam(
    alist(
        R ~ ordered_logistic(mu, cutpoints),
        # Where's intercept!
        mu <- bE * edu_norm + bA * A + bI * I + bC*C,
        c(bA, bI, bC, bE) ~ normal(0, 1),
        cutpoints ~ normal(0, 1.5)
    ), data=dat, chains=4, cores=4
)
# Much weaker effect than when split categories out!
# Could be because effect isn't linear
precis(m12.7)

# So what is the nature of the effect then if it isn't linear?
# Piecewise-linear?
# Let's look:
# NB: can't use the link function as it doesn't understand the 'append_row' Stan function
# So will have to build up link from posterior
plot_vals_gg_m126 <- function(A, C, I, E) {
    ndata <- expand_grid(A=A, C=C, I=I, E=E)
    post <- extract.samples(m12.6)
    res <- pmap_dfr(ndata, function(A, C, I, E) {
        map_dfr(1:50, function(s) {
            delta_j <- c(0, post$delta[s, ])
            phi <- post$bE[s]*sum(delta_j[1:E]) + post$bA[s]*A + post$bI[s]*I + post$bC[s]*C
            pk <- pordlogit(1:6, phi, post$kappa[s, ])
            tibble(Outcome=1:6, Prob=pk, A=A, C=C, I=I, E=E)
        }, .id="Sample") 
    })
    # Want to visualise E on the x-axis to see effect of changing E
    # So on the Y will be each outcome, somehow with the full distribution represented
    res |>
        ggplot(aes(x=as.factor(E), y=Prob, colour=as.factor(Outcome))) +
            geom_boxplot() +
            stat_smooth() +
            facet_wrap(~Outcome, scales="free")
}
# So here can see that there is indeed a relatively non-linear effect
# Noticeably in the outcome category 6, where the educational status
# effect drastically tails off
plot_vals_gg_m126(0, 1, 1, 1:8)


# Problems ----------------------------------------------------------------

# 12E1
# Ordered = results in football match (W/D/L), unordered = species of dog

# 12E2
# Ordered logistic regression uses log-cumulative-odds, so we're comparing the
# difference in odds between successive values rather than their absolute value

# 12E3
# It leads to a lower estimate of the rate, because you don't account for the fact
# that some zeros are generated separately and not as part of the process you
# are interested in

# 12E4
# The daily concentration of NO2 molecules, as it has significant variation owing
# to (amongst other things) the amount of sunlight
# Underdispersion is when there is little variance I guess? Man City's points per games?
# Annoyingly constant at 3...

# 12M1
x <- c(12, 36, 7, 41)
# Convert to proportions/probs
x_probs <- x / sum(x)
# Convert to cumulative proportions
x_cumprobs <- cumsum(x_probs)
# Convert to log-odds
logit <- function(x) log(x/(1-x))
x_cumlogodds <- logit(x_cumprobs)
x_cumlogodds

# 12M2
# Calculate likelihoods as the difference between
# successive thresholds
# Isn't this just doing the inverse of the cumulative sum?
# Pr(yi=k) = Pr(yi <= k) - Pr(yi <= k-1)
x_cumprobs - lag(x_cumprobs)

tibble(
    group=1:4,
    prob=x_probs,
    cum_prob = x_cumprobs
) |>
    mutate(cum_prob_start = lag(cum_prob, default=0)) |>
    ggplot(aes(x=group, y=cum_prob)) +
        geom_col(alpha=0.3, width=0.1) +
        geom_segment(aes(yend=cum_prob_start, xend=group), colour="blue", 
                     linewidth=2) +
        geom_point() +
        geom_line() +
        theme_minimal() +
        ylim(c(0, 1))

# 12M3
# Modify derivation of ZIP to get ZIB (binomial)
# p_zero = prob of a zero, p_process = prob of a 1 from the modelled process
# Pr(0|p_zero, p_process) = Pr(0 | p_zero) + (1 - Pr(0 | p_zero)) * Pr(0 | p_process)
# = p_zero + (1 - p_zero) * (1 - p_process)
# Pr(1|p_zero, p_process) = (1 - Pr(0 | p_zero)) * Pr(1 | p_process)
# = (1 - p_zero) * p_process
# Then define likelihood as this with the 2 parameters
# and would define both parameters in terms of co-vars using logit link

# Things to ask Richard: --------------------------------------------------
# Why is the simplex needed in m12.6, if dirichlet always returns values that sum to 1?

# 12H1
data("Hurricanes")
Hurricanes <- as_tibble(Hurricanes)
# 92 hurricans with their names, strength, damage, deaths, 
# and how feminine the names are (?!)
# Going to look at modelling death (why not damage?) as function of feminity
Hurricanes

# Very bi-modal!
# Will z-score
dens(Hurricanes$femininity)
Hurricanes$F <- as.numeric(scale(Hurricanes$femininity))

m12h1_intercept <- quap(
    alist(
        deaths ~ dpois(lambda),
        log(lambda) <- a,
        a ~ dnorm(0, 1.5)
    ), data=Hurricanes
)

m12h1 <- quap(
    alist(
        deaths ~ dpois(lambda),
        log(lambda) <- a + bF * F,
        a ~ dnorm(0, 1.5),
        bF ~ dnorm(0, 1)
    ), data=Hurricanes
)

# Adding feminity does seem to improve the model by the weight
# Although the dSE is massive!
compare(m12h1, m12h1_intercept)

# The effect is 0.24, so a 1 sd increase in feminity increases the death rate
# by exp(0.24) = 1.27
# So for a mean feminity, the death rate is exp(3) = 20, which rises to 25.4
# with 1 SD of feminity, not a massive effect, but very precisely estimated (very
# narrow PIs)
precis(m12h1)

mu <- link(m12h1)
mu_mean <- colMeans(mu)
mu_PI <- apply(mu, 2, PI)

preds <- sim(m12h1)
preds_PI <- apply(preds, 2, PI)

Hurricanes$name_year <- paste(Hurricanes$name, Hurricanes$year, sep="_")
names_ord <- Hurricanes |> arrange(deaths) |> pull(name_year)

# So we're just predicting the average case, i.e. can't handle overdispersion in the data
# And this overdispersion is either random, or isn't caused by feminity
# So it predicts well Hurricanes with a death rate ~ 20 (as this is the sample mean),
# and doesn't do well with the tails
Hurricanes |> 
    arrange(deaths) |> 
    mutate(
        mu_mean = mu_mean,
        mu_lower = mu_PI[1, ],
        mu_upper = mu_PI[2, ],
        preds_lower = preds_PI[1, ],
        preds_upper = preds_PI[2, ]
    ) |>
    mutate(name_year = factor(name_year, levels=names_ord)) |>
    ggplot(aes(y=name_year)) +
        geom_point(aes(x=deaths), colour="orange") +
        geom_point(aes(x=mu_mean), colour="steelblue") +
        geom_errorbarh(aes(xmin=mu_lower, xmax=mu_upper), colour="steelblue") +
        geom_errorbarh(aes(xmin=preds_lower, xmax=preds_upper), colour="darkgreen") +
        theme_minimal()

# Is there a relationship with feminity?
# Not really! Just that the biggest outliers also seem to have feminine names
Hurricanes |> 
    mutate(
        mu_mean = mu_mean,
        mu_lower = mu_PI[1, ],
        mu_upper = mu_PI[2, ],
        preds_lower = preds_PI[1, ],
        preds_upper = preds_PI[2, ]
    ) |>
    mutate(name_year = factor(name_year, levels=names_ord)) |>
    ggplot(aes(y=femininity)) +
        geom_point(aes(x=deaths), colour="orange") +
        geom_point(aes(x=mu_mean), colour="steelblue") +
        geom_errorbarh(aes(xmin=mu_lower, xmax=mu_upper), colour="steelblue") +
        geom_errorbarh(aes(xmin=preds_lower, xmax=preds_upper), colour="darkgreen") +
        theme_minimal()

# I.e. Diane, Camille, Agnes, Sandy. The top 4 by deaths and also high feminine names
Hurricanes |> filter(deaths > 100)

# 12H2
# Over-dispersed poisson
m12h2 <- quap(
    alist(
        deaths ~ dgampois(lambda, phi),
        log(lambda) <- a + bF * F,
        a ~ dnorm(0, 1.5),
        bF ~ dnorm(0, 1),
        phi ~ dexp(1)
    ), data=Hurricanes
)

# DONT USE WAIC/PSIS ON OVER-DISPERESED MODELS
# bF has massively reduced now to 0.03 with an SD of 0.10, so it's effectively 0
precis(m12h2)

# Can think of this as similar to having a lot of outliers in a Gaussian model
# so using Student T instead (not exactly the same as it allows for different variance
# per-observation, rather than general fatter tails, but still the same principle
# of allowing for more variance in the model)

# Visualize predictions
mu_b <- link(m12h2)
mu_mean_b <- colMeans(mu_b)
mu_PI_b <- apply(mu_b, 2, PI)

preds_b <- sim(m12h2)
preds_PI_b <- apply(preds_b, 2, PI)

# Here can see the difference is that by allowing more variance in the outcome
# it reduces the impact of the covariates (I guess unless the covariates are 
# particularly impactful)
Hurricanes |> 
    arrange(deaths) |> 
    mutate(
        mu_mean_pois = mu_mean,
        mu_lower_pois = mu_PI[1, ],
        mu_upper_pois = mu_PI[2, ],
        preds_lower_pois = preds_PI[1, ],
        preds_upper_pois = preds_PI[2, ],
        mu_mean_gampois = mu_mean_b,
        mu_lower_gampois = mu_PI_b[1, ],
        mu_upper_gampois = mu_PI_b[2, ],
        preds_lower_gampois = preds_PI_b[1, ],
        preds_upper_gampois = preds_PI_b[2, ]
    ) |>
    pivot_longer(ends_with("pois"), names_pattern = "(.+)_([gampois]+)",
                 names_to=c("metric", "model")) |> 
    pivot_wider(names_from=metric, values_from=value) |>
    mutate(name_year = factor(name_year, levels=names_ord)) |>
    ggplot(aes(y=name_year)) +
        geom_point(aes(x=deaths), colour="orange") +
        geom_point(aes(x=mu_mean), colour="steelblue") +
        geom_errorbarh(aes(xmin=mu_lower, xmax=mu_upper), colour="steelblue") +
        geom_errorbarh(aes(xmin=preds_lower, xmax=preds_upper), colour="darkgreen") +
        theme_minimal() +
        facet_wrap(~model, ncol=1)

# 12H3
dens(Hurricanes$min_pressure)
dens(Hurricanes$damage_norm)
# Will z-score and z(log) respective
Hurricanes$P <- as.numeric(scale(Hurricanes$min_pressure))
Hurricanes$D <- as.numeric(scale(log(Hurricanes$damage_norm)))

# Look good!
dens(Hurricanes$P)
dens(Hurricanes$D)

# So now want to consider interactions
# Will stick with gamma-pois as have already observed the significant inter-observation
# variance
m12h3_d <- quap(
    alist(
        deaths ~ dgampois(lambda, phi),
        log(lambda) <- a + bF * F + bD * D + bFD * F * D,
        a ~ dnorm(0, 1.5),
        c(bF, bD, bFD) ~ dnorm(0, 1),
        phi ~ dexp(1)
    ), data=Hurricanes
)

m12h3_p <- quap(
    alist(
        deaths ~ dgampois(lambda, phi),
        log(lambda) <- a + bF * F + bP * P + bFP * F * P,
        a ~ dnorm(0, 1.5),
        c(bF, bP, bFP) ~ dnorm(0, 1),
        phi ~ dexp(1)
    ), data=Hurricanes
)

m12h3_dp <- quap(
    alist(
        deaths ~ dgampois(lambda, phi),
        log(lambda) <- a + bF * F + bD * D + bP * P + bFD * D * F + bFP * F * P,
        a ~ dnorm(0, 1.5),
        c(bF, bP, bD, bFD, bFP) ~ dnorm(0, 1),
        phi ~ dexp(1)
    ), data=Hurricanes
)
# I _think_ can use WAIC to compare between gamma-poisson models
# just not between Poisson & Gamma-Poisson
# _If_ that's the case, damage has a greater impact
compare(m12h2, m12h3_d, m12h3_p, m12h3_dp)

# First compare coefs
# (Can't plot as phi stretches the scale too much)
coeftab(m12h2, m12h3_d, m12h3_p, m12h3_dp)
# All of them have the feminity effect very low
# Damage has a reasonable effect, moreso than pressure
# When looking at the full model, only damage is significant really
precis(m12h3_dp)

# What can we infer from the counter-factuals?
plot_counterfactual <- function(mod, newdata) {
    mu <- link(mod, data = newdata)
    mu_mean <- colMeans(mu)
    mu_PI <- apply(mu, 2, PI)
    preds <- sim(mod, data = newdata)
    preds_PI <- apply(preds, 2, PI)
    tibble(newdata) |>
        mutate(mu=mu_mean,
               mu_lower=mu_PI[1, ],
               mu_upper=mu_PI[2, ],
               preds_lower=preds_PI[1, ],
               preds_upper=preds_PI[2, ]
       )
}

# Let's see counter-factuals under all 3 models
ndata_D <- expand_grid(
    F = seq(min(Hurricanes$F)-0.15, max(Hurricanes$F)+0.15),
    D = seq(min(Hurricanes$D)-0.15, max(Hurricanes$D)+0.15),
    P=0
)
ndata_P <- expand_grid(
    F = seq(min(Hurricanes$F)-0.15, max(Hurricanes$F)+0.15),
    P = seq(min(Hurricanes$P)-0.15, max(Hurricanes$P)+0.15),
    D=0
)
ndata_DP <- expand_grid(
    F = seq(min(Hurricanes$F)-0.15, max(Hurricanes$F)+0.15),
    D = seq(min(Hurricanes$D)-0.15, max(Hurricanes$D)+0.15),
    P = seq(min(Hurricanes$P)-0.15, max(Hurricanes$P)+0.15)
)

# Unsure how to visualise this with 3 dependent variables, and on the outcome
# have mean, mean PI, and outcome PIs...
# Just trying for the best model by WAIC
# So this does predict that the highest death tolls are found in Hurricanes with both
# more feminine names AND more damage
# (Isn't it obvious that higher damage -> higher deaths?)
plot_counterfactual(m12h3_d, ndata_D) |>
    ggplot(aes(x=F, y=D, z=mu)) +
        geom_contour_filled()

# How strong is the effect of this interaction?
# It is 'significant'!
precis(m12h3_d)

# How does it compare to a model with D but no interaction?
m12h3_d_only <- quap(
    alist(
        deaths ~ dgampois(lambda, phi),
        log(lambda) <- a + bF * F + bD * D,
        a ~ dnorm(0, 1.5),
        c(bF, bD) ~ dnorm(0, 1),
        phi ~ dexp(1)
    ), data=Hurricanes
)
# WAIC prefers the model without the interaction, although there's not much in it
# Either way, why does this interaction exist?
# It doesn't seem plausible to answer the question!
# That Hurricanes that both do a lot of damage AND are feminine have higher death tolls
# than those that are just one or the other?
# Isn't this simply because the most damaging Hurricanes happen to be Feminine, and so 
# because it's an exponential increase in damage/deaths, adding this multiplicative
# factor in makes the model fit better?
compare(m12h3_d, m12h3_d_only)

# 12H4
# Hah I've already taken the log of damage as there is a clear skew in the data and
# it was very expected that damage is exponential
# Let's do the reverse and see what a linear relationship looks like
# Will z-scale ofc
Hurricanes$D_lin <- as.numeric(scale(Hurricanes$damage_norm))
m12h3_d_linear <- quap(
    alist(
        deaths ~ dgampois(lambda, phi),
        log(lambda) <- a + bF * F + bD * D_lin + bFD * F * D_lin,
        a ~ dnorm(0, 1.5),
        c(bF, bD, bFD) ~ dnorm(0, 1),
        phi ~ dexp(1)
    ), data=Hurricanes
)

# Yep the log model is better as I expected
compare(m12h3_d, m12h3_d_linear)

# Does the linear model have an even stronger interaction effect?
# No, it's weaker
precis(m12h3_d_linear)

# Let's look at counter-factual
ndata_D_lin <- expand_grid(
    F = seq(min(Hurricanes$F)-0.15, max(Hurricanes$F)+0.15),
    D_lin = seq(min(Hurricanes$D_lin)-0.15, max(Hurricanes$D_lin)+0.15),
    P=0
)
# Yep weaker interaction
plot_counterfactual(m12h3_d_linear, ndata_D_lin) |>
    ggplot(aes(x=F, y=D_lin, z=mu)) +
        geom_contour_filled()

# Finally then, what about linear D with no interactions?
m12h3_d_linear_only <- quap(
    alist(
        deaths ~ dgampois(lambda, phi),
        log(lambda) <- a + bF * F + bD * D_lin,
        a ~ dnorm(0, 1.5),
        c(bF, bD) ~ dnorm(0, 1),
        phi ~ dexp(1)
    ), data=Hurricanes
)

# And the 2 log models are preferred, with not much in it whether you want an interaction
# or not
compare(m12h3_d, m12h3_d_only, m12h3_d_linear, m12h3_d_linear_only)

# Let's plot predictions for the log and linear models (no interactions)
mu_log <- link(m12h3_d_only)
mu_mean_log <- colMeans(mu_log)
mu_PI_log <- apply(mu_log, 2, PI)
preds_log <- sim(m12h3_d_only)
preds_PI_log <- apply(preds_log, 2, PI)

mu_lin <- link(m12h3_d_linear_only)
mu_mean_lin <- colMeans(mu_lin)
mu_PI_lin <- apply(mu_lin, 2, PI)
preds_lin <- sim(m12h3_d_linear_only)
preds_PI_lin <- apply(preds_lin, 2, PI)

# Both models are still really bad! Struggle to pick up the number of deaths
# from the most deadly Hurricanes
# Linear does well for Camille but is subtley worse for the others
# While log does better on average, but fails at the top end
Hurricanes |> 
    arrange(deaths) |> 
    mutate(
        mu_mean_log = mu_mean_log,
        mu_lower_log = mu_PI_log[1, ],
        mu_upper_log = mu_PI_log[2, ],
        preds_lower_log = preds_PI_log[1, ],
        preds_upper_log = preds_PI_log[2, ],
        mu_mean_lin = mu_mean_lin,
        mu_lower_lin = mu_PI_lin[1, ],
        mu_upper_lin = mu_PI_lin[2, ],
        preds_lower_lin = preds_PI_lin[1, ],
        preds_upper_lin = preds_PI_lin[2, ]
    ) |>
    mutate(name_year = factor(name_year, levels=names_ord)) |>
    pivot_longer(matches("^(mu|preds)_[a-z]+_(log|lin)$"), names_pattern = "(.+)_(log|lin)",
                 names_to=c("metric", "model")) |> 
    pivot_wider(names_from=metric, values_from=value) |>
    ggplot(aes(y=name_year)) +
        geom_point(aes(x=deaths), colour="orange") +
        geom_point(aes(x=mu_mean), colour="steelblue") +
        geom_errorbarh(aes(xmin=mu_lower, xmax=mu_upper), colour="steelblue") +
        geom_errorbarh(aes(xmin=preds_lower, xmax=preds_upper), colour="darkgreen") +
        theme_minimal() +
        facet_wrap(~model, ncol=1)

# Definitely an extremely steep exponential curve that just wasn't fitted by the data
# As well as huge amount of variance
Hurricanes |>
    ggplot(aes(x=D, y=deaths, colour=F)) +
        geom_point()  +
        scale_colour_viridis_c()

# 12H5
# Do Men and Women have different average tendencies in moral reasoning?
# Hypothesis: Women are more concerned with avoiding harm, while men are more
# concerned with justice
# Using the trolley data, using contact as a proxy for harm, are women differently
# bothered by harm?
data("Trolley")
d <- as_tibble(Trolley)
edu_levels <- c(6, 1, 8, 4, 7, 2, 5, 3)
d$edu_new <- edu_levels[d$edu]
d$E <- as.integer(d$edu_new)
# reminder:
# 9,930 rows with 12 columns. data from 331 individuals, so have 30 conditions
# for each person
# Response is 1-7 ordinal with 1 being least morally permissible, and 7 being most
summary(d)
d

# So want to have a slope and intercept for men/women separately for contact and
# compare the contrasts on this
# (How would draw dag for this?)
# Are there any confounders? Or can assume randomly sampled?
# Will check the 2 obvious ones: age + edu

# Men in this study are younger...
quap(
    alist(
        male ~ dbinom(1, p),
        logit(p) ~ a + b * age,
        a ~ dnorm(0, 1),
        b ~ dnorm(0, 1)
    ), data=d |> mutate(age = age / max(age))
) |> precis()
# And men have more education...
ulam(
    alist(
        male ~ dbinom(1, p),
        logit(p) <- bE*sum(delta_j[1:E]),
        bE ~ normal(0, 1),
        vector[8]: delta_j <<- append_row(0, delta),
        simplex[7]: delta ~ dirichlet(alpha)
    ),
    data=list(E=d$E, male=d$male, alpha=rep(2, 7)), chains=4, cores=4
) |>
    precis()

# So let's control for both of those
# Should Action and Intent be included?
# They should be independent so I don't think it'll matter, might try both ways
m12h5 <- ulam(
    alist(
        R ~ ordered_logistic(phi, kappa),
        phi <- bE*sum(delta_j[1:E]) + bAge * age + bA*A + bI*I + bC[gid]*C,
        kappa ~ normal(0, 1.5),
        bC[gid] ~ normal(0, 1),
        c(bA, bI, bAge, bE) ~ normal(0, 1),
        vector[8]: delta_j <<- append_row(0, delta),
        simplex[7]: delta ~ dirichlet(alpha)
    ),
    data=list(E=d$E, male=d$male, alpha=rep(2, 7),
              age=d$age / max(d$age),
              gid=d$male + 1,
              R=d$response,
              A=d$action,
              I=d$intention,
              C=d$contact), chains=4, cores=4
)


m12h5_b <- ulam(
    alist(
        R ~ ordered_logistic(phi, kappa),
        phi <- bE*sum(delta_j[1:E]) + bAge * age + bC[gid]*C,
        kappa ~ normal(0, 1.5),
        bC[gid] ~ normal(0, 1),
        c(bAge, bE) ~ normal(0, 1),
        vector[8]: delta_j <<- append_row(0, delta),
        simplex[7]: delta ~ dirichlet(alpha)
    ),
    data=list(E=d$E, male=d$male, alpha=rep(2, 7),
              age=d$age / max(d$age),
              gid=d$male + 1,
              R=d$response,
              C=d$contact), chains=4, cores=4
)

# The question is, is there a difference between the Contact coefficient between
# Men and Women?
# Doesn't make a huge difference when adjusting for Intent and Action
precis(m12h5, pars="bC", depth=2)
precis(m12h5_b, pars="bC", depth=2)

# Calculate the contrasts
calc_contrasts <- function(mod) {
    post <- extract.samples(mod)
    tibble(women_less_men=post$bC[, 1] - post$bC[, 2]) |> precis()
}

# So women are associated with a -0.36 decrease in cum-log-probs in moral acceptability
# than men
# AND it doesn't matter whether we account for Intent or Action (which suggests they
# are independent as expectd)
# Can we convert this back to probability?
# The pipeline from Prob to cum-log-odds is:
#   - prob -> cumsum -> logit
# We can inverse logit, but cumsum isn't invertible
calc_contrasts(m12h5)
calc_contrasts(m12h5_b)

# Out of interest, what would we conclude if we'd just done a bivariate analysis?
m12h5_c <- ulam(
    alist(
        R ~ ordered_logistic(phi, kappa),
        phi <- a[gid] + bC[gid]*C,
        kappa ~ normal(0, 1.5),
        a[gid] ~ normal(0, 1),
        bC[gid] ~ normal(0, 1)
    ),
    data=list(gid=d$male + 1,
              R=d$response,
              C=d$contact), chains=4, cores=4
)
# We get a positive score, indicating that Women find it _more_ morally permissible
# when there is contact, contradictory to the hypothesis
calc_contrasts(m12h5_c)

# 12H6
data(Fish)
Fish <- as_tibble(Fish)
# Each row is a visit to a NP
# We want to know how many fish a visitor makes per hour
# But not everyone tried to fish! So have Zero-inflated data
# Also, need to adjust for the period being different (hours)
Fish
# Have livebait, camper, persons, child
# I think having livebait is a predictor of whether someone is fishing
# But also a predictor of HOW MANY fish were caught
# Yeah so 70% of livebait = 0 also had 0 fish caught
# Vs 55% of livebait = 1, so that seems a reasonable predictor
Fish |>
    count(livebait, fish_caught>0) |>
    group_by(livebait) |>
    mutate(prop = n / sum(n))

# I thought having a child might be an indicator that you're not there to fish,
# based on the average age of fishers on the river near me
# 73% of people with a child caught no fish, vs 42% of people without a child
# so that seems easonable too
Fish |>
    count(has_child=child>0, fish_caught>0) |>
    group_by(has_child) |>
    mutate(prop = n / sum(n)) |>
    arrange(desc(has_child))

# I might expect a camper is more likely to be a fisher given the early hours that 
# fishers often get up
# 67% of non-campers caught zero fish, vs 50% of campers
# So yes that seems reasonable
Fish |>
    count(camper, fish_caught>0) |>
    group_by(camper) |>
    mutate(prop = n / sum(n))

# So what about the Poisson model for number of fish caught?
# I think only bait seems reasonable here, unless maybe if there were
# more people it counts each person's fish?
# Ah yes this looks likely, more people tends to equal more fish
Fish |>
    group_by(persons) |>
    summarise(mean(fish_caught))

# Fish caught is very skewed!
# Mean is 3 but can go right up to 149
summary(Fish$fish_caught)

# So my model is:
# WasFishing ~ livebait + camper + had_child
# FishCaught ~ livebait + persons
# While also adding the offset for the number of hours
# The only transform I'll do is z-transforming Persons
m12h6 <- ulam(
    alist(
        fish_caught ~ dzipois(p, lambda),
        logit(p) <- az + bzLivebait * livebait + bzCamper * camper + bzChild * has_child,
        log(lambda) <- log_hours + af + bfLivebait * livebait + bfPersons * persons_z,
        az ~ dnorm(0, 1),  # Prior centered around 0.5 chance of fishing
        c(bzLivebait, bzCamper, bzChild) ~ dnorm(0, 1),
        af ~ dnorm(1, 15),
        c(bfLivebait, bfPersons) ~ dnorm(0, 1)
    ),
    data=Fish |> 
        mutate(
            has_child = as.integer(child > 0),
            persons_z = as.numeric(scale(persons)),
            log_hours = log(hours)
        ) |>
        select(fish_caught, livebait, camper, has_child, persons_z, log_hours),
    chains=4, cores=4
)
# MCMC went surprisingly well! Decent ess and rhat
precis(m12h6)
# Let's look at these coefficients:
# Predicting whether was fisher:
#  - Having a child has a slightly negative effect on whether someone was a fisher
#  - Being a camper only has a very small association
#  - Likewise if someone used livebait! (why would you use livebait otherwise?!)
# Predicting fish caught:
#  - Each person adds an additional 2.2x multiplier to the number of fish caught
#    which starts at 0.13
#  - Having livebait increases rate four-fold 
plot(precis(m12h6))

# Posterior predictions
mu <- link(m12h6)
mu_mean <- colMeans(mu$lambda)    # Don't want to see predictions on whether fisher 
mu_PI <- apply(mu$lambda, 2, PI)  # have ground truth to compare to
preds <- sim(m12h6)
preds_PI <- apply(preds, 2, PI)

# This model hasn't done great
# There are a bunch of groups that caught ~10 fish that the model all predicted
# > 50
# And the most succesful group (~150) were only predicted 70
# However, despite that, the model conveys the huge uncertainties in these predictions
# Although there is a blind-spot between ~20 and 50 fish where the model is
# certain these groups weren't fishing
fish_plt <- Fish |>
    mutate(
        mu_mean = mu_mean,
        mu_lower=mu_PI[1, ],
        mu_upper=mu_PI[2, ],
        preds_lower=preds_PI[1, ],
        preds_upper=preds_PI[2, ],
        group_num = row_number()
    ) |>
    arrange(fish_caught) |>
    mutate(group_num = factor(group_num, levels=group_num)) 

fish_plt |>
    ggplot(aes(x=group_num)) +
        geom_point(aes(y=fish_caught), colour="red") +
        geom_point(aes(y=mu_mean), colour="steelblue") +
        geom_errorbar(aes(ymin=preds_lower, ymax=preds_upper), colour="orange") +
        geom_errorbar(aes(ymin=mu_lower, ymax=mu_upper), colour="steelblue") +
        theme_minimal()

# Let's look at this group in more detail
# As well as the highest fishing group!
# So the group with the highest amount of fish caught were camping, were using livebait
# didn't have children, so the model should have had enough information to know
# they were likely actually fishing, but also they were visiting for 36 hours
fish_plt |>
    arrange(desc(fish_caught))

# Which wasn't the longest stay, although quite up there
# This also shows some of the comical overpredictions (these are to the left
# of the observations mentioned above where the model severely _UNDER_ predicts),
# where the model predicted a large number of fish caught, say for the 71 hours 
# event, predicting 113-129, but only 8 were caught in actuality, despite
# having the exact same predictors as the group that caught 149, AND
# having double the time, so there's not that much could be done there
fish_plt |>
    arrange(desc(hours))

# So what happened with the comical under-predictions?
# I imagine they were all not using livebait, weren't camping, and had children,
# leading to the model thinking they weren't actually fishing
# Huh not necessarily, there's a lot of instances of ~30 fish being caught
# with livebait, camping, 4 people, 0 children. But some of these were only fishing
# for a very short time, so maybe they were just extremely efficient
# I could see this making sense - people turning up for just 3 hours are going to
# be maximising their fishing time, while people fishing for 36 hours will be
# taking it easy and also needing sleep!
fish_plt |>
    filter(fish_caught > mu_upper) |>
    arrange(desc(fish_caught))
    
# 12H7
dag <- dagitty("dag { 
A -> E -> R;
A -> R;
}")
coordinates(dag) <- list(x=c(A=0, E=1, R=1), y=c(A=0, E=0, R=1))
drawdag(dag)

# Need to fit E ~ A
# And R ~ E + A
# Taking too long to fit as single model!
# Try as 2 separate
m12h7_edu <- ulam(
    alist(
        E ~ ordered_logistic(phiE, kappaE),
        phiE <- aE + bAgeE * age,
        aE ~ normal(0, 1),
        bAgeE ~ normal(0, 1),
        kappaE ~ normal(0, 1.5)
    ),
    data=list(E=d$E,
              age=d$age / max(d$age)), chains=4, cores=4
)

m12h7_response <- ulam(
    alist(
        R ~ ordered_logistic(phiR, kappaR),
        phiR <- bER*sum(delta_j[1:E]) + bAgeR * age,
        kappaR ~ normal(0, 1.5),
        c(bAgeR, bER) ~ normal(0, 1),
        vector[8]: delta_j <<- append_row(0, delta),
        simplex[7]: delta ~ dirichlet(alpha)
    ),
    data=list(E=d$E, alpha=rep(2, 7),
              age=d$age / max(d$age),
              R=d$response), chains=4, cores=4
)

# So firstly does Age have an association with Education?
# Should hope so...
# Yes, a massive effect size
precis(m12h7_edu)

# And what about the direct effects of both Age and Education on Response?
# Education's effect is far smaller now (was previously -0.31 [-0.59 - -0.07])
# But is now 0.16 with -0.2 - 0.34
# Meanwhile, age has a strong effect that is consistently negative, so it seems
# that the apparent relationship between education and moral rating is nearly all
# driven by age
precis(m12h7_response)
precis(m12.6)

# 12H8
# Is it possible that gender is a confounder too?
dag2 <- dagitty("dag { 
E <- G -> R;
A -> E -> R;
A -> R;
}")
coordinates(dag2) <- list(x=c(A=0, E=1, R=1, G=2), y=c(A=0, E=0, R=1, G=0))
drawdag(dag2)

# What would need to test for to determine the causal influence of education
# on response?
# Yes, would need to adjust for Gender too!
adjustmentSets(dag2, exposure="E", outcome="R", effect="direct")

# NB: this is very similar to a model I've already fitted in 12H5
# Where I had response ~ age + education + Contact, where contact
# had a varying slope by gender
# But I haven't actually fitted a model of just age, edu, gender before
# Will use index variable
precis(m12h5, depth=2)

m12h8 <- ulam(
    alist(
        R ~ ordered_logistic(phiR, kappaR),
        phiR <- bER*sum(delta_j[1:E]) + bAgeR * age + a[gid],
        kappaR ~ normal(0, 1.5),
        a[gid] ~ normal(0, 1),
        c(bAgeR, bER) ~ normal(0, 1),
        vector[8]: delta_j <<- append_row(0, delta),
        simplex[7]: delta ~ dirichlet(alpha)
    ),
    data=list(E=d$E, alpha=rep(2, 7),
              gid = d$male + 1,
              age=d$age / max(d$age),
              R=d$response), chains=4, cores=4
)
# This doesn't effect the findings: education still doesn't have a strong effect
precis(m12h8, depth=2)

# And if we want to see the effect of age we can use the same model
adjustmentSets(dag2, exposure="A", outcome="R", effect="direct")
# And Age still has a moderately strong association with response
