library(rethinking)
library(tidyverse)

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




# Things to ask Richard: --------------------------------------------------
# Why is the simplex needed in m12.6, if dirichlet always returns values that sum to 1?
