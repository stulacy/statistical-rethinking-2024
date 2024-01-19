
# 5.1 Spurious association ------------------------------------------------

library(rethinking)
data("WaffleDivorce")
d <- WaffleDivorce

summary(d)

# standardize
d$A <- scale(d$MedianAgeMarriage)
d$D <- scale(d$Divorce)

# Regress divorce rate on age
m5.1 <- quap(
    alist(
        D ~ dnorm(mu, sigma),
        mu <- a + bA * A,
        a ~ dnorm(0, 0.2),
        bA ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), 
    data=d
)
# Negative slope, i.e. divorce rate decreases with age at marriage
precis(m5.1)

# Prior predictive check
set.seed(10)
prior <- extract.prior(m5.1)
# extract.prior gives N (default 1e4) samples from each prior
names(prior)
length(prior$a)

# Get mean value using these values
mu <- link(m5.1, post=prior, data=list(A=c(-2, 2)))
# 2 columns as 2 data points
# I guess it only needs 2 points to form the line since the link function is deterministic
# Interestingly, the default N for link is 1e3, so it's not using all of the 1e4 prior values
dim(mu)
plot(NULL, xlim=c(-2, 2), ylim=c(-2, 2))
# Plot the first 50 lines
for (i in 1:50) lines(c(-2, 2), mu[i, ], col=col.alpha("black", 0.4))

# posterior predictions over specified range
A_seq <- seq(from=-3, to=3.2, length.out=30)
mu <- link(m5.1, data=list(A=A_seq))
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

plot(D ~ A, data=d, col=rangi2)
lines(A_seq, mu.mean, lwd=2)
shade(mu.PI, A_seq)

# Regressing marrage rate on divorce rate
d$M <- scale(d$Marriage)
m5.2 <- quap(
    alist(
        D ~ dnorm(mu, sigma),
        mu <- a + bM * M,
        a ~ dnorm(0, 0.2),
        bM ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ),
    data=d
)

# Prior predictive checks
# Look sensibleish, allows for both positive & negative relationships
set.seed(10)
prior <- extract.prior(m5.2)
mu <- link(m5.2, post=prior, data=list(M=c(-2, 2)))
plot(NULL, xlim=c(-2, 2), ylim=c(-2, 2))
# Plot the first 50 lines
for (i in 1:50) lines(c(-2, 2), mu[i, ], col=col.alpha("black", 0.4))

# Posterior plot
M_seq <- seq(from=-3, to=3.2, length.out=30)
mu <- link(m5.2, data=list(M=M_seq))
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)
plot(D ~ M, data=d, col=rangi2)
lines(M_seq, mu.mean, lwd=2)
shade(mu.PI, A_seq)

# Weaker relationship?
# 0.35 vs -0.57
precis(m5.1)
precis(m5.2)


# 5.1.2 Testable implications ---------------------------------------------

# Drawing DAGs
library(dagitty)
dag5.1 <- dagitty(" dag {
    A -> D
    A -> M
    M -> D
}")
# Don't need to provide coordinates, only for controlling plot layout
coordinates(dag5.1) <- list(x=c(A=0, D=1, M=2), y=c(A=0, D=1, M=0))
drawdag(dag5.1)

DMA_dag2 <- dagitty("dag { D <- A -> M }")
drawdag(DMA_dag2)

# Can infer implied conditional independencies
# D is independent of M when controlling for A
impliedConditionalIndependencies(DMA_dag2)
# No conditional independencies in first DAG
impliedConditionalIndependencies(dag5.1)


# 5.1.4 Approximating the posterior ---------------------------------------

# Fitting multiple regression model 5.10
m5.3 <- quap(
    alist(
        D ~ dnorm(mu, sigma),
        mu <- a + bM*M + bA*A,
        a ~ dnorm(0, 0.2),
        bM ~ dnorm(0, 0.5),
        bA ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ),
    data=d
)
# Now no effect of bM!
precis(m5.3)

# Really nice function from rethinking package to see how coefs change over models
# So marriage rate was correlated with divorce rate in univariate model, but once control for age
# this effect disappears. So the first dag looks wrong as it had direct M-D connection
# And we've tested the implied conditional independency from dag2 and the results support this
plot(coeftab(m5.1, m5.2, m5.3), par=c("bA", "bM"))

# Coeftab on its own just displays a summary of estimates
coeftab(m5.1, m5.2, m5.3)

# 5.1.5.1 Predictor residual plots ----------------------------------------
# Have D ~ A + M
# Now want to see individual effect of a predictor holding
# the others constant
# So firstly model M ~ A
m5.4 <- quap(
    alist(
        M ~ dnorm(mu, sigma),
        mu <- a + bAM * A,
        a ~ dnorm(0, 0.2),
        bAM ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ),
    data = d
)
# Now compute RESIDUALS from MAP/mean mu
mu <- link(m5.4)
mu_mean <- apply(mu, 2, mean)
mu_resid <- d$M - mu_mean

# To use them, plot against D
# "... displaying the linear relationship between divorce and marrage rates,
# having statistically controlled for median age of marriage"
# "...brings home the message that regression models measure the
# remaining association of each predictor with the outcome,
# after already knowing the other predictors."
# Haven't added code to clean up this plot to replicate
# Fig 5.4. Looks like then fits another regression
# on D ~ resids to show that M has little effect
# on D after taking A into account. Then repeats for
# A ~ M (right hand panel Fig 5.4)
plot(mu_resid, d$D)

# 5.1.5.2 Posterior prediction plots --------------------------------------
# Model of D ~ M + A
m5.3

# Firstly get mean predictions **for the original data**
mu <- link(m5.3)
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI)

# Now get posterior predictive **for the original data**
D_sim <- sim(m5.3)
D_PI <- apply(D_sim, 2, PI)

# Plot! This is the type of diagnostic plot I frequently do
# Although NB: This plot **doesn't show the posterior predictive**
# Just the posterior mean divorce rate
# Rather than running a sweep of x-values and plotting that line,
# Instead it's just plotting the predicted points along with their 89% PI
# Is this more useful or less useful than having the regression line
# with it's PI shaded region?
# Why not both?! Could easily plot that underneath this
plot(mu_mean ~ d$D, col=rangi2, ylim=range(mu_PI),
     xlab="Observed divorce", ylab="Predicted divorce")
abline(a=0, b=1, lty=2)
for (i in 1:nrow(d)) {
    lines(rep(d$D[i], 2), mu_PI[, i], col=rangi2)
}

# Genuinely had never head of the 'identify' function before!
# I'd rather use plotly for this but still good to know for base plots
identify(x=d$D, y=mu_mean, labels=d$Loc)

# 5.18 Overthinking -------------------------------------------------------

N <- 100
x_real <- rnorm(N)
x_spur <- rnorm(N, x_real)
y <- rnorm(N, x_real)
d <- data.frame(y, x_real, x_spur)
# Can see correlation between x_spur and y, even though it isn't 
# explicitly correlated
pairs(d)
# But when model it the model picks out
# that x_real is the real driver
summary(lm(y ~ x_real + x_spur))
d_5.18 <- dagitty("DAG { y <- x_real -> x_spur }")
# This is the DAG for this dataset
drawdag(d_5.18)

# 5.1.5.3 Counterfactual plots --------------------------------------------
data("WaffleDivorce")
d <- list()
# Is this basically a wrapper for scale?
# But it returns a vector rather than a column vector?
d$A <- standardize(WaffleDivorce$MedianAgeMarriage)
d$D <- standardize(WaffleDivorce$Divorce)
d$M <- standardize(WaffleDivorce$Marriage)

# Using this original dag
drawdag(dag5.1)
# To model this will use 2 models:
# One for A and M on D, and one for A on M
# Or could think of this as one single joint model

m5.3_A <- quap(
    alist(
        ## A -> D <- M
        D ~ dnorm(mu, sigma),
        mu <- a + bM*M + bA*A,
        a ~ dnorm(0, 0.2),
        bM ~ dnorm(0, 0.5),
        bA ~ dnorm(0, 0.5),
        sigma ~ dexp(1),
        ## A -> M
        M ~ dnorm(mu_M, sigma_M),
        mu_M <- aM + bAM*A,
        aM ~ dnorm(0, 0.2),
        bAM ~ dnorm(0, 0.5),
        sigma_M ~ dexp(1)
    ),
    data=d
)
# Get the strong bAM effect, which suggests
# A reduces M.
# Also get the low M-D effect since adjusting for A
precis(m5.3_A)

# So want to simulate what happens if A changes
# So need a sweep over A, which we'll use to then generate M,
# and then D
A_seq <- seq(from=-2, to=2, length.out=30)

# Use 'vars' argument to tell it to A on M first,
# then the D model
# Not hard to implement ourselves but nice that can do it easily with
# this function. I wish in Stan it were so easy to generate predictions - 
# Either have to do it yourself with the posterior, or be aware of all
# predictions will want to run when writing the model
sim_dat <- data.frame(A = A_seq)
s <- sim(m5.3_A, data=sim_dat, vars=c("M", "D"))
# Get 1 entry for each likelihood
names(s)
# With the usual n_sims x n_data_points matrix
dim(s$M)

# Plot counter-factual
# D
plot(sim_dat$A, colMeans(s$D), ylim=c(-2, 2), type='l',
     xlab="manipulated A", ylab="counterfactual D")
shade(apply(s$D, 2, PI), sim_dat$A)
mtext("Total counterfactual effect of A on D")

# M
plot(sim_dat$A, colMeans(s$M), ylim=c(-2, 2), type='l',
     xlab="manipulated A", ylab="counterfactual M")
shade(apply(s$M, 2, PI), sim_dat$A)
mtext("Total counterfactual effect of A on M")

# Now want to see effect of M on D, so hold A constant here
# Which breaks the A->M edge, then sim
sim_dat <- data.frame(M=seq(-2, to=2, length.out=30), A=0)
s <- sim(m5.3_A, data=sim_dat, vars="D")

# No strong effect!
plot(sim_dat$M, colMeans(s), ylim=c(-2, 2), type='l',
     xlab="manipulated M", ylab="counterfactual D")
shade(apply(s, 2, PI), sim_dat$M)
mtext("Total counterfactual effect of M on D")

# However, couldn't we just have got this effect from the coefficient of M
# from the multivariate model with A as a predictor (same as holding age constant)?
# So the slope from the counter factual plot is -0.067571
lm(colMeans(s) ~ sim_dat$M)
# Which seems to be the same as the -0.07 value from the model, so yes there's no need to do this
# counterfactual business for this effect!
precis(m5.3_A)

# And 
quap(
    alist(
        ## A -> D <- M
        D ~ dnorm(mu, sigma),
        mu <- a + bM*M + bA*A,
        a ~ dnorm(0, 0.2),
        bM ~ dnorm(0, 0.5),
        bA ~ dnorm(0, 0.5),
        sigma ~ dexp(1),
    ),
    data=d
)

# NB: often not possible to produce un-confounded causal effect
# But this technique is still useful

# 5.2 Masked relationship -------------------------------------------------
data(milk)
d <- milk
str(d)

d$K <- scale(d$kcal.per.g)
d$N <- scale(d$neocortex.perc)
d$M <- scale(log(d$mass))

# Have not seen a z-transform of a logged var before!
# Raw var has long right tail
dens(d$mass)
# Log transform makes it more normal
dens(log(d$mass))
# z-scale just places density between -2, 2 but keeps shape the same
dens(scale(log(d$mass)))

# Univariate model of calories ~ neocortex
# Doesn't work! Have missing values
#m5.5_draft <- quap(
#    alist(
#        K ~ dnorm(mu, sigma),
#        mu <- a + bN*N,
#        a ~ dnorm(0, 1),
#        bN ~ dnorm(0, 1),
#        sigma ~ dexp(1)
#    ),
#    data=d
#)
summary(d)

dcc <- d[complete.cases(d$K, d$N, d$M), ]
summary(dcc)

# Fit model using complete data set
m5.5_draft <- quap(
    alist(
        K ~ dnorm(mu, sigma),
        mu <- a + bN*N,
        a ~ dnorm(0, 1),
        bN ~ dnorm(0, 1),
        sigma ~ dexp(1)
    ),
    data=dcc
)

# Prior predictive
# Not great priors, would expect
# them to pass through 0,0
prior <- extract.prior(m5.5_draft)
xseq <- c(-2, 2)
mu <- link(m5.5_draft, post=prior, data=list(N=xseq))
plot(NULL, xlim=xseq, ylim=xseq)
for (i in 1:50) {
    lines(xseq, mu[i,], col=col.alpha("black", 0.3))
}

# Richard suggests using a tighter prior on alpha
# to make it closer to zero, since with 2 standardized
# variables, the expected value of the outcome should be zero
# And makes the slope a bit tighter as well
m5.5 <- quap(
    alist(
        K ~ dnorm(mu, sigma),
        mu <- a + bN*N,
        a ~ dnorm(0, 0.2),   # Was 1 before
        bN ~ dnorm(0, 0.5),  # Was 1 before
        sigma ~ dexp(1)
    ),
    data=dcc
)
precis(m5.5)

# Much better!
prior <- extract.prior(m5.5)
xseq <- c(-2, 2)
mu <- link(m5.5, post=prior, data=list(N=xseq))
plot(NULL, xlim=xseq, ylim=xseq)
for (i in 1:50) {
    lines(xseq, mu[i,], col=col.alpha("black", 0.3))
}

# Still not significant slope
precis(m5.5)

# Plot posterior mean (difference from prior predictive is that
# the parameter estimates are from the posterior not the prior)
# "Weakly positive" trend
xseq <- seq(from=min(dcc$N)-0.15, to=max(dcc$N)+0.15, length.out=13)
mu <- link(m5.5, data=list(N=xseq))
mu_mean <- colMeans(mu)
mu_PI <- apply(mu, 2, PI)
plot(K~N, data=dcc)
lines(xseq, mu_mean, lwd=2)
shade(mu_PI, xseq)

# Now model calories (K) ~ mass (M)
# NB: remember we log transformed mass
# Looks weakly negative trend here!
m5.6 <- quap(
    alist(
        K ~ dnorm(mu, sigma),
        mu <- a + bM*M,
        a ~ dnorm(0, 0.2),   # Was 1 before
        bM ~ dnorm(0, 0.5),  # Was 1 before
        sigma ~ dexp(1)
    ),
    data=dcc
)
precis(m5.6)

# Looks better but still a lot of uncertainty
xseq <- seq(from=min(dcc$M)-0.15, to=max(dcc$M)+0.15, length.out=13)
mu <- link(m5.6, data=list(M=xseq))
mu_mean <- colMeans(mu)
mu_PI <- apply(mu, 2, PI)
plot(K~M, data=dcc)
lines(xseq, mu_mean, lwd=2)
shade(mu_PI, xseq)

# Now to include both variables in a model
m5.7 <- quap(
    alist(
        K ~ dnorm(mu, sigma),
        mu <- a + bN*N + bM*M,
        a ~ dnorm(0, 0.2),   # Was 1 before
        bN ~ dnorm(0, 0.5),  # Was 1 before
        bM ~ dnorm(0, 0.5),  # Was 1 before
        sigma ~ dexp(1)
    ),
    data=dcc
)
# Now both coefficients seem to be both stronger
# and less variance
precis(m5.7)
plot(coeftab(m5.5, m5.6, m5.7))

# How did this happen?
# Here we know we have 2 variables correlated with outcome (one +ve, other -ve),
# but both predictors are correlated with each other
pairs(~K+M+N, dcc)

# 3 possible dags
# M causes K & N, N also has direct cause with K
# This direct cause doesn't appear in the univariate model
# why? Because of confounding from M?
dag_1 <- dagitty("dag { K <- M -> N -> K }")
drawdag(dag_1)

# Same as above but N causes M this time
dag_2 <- dagitty("dag { K <- N -> M -> K }")
drawdag(dag_2)

# Here U is a common factor that causes both M and N
dag_3 <- dagitty("dag {K[o]M[e]U[u]N[e] K <- M <- U -> N -> K }")
drawdag(dag_3)

# NB: CANNOT DIFFERENTIATE BETWEEN THESE CANDIDATE DAGS
# FROM THE DATA ALONE! 
# This is because all have the same conditional independencies
impliedConditionalIndependencies(dag_1)
impliedConditionalIndependencies(dag_2)
impliedConditionalIndependencies(dag_3)

# Counter factual plots
# First vary M, holding N constant, i.e. breaking the N->M edge
xseq <- seq(from=min(dcc$M)-0.15, to=max(dcc$M)+0.15, length.out=30)
mu <- link(m5.7, data=data.frame(M=xseq, N=0))
mu_mean <- colMeans(mu)
mu_PI <- apply(mu, 2, PI)
plot(NULL, xlim=range(dcc$M), ylim=range(dcc$K))
lines(xseq, mu_mean, lwd=2)
shade(mu_PI, xseq)

# And likewise for N holding M constant
xseq <- seq(from=min(dcc$N)-0.15, to=max(dcc$N)+0.15, length.out=30)
mu <- link(m5.7, data=data.frame(N=xseq, M=0))
mu_mean <- colMeans(mu)
mu_PI <- apply(mu, 2, PI)
plot(NULL, xlim=range(dcc$N), ylim=range(dcc$K))
lines(xseq, mu_mean, lwd=2)
shade(mu_PI, xseq)

# Both of these look significant, so the N and M direct effects on K seem realistic
# (But could tell this already from the multiple regression coefs)

# Simulating these relationships
# Dag 1: M -> K <- N
# M -> N
n <- 100
M <- rnorm(n)
N <- rnorm(n, M)
K <- rnorm(n, N-M)
d_sim <- data.frame(K=K, N=N, M=M)
# Then if use this as data frame in models 5.5, 5.6, 5.7 we see the same pattern

m5.5_sim <- quap(
    alist(
        K ~ dnorm(mu, sigma),
        mu <- a + bN*N,
        a ~ dnorm(0, 0.2),   # Was 1 before
        bN ~ dnorm(0, 0.5),  # Was 1 before
        sigma ~ dexp(1)
    ),
    data=d_sim
)
m5.6_sim <- quap(
    alist(
        K ~ dnorm(mu, sigma),
        mu <- a + bM*M,
        a ~ dnorm(0, 0.2),   # Was 1 before
        bM ~ dnorm(0, 0.5),  # Was 1 before
        sigma ~ dexp(1)
    ),
    data=d_sim
)
m5.7_sim <- quap(
    alist(
        K ~ dnorm(mu, sigma),
        mu <- a + bN*N + bM*M,
        a ~ dnorm(0, 0.2),   # Was 1 before
        bN ~ dnorm(0, 0.5),  # Was 1 before
        bM ~ dnorm(0, 0.5),  # Was 1 before
        sigma ~ dexp(1)
    ),
    data=d_sim
)
# Can see the effects getting stronger in the multiple
# regression model
plot(coeftab(m5.5_sim, m5.6_sim, m5.7_sim))

# Could also do this for the other 2 dags
# Dag 2 is the same as Dag 1 but swapping M & N
# Dag 3 is more interesting, would look like this
# M -> K <- N
# M <- U -> N
# This would also produce results consistent with the actual results
n <- 100
U <- rnorm(n)
M <- rnorm(n, U)
N <- rnorm(n, U)
K <- rnorm(n, N-M)

# Categorical variables ---------------------------------------------------
data("Howell1")
d <- Howell1
# Going to be looking at sex, an indicator variable [0, 1]
str(d)

# One motivation for using index variable rather than indicator
# is that the intercept is just the average height
# for the group with index 0, rather than average overall
# It also means there is more uncertainty about the '1' case,
# because it has uncertanty in both intercept and slope, where the '0'
# case just has uncertainty in the intercept
# Instead use index variable
d$sex <- ifelse(d$male == 1, 2, 1)

# Fit model of height ~ sex using indicator
# Here assigning same prior to both sexes but could use different
m5.8 <- quap(
    alist(
        height ~ dnorm(mu, sigma),
        mu <- a[sex],
        a[sex] ~ dnorm(178, 20),
        sigma ~ dunif(0, 50)
    ),
    data=d
)
# Can see women are shorter on average
precis(m5.8, depth=2)

# To get the expected difference between males and females
# we must find the expected posterior difference, NOT the difference
# of the 2 expected posteriors!

# GOOD: Contrasts method
# here the mean difference is -7.65, i.e. women are 7.65 shorter on average
post <- extract.samples(m5.8)
post$diff_fm <- post$a[, 1] - post$a[, 2]
precis(post, depth=2)

# BAD: Difference of expected posteriors
# This does give the same answer but just as a point estimate! Should always work in full posterior
mus <- colMeans(post$a)
mus[2] - mus[1]

# 5.3.2 Many categories ---------------------------------------------------
data(milk)
d <- milk
# 4 categories!
unique(d$clade)

# Easiest way to get index variables if the vector is already factor-encoded
d$clade_id <- as.integer(d$clade)

# Model average milk energy per group
d$K <- scale(d$kcal.per.g)

m5.9 <- quap(
    alist(
        K ~ dnorm(mu, sigma),
        mu <- a[clade_id],
        a[clade_id] ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ),
    data=d
)
labels <- paste("a[", 1:4, "]:", levels(d$clade), sep="")
# Hadn't realised precis had a plot method!
plot(precis(m5.9, depth=2, pars="a"), labels=labels, 
      xlab="expected kcal (std)")

# Problems ----------------------------------------------------------------

# 5E1
# 2 and 4 both have multiple predictors

# 5E2
# latitude ~ alpha + b_1 * animal + b_2 * plant diversity

# 5E3
# time ~ alpha + b_1 * funding + b_2 * size
# b_1, b_2 > 0

# 5E4
# 1, 3, 4, 5 are all the same, the intercept will take on A's value 
# in 1 & 3, and A is explicit in 4 and 5.
# 2 is different as it has a separate intercept

# 5M1
# The example from the chapter was on marriage rates
# So in the example in the book, both Age and Marriage rate
# were associated with Divorce rates in their separate regressions
# (m5.1 and m5.2 respectively), but when included in the same model
# the marriage effect shrinks (m5.3)
# This implies the DAG D <- A -> M is plausible, as we tested the
# implied conditional independencies of D _||_ M | A, which is what
# m5.3 was doing, and the low effect of M -> D implies this DAG
# is feasible.
# The first DAG mplies that all pairs of variables should be associatd 
# regardless of what we condition on, but the test rules this out
plot(coeftab(m5.1, m5.2, m5.3), pars=c("bM", "bA"))

# So then we want a model like this:
dag <- dagitty("dag { B <- A -> Y }")
drawdag(dag)

N <- 100
A <- rnorm(N)
B <- rnorm(N, A)
Y <- rnorm(N, A)

# Just using lm to save time... Would expect same results
# from Bayesian
# A is associated with Y
summary(lm(Y ~ A))
# B is associated with Y
summary(lm(Y ~ B))
# But when we control for A, B is no longer significant
summary(lm(Y ~ A + B))

# 5M2
# Masked relationship:
# This is the milk - neocortex - mass example from 5.2
# Have the 3 possible dags
# But can't resolve between them with the data we have
drawdag(dag_1)
drawdag(dag_2)
drawdag(dag_3)
# The data showed that in the Univariate models there was little to noeffect
# but after controlling it was significant!
# This is in constrast to the confounder above, where there was an effect
# in the Univariate, that disappeared in the multivariate
plot(coeftab(m5.5, m5.6, m5.7), pars=c("bN", "bM"))

# Can simulate, say using the dag_3 setup with an unobserved confounder
N <- 1e2
U <- rnorm(N)
A <- rnorm(N, U)
B <- rnorm(N, U)
Y <- rnorm(N, A-B)

summary(lm(Y ~ A)) # 0.4 +/- 0.1
summary(lm(Y ~ B)) #  -0.45 +/- 0.1
summary(lm(Y ~ A + B))  # 1 +/- 0.1 and -1 +/- 0.1

# What happens if both have a positive association with the outcome?
N <- 1e2
U <- rnorm(N)
A <- rnorm(N, U)
B <- rnorm(N, U)
Y <- rnorm(N, A+B)

# So here we have the opposite, in that we have a strong effect
# in the Univariate that gets slightly dampened in the multivariate
# So having a positive & negative effect are crucial to masking
summary(lm(Y ~ A)) # 1.3 +/- 0.1
summary(lm(Y ~ B)) #  1.5 +/- 0.1
summary(lm(Y ~ A + B))  # 0.9 +/- 0.08 and 1 +/- 0.08

# 5M3
# More divorces meaning more dating once they become single = more marriages!
# I can't think of a confounding factor or a causal diagram to draw and test
# Unless we just use the age example from the chapter

# 5M4
# Add Mormom state data to divorce model to see if reduces Utah outlier
# Mormon data obtained from: https://worldpopulationreview.com/state-rankings/mormon-population-by-state
# State data obtained from: https://worldpopulationreview.com/state-rankings/mormon-population-by-state
library(tidyverse)
d <- WaffleDivorce
d$A <- scale(d$MedianAgeMarriage)
d$D <- scale(d$Divorce)
d$M <- scale(d$Marriage)

# Mormon data is from 2022 despite the filename!
df_mormon <- read_csv("data/mormon-population-by-state-2024.csv")
# Will use 2023 population data as is closest to 2022
df_pop <- read_csv("data/state-population-table.csv")
# Calculate % mormon and join with Waffle data
d2 <- df_mormon |>
    inner_join(df_pop, by="state") |>
    mutate(mormon_pct = MormonPopulation2022 / pop2023) |>
    select(state, mormon_pct) |>
    left_join(d, by=c("state"="Location")) |>
    select(state, mormon_pct, A, D, M)  
# Nevada is missing divorce data!
d2[!complete.cases(d2), ]
# And the Waffle data has District of Columbia which is missing from the mormon data
d |> anti_join(d2, by=c("Location"="state"))

# Move forward with the 49 states do have!
d2 <- d2[complete.cases(d2), ]

# Question hints we might want to standardize this
# It is very skewed...
hist(d2$mormon_pct)
# 64% of Utah is mormon, that's ridiculous!
d2 |> arrange(desc(mormon_pct))
# A log then z-score should be enough
d2$LDS <- scale(log(d2$mormon_pct))
# Yep that looks better
hist(d2$LDS)

# Now for model
summary(lm(D ~ A + M + LDS, data=d2))
m5m4 <- quap(
    alist(
        D ~ dnorm(mu, sigma),
        mu <- a + bA * A + bM * M + bL * LDS,
        a ~ dnorm(0, 0.2),
        bA ~ dnorm(0, 0.5),
        bM ~ dnorm(0, 0.5),
        bL ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ),
    data = d2
)
precis(m5m4)
plot(coeftab(m5m4))

# Does this give better predictions?
# It has fewer outliers that's for sure!
plot_post_waffle <- function(mod, data) {
    mu <- link(mod)
    mu_mean <- apply(mu, 2, mean)
    mu_PI <- apply(mu, 2, PI)
    
    # Now get posterior predictive **for the original data**
    D_sim <- sim(mod)
    D_PI <- apply(D_sim, 2, PI)
    
    plot(mu_mean ~ data$D, col=rangi2, ylim=range(mu_PI),
         xlab="Observed divorce", ylab="Predicted divorce")
    abline(a=0, b=1, lty=2)
    for (i in 1:nrow(data)) {
        lines(rep(data$D[i], 2), mu_PI[, i], col=rangi2)
    }
}

plot_post_waffle(m5.3, WaffleDivorce |> mutate(D = scale(Divorce)))
plot_post_waffle(m5m4, d2)

# 5M5
# The question implies the following dag
dag_5m5 <- dagitty("dag {
  G[e]
  D[e]
  R[e]
  O[o]
  G -> O
  G -> D
  D -> O
  D -> R
  R -> O
}")
drawdag(dag_5m5)
# This DAG implies G _||_ R | D
# So would model G ~ R + D and see if R slope is far from zero
# If it isn't, then this DAG is plausible. 
# If it is far from zero, then this DAG isn't plausible
impliedConditionalIndependencies(dag_5m5)

# We want to know if G -> O is a 'real' association,
# what needs to be controlled for?
# Need to control for D otherwise 

# 5H1
data(foxes)
# Hasn't been scaled!
summary(foxes)
# Do that now for numerical data
foxes_std <- foxes |>
            mutate(
                avgfood = scale(avgfood),
                groupsize = scale(groupsize),
                area = scale(area),
                weight = scale(weight)
            )
m5h1_a <- quap(
    alist(
        weight ~ dnorm(mu, sigma),
        mu <- a + bA * area,
        a ~ dnorm(0, 0.5),
        bA ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ),
    data=foxes_std
)
m5h1_b <- quap(
    alist(
        weight ~ dnorm(mu, sigma),
        mu <- a + bG * groupsize,
        a ~ dnorm(0, 0.5),
        bG ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ),
    data=foxes_std
)

# Now choose values to iterate over
summary(foxes_std)
# Will just use Richard's way of limits +/- 0.15

area_seq <- seq(min(foxes_std$area)-0.15,
                max(foxes_std$area)+0.15,
                length.out=100)
mu_a <- link(m5h1_a, data=list(area=area_seq))
mu_a_mean <- colMeans(mu_a)
mu_a_PI <- apply(mu_a, 2, PI, 0.95)

groupsize_seq <- seq(min(foxes_std$groupsize)-0.15,
                max(foxes_std$groupsize)+0.15,
                length.out=100)
mu_b <- link(m5h1_b, data=list(groupsize=groupsize_seq))
mu_b_mean <- colMeans(mu_b)
mu_b_PI <- apply(mu_b, 2, PI, 0.95)

# Group size looks correlated with weight
# But not area
plot(coeftab(m5h1_a, m5h1_b))

# No effect at all for area!
plot(weight ~ area, data=foxes_std, col=rangi2)
lines(area_seq, mu_a_mean, lwd=2)
shade(mu_a_PI, area_seq)

# Still a weak effect but stronger!
plot(weight ~ groupsize, data=foxes_std, col=rangi2)
lines(groupsize_seq, mu_b_mean, lwd=2)
shade(mu_b_PI, groupsize_seq)

# 5H2
m5h1_c <- quap(
    alist(
        weight ~ dnorm(mu, sigma),
        mu <- a + bA * area + bG * groupsize,
        a ~ dnorm(0, 0.5),
        bA ~ dnorm(0, 0.5),
        bG ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ),
    data=foxes_std
)

# Both are now more significant!
# This looks like masking!
# NB: I looked at coefficients rather than 'counter factual' plots
# Since the main purpose of the counter factual is to see the slope,
# which this plot shows already
plot(coeftab(m5h1_a, m5h1_b, m5h1_c), pars=c("bA", "bG"))

# 5H3
m5h3_a <- quap(
    alist(
        weight ~ dnorm(mu, sigma),
        mu <- a + bF * avgfood + bG * groupsize,
        a ~ dnorm(0, 0.5),
        bF ~ dnorm(0, 0.5),
        bG ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ),
    data=foxes_std
)

m5h3_b <- quap(
    alist(
        weight ~ dnorm(mu, sigma),
        mu <- a + bA * area + bF * avgfood + bG * groupsize,
        a ~ dnorm(0, 0.5),
        bA ~ dnorm(0, 0.5),
        bF ~ dnorm(0, 0.5),
        bG ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ),
    data=foxes_std
)

plot(coeftab(m5h1_a, m5h1_b, m5h1_c, m5h3_a, m5h3_b), pars=c("bA", "bG", "bF"))

# Questions: 
# a) Is avgfood or area a better predictor of body weight?
# b) When avgfood and area are in same model their effects are reduced

# a) So this is between m5h1_c (area & groupsize) vs m5h3_a (avgfood & groupsize)
# In these, area had 0.41, while avgfood had 0.48
# However, the uncertainty is larger on avgfood
# I'd go with avgfood
coeftab(m5h1_a, m5h1_b, m5h1_c, m5h3_a, m5h3_b)

# b) Looks like confounding again! But in what way?
# 
