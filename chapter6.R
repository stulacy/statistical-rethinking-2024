library(rethinking)
library(dagitty)
library(tidyverse)

# 
set.seed(1914)
N <- 200 # num grant proposals
p <- 0.1 # proportion to select
nw <- rnorm(N) # newsworthiness - uncorrelated
tw <- rnorm(N) # trustworthiness - uncorrelated
s <- nw + tw # total score is combination
q <- quantile(s, 1-p) # top 10% threshold
selected <- ifelse(s >= q, TRUE, FALSE)
# And yet, get correlation!
# This _isn't_ confounding as we don't have a common factor!
# So what is the dag?
cor(tw[selected], nw[selected]) 

# We're conditioning on selected, which makes newsworthiness and trustworthiness appear correlated
# COLLIDER BIAS!
dag_1 <- dagitty("dag { 
    newsworthiness -> score
    trustworthiness -> score
    score -> selected
}")
drawdag(dag_1)


# 6.1 Multicollinearity ---------------------------------------------------
# very strong correlation between 1 or more predictor variables
# helpful for prediction, less so for inference!
# Will suggest that none of variables are associated with outcome even though they are

# Predicting height using length of legs as predictors
N <- 100
set.seed(909)
height <- rnorm(N, 10, 2)
leg_prop <- runif(N, 0.4, 0.5) # leg as proportion of height
leg_left <- leg_prop*height + rnorm(N, 0, 0.02) # add error so legs aren't same size exactly
leg_right <- leg_prop*height + rnorm(N, 0, 0.02)
d <- tibble(height, leg_left, leg_right)
# Beta is average height (10) / average leg height (4.5) = 2.2

# Very vague priors!
m6.1 <- quap(
    alist(
        height ~ dnorm(mu, sigma),
        mu <- a + bl*leg_left + br*leg_right,
        a ~ dnorm(10, 100),
        bl ~ dnorm(2, 10),
        br ~ dnorm(2, 10),
        sigma ~ dexp(1)
    ),
    data=d
)
# note how both bl and br aren't equal to 2.2!
precis(m6.1)
# And huge uncertainty!
plot(precis(m6.1))

post <- extract.samples(m6.1)
# Note that the posterior are perfectly anti-correlated!
plot(bl ~ br, post, col=col.alpha(rangi2, 0.1), pch=16)

# The model is: height ~ N(a + bl*xl + br*xr)
# But xl ~= xr
# Which is effectively height ~ N(a + (bl + br) * x) where x ~= xl ~= xr
# So their combined density is around the expected value of 2.2
dens(post$bl + post$br)

# If just use one leg get the expected value
# And much smaller uncertainty!
m6.2 <- quap(
    alist(
        height ~ dnorm(mu, sigma),
        mu <- a + bl*leg_left,
        a ~ dnorm(10, 100),
        bl ~ dnorm(2, 10),
        sigma ~ dexp(1)
    ),
    data=d
)
plot(precis(m6.2))

# 6.1.2 Multicollinear milk -----------------------------------------------
# Predicting calorific content of milk from % fat and % lactose
data(milk)
d <- milk
d$K <- scale(d$kcal.per.g)
d$F <- scale(d$perc.fat)
d$L <- scale(d$perc.lactose)

# Firstly fit bivariate (I would call these univariate as there's a single 
# predictor but I appreciate that's technically incorrect) models
m6.3 <- quap(
    alist(
        K ~ dnorm(mu, sigma),
        mu <- a + bF*F,
        a ~ dnorm(0, 0.2),  # Usual priors! seems a shame to use the same weakly informative priors after going on about how useful it is to impart scientific knowledge
        bF ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), data=d
)

m6.4 <- quap(
    alist(
        K ~ dnorm(mu, sigma),
        mu <- a + bL*L,
        a ~ dnorm(0, 0.2),  # Usual priors! seems a shame to use the same weakly informative priors after going on about how useful it is to impart scientific knowledge
        bL ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), data=d
)

# Fat is positively correlated with calories, lactose negative
coeftab(m6.3, m6.4)
# Looking at the pairs plot and fat is negatively correlated with lactose!
pairs( ~ kcal.per.g + perc.fat + perc.lactose, data=d, col=rangi2)

# So what happens when we include them both in a model?

m6.5 <- quap(
    alist(
        K ~ dnorm(mu, sigma),
        mu <- a + bL*L + bF * F,
        a ~ dnorm(0, 0.2),
        bL ~ dnorm(0, 0.5),
        bF ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), data=d
)
# The effects get weaker! And closer to zero!
# Is this masking again?
# It's got the same DAG as masking, but in masking the 2 direct effects had different
# directions and so cancelled out (which makes sense actually if the overall effect is the sum of both)
# Or is it more because the 2 variables are strongly correlated?
plot(coeftab(m6.3, m6.4, m6.5), pars=c('bF', 'bL'))

dag_m65 <- dagitty("dag {
  L[o]
  D[u]
  F[e]
  K[o]
  D -> L
  D -> F
  L -> K
  F -> K
}")
drawdag(dag_m65)

# Overthinking 6.12
# Fits a multivariate model of calories as a function of fat + a correlated parameter
# The correlated parameter is swept over values from r=0, ..., 0.99
# And for each r value we fit 100 models and calculate the average standard error
# Just check what the sqrt(diag(vcov(m)))[2] line is doing:
foo <- lm(kcal.per.g ~ perc.fat + perc.lactose, data=milk)
vcov(foo)  # Gets the covariance matrix
diag(vcov(foo)) # Gets the standard errors (as variance)
sqrt(diag(vcov(foo))) # Gets standard errors (as stddev)
sqrt(diag(vcov(foo)))[2] # Gets the corresponding value for the 2nd parameter, perc.fat
summary(foo) # This matches up with what we see using summary

# So the simulation shows that uncertainty in a parameter increases with collinearity
d <- milk
sim.coll <- function(r=0.9) {
    d$x <- rnorm(nrow(d), mean=r*d$perc.fat,
                 sd=sqrt((1-r**2)*var(d$perc.fat)))
    m <- lm(kcal.per.g ~ perc.fat + x, data=d)
    sqrt(diag(vcov(m)))[2] # stddev of perc.fat
}
rep.sim.coll <- function(r=0.9, n=100) {
    stddev <- replicate(n, sim.coll(r))
    mean(stddev)
}
r.seq <- seq(0, 0.99, by=0.01)
sttdev <- sapply(r.seq, function(z) rep.sim.coll(r=z, n=100))
plot(sttdev ~ r.seq, type='l', col=rangi2, lwd=2, xlab="correlation")


# 6.2 Post-treatment bias -------------------------------------------------
set.seed(71)
# Simulating plants that are treated and the effect on fungus and their final height
N <- 100  # Number of plants
h0 <- rnorm(N, 10, 2)  # Initial heights
# assign treatments and simulate fungus and growth
treatment <- rep(0:1, each=N/2)
table(treatment)  # evenly split
fungus <- rbinom(N, size=1, prob=0.5 - treatment*0.4)
table(fungus) # fungus is more likely in untreated plants
h1 <- h0 + rnorm(N, 5-3*fungus) # final height is a function of fungus, NOT treatment
hist(h1)  
# But if include fungus in the model get biased results because it is a POST-TREATMENT effect
dag_m62 <- dagitty("dag {
  T -> F
  F -> H1
  H0 -> H1
}")
# And this is the DAG
drawdag(dag_m62)

d <- tibble(h0=h0, h1=h1, treatment=treatment, fungus=fungus)
precis(d)

# Choosing priors
# model h1 ~ N(h0*p, sigma)
# where h1 = final height, h0=start height, p = h1/h0
# p = 3 = tripled in height
# want p > 0, but can be less than 1 (reduced in height)
# Use lognormal(0, 0.25)
dens(rlnorm(1e4, 0, 0.25))

# Fit model!
m6.6 <- quap(
    alist(
        h1 ~ dnorm(mu, sigma),
        mu <- h0*p,
        p ~ dlnorm(0, 0.25),
        sigma ~ dexp(1)
    ), data=d
)
# 40% growth on average
precis(m6.6)

# Now want to include treatment and fungus in the model, as obviously want to
# know the treatment effect
# p = a + bT * T + bF * F
# Where T = Treatment and F = Fungus
# alpha still has the log-normal prior, but betas can be in any direction
m6.7 <- quap(
    alist(
        h1 ~ dnorm(mu, sigma),
        mu <- h0*p,
        p <- a + bT*treatment + bF*fungus,
        a ~ dlnorm(0, 0.25),
        bT ~ dnorm(0, 0.5),
        bF ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), data=d
)
# a is the same, as this is when there is no treatment and no fungus
# 0 treatment effect!
# But strong negative fungus effect
precis(m6.7)
# This because this model answers the question: "when we know if a plant developed fungus,
# does knowing the soil treatment change our guess of its final height?"
# And the answer is NO because we expect treatment to effect height THROUGH fungus

# So to actually see the treatment effect we can just remove fungus from the model
m6.8 <- quap(
    alist(
        h1 ~ dnorm(mu, sigma),
        mu <- h0*p,
        p <- a + bT*treatment,
        a ~ dlnorm(0, 0.25),
        bT ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), data=d
)
# Now we see that treatment has a positive effect on the final height
precis(m6.8)

# Going back to the dag
# And we can say that if we include F (post-treatment effect) in the model
# Then we block the path from treatment to outcome
# Or "conditioning on F induces D-Separation
drawdag(dag_m62)

# The third one here is the problem, T is independent of H1
# when we condition on F. This isn't a problem with the model or dag,
# these are correct! The problem only occurs when we naively add F 
# into the model
# The other 2 are useful diagnostics, i.e. are H0 and F independent,
# and are H0 and T independent? Would quickly test in a diagnostic bivariate
# model. If these fail then the DAG is wrong!
# The third one we can test as well, which we did in model m6.7 and sure enough
# the test passed so the DAG is accurate
impliedConditionalIndependencies(dag_m62)

# Post-treatment conditioning can also thinking a treatment DOES work when it
# DOESNT
# Below, fungus doesn't effect growth, but M effects growth and it also effects
# fungus. This could be moisture. 
# H1 ~ T shows no effect, but if you condition on T it will show an association
# (I think this was similar to a previous chapter example)
dag_m63 <- dagitty("dag {
  M[u]
  H0 -> H1
  M -> H1
  M -> F
  T -> F
}")
drawdag(dag_m63)

set.seed(71)
N <- 1e3
h0 <- rnorm(N, 10, 2)
treatment <- rep(0:1, each=N/2)
M <- rbern(N)
fungus <- rbinom(N, size=1, prob=0.5 - treatment*0.4 + 0.4*M)
h1 <- h0 + rnorm(N, 5+3*M)
d2 <- tibble(h0, h1, treatment, fungus)

m6.7_b <- quap(
    alist(
        h1 ~ dnorm(mu, sigma),
        mu <- h0*p,
        p <- a + bT*treatment + bF*fungus,
        a ~ dlnorm(0, 0.25),
        bT ~ dnorm(0, 0.5),
        bF ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), data=d2
)
m6.8_b <- quap(
    alist(
        h1 ~ dnorm(mu, sigma),
        mu <- h0*p,
        p <- a + bT*treatment,
        a ~ dlnorm(0, 0.25),
        bT ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), data=d2
)
# So just treatment (6.8_b) had no effect
# But adding fungus (6.7_b) means it looks like it did have an effect!
plot(coeftab(m6.7_b, m6.8_b), pars=c("bT", "bF"))

# 6.3 Collider bias -------------------------------------------------------
# "When conditionon a collider it creates statistical associations amongst
# its cause"
drawdag(dagitty("dag { T -> S <- N }"))

# 6.3.1 Collider of false sorrow ------------------------------------------
# Example of marriage, age, and happiness
# Marriage is a collider
# happiness causes marriage
# But increased age also causes marriage
# If we condition on marriage, it will look like age and happiness
# are correlated
drawdag(dagitty("dag { H -> M <- A }"))

# Richard coded an agent based simulation in the package
d <- sim_happiness(seed=1977, N_years=1000)
precis(d)

# mu = alpha[married] + bA * Age
# restrict to adults and standardize to [0, 1]
# This is done because we imagine a scenario where happiness is at a maximum at
# age 18 and minimum at age 65 EVEN THOUGH we simulated a model where age
# doesn't influece happiness directly
# So we want age 0 to be 18, and then useful to have 0-1 rather than 0-(65-18)
d2 <- d[d$age > 17, ]
d2$A <- (d2$age - 18) / (65 - 18)
precis(d2)

# Happiness is [-2, 2]
# So strongest relationship is 4 ((2- -2) / 1)
# so set prior to half of 4
# for intercept is happiness at age 0, so give a prior that covers
# most of the range 
# Why doesn't this work?! Using same seed as book, and won't work with any other seeds!
d2$mid <- d2$married + 1
m6.9 <- quap(
    alist(
        happiness ~ dnorm(mu, sigma),
        mu <- a[mid] + bA*A,
        a[mid] ~ dnorm(0, 1),
        bA ~ dnorm(0, 2),
        sigma ~ dexp(1)
    ), data=d2
)
precis(m6.9, depth=2)
# with happiness

# But when remove marriage status there is no age effect
# The association was because we conditioned on a collider!
m6.10 <- quap(
    alist(
        happiness ~ dnorm(mu, sigma),
        mu <- a + bA*A,
        a ~ dnorm(0, 1),
        bA ~ dnorm(0, 2),
        sigma ~ dexp(1)
    ), data=d2
)
precis(m6.10)

# 6.3.2 The haunted DAG ---------------------------------------------------
# Resume from 6.25 on p185
# Looking at effect on children's education from grandparents, parents, and common
# unobserved factors
# This diagram in particular!
dag_632 <- dagitty("dag {
    U[u]
    G -> P
    G -> C
    P -> C
    U -> P
    U -> C
}")
# We want to identify G->C and P->C
# BUT! If we condition of P, it will bias inference G->C (how? this isn't a collider)
drawdag(dag_632)

# Simulate this model
N <- 200
b_GP <- 1 # G->P
b_GC <- 0 # No real effect of G on C
b_PC <- 1 # P->C
b_U <- 2  # C <- U -> P

set.seed(1)
# Describe U as a categorical variable
U <- 2*rbern(N, 0.5) - 1
G <- rnorm(N)
P <- rnorm(N, b_GP*G + b_U * U)
C <- rnorm(N, b_PC * P + b_GC*G + b_U *U)
d <- tibble(C, P, G, U)

m6.11 <- quap(
    alist(
        C ~ dnorm(mu, sigma),
        mu <- a + b_PC * P + b_GC*G,
        a ~ dnorm(0, 1),
        c(b_PC, b_GC) ~ dnorm(0, 1),
        sigma ~ dexp(1)
    ), 
    data=d
)
# Parental effect is far too high (defined as 1!)
# There's also negative grandparent effect!
precis(m6.11)

# It is collider bias again, with P the collider, opening up the path from G to U to C
# Conditioning on parent infers something about U, i.e. if know that parent
# has had positive effect, then likely that U is the opposite which leads to
# negative children effect associated with grandparent

# THE ONLY WAY TO PREVENT THIS IS TO CONDITION ON U
m6.12 <- quap(
    alist(
        C ~ dnorm(mu, sigma),
        mu <- a + b_PC * P + b_GC*G + b_U*U,
        a ~ dnorm(0, 1),
        c(b_PC, b_GC, b_U) ~ dnorm(0, 1),
        sigma ~ dexp(1)
    ), 
    data=d
)
# Now direct effect of G->C is 0
precis(m6.12)

# 6.4 Confronting confounding ---------------------------------------------
# "...Confounding [is] any context in which the association between an outcome Y
# and a predictor of interest X is not as it would be, if we had experimentally
# determined the values of X.
dag_64a <- dagitty("dag { W <- E <- U -> W} ")
# Here, there are 2 paths connecting E and W
# NB THE DIRECTION OF THE ARROWS DOESN'T MATTER
drawdag(dag_64a)
# So if we model E ~ W, we'll get a biased result
# as it will be _confounded_ by U
# However! only W->E is causal; E <- U -> W is a statistical
# association only

# To highlight the E->W effect we can either control U in an experiment,
# or control for it in a statistical model

# 6.4.1 Shutting the backdoor ---------------------------------------------
# Want to block all confounding paths (again regardless of direction of arrow)

# The Fork
# Z causes a correlation between X and Y. Control for Z to see that X and Y are
# independent
dag_1 <- dagitty("dag { X <- Z -> Y }")
coordinates(dag_1) <- list(x=c(X=0, Z=1, Y=2), y=c(X=0, Z=1, Y=0))
drawdag(dag_1)

# The Pipe
# post-treatment bias example. If condition on Z then we block the path between
# X and Y
dag_2 <- dagitty("dag { X -> Z -> Y }")
coordinates(dag_2) <- list(x=c(X=0, Z=1, Y=2), y=c(X=0, Z=1, Y=2))
drawdag(dag_2)

# The Collider
# Opposite of the fork: conditioning here INTRODUCES a correlation between X and Y
dag_3 <- dagitty("dag { X -> Z <- Y }")
coordinates(dag_3) <- list(x=c(X=0, Z=1, Y=2), y=c(X=1, Z=0, Y=1))
drawdag(dag_3)

# The Descendant
# Conditioning on a descendant is like conditioning on the parent variable but weaker
# So in this example controlling for D will partially open the path X-Y
# And likewise for Pipes and Forks
dag_4 <- dagitty("dag { X -> Z <- Y 
                 Z -> D}")
coordinates(dag_4) <- list(x=c(X=0, Z=1, Y=2, D=1), y=c(X=5, Z=0, Y=5, D=4))
drawdag(dag_4)

# Fortunately can identify which variables need to control for to analyse
# specified paths
# Algorithm:
# 1) List all paths from X->Y
# 2) Classify the paths as open or closed (path is open UNLESS contains collider)
# 3) Classify the paths as backdoor or not (backdoor has arrow into X)
# 4) All backdoor paths that are open need closing, identify which variable to condition on

# dagitty can do this for us!
library(dagitty)
dag_6.1 <- dagitty( "dag {
U [unobserved]
X -> Y
X <- U <- A -> C -> Y
U -> B <- C
}")
coordinates(dag_6.1) <- list(x=c(X=0, U=0, A=1, B=1, C=2, Y=2),
                             y=c(X=3, U=1, B=2, A=0, C=1, Y=3))
drawdag(dag_6.1)
# I.e. can condition on A OR C (we told dagitty U is unobserved)
# C better because it "could help precision of X-Y" (how?!)
adjustmentSets( dag_6.1 , exposure="X" , outcome="Y" )

# The data alone can never tell us when a DAG is right, only when a DAG is WRONG

# Problems ----------------------------------------------------------------

# 6E1: 
#   - Multicollinearity (including correlated predictors in the same model dampening direct effects)
#   - Post treatment bias (controlling for an effect between T and outcome nullifies T effect)
#   - Collider bias (controlling for a collider INTRODUCES an association between 2 exposures that isn't causal)

# 6E2: Multi-collinearity: if include both NO2 and O3 in a model on health outcomes will get biased results
# because NO2 and O3 are anticorrelated

# 6E3:
# Fork: X _||_ Y | Z
impliedConditionalIndependencies(dagitty("dag { X <- Z -> Y }"))

# Pipe: X _||_ Y | Z
impliedConditionalIndependencies(dagitty("dag { X -> Z -> Y }"))

# The Collider
# X _||_ Y
impliedConditionalIndependencies(dagitty("dag { X -> Z <- Y }"))

# The Descendant
# Multiple!
# This is a descendant of a collider so X _||_ Y
# And 2 ones on D:
# D _||_ X | Z
# D _||_ Y | Z
impliedConditionalIndependencies(dagitty("dag { X -> Z <- Y 
                 Z -> D}"))

# 6E4
# Because there should be no association between 2 exposures, but conditioning
# on something has introduced one (i.e. people who answer phone surveys)

# 6M1
dag_6.1_b <- dagitty( "dag {
U [u]
V [u]
V -> C
V -> Y
X -> Y
X <- U <- A -> C -> Y
U -> B <- C
}")
coordinates(dag_6.1_b) <- list(x=c(X=0, U=0, A=1, B=1, C=2, Y=2, V=3),
                             y=c(X=3, U=1, B=2, A=0, C=1, Y=3, V=2))
drawdag(dag_6.1_b)

# Backdoor Paths from X to Y aside from direct:
# X->U->B->C->Y (SHUT, B is a colliders)
# X->U->B->C->V->Y (SHUT, B and C are colliders)
# X->U->A->C->Y (OPEN)
# X->U->A->C->V->Y (SHUT, C is a collider)

# Could condition on A only, as C is a collider through V
# And dagitty agrees!
adjustmentSets(dag_6.1_b, exposure="X", outcome="Y")

# 6M2
foo <- dagitty("dag { X->Z->Y}")
coordinates(foo) <- list(x=c(X=0, Z=1, Y=2), y=c(X=0, Z=1, Y=2))
drawdag(foo)

N <- 100
X <- rnorm(N)
Z <- rnorm(N, X)
Y <- rnorm(N, Z)

# High correlation between the direct effects
cor(X, Z)
cor(Z, Y)

# Relatively high correlation on X->Y (but we know this is through Z)
cor(X, Y)

m6m2 <- quap(
    alist(
        Y ~ dnorm(mu, sigma),
        mu <- a + bX * X + bZ*Z,
        a ~ dnorm(0, 0.2),
        c(bX, bZ) ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), 
    data=tibble(X, Y, Z)
)
# Model has the direct effect of X on Y as very low!
precis(m6m2)

# Which fits with the implied conditional probabilities
impliedConditionalIndependencies(foo)


# The difference here is that there is a direct causal effect of X on Z, rather
# than being mediated through a confounder

# 6M3
# a) Backdoor paths:
# X->Z->Y (open)
# X->Z->A->Y (open)
# So condition on Z to shut both paths
# Check answer with dagitty: correct!
adjustmentSets(dagitty("dag { 
  X <- Z <- A -> Y
  Z -> Y
  X -> Y
}"), exposure="X", outcome="Y")

# b) Backdoor paths:
# None!
# Correct!
adjustmentSets(dagitty("dag { 
  X -> Z <- A -> Y
  Z -> Y
  X -> Y
}"), exposure="X", outcome="Y")

# c) Backdoor paths:
# X->A->Z->Y (shut as Z is a collider)
# No need to condition on anything
# Correct!
adjustmentSets(dagitty("dag { 
  X <- A -> Z <- Y
  X -> Z
  X -> Y
}"), exposure="X", outcome="Y")

# d) Backdoor paths:
# X->A->Z->Y (open)
# So condition on A or Z
# Oh, just A, I wonder why
# Oh, because Z is a sneaky collider
adjustmentSets(dagitty("dag { 
  X <- A -> Z -> Y
  X -> Z
  X -> Y
}"), exposure="X", outcome="Y")

# 6H1
data("WaffleDivorce")
d <- WaffleDivorce
d$S <- as.numeric(d$South) + 1
d$A <- scale(d$MedianAgeMarriage)
d$M <- scale(d$Marriage)
d$W <- scale(d$WaffleHouses)
d$D <- scale(d$Divorce)
# Firstly draw the dag from the book
dag_6h1 <- dagitty("dag {
  D <- A <- S -> W -> D
  S -> M <- A
  M -> D
}")
coordinates(dag_6h1) <- list(x=c(A=0, S=0, M=1, W=2, D=2),
                             y=c(A=2, S=0, M=1, W=0, D=2))
drawdag(dag_6h1)

# What are the implied conditional probs?
impliedConditionalIndependencies(dag_6h1)
# Ok let's test all these!
# A _||_ W | S : Marriage age is independent of waffle houses given Southernness
# Yep that's independent all right
quap(
    alist(
        A ~ dnorm(mu, sigma),
        mu <- a[S] + bW * W,
        a[S] ~ dnorm(0, 0.2),
        c(bW) ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), 
    data=d
) |> precis(depth=2)

# M _||_ W | S : Marriage rate is independent of waffle houses given Southernness
# Also independent!
quap(
    alist(
        M ~ dnorm(mu, sigma),
        mu <- a[S] + bW * W,
        a[S] ~ dnorm(0, 0.2),
        c(bW) ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), 
    data=d
) |> precis(depth=2)

# The big one:
# D _||_ S | A, M, W : Divorce rate is independent of Southernness given
# age @ marriage, marriage rate, and wafflehouses
# Yep both southernesses are independent of divorce rate with these controls!
quap(
    alist(
        D ~ dnorm(mu, sigma),
        mu <- a[S] + bW * W + bM * M + bA * A,
        a[S] ~ dnorm(0, 0.2),
        c(bW, bM, bA) ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), 
    data=d
) |> precis(depth=2)

# In which case to estimate W->D
# Can either condition on A & M, or just W
adjustmentSets(dag_6h1, exposure="W", outcome="D")
# Try both ways and compare estimates

m6h1_a <- quap(
    alist(
        D ~ dnorm(mu, sigma),
        mu <- a[S] + bW * W,
        a[S] ~ dnorm(0, 0.2),
        c(bW) ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), 
    data=d
)
m6h1_b <- quap(
    alist(
        D ~ dnorm(mu, sigma),
        mu <- a + bW * W + bM * M + bA * A,
        a ~ dnorm(0, 0.2),
        c(bW, bM, bA) ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), 
    data=d
)

# Not significant in either model, although subtly different estimates 
# (stronger effect when adjusting for M & A rather than just S)
plot(coeftab(m6h1_a, m6h1_b))

# 6H2
# Did this as part of 6H1!
# None of the tests failed for me (should they have?!)
# I would say that having Southerness influence many of these factors but not 
# divorce rate directly seems a bit strange. If the south is really that different
# in attitudes to marriage (reflected in different marriage rate and age at marriage),
# then I wouldn't be surprised if it had a direct association to D

# 6H3
dag_6h3 <- dagitty("dag {
    A -> F -> G -> W <- F
                   }")
coordinates(dag_6h3) <- list(x=c(F=0, A=1, W=1, G=2),
                               y=c(F=1, A=0, W=2, G=1))
drawdag(dag_6h3)

# No adjustment needed, can infer directly
# If adjusted for F would be getting into post-treatment bias
adjustmentSets(dag_6h3, exposure = "A", outcome="W")

data(foxes)
d <- foxes
d$F <- scale(d$avgfood)
d$A <- scale(d$area)
d$G <- scale(d$groupsize)
d$W <- scale(d$weight)
summary(d)

m6h3 <- quap(
    alist(
        W ~ dnorm(mu, sigma),
        mu <- a + bA * A,
        a ~ dnorm(0, 0.2),
        c(bA) ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), 
    data=d
)
# Zero effect! 1 stdev change in area leads to 0.02 stdevs change in weight
precis(m6h3) |> plot()

# Prior predictive
prior <- extract.prior(m6h3)
mu <- link(m6h3, post=prior, data=list(A=c(min(d$A)-0.15, max(d$A)+0.15)))
plot(NULL, xlim=c(-2, 2), ylim=c(-2, 2))
# Plot the first 50 lines, looks sensible!
# Just using the standard a ~ N(0, 0.2), b ~ N(0, 0.5) priors been using
# throughout
for (i in 1:50) lines(c(-2, 2), mu[i, ], col=col.alpha("black", 0.4))

# 6H4
# Still no need to adjust for Food either
# Am I missing something in not needing to adjust for either of these 2 Qs?
# But there's no backdoor path so it makes sense
adjustmentSets(dag_6h3, exposure = "F", outcome="W")

# Also see no real effect here, just a slight decrease again of by 0.2
# Which doesn't seem right at all, surely adding food can only increase weight?!
quap(
    alist(
        W ~ dnorm(mu, sigma),
        mu <- a + bF * F,
        a ~ dnorm(0, 0.2),
        c(bF) ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), 
    data=d
) |> precis() |> plot()

# 6H5
# Just need to adjust for the backdoor path through G
# Would expect increasing group size to decrease weight!
adjustmentSets(dag_6h3, exposure = "G", outcome="W")

# Yep finally see a significant effect, only just though!
# And still far lower effect size than might expect (-0.16 stdevs in weight
# for 1 stdev increase in group size)
quap(
    alist(
        W ~ dnorm(mu, sigma),
        mu <- a + bG * G,
        a ~ dnorm(0, 0.2),
        c(bG) ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), 
    data=d
) |> precis() |> plot()

# The question then is how can we explain this?
# Is the DAG correct?
# Let's look at implied probs:
impliedConditionalIndependencies(dag_6h3)

# Ah, this fails, G does have an association (albeit small)
# A _||_ G | F
quap(
    alist(
        A ~ dnorm(mu, sigma),
        mu <- a + bG * G + bF * F,
        a ~ dnorm(0, 0.2),
        c(bG, bF) ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), 
    data=d
) |> precis() |> plot()

# A _||_ W | F
# This doesn't fail!
quap(
    alist(
        A ~ dnorm(mu, sigma),
        mu <- a + bW * W + bF * F,
        a ~ dnorm(0, 0.2),
        c(bW, bF) ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), 
    data=d
) |> precis() |> plot()

# So the problem is that A and G aren't independent like the DAG suggests they are 
drawdag(dag_6h3)

# It doesn't seem unreasonable to believe that the area affects the group size
# But that doesn't explain everything :
# Why doesn't increasing the amount of food available increase the wolves weight?
# Why doesn't increasing the area increase their weight?
# Why is group sizes effect so low?
# Why are area and group size effected?
# Can this be explained by:
#  - redrawing diagram (adding path from A to G)
#  - adding group number in some way?
#  - adding unobserved latent variable?
