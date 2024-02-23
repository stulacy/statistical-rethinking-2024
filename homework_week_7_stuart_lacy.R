library(rethinking)
library(tidyverse)
library(dagitty)


# Question 1 --------------------------------------------------------------
data("bangladesh")
# But district is non-contiguous
max(bangladesh$district)
length(unique(bangladesh$district))

# So make it so!
bangladesh$district_id <- as.integer(as.factor(bangladesh$district))

# Use the dag from lecture 13
dag <- dagitty("dag { 
    C <- A -> K;
    K -> C;
    K <- U -> C;
    U <- D -> C;
}")
coordinates(dag) <- list(
    x=c(C=3, A=0, D=6, K=1, U=5),
    y=c(C=0, A=1, D=1, K=3, U=3)
)
drawdag(dag)

# We want to estimate the direct effect of urban residence on contraceptive usage
# The lecture estimated the TOTAL effect by conditioning on G
# As demonstrated by the adjustment set
adjustmentSets(dag, exposure="U", outcome="C", effect = "total")

# For the DIRECT effect we need to condition on Age and Kids as well
# Is this as straightforward as just including them in the model as covariates?
# In the lecture this allowed each district to have a separate urban effect, but with a multi-level prior
# However, how to accomplish this for continuous Age and Kids?
# Is just having separate covariates enough? Could include a multi-level prior on Age and Kids
# by district again? In which case could make a correlated multi-variate model with intercept, slope urban, slope age, slope kids
adjustmentSets(dag, exposure="U", outcome="C", effect = "direct")

# Children is discrete rather than continuous, but only has values 1-4 (NB: biased by not having women with no children?!)
# Ideally would use ordinal predictor... will get working as continuous first
dens(bangladesh$living.children)
# Age is continuous and already centered (although will want to scale to get SD=1)
dens(bangladesh$age.centered)

# This is the original model
m_urban <- ulam(
    alist(
        contraception ~ bernoulli(p),
        logit(p) <- a[district_id] + b[district_id]*urban,
        
        # Get parameters back on original scale
        transpars>vector[60]:a <<- abar[1] + v[, 1],
        transpars>vector[60]:b <<- abar[2] + v[, 2],
        transpars> matrix[60, 2]:v <- compose_noncentered(sigma, L_Rho, Z),
        
        # non-centered priors
        matrix[2, 60]:Z ~ normal(0, 1),
        vector[2]: abar ~ normal(0, 1),
        cholesky_factor_corr[2]:L_Rho ~ lkj_corr_cholesky(4),
        vector[2]: sigma ~ exponential(1),
        
        # Convert Cholesky to Corr matrix
        gq> matrix[2, 2]:Rho <<- Chol_to_Corr(L_Rho)
    ), data=bangladesh |> select(contraception=use.contraception, district_id, urban),
    cores=4, chains=4, log_lik = TRUE
)
# The average Urban effect is the second slope, i.e. abar[2] = 0.67 with a strong positive effect
# This is a log-odds increase
precis(m_urban, depth=2)
# Or a 95% relative increase in chance of contraception, which is very high
exp(0.67)

# What happens when adjusting for age and kids (both scaled)?
m_urban_direct <- ulam(
    alist(
        contraception ~ bernoulli(p),
        logit(p) <- a[district_id] + b[district_id]*urban + d[district_id]*age + e[district_id]*kids,
        
        # Get parameters back on original scale
        transpars>vector[60]:a <<- abar[1] + v[, 1],
        transpars>vector[60]:b <<- abar[2] + v[, 2],
        transpars>vector[60]:d <<- abar[3] + v[, 3],
        transpars>vector[60]:e <<- abar[4] + v[, 4],
        transpars> matrix[60, 4]:v <- compose_noncentered(sigma, L_Rho, Z),
        
        # non-centered priors
        matrix[4, 60]:Z ~ normal(0, 1),
        vector[4]: abar ~ normal(0, 1),
        cholesky_factor_corr[4]:L_Rho ~ lkj_corr_cholesky(4),
        vector[4]: sigma ~ exponential(1),
        
        # Convert Cholesky to Corr matrix
        gq> matrix[4, 4]:Rho <<- Chol_to_Corr(L_Rho)
    ), data=bangladesh |> 
        mutate(age = as.numeric(scale(age.centered)), kids=as.numeric(scale(living.children))) |>
        select(contraception=use.contraception, district_id, urban, age, kids),
    cores=4, chains=4, log_lik = TRUE
)
# Now the urban effect is 0.75, so it's actually got stronger!
precis(m_urban_direct, depth=2)
# This is a 2.1 relative effect, which is again quite high!
# In these terms it's definitely stronger than before, but it hasn't really changed the interpretation, which 
# I guess means that the effect of urban living isn't mediated through children or age
# HOWEVER: how is it possible to have a direct effect > total ?!
exp(0.75)

# Let's look at the overlap
# Yep this also shows the direct effect as being higher...
post_urban$abar[, 2]
post_direct <- extract.samples(m_urban_direct, pars=c("abar"))
tibble(
    sample=1:2000,
    total=post_urban$abar[, 2],
    direct=post_direct$abar[, 2]
) |>
    pivot_longer(-sample) |>
    ggplot(aes(x=value, fill=name)) +
        geom_density(alpha=0.4) +
        theme_classic() +
        theme(legend.position = "bottom") +
        labs(x="Log-odds of urban living on contraception", y="Density")

# Question 2 --------------------------------------------------------------
# Using the same dag, what is the (presumed direct?) effect of number of surviving children on contraceptive use?
drawdag(dag)
# Both total and direct have the same adjustment sets!
# Which incidentally neglect D...
# How to handle this? My entire model is centered around clustering on D
adjustmentSets(dag, exposure="K", outcome="C", effect = "total")
adjustmentSets(dag, exposure="K", outcome="C", effect = "direct")

# The solutions say that the adjustment set is A and U, but we want to include D as a competing cause...
# abar[4] has the average kid effect over districts at 0.53 with a tight CI
precis(m_urban_direct, depth=2, pars="abar")
# Each additional kid increases the relative log-odds by 1.7, which is a moderate effect
# NB: the solutions uses an ordinal predictor here which gives nice non-linear results, as the effect
# from 1-2 is massive, but dramatically drops off after that
exp(0.53)

# Question 3 --------------------------------------------------------------
# Write a synthetic data simulation for this problem set
# Simulation should include:
#  - district
#  - urban
#  - children
#  - age
#  - contraceptive use
# Best option is to make a dynamic population simulation where individuals age and produce children...
# Could give a decent survival model, but would also need a model for producing children...
# The solution states that if take an agent based approach then can have the feedback from using contraception
# on odds of having more children
# I actually wrote a Discrete Event Simulation Rcpp code that could _potentially_ run this
# If I think of it as states of infancy -> motherhood -> menopause
# I'm not sure if it could discrete predictors such as the probability of being urban / rural or contraceptive use
# It would be a much easier mental model to code this as an agent based simulation, although that could
# be much slower to run



