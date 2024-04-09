library(rethinking)
library(tidyverse)

# 16.1 Geometric people ---------------------------------------------------
data(Howell1)
d <- Howell1 |> 
        as_tibble() |>
        mutate(
            w = weight / mean(weight),
            h = height / mean(height)
        )

m16.1 <- ulam(
    alist(
        w ~ dlnorm(mu, sigma),
        exp(mu) <- 3.141593 * k * p^2 * h^3,
        k ~ exponential(0.5),
        p ~ beta(2, 18),
        sigma ~ exponential(1)
    ), data=d, chains=4, cores=4, log_lik = TRUE
)
precis(m16.1)
post <- extract.samples(m16.1)

# Posterior for beta has shifted along
dens(rbeta(1e4, 2, 18))
dens(post$p, col="red", add=TRUE)

# So has k. We had 544 observations so a decent chunk
dens(rexp(1e4, 0.5))
dens(post$k, col="red", add=TRUE)

# k is the density i.e. 5 mean units of mass from an increase in p^2 or h^3 by 1...
# Easiest to just plot it rather than try to solve analytically (at least for me)
ndata <- tibble(h = seq(min(d$h), max(d$h), length.out=100))
mu <- exp(link(m16.1, data=ndata))
# Extremely small uncertainty around mu!
# There's 2 non-linearities at play here: height^3 and exp(), which is why I was getting
# lost trying to analytically identify effect from coefficients
# Also remember the non-identifiability of k & p
# Also doesn't work well at lower heights!
# Probably p changes, or even the full functional form
ndata |>
    mutate(
        height_cm = h * mean(d$height),
        mu_mean = colMeans(mu) * mean(d$weight),
        mu_lower = apply(mu, 2, PI)[1, ] * mean(d$weight),
        mu_upper = apply(mu, 2, PI)[2, ] * mean(d$weight)
    ) |>
    ggplot(aes(x=height_cm, y=mu_mean)) +
        geom_ribbon(aes(ymin=mu_lower, ymax=mu_upper), alpha=0.4, fill="orange") +
        geom_point(data=d |> rename(height_cm=height, mu_mean=weight)) +
        geom_line() +
        theme_classic() + 
        labs(x="Height (cm)", y="Weight (kg)")

# Can see massive correlation between k and p!
# Because we used informative priors we could actually fit these, otherwise it might have struggled
#pairs(m16.1)

# Can test that by using wider priors
m16.1_b <- ulam(
    alist(
        w ~ dlnorm(mu, sigma),
        exp(mu) <- 3.141593 * k * p^2 * h^3,
        p ~ beta(20, 20),
        k ~ exponential(0.1),
        sigma ~ exponential(1)
    ), data=d, chains=4, cores=4
)
# Huh, still fit it ok
#pairs(m16.1_b)

# If plot on the outcome scale with the additional uncertainty from sigma we can see the far
# greater uncertainty
mu <- sim(m16.1, data=ndata)
ndata |>
    mutate(
        height_cm = h * mean(d$height),
        mu_mean = colMeans(mu) * mean(d$weight),
        mu_lower = apply(mu, 2, PI)[1, ] * mean(d$weight),
        mu_upper = apply(mu, 2, PI)[2, ] * mean(d$weight)
    ) |>
    ggplot(aes(x=height_cm, y=mu_mean)) +
        geom_ribbon(aes(ymin=mu_lower, ymax=mu_upper), alpha=0.4, fill="orange") +
        geom_point(data=d |> rename(height_cm=height, mu_mean=weight)) +
        geom_line() +
        theme_classic() + 
        labs(x="Height (cm)", y="Weight (kg)")

# 16.2 Hidden minds and observed behaviour --------------------------------
data("Boxes")
precis(Boxes)
boxes <- as_tibble(Boxes) |>
    mutate(gender = as.factor(gender), majority_first=as.factor(majority_first),
           culture = as.factor(culture), y=as.factor(y))
summary(boxes)
# y is 1 for unchosen, 2 for majority, 3 for minority
# Would expect more unchosen than minority!
# majority first means were the 3 majority colours demonstrated before the 1 minority colour

# Simulate a scenario where half children choose randomly and half follow majority
set.seed(7)
N <- 30 # children
y1 <- sample(1:3, size=N/2, replace=TRUE)
y2 <- rep(2, N/2)

# Combine and shuffle
y <- sample(c(y1, y2))

# 73% for majority despite only half deliberately following it 
# Can imagine an infinite number of such scenarios in which we get majority chosen 45% of the time
# as in the study
sum(y==2)/N

# Going to consider 5 strategies that seem relatively plausible (could imagine others)
#  1. Follow majority
#  2. Follow minority
#  3. Choose the unchosen colour
#  4. Choose randomly
#  5. Choose colour that was demonstrated first

# Of course since we're Bayesian we'll assume children choose scenarios randomly
# So use a Dirichlet prior, weak-uniform for scenario probabilities
# For each strategy we have a prob of outcome: Pr(yi | sj) for j in 1-5
# Average over them (marginalize) with Pr(yi) = sum(pj * Pr(yi | sj))
# Is this the same as a catgorical likelihood?
# Doesn't look like it to me, as that is the p1*p2*p3 for each of the 3 outcomes
# Whereas here we are just calculating the probability of the outcome we observed...
# This is what we did with the missing data cat example
# Ah maybe it's because there are different forms of the Categorical likelihood:
# On wikipedia it shows 3 forms: https://en.wikipedia.org/wiki/Categorical_distribution
# The first is what we're using, just the p(x == i) = pi
# Whereas the other 2 calculate the probability for each outcome

# Here's the model in Stan!
data("Boxes_model")
cat(Boxes_model)

dat_list <- list(
    N=nrow(Boxes),
    y=Boxes$y,
    majority_first=Boxes$majority_first
)

m16.2 <- stan(model_code=Boxes_model, data=dat_list, chains=3, cores=3)

# Plot marginal posterior for p
# NB: what does marginal mean here?
p_labels <- c("1 Majority", "2 Minority", "3 Unchosen", "4 Random", "5 Follow first")
plot(precis(m16.2, depth=2), labels=p_labels)

# 16.3 Ordinary differential nut cracking ---------------------------------
data(Panda_nuts)
Panda_nuts <- as_tibble(Panda_nuts)
# Interesting in: nuts_opened, seconds (time taken), and age
# Interested in how this skill develops with age
summary(Panda_nuts)

# Prior predictive
N <- 1e4
phi <- rlnorm(N, log(1), 0.1)
k <- rlnorm(N, log(2), 0.25)
theta <- rlnorm(N, log(5), 0.25)

# Relative growth curve
map_dfr(1:30, function(i) {
    tibble(
        sample=i,
        age=seq(0, 1.5, length.out=100),
        mass=1-exp(-k[i]*age)
    )
}) |>
    ggplot(aes(x=age, y=mass, group=sample)) +
        geom_line() +
        theme_classic()

# nuts per second
map_dfr(1:30, function(i) {
    tibble(
        sample=i,
        age=seq(0, 1.5, length.out=100),
        mass=1-exp(-k[i]*age),
        nuts=phi[i] * (mass)**(theta[i])
    )
}) |>
    ggplot(aes(x=age, y=nuts, group=sample)) +
        geom_line() +
        theme_classic()

# Fit model
dat_list <- list(
    n = as.integer(Panda_nuts$nuts_opened),
    age=Panda_nuts$age/max(Panda_nuts$age),
    seconds=Panda_nuts$seconds
)

m16.4 <- ulam(
    alist(
        n ~ poisson(lambda),
        lambda <- seconds*phi*(1-exp(-k*age))^theta,
        phi ~ lognormal(log(1), 0.1),
        k ~ lognormal(log(2), 0.25),
        theta ~ lognormal(log(5), 0.25)
    ), data=dat_list, chains=4, cores=4
)

precis(m16.4)
# Decent mixing
traceplot(m16.4)

# Posterior developmental curve
ndata <- expand_grid(seconds=1, age=seq(0, 1, length.out=100))
mu <- link(m16.4, data=ndata)
ndata |>
    mutate(
        mu_mean = colMeans(mu),
        mu_lower = apply(mu, 2, PI)[1, ],
        mu_upper = apply(mu, 2, PI)[2, ],
        age = age * max(Panda_nuts$age),
        mu_mean = mu_mean,
        mu_lower = mu_lower,
        mu_upper = mu_upper
   ) |>
    ggplot(aes(x=age, y=mu_mean)) +
        geom_ribbon(aes(ymin=mu_lower, ymax=mu_upper), alpha=0.3) +
        geom_line(colour="steelblue", linewidth=1) +
        theme_classic() +
        geom_point(aes(y=nuts_per_second), data=Panda_nuts |> mutate(nuts_per_second = nuts_opened / seconds))


# 16.4 Population Dynamics ------------------------------------------------
data("Lynx_Hare")
Lynx_Hare <- as_tibble(Lynx_Hare)
Lynx_Hare |>
    pivot_longer(-Year) |>
    ggplot(aes(x=Year, y=value, colour=name)) +
        geom_point() +
        geom_line() +
        theme_classic() +
        scale_colour_brewer("", palette = "Dark2") +
        theme(
            legend.position = "bottom"
        )

# Simulate these dynamics in discrete time
sim_lynx_hare <- function(n_steps, init, theta, dt=0.002) {
    L <- rep(NA, n_steps)
    H <- rep(NA, n_steps)
    L[1] <- init[1]
    H[1] <- init[2]
    for (i in 2:n_steps) {
        H[i] <- H[i-1] + dt*H[i-1]*(theta[1] - theta[2] * L[i-1])
        L[i] <- L[i-1] + dt*L[i-1]*(theta[3] * H[i-1] - theta[4])
    }
    tibble(
        Lynx=L,
        Hare=H,
        step=1:n_steps
    )
}

# Run with some parameters
# Produces osciallating relationships
theta <- c(0.5, 0.05, 0.025, 0.5)
z <- sim_lynx_hare(1e4, as.numeric(Lynx_Hare[1, 2:3]), theta)
z |>
    pivot_longer(-step) |>
    ggplot(aes(x=step, y=value, colour=name)) +
        geom_point() +
        geom_line() +
        theme_classic() +
        scale_colour_brewer("", palette = "Dark2") +
        theme(
            legend.position = "bottom"
        )

# Prior predictives
N <- 1e4
Ht <- 1e4
p <- rbeta(1e4, 2, 18)
h <- rbinom(N, Ht, p)
h <- round(h/1000, 2)
dens(h, xlab="Thousands of pelts")

# Stan model with ODE
data("Lynx_Hare_model")
cat(Lynx_Hare_model)

dat_list <- list(
    N = nrow(Lynx_Hare),
    pelts=Lynx_Hare[, 2:3]
)

# One chain failed! Is that expected?
m16.5 <- stan(model_code=Lynx_Hare_model, data=dat_list, chains=3, cores=3, control=list(adapt_delta=0.95))

# The link and sim functions don't work here since these models weren't fitted with ulam
link(m16.5)
#sim(m16.5)

# Fortunately Richard saved the posterior predictions in GQ
post <- extract.samples(m16.5)
res <- Lynx_Hare |>
    mutate(
        lynx_mu_mean = colMeans(post$pelts_pred[, , 1]),
        lynx_mu_lower = apply(post$pelts_pred[, , 1], 2, PI)[1, ],
        lynx_mu_upper = apply(post$pelts_pred[, , 1], 2, PI)[2, ],
        hare_mu_mean = colMeans(post$pelts_pred[, , 2]),
        hare_mu_lower = apply(post$pelts_pred[, , 2], 2, PI)[1, ],
        hare_mu_upper = apply(post$pelts_pred[, , 2], 2, PI)[2, ],
        lynx_pop_mean = colMeans(post$pop[, , 1]),
        lynx_pop_lower = apply(post$pop[, , 1], 2, PI)[1, ],
        lynx_pop_upper = apply(post$pop[, , 1], 2, PI)[2, ],
        hare_pop_mean = colMeans(post$pop[, , 2]),
        hare_pop_lower = apply(post$pop[, , 2], 2, PI)[1, ],
        hare_pop_upper = apply(post$pop[, , 2], 2, PI)[2, ]
    ) |>
    rename(lynx_mu_actual=Lynx, hare_mu_actual=Hare) |>
    pivot_longer(-c(Year), names_pattern="(.+)_(.+)_(.+)", names_to=c("animal", "outcome", "statistic"))  |>
    pivot_wider(names_from=statistic, values_from=value)

# On the observed pelts scale
# TODO Why doesn't this work? need it for 16H3
res |>
    filter(outcome == 'mu') |>
    ggplot(aes(x=Year, colour=animal, fill=animal)) +
        geom_ribbon(aes(ymin=lower, ymax=upper, y=mean), alpha=0.3) +
        geom_point(aes(y=actual)) +
        geom_line(aes(y=mean)) +
        theme_classic() +
        scale_colour_brewer("", palette="Dark2") +
        scale_fill_brewer("", palette="Dark2") +
        labs(x="Year", y="Thousands of pelts") +
        theme(
            legend.position = "bottom"
        )

# Far smoother for populations where don't have measurement error 
res |>
    filter(outcome == 'pop') |>
    ggplot(aes(x=Year, colour=animal, fill=animal)) +
        geom_ribbon(aes(ymin=lower, ymax=upper, y=mean), alpha=0.3) +
        geom_line(aes(y=mean)) +
        theme_classic() +
        scale_colour_brewer("", palette="Dark2") +
        scale_fill_brewer("", palette="Dark2") +
        labs(x="Year", y="Thousands of animals") +
        theme(
            legend.position = "bottom"
        )


# Problems ----------------------------------------------------------------
# 16E1: Disadvantages of GLMs are they are linear by default, can tunnel-vision
# you into thinking solely in terms of GLMs, parameters don't have natural scientific
# interpretations

# 16E2:
# Famous scientific models that aren't GLMs:
#  - time-series models
#  - E=MC^2
#  - Logistic growth

# 16E3: 
# What scientific models can be transformed through logarithms into GLMs?
# log(E) = log(M) + 2 * log(C)
# i.e. intercept is 0, bM is 1, bC is 2

# 16M1: 
# Refit m16.1 with exponent as free-parameter - is the correct value of 3 identified?
d <- Howell1 |> 
        as_tibble() |>
        mutate(
            w = weight / mean(weight),
            h = height / mean(height)
        )

m16.1_2 <- ulam(
    alist(
        w ~ dlnorm(mu, sigma),
        exp(mu) <- 3.141593 * k * p^2 * h^x,
        p ~ beta(2, 18),
        k ~ exponential(0.5),
        x ~ exponential(1),
        sigma ~ exponential(1)
    ), data=d, chains=4, cores=4, log_lik=TRUE
)
# Does it identify 3?, slightly lower at 2.32
precis(m16.1_2)

# This has a shallower curve at first so it fits the younger population better
# But fares less well with the older
ndata <- tibble(h = seq(min(d$h), max(d$h), length.out=100))
mu_1 <- exp(link(m16.1, data=ndata))
mu_2 <- exp(link(m16.1_2, data=ndata))
ndata |>
    mutate(
        height_cm = h * mean(d$height),
        mu_fixed_mean = colMeans(mu_1) * mean(d$weight),
        mu_fixed_lower = apply(mu_1, 2, PI)[1, ] * mean(d$weight),
        mu_fixed_upper = apply(mu_1, 2, PI)[2, ] * mean(d$weight),
        mu_free_mean = colMeans(mu_2) * mean(d$weight),
        mu_free_lower = apply(mu_2, 2, PI)[1, ] * mean(d$weight),
        mu_free_upper = apply(mu_2, 2, PI)[2, ] * mean(d$weight)
    ) |>
    pivot_longer(starts_with("mu_"), names_pattern = "mu_(.+)_(.+)",
                 names_to=c("model", "metric")) |>
    pivot_wider(names_from=metric, values_from=value) |>
    ggplot(aes(x=height_cm, y=mean)) +
        geom_ribbon(aes(ymin=lower, ymax=upper, fill=model), alpha=0.4) +
        geom_point(data=d |> rename(height_cm=height, mean=weight)) +
        geom_line(aes(colour=model)) +
        theme_classic() + 
        labs(x="Height (cm)", y="Weight (kg)")

# How does this impact predictive accuracy?
# The model with the free parameter is much preferred!
compare(m16.1, m16.1_2)

# 16M2:
# Predictive prior simulation for this dataset
prior_predictive <- function(h, k_func, p_func, N=1e4, N_plot=20) {
    ks <- k_func(N)
    ps <- p_func(N)
     out <- map_dfr(1:N, function(i) {
         out = 3.141593 * ks[i] * ps[i]^2 * h^3
         tibble(
             height=h,
             weight=out,
             sample=i
         )
     })
     out |>
         mutate(
             height = height * mean(d$height),
             weight = weight * mean(d$weight)
             
         ) |>
         filter(sample <= N_plot) |>
         ggplot(aes(x=height, y=weight, group=sample)) +
           geom_line(alpha=0.5) +
           theme_classic() + 
            facet_wrap(~sample)
}

# With the default priors it looks ok, although likely too flat in most
# situations
# Want to at least height ~60kg by 160cm in nearly every simulation
h_seq <- seq(min(d$h)-0.2, max(d$h)+0.2, length.out=100)
prior_predictive(
    h_seq,
    function(N) rexp(N, 0.5),
    function(N) rbeta(N, 2, 18),
    N_plot=20
)

# How can we fix this?
# Function is: 3.141593 * k * p^2 * h^3
# Only have 2 dials, k and p, both of which need to be positive
# and p is in [0, 1]
# I think it's more likely that k is what needs tweaking
# Will firstly make k very flat to see if this makes much difference?
# Nope not really
prior_predictive(
    h_seq,
    function(N) rexp(N, 0.5),
    function(N) rbeta(N, 1, 1),
    N_plot=20
)

# In which case let's try increasing k
# I will do this with a half normal centered around 2 for starters
# Thta's very flat still!
prior_predictive(
    h_seq,
    function(N) abs(rnorm(N, 2, 1)),
    function(N) rbeta(N, 2, 18),
    N_plot=20
)

# The median draw from the original exp(0.5) prior is 1.4
# and this produced very flat curves as saw above, so maybe
# should be more ambitions
summary(rexp(1e4, 0.5))

# N(4, 1) on k does a lot better actually
# Still not quite enough curve - i.e. would want to see at least 
# 60kg being reached in 90% of simulations
prior_predictive(
    h_seq,
    function(N) abs(rnorm(N, 4, 1)),
    function(N) rbeta(N, 2, 18),
    N_plot=20
)

# So can try even higher at N(7, 1) and that looks about right
# for the magnitude, although there are still some flat curves
prior_predictive(
    h_seq,
    function(N) abs(rnorm(N, 7, 1)),
    function(N) rbeta(N, 2, 18),
    N_plot=20
)

# In which case can bump up the beta prior on p to 2/10
# And that has helped a lot, some ludicrous curves now!
# Will leave it at this, slightly better than what had before
prior_predictive(
    h_seq,
    function(N) abs(rnorm(N, 7, 1)),
    function(N) rbeta(N, 2, 10),
    N_plot=20
)

# 16M3:
# Priors are N(1, 0.5) for thetas 1 & 3
# N(0.05, 0.05) for thetas 2 & 4
# Can use the sim_lynx_hare function for the generating values

# Let's try a handful of runs using the default priors
# and it fails miserably! no oscillating behaviour, lynx populations
# blow up, I don't see how they would blow up because the birth rate is dependent
# upon hares, but if that is 0 then would expect to have no Lynxes
# which the system is settling towards but it's taking a while
N <- 1e4
thetas <- cbind(
    rnorm(N, 1, 0.5),
    rnorm(N, 0.05, 0.05),
    rnorm(N, 1, 0.5),
    rnorm(N, 0.05, 0.05)
)
map_dfr(1:5, function(i) {
    out = sim_lynx_hare(1e4, as.numeric(Lynx_Hare[1, 2:3]), thetas[i, ]) |>
        mutate(sim=i)
}) |>
    pivot_longer(-c(step, sim)) |>
    ggplot(aes(x=step, y=value, colour=name)) +
        geom_point() +
        geom_line() +
        facet_wrap(~sim) +
        theme_classic() +
        scale_colour_brewer("", palette = "Dark2") +
        theme(
            legend.position = "bottom"
        )

# Thetas 1&3 are birth rates, while thetas 2&4 are death rates
# of hares and lynxes respectively
# So firstly want to have more hares born than lynxes
# So will adjust priors to ensure this
# Will half lynxes birth rate
# That has helped a little to stop the Lynxes blowing up, but they are still
# much more populous than hares over the entire time-history
thetas <- cbind(
    rnorm(N, 1, 0.5),
    rnorm(N, 0.05, 0.05),
    rnorm(N, 0.5, 0.5),
    rnorm(N, 0.05, 0.05)
)
map_dfr(1:5, function(i) {
    out = sim_lynx_hare(1e4, as.numeric(Lynx_Hare[1, 2:3]), thetas[i, ]) |>
        mutate(sim=i)
}) |>
    pivot_longer(-c(step, sim)) |>
    ggplot(aes(x=step, y=value, colour=name)) +
        geom_point() +
        geom_line() +
        facet_wrap(~sim) +
        theme_classic() +
        scale_colour_brewer("", palette = "Dark2") +
        theme(
            legend.position = "bottom"
        )


# Let's whack up the Lynxes death rate
# Oh no, that's got stupid values
thetas <- cbind(
    rnorm(N, 1, 0.5),
    rnorm(N, 0.05, 0.05),
    rnorm(N, 0.5, 0.5),
    rnorm(N, 0.1, 0.05)
)
map_dfr(1:5, function(i) {
    out = sim_lynx_hare(1e4, as.numeric(Lynx_Hare[1, 2:3]), thetas[i, ]) |>
        mutate(sim=i)
}) |>
    pivot_longer(-c(step, sim)) |>
    ggplot(aes(x=step, y=value, colour=name)) +
        geom_point() +
        geom_line() +
        facet_wrap(~sim) +
        theme_classic() +
        scale_colour_brewer("", palette = "Dark2") +
        theme(
            legend.position = "bottom"
        )

# Let's whack up the hare birth rate instead then in an effort to build up
# the hare population, and that's also gone to infinity in a couple of sims
thetas <- cbind(
    rnorm(N, 5, 0.5),
    rnorm(N, 0.05, 0.05),
    rnorm(N, 0.5, 0.5),
    rnorm(N, 0.05, 0.05)
)
map_dfr(1:5, function(i) {
    out = sim_lynx_hare(1e4, as.numeric(Lynx_Hare[1, 2:3]), thetas[i, ]) |>
        mutate(sim=i)
}) |>
    pivot_longer(-c(step, sim)) |>
    ggplot(aes(x=step, y=value, colour=name)) +
        geom_point() +
        geom_line() +
        facet_wrap(~sim) +
        theme_classic() +
        scale_colour_brewer("", palette = "Dark2") +
        theme(
            legend.position = "bottom"
        )

# Instead let's go back to the default values and rather than *increasing* anything
# we will just *decrease* certain things. This should avoid blowing up.
# So we want to decrease Lynx birth rate, and decrease Hare death rate
# Will halve these
# And it still blows up in one and the others don't look that good, although lynx
# are dying off fast
thetas <- cbind(
    rnorm(N, 1, 0.5),
    rnorm(N, 0.025, 0.05),
    rnorm(N, 0.5, 0.5),
    rnorm(N, 0.05, 0.05)
)
map_dfr(1:5, function(i) {
    out = sim_lynx_hare(1e4, as.numeric(Lynx_Hare[1, 2:3]), thetas[i, ]) |>
        mutate(sim=i)
}) |>
    pivot_longer(-c(step, sim)) |>
    ggplot(aes(x=step, y=value, colour=name)) +
        geom_point() +
        geom_line() +
        facet_wrap(~sim, scales="free") +
        theme_classic() +
        scale_colour_brewer("", palette = "Dark2") +
        theme(
            legend.position = "bottom"
        )

# Let's speed up the lynx dying process
# Had to whack it up a fair bit to get the numbers of lynx to come down
# Although now 3/5 are infinite
thetas <- cbind(
    rnorm(N, 1, 0.5),
    rnorm(N, 0.025, 0.05),
    rnorm(N, 0.5, 0.5),
    rnorm(N, 0.5, 0.05)
)
map_dfr(1:5, function(i) {
    out = sim_lynx_hare(1e4, as.numeric(Lynx_Hare[1, 2:3]), thetas[i, ]) |>
        mutate(sim=i)
}) |>
    pivot_longer(-c(step, sim)) |>
    ggplot(aes(x=step, y=value, colour=name)) +
        geom_point() +
        geom_line() +
        facet_wrap(~sim, scales="free") +
        theme_classic() +
        scale_colour_brewer("", palette = "Dark2") +
        theme(
            legend.position = "bottom"
        )

# Can also try reducing the lynx birth rate
# And finally we have an oscillation! But alas only in 1 scenarios whereas
# in the others we get infinite
thetas <- cbind(
    rnorm(N, 1, 0.5),
    rnorm(N, 0.025, 0.05),
    rnorm(N, 0.25, 0.1),
    rnorm(N, 0.5, 0.05)
)
map_dfr(1:5, function(i) {
    out = sim_lynx_hare(1e4, as.numeric(Lynx_Hare[1, 2:3]), thetas[i, ]) |>
        mutate(sim=i)
}) |>
    pivot_longer(-c(step, sim)) |>
    ggplot(aes(x=step, y=value, colour=name)) +
        geom_point() +
        geom_line() +
        facet_wrap(~sim, scales="free") +
        theme_classic() +
        scale_colour_brewer("", palette = "Dark2") +
        theme(
            legend.position = "bottom"
        )

# High variance is often the cause of things blowing up so
# will try reducing the variance on the Hare birth rate 
# That's helped a bit (3/5 oscillate now)
thetas <- cbind(
    rnorm(N, 1, 0.1),
    rnorm(N, 0.025, 0.05),
    rnorm(N, 0.25, 0.1),
    rnorm(N, 0.5, 0.05)
)
map_dfr(1:5, function(i) {
    out = sim_lynx_hare(1e4, as.numeric(Lynx_Hare[1, 2:3]), thetas[i, ]) |>
        mutate(sim=i)
}) |>
    pivot_longer(-c(step, sim)) |>
    ggplot(aes(x=step, y=value, colour=name)) +
        geom_point() +
        geom_line() +
        facet_wrap(~sim, scales="free") +
        theme_classic() +
        scale_colour_brewer("", palette = "Dark2") +
        theme(
            legend.position = "bottom"
        )

# But can go further. Let's increase hare births again
thetas <- cbind(
    rnorm(N, 2, 0.1),
    rnorm(N, 0.025, 0.05),
    rnorm(N, 0.25, 0.1),
    rnorm(N, 0.5, 0.05)
)
map_dfr(1:5, function(i) {
    out = sim_lynx_hare(1e4, as.numeric(Lynx_Hare[1, 2:3]), thetas[i, ]) |>
        mutate(sim=i)
}) |>
    pivot_longer(-c(step, sim)) |>
    ggplot(aes(x=step, y=value, colour=name)) +
        geom_point() +
        geom_line() +
        facet_wrap(~sim, scales="free") +
        theme_classic() +
        scale_colour_brewer("", palette = "Dark2") +
        theme(
            legend.position = "bottom"
        )

# And getting back to blowing up!
# Let's curtail the lynx again instead
# The death rate isn't changing that much so will reduce their birth rate some more
# And there we go! (NB: the key here was also to reduce the variance. I think
# having negative values too often is what causes the blowing up)
thetas <- cbind(
    rnorm(N, 2, 0.1),
    rnorm(N, 0.025, 0.05),
    rnorm(N, 0.1, 0.05),
    rnorm(N, 0.5, 0.05)
)
map_dfr(1:5, function(i) {
    out = sim_lynx_hare(1e4, as.numeric(Lynx_Hare[1, 2:3]), thetas[i, ]) |>
        mutate(sim=i)
}) |>
    pivot_longer(-c(step, sim)) |>
    ggplot(aes(x=step, y=value, colour=name)) +
        geom_point() +
        geom_line() +
        facet_wrap(~sim, scales="free") +
        theme_classic() +
        scale_colour_brewer("", palette = "Dark2") +
        theme(
            legend.position = "bottom"
        )

# We've got oscillations now, but the Lynx populations are still very high
# With another halving and it looks ok
# Will call it a day there
thetas <- cbind(
    rnorm(N, 2, 0.1),
    rnorm(N, 0.025, 0.05),
    rnorm(N, 0.05, 0.025),
    rnorm(N, 0.5, 0.05)
)
map_dfr(1:5, function(i) {
    out = sim_lynx_hare(1e4, as.numeric(Lynx_Hare[1, 2:3]), thetas[i, ]) |>
        mutate(sim=i)
}) |>
    pivot_longer(-c(step, sim)) |>
    ggplot(aes(x=step, y=value, colour=name)) +
        geom_point() +
        geom_line() +
        facet_wrap(~sim, scales="free") +
        theme_classic() +
        scale_colour_brewer("", palette = "Dark2") +
        theme(
            legend.position = "bottom"
        )

# 16M4
# Using sphere model instead of a cylinder
# I.e. volume = 4/3*pi*r^3
# Now don't even have height encoded in the model but height is the diameter
# as we're imagining a person filling the sphere
# So r = h/2
# So V = 4/3 * pi * (h/2)^3  = 1/6 * pi * h^3
# And also need the parameter k to express translation between volume to weight
# W = k * 1/6 * pi * h^3
# But, the old model was k*pi*p^2 * h^3
# So it's just now going to fit k/6 = k*p^2
# So the functional form hasn't changed but the parameterisation has
# I expect it will have very similar predictions to the cylinder model
m16m4 <- ulam(
    alist(
        w ~ dlnorm(mu, sigma),
        exp(mu) <- k * 1/6 * 3.141593 * h^3,
        k ~ exponential(0.5),
        sigma ~ exponential(1)
    ), data=d, chains=4, cores=4, log_lik = TRUE
)
# Now k has reduced massively to 1.81 from 5.53
precis(m16m4)
precis(m16.1)

# But remember that k and p are non-identifiable 
# So just using the mean parameters, does k/6 = k*p^2?
# Near enough
1.81 / 6
5.53 * 0.25^2

# In terms of WAIC there is minimal difference as expected
compare(m16m4, m16.1)

# Do the posterior predictions look the same?
# Identical!
ndata <- tibble(h = seq(min(d$h), max(d$h), length.out=100))
mu_1 <- exp(link(m16.1, data=ndata))
mu_2 <- exp(link(m16m4, data=ndata))
ndata |>
    mutate(
        height_cm = h * mean(d$height),
        mu_cylinder_mean = colMeans(mu_1) * mean(d$weight),
        mu_cylinder_lower = apply(mu_1, 2, PI)[1, ] * mean(d$weight),
        mu_cylinder_upper = apply(mu_1, 2, PI)[2, ] * mean(d$weight),
        mu_sphere_mean = colMeans(mu_2) * mean(d$weight),
        mu_sphere_lower = apply(mu_2, 2, PI)[1, ] * mean(d$weight),
        mu_sphere_upper = apply(mu_2, 2, PI)[2, ] * mean(d$weight)
    ) |>
    pivot_longer(starts_with("mu_"), names_pattern = "mu_(.+)_(.+)",
                 names_to=c("model", "metric")) |>
    pivot_wider(names_from=metric, values_from=value) |>
    ggplot(aes(x=height_cm, y=mean)) +
        geom_ribbon(aes(ymin=lower, ymax=upper, fill=model), alpha=0.4) +
        geom_point(data=d |> rename(height_cm=height, mean=weight)) +
        geom_line(aes(colour=model)) +
        theme_classic() + 
        labs(x="Height (cm)", y="Weight (kg)")

# 16H1
# Phi is the parameter that simplifies alpha * beta^theta * Mmax
# Mmax we are told can vary by sex with males on average being higher
# To account for this we can add an additional multiplicative parameter that is the
# additional max mass of being a male
# This is easier than using the non-simplified model and adding an additive offset to MMax.
dat_list <- list(
    n = as.integer(Panda_nuts$nuts_opened),
    age=Panda_nuts$age/max(Panda_nuts$age),
    seconds=Panda_nuts$seconds,
    sex=as.integer(Panda_nuts$sex)-1  # Male is 1, Female 0
)

# Rather than modelling the additional male effect directly, we'll subtract one
# i.e. so that for an additional 50% max mass, the coefficient will be 0.5, rather than
# 1.5. This is achieved by adding 1 to the male multiplier
# Rather than assigning a constrained parameter, we'll use a normal distribution that allows
# for negative effects but with the bulk of the mass positive
m16h1 <- ulam(
    alist(
        n ~ poisson(lambda),
        lambda <- seconds*(1+male*sex)*phi*(1-exp(-k*age))^theta,
        male ~ normal(0.2, 0.1),
        phi ~ lognormal(log(1), 0.1),
        k ~ lognormal(log(2), 0.25),
        theta ~ lognormal(log(5), 0.25)
    ), data=dat_list, chains=4, cores=4
)

# Male has a 40% difference (technically in either max mass, or mass to strength ratio,
# or strength to nut opening ability)
precis(m16h1)

# 16H2
# Allow either phi or k to vary by indiviudal
# Phi does multiple things: mass -> strength, strength -> nut opening, max mass
# k is rate of increasing mass with age
# Makes more sense to allow phi to vary by individual since it does a lot more than k
# And so would expect to have greater variability

# Just check have continuous chimpanzee ids
# And we do!
Panda_nuts |>
    distinct(chimpanzee) |>
    arrange(chimpanzee) |>
    print(n=Inf)
dat_list$id <- Panda_nuts$chimpanzee
    
m16h2 <- ulam(
    alist(
        n ~ poisson(lambda),
        lambda <- seconds*(1+male*sex)*phi[id]*(1-exp(-k*age))^theta,
        male ~ normal(0.2, 0.1),
        phi[id] ~ lognormal(phi_mu, phi_sigma),
        phi_mu ~ lognormal(log(1), 0.1),
        phi_sigma ~ exponential(1),
        k ~ lognormal(log(2), 0.25),
        theta ~ lognormal(log(5), 0.25)
    ), data=dat_list, chains=4, cores=4, log_lik=TRUE
)
# A lot of variance between individuals! (phi sigma)
# Pretty bad Rhat, should really non-center but unsure how to do that for log-normal
precis(m16h2, depth=2)

# Do predictions look reasonable?
ndata <- expand_grid(
    seconds=1,
    age=seq(0, 1, length.out=100),
    distinct(Panda_nuts, sex, id=chimpanzee)
) |>
    mutate(sex = as.integer(sex)-1)
mu <- link(m16h2, data=ndata)
# Yep lots of very different growth curves!
ndata |>
    mutate(
        mu_mean = colMeans(mu),
        mu_lower = apply(mu, 2, PI)[1, ],
        mu_upper = apply(mu, 2, PI)[2, ]
    ) |>
    ggplot(aes(x=age, group=id)) +
        geom_ribbon(aes(y=mu_mean, ymin=mu_lower, ymax=mu_upper), alpha=0.1) +
        geom_line(aes(y=mu_mean)) +
        theme_classic()

# And while we're at it, the sex difference from the previous model
ndata <- expand_grid(
    seconds=1,
    age=seq(0, 1, length.out=100),
    sex = c(0, 1)
)
mu <- link(m16h1, data=ndata)
# Yep males have higher end mass
ndata |>
    mutate(
        mu_mean = colMeans(mu),
        mu_lower = apply(mu, 2, PI)[1, ],
        mu_upper = apply(mu, 2, PI)[2, ],
        sex=as.factor(sex)
    ) |>
    ggplot(aes(x=age, colour=sex, fill=sex)) +
        geom_ribbon(aes(y=mu_mean, ymin=mu_lower, ymax=mu_upper), alpha=0.1) +
        geom_line(aes(y=mu_mean)) +
        theme_classic()

# What Richard did to solve the non-centered issue, and also to ensure that the lambda is constrained
# to be positive (I used a log-normal), was to use an Exponential varying effect
# I.e. each effect is *multiplicative*
# He also uses a non-centered exponential in particular
# I think what this is doing is having the original phi
# on its original units, and we're just having individuals
# with different multiplicative effects, rather than
# completely different phis per person
# Then each individuals effect is exp(tau)
# But this is non-centered so draw individual from exp(1)
# and multiply this by group scale tau
m16h2_exp <- ulam(
    alist(
        n ~ poisson(lambda),
        lambda <- seconds*(1+male*sex)*(phi * z[id] * tau)*(1-exp(-k*age))^theta,
        male ~ normal(0.2, 0.1),
        phi ~ lognormal(log(1), 0.1),
        z[id] ~ exponential(1),
        tau ~ exponential(1),
        k ~ lognormal(log(2), 0.25),
        theta ~ lognormal(log(5), 0.25),
        gq> vector[id]:zz <<- z*tau
    ), data=dat_list, chains=4, cores=4, log_lik = TRUE
)
# Decent Rhat!
precis(m16h2_exp, depth=2)
# tau is relatively high, indicating a decent amount of
# variance

# How do predictions look compared to my log-normal model?
# The exponential model looks far more sensible!
ndata <- expand_grid(
    seconds=1,
    age=seq(0, 1, length.out=100),
    distinct(Panda_nuts, sex, id=chimpanzee)
) |>
    mutate(sex = as.integer(sex)-1)
mu_ln <- link(m16h2, data=ndata)
mu_exp <- link(m16h2_exp, data=ndata)
ndata |>
    mutate(
        mu_ln_mean = colMeans(mu_ln),
        mu_ln_lower = apply(mu_ln, 2, PI)[1, ],
        mu_ln_upper = apply(mu_ln, 2, PI)[2, ],
        mu_exp_mean = colMeans(mu_exp),
        mu_exp_lower = apply(mu_exp, 2, PI)[1, ],
        mu_exp_upper = apply(mu_exp, 2, PI)[2, ]
    ) |>
    pivot_longer(starts_with("mu"),
                 names_pattern="mu_(.+)_(.+)",
                 names_to=c("model", "stat")) |>
    pivot_wider(names_from=stat) |>
    ggplot(aes(x=age, group=id)) +
        geom_ribbon(aes(y=mean, ymin=lower, ymax=upper), alpha=0.1) +
        geom_line(aes(y=mean)) +
        facet_wrap(~model) +
        theme_classic()

# The exponential model is prefered by WAIC
compare(m16h2, m16h2_exp)

# 16H3
# Fitting lagged variables and comparing to ODE model
dat_ar1 <- list(
    L = Lynx_Hare$Lynx[2:21],
    L_lag1 = Lynx_Hare$Lynx[1:20],
    H = Lynx_Hare$Hare[2:21],
    H_lag1 = Lynx_Hare$Hare[1:20]
)

m16h3 <- ulam(
    alist(
        L ~ lognormal(log(mu_l), sigma_l),
        H ~ lognormal(log(mu_h), sigma_h),
        mu_l <- alpha_l + beta_ll * L_lag1 + beta_lh*H_lag1,
        mu_h <- alpha_h + beta_hh * H_lag1 + beta_hl*L_lag1,
        c(alpha_l, alpha_h) ~ normal(0, 1),
        c(beta_ll, beta_hh) ~ normal(0.5, 1),
        beta_lh ~ normal(0.5, 1),
        beta_hl ~ normal(-0.5, 1), # NB: lynxes leads to decreasing hares
        sigma_l ~ exponential(1),
        sigma_h ~ exponential(1)
    ),
    data=dat_ar1, chains=4, cores=4, log_lik=TRUE
)
# Diagnostics are ok
precis(m16h3)

# Let's look at predictions - can't compare to the ODE model
# since the predictions for pelts isn't working (see m16.5)
# Looks ok but very weak effect between them
# I.e. lynx population isn't dropping as quickly as would expect
mu <- link(m16h3)
Lynx_Hare |>
    filter(Year > 1900) |>
    mutate(
        mu_lynx_mean = colMeans(mu$mu_l),
        mu_lynx_lower = apply(mu$mu_l, 2, PI)[1, ],
        mu_lynx_upper = apply(mu$mu_l, 2, PI)[2, ],
        mu_hare_mean = colMeans(mu$mu_h),
        mu_hare_lower = apply(mu$mu_h, 2, PI)[1, ],
        mu_hare_upper = apply(mu$mu_h, 2, PI)[2, ]
    ) |>
    pivot_longer(
        starts_with("mu"),
        names_pattern = "mu_(.+)_(.+)",
        names_to=c("animal", "stat") 
    ) |>
    pivot_wider(names_from=stat) |>
    ggplot(aes(x=Year, y=mean, colour=animal, fill=animal)) +
        geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
        geom_line() +
        theme_classic() +
        scale_fill_brewer("", palette="Dark2") +
        scale_colour_brewer("", palette="Dark2")

# Richard states this is because the ODE allows for interactions
# I.e. lynx and hare populations effect each other *multiplicatively*
# whereas here it is just by addition
m16h3_interaction <- ulam(
    alist(
        L ~ lognormal(log(mu_l), sigma_l),
        H ~ lognormal(log(mu_h), sigma_h),
        mu_l <- alpha_l + beta_ll * L_lag1 + beta_lh*H_lag1*L_lag1,
        mu_h <- alpha_h + beta_hh * H_lag1 + beta_hl*L_lag1*H_lag1,
        c(alpha_l, alpha_h) ~ normal(0, 1),
        c(beta_ll, beta_hh) ~ normal(0.5, 1),
        beta_lh ~ normal(0.5, 1),
        beta_hl ~ normal(-0.5, 1), # NB: lynxes leads to decreasing hares
        sigma_l ~ exponential(1),
        sigma_h ~ exponential(1)
    ),
    data=dat_ar1, chains=4, cores=4, log_lik=TRUE
)
# Diagnostics are ok
precis(m16h3_interaction)

mu_int <- link(m16h3_interaction)
# Much better fit now
Lynx_Hare |>
    filter(Year > 1900) |>
    mutate(
        mu_lynx_mean = colMeans(mu_int$mu_l),
        mu_lynx_lower = apply(mu_int$mu_l, 2, PI)[1, ],
        mu_lynx_upper = apply(mu_int$mu_l, 2, PI)[2, ],
        mu_hare_mean = colMeans(mu_int$mu_h),
        mu_hare_lower = apply(mu_int$mu_h, 2, PI)[1, ],
        mu_hare_upper = apply(mu_int$mu_h, 2, PI)[2, ]
    ) |>
    pivot_longer(
        starts_with("mu"),
        names_pattern = "mu_(.+)_(.+)",
        names_to=c("animal", "stat") 
    ) |>
    pivot_wider(names_from=stat) |>
    ggplot(aes(x=Year, y=mean, colour=animal, fill=animal)) +
        geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
        geom_line() +
        theme_classic() +
        scale_fill_brewer("", palette="Dark2") +
        scale_colour_brewer("", palette="Dark2")

# Can we compare models?
# The non-interaction model is much prefered actually
# By a crazy difference too!
# Why is that the case?
# How as adding the interactions made such a huge difference
# to the number of effective parameters?
compare(m16h3, m16h3_interaction)

# Here's the actual data and this looks to me like the
# interaction model fits better
# Shame can't compare the ODE model as it doesn't
# have likelihoods
# It also struggles to fit now, although it was definitely
# fine at first. I'm not sure why this is the case...
Lynx_Hare |>
    pivot_longer(-Year) |>
    ggplot(aes(x=Year, y=value, colour=name)) +
        geom_line() +
        theme_classic() +
        scale_colour_brewer("", palette="Dark2")

# 16H4
# 2-step lag
# Will use the interaction model
#dat_dar2 <- Lynx_Hare |>
foo <- Lynx_Hare |>
    rename(L=Lynx, H=Hare) |>
    mutate(
        L_lag1 = lag(L),
        L_lag2 = lag(L_lag1),
        H_lag1 = lag(H),
        H_lag2 = lag(H_lag1)
    ) |>
    filter(Year >= 1902)
m16h4 <- ulam(
    alist(
        L ~ lognormal(log(mu_l), sigma_l),
        H ~ lognormal(log(mu_h), sigma_h),
        mu_l <- alpha_l + beta_ll * L_lag1 + beta_ll2*L_lag2 + beta_lh*H_lag1*L_lag1 + beta_lh2*H_lag2*L_lag2,
        mu_h <- alpha_h + beta_hh * H_lag1 + beta_hh2*L_lag2 + beta_hl*L_lag1*H_lag1 + beta_hl2*H_lag2*L_lag2,
        c(alpha_l, alpha_h) ~ normal(0, 1),
        c(beta_ll, beta_hh, beta_ll2, beta_hh2) ~ normal(0.5, 1),
        c(beta_lh, beta_lh2) ~ normal(0.5, 1),
        c(beta_hl, beta_hl2) ~ normal(-0.5, 1), # NB: lynxes leads to decreasing hares
        sigma_l ~ exponential(1),
        sigma_h ~ exponential(1)
    ),
    data=foo, chains=4, cores=4, log_lik=TRUE
)
# Diagnostics are ok
precis(m16h4)

# Any effect on predictions?
mu_2 <- link(m16h4)
# Lynx overshoot!
# Doesn't look much different to the 1-lag model
Lynx_Hare |>
    filter(Year > 1901) |>
    mutate(
        mu_lynx_mean = colMeans(mu_2$mu_l),
        mu_lynx_lower = apply(mu_2$mu_l, 2, PI)[1, ],
        mu_lynx_upper = apply(mu_2$mu_l, 2, PI)[2, ],
        mu_hare_mean = colMeans(mu_2$mu_h),
        mu_hare_lower = apply(mu_2$mu_h, 2, PI)[1, ],
        mu_hare_upper = apply(mu_2$mu_h, 2, PI)[2, ]
    ) |>
    pivot_longer(
        starts_with("mu"),
        names_pattern = "mu_(.+)_(.+)",
        names_to=c("animal", "stat") 
    ) |>
    pivot_wider(names_from=stat) |>
    ggplot(aes(x=Year, y=mean, colour=animal, fill=animal)) +
        geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
        geom_line() +
        theme_classic() +
        scale_fill_brewer("", palette="Dark2") +
        scale_colour_brewer("", palette="Dark2")

# What about WAIC?
# Given the far more parameters now I'd expect the 2-lag model
# to do even worse
# Hah yep it does even worse with 1 fewer observation as well!
# To be fair there are only 21 observations and 12 parameters now
# And not much difference from 1-lag model
compare(m16h3, m16h3_interaction, m16h4)

# 16H5
data(Mites)
Mites <- as_tibble(Mites)
# 3 cols:
#   - day
#   - # prey
#   - # predator
# No measurement error unlike pelts
# Use the same Lotka-Volterra ODE from the Lynx-Hare model
# Adapt the stancode, taking care to rescale priors
# This means will also need to rescale parameters
Mites

# Let's look at the data first
# And we see a similar trend as Lynx Hares but less pronounced
# I.e. on the 3rd peak of the prey there isn't a corresponding predator peak
Mites |>
    pivot_longer(-day) |>
    ggplot(aes(x=day, y=value, colour=name)) +
        geom_line() +
        theme_classic() +
        scale_colour_brewer("", palette="Dark2") +
        theme(legend.position = "bottom")

# Firstly want to look at the code and modify it:
cat(Lynx_Hare_model)

# The measurement errors stuff want to remove is just the p multiplier
# corresponding to trap success
# still want to model a latent population with our observations of it

# Let's just try with default priors first
# I.e. the same as Lynx Hare
m16h5_code <- "
functions {
  real[] dpop_dt( real t,                 // time
                real[] pop_init,          // initial state {lynx, hares}
                real[] theta,             // parameters
                real[] x_r, int[] x_i) {  // unused
    real L = pop_init[1];
    real H = pop_init[2];
    real bh = theta[1];
    real mh = theta[2];
    real ml = theta[3];
    real bl = theta[4];
    // differential equations
    real dH_dt = (bh - mh * L) * H;
    real dL_dt = (bl * H - ml) * L;
    return { dL_dt , dH_dt };
  }
}
data {
  int<lower=0> N;              // number of measurement times
  real<lower=0> mites[N,2];    // measured populations
  real<lower=0> times_measured[N-1];
}
parameters {
  real<lower=0> theta[4];      // { bh, mh, ml, bl }
  real<lower=0> pop_init[2];   // initial population state
  real<lower=0> sigma[2];      // measurement errors
}
transformed parameters {
  real pop[N, 2];
  pop[1,1] = pop_init[1];
  pop[1,2] = pop_init[2];
  pop[2:N,1:2] = integrate_ode_rk45(
    dpop_dt, pop_init, 0, times_measured, theta,
    rep_array(0.0, 0), rep_array(0, 0),
    1e-5, 1e-3, 5e2);
}
model {
  // priors
  theta[{1,3}] ~ normal( 1 , 0.5 );    // bh,ml
  theta[{2,4}] ~ normal( 0.05, 0.05 ); // mh,bl
  sigma ~ exponential( 1 );
  pop_init ~ lognormal( log(10) , 1 );
  // observation model
  // connect latent population state to observed pelts
  for ( t in 1:N )
    for ( k in 1:2 )
      mites[t,k] ~ lognormal( log(pop[t,k]) , sigma[k] );
}
generated quantities {
  real mites_pred[N,2];
  for ( t in 1:N )
    for ( k in 1:2 )
      mites_pred[t,k] = lognormal_rng( log(pop[t,k]) , sigma[k] );
}
"

dat_list <- list(
    N=nrow(Mites),
    mites=Mites |> select(predator, prey) |> as.matrix(),
    times_measured=Mites$day[2:nrow(Mites)]
)

m16h5 <- stan(model_code=m16h5_code, 
              data=dat_list,
              chains=3,
              cores=3,
              control=list(adapt_delta=0.95))
# Very poor mixing - probably from poor priors
# Theta 2 & 4 having minimal effect (birth rate of prey and death of predator
# respectively)
# I.e. this is suggesting to me that the birth rate of the predator seems to lack
# both auto-correlation and correlation with the prey death rate
precis(m16h5, depth=2)

# Very different scales from the observed mites (left) and the latent population
# (right) - due to the observed mites being a draw from log-normal of latent pop
post <- extract.samples(m16h5)
res <- Mites |>
    mutate(
        predator_mu_mean = colMeans(post$mites_pred[, , 1]),
        predator_mu_lower = apply(post$mites_pred[, , 1], 2, PI)[1, ],
        predator_mu_upper = apply(post$mites_pred[, , 1], 2, PI)[2, ],
        prey_mu_mean = colMeans(post$mites_pred[, , 2]),
        prey_mu_lower = apply(post$mites_pred[, , 2], 2, PI)[1, ],
        prey_mu_upper = apply(post$mites_pred[, , 2], 2, PI)[2, ],
        predator_pop_mean = colMeans(post$pop[, , 1]),
        predator_pop_lower = apply(post$pop[, , 1], 2, PI)[1, ],
        predator_pop_upper = apply(post$pop[, , 1], 2, PI)[2, ],
        prey_pop_mean = colMeans(post$pop[, , 2]),
        prey_pop_lower = apply(post$pop[, , 2], 2, PI)[1, ],
        prey_pop_upper = apply(post$pop[, , 2], 2, PI)[2, ]
    ) |>
    rename(predator_mu_actual=predator, prey_mu_actual=prey) |>
    pivot_longer(-c(day), names_pattern="(.+)_(.+)_(.+)", names_to=c("animal", "outcome", "statistic"))  |>
    pivot_wider(names_from=statistic, values_from=value)

# So now the question is how to improve this model?
# Richard shows some prior predictive checks that I didn't do myself (basically
# what priors on theta make sense)
# Also make the initial population have a prior on the observed data
# So actually this is including measurement errors!
# Unsure how to not do this, as when I tried it failed - 
# it seemed to be expecting a parameter rather than data for the
# output of integrate_ode_rk45

# I also don't see how thetas match up, as the priors describe
# theta[3] as bl, but the function dpop_dt shows that
# theta[3] is actually ml
m16h5_code_2 <- "
functions {
  real[] dpop_dt( real t,                 // time
                real[] pop_init,          // initial state {lynx, hares}
                real[] theta,             // parameters
                real[] x_r, int[] x_i) {  // unused
    real L = pop_init[1];
    real H = pop_init[2];
    real bh = theta[1];
    real mh = theta[2];
    real ml = theta[3];
    real bl = theta[4];
    // differential equations
    real dH_dt = (bh - mh * L) * H;
    real dL_dt = (bl * H - ml) * L;
    return { dL_dt , dH_dt };
  }
}
data {
  int<lower=0> N;              // number of measurement times
  real<lower=0> mites[N,2];    // measured populations
  real<lower=0> times_measured[N-1];
}
parameters {
  real<lower=0> theta[4];      // { bh, mh, ml, bl }
  real<lower=0> pop_init[2];   // initial population state
  real<lower=0> sigma[2];      // measurement errors
}
transformed parameters {
  real pop[N, 2];
  pop[1,1] = pop_init[1];
  pop[1,2] = pop_init[2];
  pop[2:N,1:2] = integrate_ode_rk45(
    dpop_dt, pop_init, 0, times_measured, theta,
    rep_array(0.0, 0), rep_array(0, 0),
    1e-5, 1e-3, 5e2);
}
model {
  // priors
  theta[1] ~ normal(3*0.5, 1);         // bh
  theta[2] ~ normal(0.01*0.5, 0.1);    // mh
  theta[3] ~ normal(0.001*0.5, 0.1);   // bl
  theta[4] ~ normal(1*0.5, 1);         // ml
  sigma ~ exponential( 1 );
  pop_init[1] ~ normal( mites[1, 1] , 50 );
  pop_init[2] ~ normal( mites[1, 2] , 50 );
  // observation model
  // connect latent population state to observed pelts
  for ( t in 1:N )
    for ( k in 1:2 )
      mites[t,k] ~ lognormal( log(pop[t,k]) , sigma[k] );
}
generated quantities {
  real mites_pred[N,2];
  for ( t in 1:N )
    for ( k in 1:2 )
      mites_pred[t,k] = lognormal_rng( log(pop[t,k]) , sigma[k] );
}
"

m16h5_2 <- stan(model_code=m16h5_code_2, 
              data=dat_list,
              chains=3,
              cores=3,
              control=list(adapt_delta=0.95))
# Adding the more constrained priors, and in particular the prior on initial population,
# have made a massive effect on Rhat, reducing it quite significantly for 2 of the theta values.
precis(m16h5, depth=2)
precis(m16h5_2, depth=2)

# Very different scales from the observed mites (left) and the latent population
# (right) - due to the observed mites being a draw from log-normal of latent pop
post_2 <- extract.samples(m16h5_2)
res2 <- Mites |>
    mutate(
        predator_mu_mean = colMeans(post_2$mites_pred[, , 1]),
        predator_mu_lower = apply(post_2$mites_pred[, , 1], 2, PI)[1, ],
        predator_mu_upper = apply(post_2$mites_pred[, , 1], 2, PI)[2, ],
        prey_mu_mean = colMeans(post_2$mites_pred[, , 2]),
        prey_mu_lower = apply(post_2$mites_pred[, , 2], 2, PI)[1, ],
        prey_mu_upper = apply(post_2$mites_pred[, , 2], 2, PI)[2, ],
        predator_pop_mean = colMeans(post_2$pop[, , 1]),
        predator_pop_lower = apply(post_2$pop[, , 1], 2, PI)[1, ],
        predator_pop_upper = apply(post_2$pop[, , 1], 2, PI)[2, ],
        prey_pop_mean = colMeans(post_2$pop[, , 2]),
        prey_pop_lower = apply(post_2$pop[, , 2], 2, PI)[1, ],
        prey_pop_upper = apply(post_2$pop[, , 2], 2, PI)[2, ]
    ) |>
    rename(predator_mu_actual=predator, prey_mu_actual=prey) |>
    pivot_longer(-c(day), names_pattern="(.+)_(.+)_(.+)", names_to=c("animal", "outcome", "statistic"))  |>
    pivot_wider(names_from=statistic, values_from=value)

# The first model has massive outliers in the predators, presumably from too vague priors
res |>
    ggplot(aes(x=day, colour=animal, fill=animal)) +
        geom_ribbon(aes(y=mean, ymin=lower, ymax=upper), alpha=0.3) +
        geom_line(aes(y=mean), na.rm=T) +
        theme_classic() +
        facet_wrap(~outcome) +
        scale_colour_brewer("", palette="Dark2") +
        scale_fill_brewer("", palette="Dark2")

# The second model has more realistic values, but no cycles
# Much lower values in the population, i.e. the modelled measurement error is high
res2 |>
    ggplot(aes(x=day, colour=animal, fill=animal)) +
        geom_ribbon(aes(y=mean, ymin=lower, ymax=upper), alpha=0.3) +
        geom_line(aes(y=mean), na.rm=T) +
        theme_classic() +
        facet_wrap(~outcome) +
        scale_colour_brewer("", palette="Dark2") +
        scale_fill_brewer("", palette="Dark2")

# How does it compare to the actual data?
# In the actual data there are measurements of up to 2k and clear cycles
res2 |>
    filter(!is.na(actual)) |>
    ggplot(aes(x=day, y=actual, colour=animal, fill=animal)) +
        geom_line(na.rm=T) +
        theme_classic() +
        scale_colour_brewer("", palette="Dark2") +
        scale_fill_brewer("", palette="Dark2")

# Questions for Richard ---------------------------------------------------
# What exactly does 'marginal' mean? It's used when we average over a discrete parameter (marginalizes out),
# i.e. m16.2, but also used for the 'marginal posterior', when to me these are the posterior parameters?
# Or does it mean marginal as in after the scenario has been taken into account?

# Why does m16.2 seem to not use the Categorical likelihood?
