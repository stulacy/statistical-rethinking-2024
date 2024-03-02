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
        p ~ beta(2, 18),
        k ~ exponential(0.5),
        sigma ~ exponential(1)
    ), data=d, chains=4, cores=4
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
pairs(m16.1)

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
pairs(m16.1_b)

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
sim(m16.5)

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
# Fits quite well!
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

# Questions for Richard ---------------------------------------------------
# What exactly does 'marginal' mean? It's used when we average over a discrete parameter (marginalizes out),
# i.e. m16.2, but also used for the 'marginal posterior', when to me these are the posterior parameters?
# Or does it mean marginal as in after the scenario has been taken into account?

# Why does m16.2 seem to not use the Categorical likelihood?

