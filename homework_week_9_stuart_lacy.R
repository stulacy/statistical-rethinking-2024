library(tidyverse)
library(rethinking)
library(patchwork)

# Question 1 --------------------------------------------------------------
data("Achehunting")
d <- Achehunting |> as_tibble()
# Each row is a hunting trip by 1 of 147 men (id)
# Also have time of trip (month/day/year), age of hunter (age),
# hours of trip (hours), how much meat they returned (kg.meat)
d

# 147 men
length(unique(d$id))

# It doesn't say what 'datatype' is, but it's either a 1, 1.1, or 3
d |> count(datatype)

# Have a lot of missing durations (hours)
summary(d)

# Estimate effect of age on probability of 'success' (i.e. kg.meat > 0)
# I.e. just success ~ age
# Have ages from 11 to 75, so a linear model doesn't quite seem right here
# But can start off with it
d$success <- as.integer(d$kg.meat > 0)
m1 <- ulam(
    alist(
        success ~ bernoulli(p),
        logit(p) <- a + b * age,
        a ~ normal(0, 1),
        b ~ normal(0, 0.5)
    ), data = d |> mutate(age = as.numeric(scale(age))) |> select(age, success),
    chains=4, cores=4
)
# Model converges well!
precis(m1)

# Let's look at predictions
age_z <- as.numeric(scale(d$age))
x_seq <- seq(min(age_z)-0.15, max(age_z)+0.15, length.out=100)
mu <- link(m1, data=tibble(age=x_seq))
# This approach doesn't work for Binary outcomes
tibble(
    age=x_seq * sd(d$age) + mean(d$age),
    mean = colMeans(mu),
    lower=apply(mu, 2, PI)[1, ],
    upper=apply(mu, 2, PI)[2, ]
) |>
    ggplot(aes(x=age, y=mean)) +
        geom_ribbon(aes(ymin=lower, ymax=upper),alpha=0.4) +
        geom_line() + 
        geom_point(aes(y=success), data=d, colour="steelblue") +
        theme_classic()

# Instead can plot the outcome by bins
# _Should_ be enough data to bin by 5
plot(d$age)
# This is fine on the whole, only 9 data points for age in [75, 79] but that's ok
d |>
    mutate(age_5 = floor(age/5)*5) |>
    count(age_5)

# Ok that clearly shows an inverted bathtub shape!
# So success is low at first (pre adult), at its peak during adult, then decreases in old age
tibble(
    age=x_seq * sd(d$age) + mean(d$age),
    mean = colMeans(mu),
    lower=apply(mu, 2, PI)[1, ],
    upper=apply(mu, 2, PI)[2, ]
) |>
    ggplot(aes(x=age, y=mean)) +
        geom_ribbon(aes(ymin=lower, ymax=upper),alpha=0.4) +
        geom_line() + 
        geom_point(aes(y=success), size=4, shape=1, 
                   data=d |>
                       mutate(age_5 = floor(age/5)*5) |>
                       group_by(age=age_5) |>
                       summarize(success = mean(success)),
                   colour="steelblue") +
        theme_classic()

# I tried several different growth curve models.
# I started with a Weibull but it was too easy to get to 
# infeasible predictions. I then tried a Generalised Pareto
# and while it stuck to the 0 lower bound, it gave outputs > 1
generalised_pareto <- function(age, a, b, c) {
    # https://www.degruyter.com/document/doi/10.1515/phys-2022-0047/html
    # Isn't constrained to [0,1] - has no upper bound
    # Works well for the pooled model, but doesn't handle
    # when adding varying effects, I think because of
    # the difficulty in ensuring that the output is within [0,1]
    1/c * exp(-b/age) * (1 + b/age) * (1 + a/c * age * exp(-b/age))^(-1/a - 1)
}

# I tried a couple of growth-specific models next, where the benefit
# is that they are parameterised to provide horizontal asymptotes
# and meaningful parameters
# But they don't have enough flexibility to provide the inverted
# bathtub shape, only S-curve
janoschek <- function(age, delta, k) {
    # http://www.pisces-conservation.com/growthhelp/index.html?janoschek.htm
    # Doesn't have enough flexibility, just plots S curve
    1 - exp(-k*age**delta)
}

# This 5 parameter Richards (counting the 2 asymptotes) only produced
# a horizontal line!
richards <- function(age, tm, k, T2) {
    # http://www.pisces-conservation.com/growthhelp/richards_curve.htm
    # Only gives straight line!
   1 / (1 + T2 * exp(-k * (age-tm)))^(1/T2);
}

# I came across a blog post from Andrew Gelman stating that his preferred
# solution is to make the output the sum of 2 separate growth functions
# with the first handling the upwards trajectory and the second the downwards
# So I'll use a mirrored janoschek for the downwards portion
# HOWEVER: I'll multiply the 2 curves together to ensure that the output
# is limited to [0, 1]
x_seq <- seq(0.0000001, 1.1, length.out=100)
# Testing parameters to see if this is possible and sure enough we can get
# something that looks reasonable - the question is will it converge?
# Seem to need all parameters to be positive
janoschek_falling <- function(age, delta, k) {
    # http://www.pisces-conservation.com/growthhelp/index.html?janoschek.htm
    # Doesn't have enough flexibility, just plots S curve
    exp(-k*age**delta)
}

# This looks brilliant, but when I came to fit it it really struggled
# with having 4 parameters to fit, giving crazy high Rhat and low n_eff
tibble(
    x=x_seq,
    rising=janoschek(x, 5, 9),
    falling=janoschek_falling(x, 3, 2),
    combined = rising * falling
) |>
    pivot_longer(-x) |>
    mutate(name = factor(name, levels=c("rising", "falling", "combined"))) |>
    ggplot(aes(x=x, y=value, colour=name)) +
        geom_line() +
        theme_classic() +
        theme(
            legend.position = "bottom"
        )

# Instead I'll remove the exponent from the falling function so there's just
# an exponential decay
janoschek_falling_1param <- function(age, k) {
    # Simplified version of the Janoschek without the inflexion point
    # i.e. simple linear decay
    exp(-k*age)
}

# This doesn't look as good but let's see what the model finds
tibble(
    x=x_seq,
    rising=janoschek(x, 3, 5),
    falling=janoschek_falling_1param(x, 1),
    combined = rising * falling
) |>
    pivot_longer(-x) |>
    mutate(name = factor(name, levels=c("rising", "falling", "combined"))) |>
    ggplot(aes(x=x, y=value, colour=name)) +
        geom_line() +
        theme_classic() +
        theme(
            legend.position = "bottom"
        )

# I'll do some prior predictive checks as well to really ensure that we end up
# in the right neighbourhood
N <- 1e4
d1 <- rnorm(N, 2, 0.5)
k1 <- rnorm(N, 2, 0.5)
k2 <- rnorm(N, 2, 0.5)
res <- sapply(x_seq, function(x) {
    janoschek(x, d1, k1) * janoschek_falling_1param(x, k2)
})
# Plot 50 samples: this looks ok actually!
# Even with some slightly negative values we aren't getting anything too
# unreasonable
# Will try and push forward without a constraint, just with this
# prior to keep the curves in the right direction
colnames(res) <- paste0("x", x_seq)
res |>
    as_tibble() |>
    mutate(sample = 1:N) |>
    pivot_longer(-sample, names_pattern="x(.+)", names_to="x", 
                 values_to="rate") |>
    mutate(age = as.numeric(x) * max(d$age)) |>
    filter(sample < 50) |>
    ggplot(aes(x=age, y=rate, group=sample)) +
        geom_line(alpha=0.3) +
        theme_classic()

m2_code <- "
data{
    int N;
    array[N] int success;
    vector[N] age;
}
parameters{
     real k_1;
     real k_2;
     real delta;
}
model{
    vector[N] p;
    k_1 ~ normal( 2, 0.5 );
    delta ~ normal( 2 , 0.5 );
    k_2 ~ normal( 2, 0.5 );
    for ( i in 1:N ) {
        p[i] = (1 - exp(-k_1*age[i]^delta)) * exp(-k_2*age[i]);
    }
    success ~ bernoulli( p );
}
"
d_list <- d |> 
            mutate(
                age = age / max(age)
            ) |> 
            select(age, success) |> 
            as.list()
d_list$N <- nrow(d)
m2 <- stan(model_code=m2_code, chains=4, cores=4, data=d_list)
# It has fitted without issues, even if n_eff is relatively low
# That's much better than the 4-parameter model
precis(m2)

# Now to plot
age_norm <- d$age / max(d$age)
x_seq <- seq(0.0000001, 1.1, length.out=100)
post <- extract.samples(m2)
mu <- t(sapply(1:2000, function(i) {
    janoschek(x_seq, post$delta[i], post$k_1[i]) * janoschek_falling_1param(x_seq, post$k_2[i])
}))

# That looks pretty decent!
# Got a better fit when had the second exponent but ah well
tibble(
    age=x_seq * max(d$age),
    mean = colMeans(mu),
    lower=apply(mu, 2, PI)[1, ],
    upper=apply(mu, 2, PI)[2, ]
) |>
    ggplot(aes(x=age, y=mean)) +
        geom_ribbon(aes(ymin=lower, ymax=upper),alpha=0.4) +
        geom_line() + 
        geom_point(aes(y=success), size=4, shape=1, 
                   data=d |>
                       mutate(age_5 = floor(age/5)*5) |>
                       group_by(age=age_5) |>
                       summarize(success = mean(success)),
                   colour="steelblue") +
        theme_classic()


# Question 2 --------------------------------------------------------------
# Allow each individual to have a different effect of age
# I'm taking the tip from the solutions here to think about WHICH parameters
# need to vary by age. Obviously it would be nice for everything to do so, but
# that makes sampling far more challenging.
# Instead can think that the age at which a hunter peaks should be universal (i.e.
# delta), but the growth rates should be unique to each hunter.
# Richard also provides the hint of sampling varying effects on the log-scale before
# exponentiating to ensure positivity, which is what has been hindering me
# for a bit
m3_code <- "
data{
    int N;
    int N_hunters;
    array[N] int success;
    array[N] int hunter;
    vector[N] age;
}
parameters{
    // Individual level
    vector[N_hunters] k_2_z;
    vector[N_hunters] k_1_z;
    real delta;
    
    // Hyper-parameters
    real k1_mu;
    real k2_mu;
    real<lower=0> k1_sigma;
    real<lower=0> k2_sigma;
}
transformed parameters {
    vector[N_hunters] k_2;
    vector[N_hunters] k_1;
    k_2 = k_2_z * k2_sigma + k2_mu;
    k_1 = k_1_z * k1_sigma + k1_mu;
}
model{
    vector[N] p;
    delta ~ normal( 2 , 0.5);
    k_1_z ~ normal(0, 1);
    k_2_z ~ normal(0, 1);
    
    // Hyper parameters
    k1_mu ~ normal(0, 0.5);
    k2_mu ~ normal(0, 0.5);
    k1_sigma ~ exponential(1);
    k2_sigma ~ exponential(1);
    
    for ( i in 1:N ) {
        p[i] = (1 - exp(-exp(k_1[hunter[i]])*age[i]^delta)) * exp(-exp(k_2[hunter[i]])*age[i]);
    }
    success ~ bernoulli( p );
} 
"
d_list <- d |> 
            mutate(
                age = age / max(age),
                id = as.integer(as.factor(id))
            ) |> 
            select(age, success, hunter=id) |> 
            as.list()
d_list$N <- nrow(d)
d_list$N_hunters <- length(unique(d$id))
# Nope won't work because it is getting negative probabilities
# I could use a min(0, p) function but that makes sampling problematic
# as it's a discrete step with no smooth gradient
# Ideally would parameterise the model in such a way that wasn't possible
# But I'm not sure how to do that?
m3 <- stan(model_code=m3_code, chains=4, cores=4, data=d_list)
# That's worked!
precis(m3, depth=2, pars=c("delta", "k1_mu", "k2_mu", "k1_sigma", "k2_sigma"))

# So the question is does age or hunter have a bigger effect?
# Let's look at the effect of hunter
post <- extract.samples(m3)
n_hunter_plot <- 16
n_draws <- 20
hunter_ids <- sample(1:d_list$N_hunters, size=n_hunter_plot)
res <- map_dfr(hunter_ids, function(hid) {
    map_dfr(1:n_draws, function(i) {
        tibble(
            sample=i,
            hunter=hid,
            age=x_seq * max(d$age),
            pred = janoschek(x_seq, post$delta[i], exp(post$k_1[i, hid])) * janoschek_falling_1param(x_seq, exp(post$k_2[i, hid]))
        )
    })
}) 

# The plot below shows predicted growth curves for different hunters
# along with the full dataset (blue) and the hunter-specific data (orange)
# It shows a fair amount of inter-hunter variability, i.e. Hunter 41
# is clearly a very successful hunter, having 100% success in one age bin
res |>
    ggplot(aes(x=age, y=pred)) +
        geom_line(aes(group=sample)) + 
        facet_wrap(~hunter) +
        geom_point(aes(y=success), size=3, shape=21, 
                   data=d |>
                       mutate(age_5 = floor(age/5)*5,
                              hunter = as.integer(as.factor(id))) |>
                       filter(hunter %in% hunter_ids) |>
                       group_by(hunter, age=age_5) |>
                       summarize(success = mean(success)),
                   fill="orange", colour="black") +
        geom_point(aes(y=success), size=2, shape=21, 
                   data=d |>
                       mutate(age_5 = floor(age/5)*5) |>
                       group_by(age=age_5) |>
                       summarize(success = mean(success)),
                   fill="steelblue", colour="white") +
        theme_classic() +
        labs(x="Age", y="Hunting success rate")

# But how to formally answer whether age or hunter has a bigger effect?
# We want to calculate "variation in success across age, averaged over hunters, 
# and compare it to the variation across hunters, averaged over age"

# compute variation across age averaged over hunters
# For each hunter
#   For each sample
#     Calculate posterior of success rate at all ages
#     Calculate variation over age
#   Calculate mean variation over samples
# Calculate mean over hunters

# So basically want a large dataframe of all posteriors
# Although since this runs out of memory will need to calculate it on the fly
v_HA <- sapply(1:d_list$N_hunters, function(hid) {
    # This is n_samples x n_ages
    # So now have posterior over all ages / sample / hunter
    v <- sapply( x_seq , function(x) {
        janoschek(x, post$delta, exp(post$k_1[, hid])) * janoschek_falling_1param(x, exp(post$k_2[, hid]))
    })
    # Immediate objective is the VARIATION ACROSS AGE
    vars <- apply(v, 1, var)
    # Then average this over samples
    sample_mean <- mean(vars)
})
# Then finally average over hunters
v_HA_mean <- mean(v_HA)

# compute variation across hunters averaged over ages
# For each age
#   For each sample
#     Calculate posterior of success rate for all hunters
#     Calculate variation over hunters
#   Calculate mean variation over samples
# Calculate mean over ages
v_AH <- sapply(x_seq, function(x) {
    # This is n_samples x n_hunters
    # So now have posterior over all ages / sample / hunter
    v <- sapply( 1:d_list$N_hunters , function(hid) {
        janoschek(x, post$delta, exp(post$k_1[, hid])) * janoschek_falling_1param(x, exp(post$k_2[, hid]))
    })
    # Immediate objective is the VARIATION ACROSS HUNTERS
    vars <- apply(v, 1, var)
    # Then average this over samples
    sample_mean <- mean(vars)
})
# Then finally average over age
v_AH_mean <- mean(v_AH)

# More variation in hunters than age
c("variation over hunters"=v_HA_mean, "variation over age"=v_AH_mean)

# Plotting these components
p1 <- tibble(
    age=x_seq*max(d$age),
    variance=v_AH
) |>
    ggplot(aes(x=age, y=variance)) +
        geom_line() +
        theme_classic() +
        labs(x="Age", y="Variation across hunters")

d2 <- tibble(
    hunter=1:d_list$N_hunters,
    variance=v_HA
) 
d2_ord <- d2 |> arrange(variance) |> pull(hunter)
p2 <- d2 |>
    mutate(hunter = factor(hunter, levels=d2_ord)) |>
    ggplot(aes(x=hunter, y=variance)) +
        geom_point() +
        theme_classic() +
        theme(
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
        ) +
        labs(x="Hunter", y="Variation across age")
# There is more variance between hunters at the highest ages (left)
# On average there is only moderate variance over age for hunters (right)
# although some of the hunters exhibit huge variance over age (top 4 or so points)
p1 | p2

# Plotting variance against success rate shows that it is indeed the best hunters
# who vary more, i.e. they develop over their lifetime
d2 |>
    inner_join(d |> mutate(hunter=as.integer(as.factor(id))),
               by="hunter") |>
    mutate(hunter = factor(hunter, levels=d2_ord)) |>
    group_by(hunter, variance) |>
    summarise(success = mean(success)) |>
    ungroup() |>
    ggplot(aes(x=variance, y=success)) +
        geom_point() +
        theme_classic() +
        labs(x="Variance of hunter over ages", y="Average success rate") 

# There are some individuals who get a perfect success rate
# But these only had a handful of hunts!
d |>
    mutate(hunter=as.integer(as.factor(id))) |>
    group_by(hunter) |>
    summarise(success = mean(success), n_hunts = n()) |>
    arrange(desc(success))

# The median number of hunts is 26 so can ignore these points
d |>
    count(id) |>
    summary()
    
# Question 3 --------------------------------------------------------------
# Not even going to attempt - has taken me too long to just find
# a growth model that I can fit with varying effects


