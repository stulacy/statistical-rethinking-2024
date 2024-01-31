library(rethinking)
library(tidyverse)

data(rugged)
d <- rugged
d$log_gdp <- log(d$rgdppc_2000)
dd <- d[complete.cases(d$rgdppc_2000), ]
# Why the choice of these particular rescaling?
# Have never seen dividing by the mean before
dd$log_gdp_std <- dd$log_gdp / mean(dd$log_gdp)
# And why just rescale 0-1?
dd$rugged_std <- dd$rugged / max(dd$rugged)

# Raw log GDP, normal ish, but on range 6-11
dens(dd$log_gdp)
# Divided by mean means centered around 1, with 0.7-1.3
dens(dd$log_gdp_std)
# Z-transform makes it from -2 - 2
dens(scale(dd$log_gdp))
# What about dividing by sd rather than mean?
# Ah has barely any effect, ranges changes from 6-11 to 5-10
dens(dd$log_gdp / sd(dd$log_gdp))
# What about scaling by max? Would give 0.5-1, why not use?
# Maybe just as can't interpret so easily, what is the maximum magnitude?
# Whereas for ruggedness it's more straight forward
dens(dd$log_gdp / max(dd$log_gdp))

# As for rugged:
# By default has skew
dens(dd$rugged)
# Scaling that to 0-1 just changes scale
dens(dd$rugged_std)
# Log makes it -6-2
dens(log(dd$rugged))

# It seems that the criteria are:
#   - For positive values want them to be positive (even though can enforce this in the model)
#   - Want to be as close to 0-1 range as possible
#   - Ideally want 0 or 1 as the mean, for the predictors want 0 as mean to make intercept easier to fit
#     for the output want some easy value for mean, either 0 (z-transform) or 1 (divide by mean)
# Ah he actually explains that because zero ruggedness is meaningful (i.e.
# it's a possible value on the original scale), want to keep that, hence scaling
# For LogGDP the new values mean that 1 is average, 0.8 means 80% of average etc...
# Maybe this is more useful than 0-1 when working with magnitudes?

# Ok carrying on:
# Will model log(yi) ~ Normal(ui, sigma)
# ui = a + b(ri - rbar)
# Have already logged yi (GDP)
# But ri is 0-1, and are now going to center this, why not center it earlier?
# He says that doing this centering just makes it easier to think of priors
# But why didn't we add this centering into the 0-1 transformation as well?
# Now we lose the nice interpretation of 0 ruggedness

# So this is raw rugged
dens(dd$rugged)
# Scaling to 0-1
dens(dd$rugged_std)
# Centering after scaling gives -0.2 - 0.8
dens(dd$rugged_std - mean(dd$rugged_std))
# Centering without scaling just gives -2 - 4
dens(dd$rugged - mean(dd$rugged))

# ANYWAY moving on with the model
# intercept is logGDP when ruggedness is at its mean
# So this should also be at its mean, which due to its scaling we've defined as
# 1
# So a prior centered around 1 makes sense, try with variance 1
# For slope center on zero as don't want to bias the sign, and start with 1

# Now to predictive prior (try dexp(1) for sigma)
rugged_bar <- mean(dd$rugged_std)
m8.1 <- quap(
    alist(
        log_gdp_std ~ dnorm(mu, sigma),
        mu <- a + b*(rugged_std - rugged_bar),
        a ~ dnorm(1, 1),
        b ~ dnorm(0, 1),
        sigma ~ dexp(1)
    ),
    data=dd
)
set.seed(7)
prior <- extract.prior(m8.1)

# set up the plot dimensions
plot( NULL , xlim=c(0,1) , ylim=c(0.5,1.5) ,
      xlab="ruggedness" , ylab="log GDP" )
abline( h=min(dd$log_gdp_std) , lty=2 )
abline( h=max(dd$log_gdp_std) , lty=2 )
# draw 50 lines from the prior - lots of unfeasible priors!
# Want most lines to pass through (0.215 (rugged_bar), 1) (need tighter alpha)
# Lots of lines passing beyond the max and min logGDP values (need tighter beta)
rugged_seq <- seq( from=-0.1 , to=1.1 , length.out=30 )
mu <- link( m8.1 , post=prior , data=data.frame(rugged_std=rugged_seq) )
for ( i in 1:50 ) lines( rugged_seq , mu[i,] , col=col.alpha("black",0.3) )
# The maximum slope is (1.3 - 0.7) / (1 - 0) = 0.6
# But under N(0, 1) prior more than 50% of slopes will have greater slopes than this
sum(abs(prior$b) > 0.6) / length(prior$b)

m8.2 <- quap(
    alist(
        log_gdp_std ~ dnorm(mu, sigma),
        mu <- a + b*(rugged_std - rugged_bar),
        a ~ dnorm(1, 0.1),
        b ~ dnorm(0, 0.3),
        sigma ~ dexp(1)
    ),
    data=dd
)

# Looks a lot better!
# Still some very implausible lines but better now at least
set.seed(7)
prior <- extract.prior(m8.2)
plot( NULL , xlim=c(0,1) , ylim=c(0.5,1.5) ,
      xlab="ruggedness" , ylab="log GDP" )
abline( h=min(dd$log_gdp_std) , lty=2 )
abline( h=max(dd$log_gdp_std) , lty=2 )
rugged_seq <- seq( from=-0.1 , to=1.1 , length.out=30 )
mu <- link( m8.2 , post=prior , data=data.frame(rugged_std=rugged_seq) )
for ( i in 1:50 ) lines( rugged_seq , mu[i,] , col=col.alpha("black",0.3) )

# Still no real association between ruggedness & logGDP however regardless of priors
nrow(dd)
precis(m8.2)

# 8.1.2 Adding an indicator variable isn't enough -------------------------
# An indicator variable is effectively a change in intercept, NOT slope
# But for the reasons went over before, will prefer an *index* variable
# rather than an *indicator*
# I.e. if we have a prior for an indicator variable then it implies greater prior
# uncertainty for the category with indicator = 1 than =0 as these samples
# are affected by 2 priors rather than 1
# Create index variable
dd$cid <- ifelse(dd$cont_africa == 1, 1, 2)

m8.3 <- quap(
    alist(
        log_gdp_std ~ dnorm(mu, sigma),
        mu <- a[cid] + b*(rugged_std - rugged_bar),
        a[cid] ~ dnorm(1, 0.1),
        b ~ dnorm(0, 0.3),
        sigma ~ dexp(1)
    ),
    data=dd
)
# This new model with the continent-specific intercept is much better by WAIC
compare(m8.2, m8.3)

# So here the intercept is much lower for Africa, so Africa has a 
# below average logGDP compared to other continents
precis(m8.3, depth=2)

# Test contrasts:
# Difference is reliably below 0
post <- extract.samples(m8.3)
diff_a1_a2 <- post$a[, 1] - post$a[, 2]
PI(diff_a1_a2)

# Now look at posterior predictions
rugged.seq <- seq(from=-0.1, to=1.1, length.out=50)
mu.NotAfrica <- link(m8.3, data=data.frame(cid=2, rugged_std=rugged.seq))
mu.Africa <- link(m8.3, data=data.frame(cid=1, rugged_std=rugged.seq))

mu.NotAfrica_mu <- colMeans(mu.NotAfrica)
mu.NotAfrica_ci <- apply(mu.NotAfrica, 2, PI, prob=0.97)
mu.Africa_mu <- colMeans(mu.Africa)
mu.Africa_ci <- apply(mu.Africa, 2, PI, prob=0.97)

stats <- tibble(
    mean = c(mu.NotAfrica_mu, mu.Africa_mu),
    lower = c(mu.NotAfrica_ci[1, ], mu.Africa_ci[1, ]),
    upper = c(mu.NotAfrica_ci[2, ], mu.Africa_ci[2, ]),
    cid = factor(rep(c(2, 1), each=50), levels=c(1, 2), labels=c("Africa", "Not Africa")),
    rugged_std = rep(rugged.seq, 2)
)

# NB: Africa has lower GDP through intercept, but has same slope (regression line is parallel)
dd |>
    mutate(cid = factor(cid, levels=c(1, 2), labels=c("Africa", "Not Africa"))) |>
    select(country, cid, rugged_std, log_gdp_std) |>
    ggplot(aes(x=rugged_std, colour=cid, fill=cid)) +
        geom_ribbon(aes(ymin=lower, ymax=upper), data=stats, alpha=0.3) +
        geom_point(aes(y=log_gdp_std)) +
        geom_line(aes(y=mean), data=stats) +
        scale_colour_discrete("") +
        scale_fill_discrete("") +
        labs(title="Index variable on intercept") +
        theme_minimal() +
        theme(legend.position = "bottom")

# Instead, want an interaction between rugged slope and continent
# Can achieve this by giving a separate slope to each continent
m8.4 <- quap(
    alist(
        log_gdp_std ~ dnorm(mu, sigma),
        mu <- a[cid] + b[cid]*(rugged_std - rugged_bar),
        a[cid] ~ dnorm(1, 0.1),
        b[cid] ~ dnorm(0, 0.3),
        sigma ~ dexp(1)
    ),
    data=dd
)
# Adding this slope in doesn't make such a _huge_ difference as before
# the dWAIC is 7.2, but the dSE is 6.69, so can't reliably say it makes the model's *predictions*
# a whole lot better
# HOWEVER, it does still take most of the weight, and is more consistent with our domain knowledge
compare(m8.2, m8.3, m8.4)

# Note now that we have a slope in the opposite direction for Africa! (+ve vs -ve in non-Africa)
precis(m8.4, depth=2)

# If compare using PSIS then we get warnings about outliers (influential points)
# But otherwise still the same general interpretation
compare(m8.2, m8.3, m8.4, func=PSIS)

# Can see one outlier there with k > 0.5 and a few high up
plot(PSIS(m8.4, pointwise=TRUE)$k)
# The first outlier is Lesotho, which has a very high rugged index (its the most rugged) but low GDP
# Seychelles is also an interesting example as it has high ruggedness and high GDP, the highest African GDP
# Lesotho doesn't fit the trend of high ruggedness -> high GDP in Africa, and I guess Seychelles is just an 
# outlier due to its high GDP
dd |>
    mutate(k = PSIS(m8.4, pointwise=TRUE)$k) |>
    select(country, rugged_std, log_gdp_std, cid, k) |>
    arrange(desc(k)) |>
    head()

# Plot the new posterior predictive so can see the Seychelles in context
mu.NotAfrica <- link(m8.4, data=data.frame(cid=2, rugged_std=rugged.seq))
mu.Africa <- link(m8.4, data=data.frame(cid=1, rugged_std=rugged.seq))

mu.NotAfrica_mu <- colMeans(mu.NotAfrica)
mu.NotAfrica_ci <- apply(mu.NotAfrica, 2, PI, prob=0.97)
mu.Africa_mu <- colMeans(mu.Africa)
mu.Africa_ci <- apply(mu.Africa, 2, PI, prob=0.97)

stats <- tibble(
    mean = c(mu.NotAfrica_mu, mu.Africa_mu),
    lower = c(mu.NotAfrica_ci[1, ], mu.Africa_ci[1, ]),
    upper = c(mu.NotAfrica_ci[2, ], mu.Africa_ci[2, ]),
    cid = factor(rep(c(2, 1), each=50), levels=c(1, 2), labels=c("Africa", "Not Africa")),
    rugged_std = rep(rugged.seq, 2)
)

dd <- dd |> mutate(cid = factor(cid, levels=c(1, 2), labels=c("Africa", "Not Africa"))) 

# NB: Africa has lower GDP through intercept, but has same slope (regression line is parallel)
# NB: Sometimes when running this Lesotho has a high k, othertimes not at all
dd |>
    select(country, cid, rugged_std, log_gdp_std) |>
    ggplot(aes(x=rugged_std, colour=cid, fill=cid)) +
        geom_ribbon(aes(ymin=lower, ymax=upper), data=stats, alpha=0.3) +
        geom_point(aes(y=log_gdp_std)) +
        geom_line(aes(y=mean), data=stats) +
        geom_text(aes(y=log_gdp_std, label=country), data=dd |> mutate(k = PSIS(m8.4, pointwise = TRUE)$k) |> filter(k > 0.4)) +
        scale_colour_discrete("") +
        scale_fill_discrete("") +
        labs(title="Index variable on intercept") +
        theme_minimal() +
        theme(legend.position = "bottom")

# Interpreting this model:
#   - How much does the association between ruggedness and logGDP depend on whether the nation is in Africa? (what m8.4 analysed)
#   - How much does the association of Africa with logGDP depend upon ruggedness?

# To look at the second phrasing can compute the difference between a country's GDP whether it is in Africa or not,
# conditioning on ruggedness (holding constant)
# COUNTER-FACTUAL PLOT!
# Useful to consider this perspective, as we often tend to favour the intepretation where we can imagine changing one of the predictors
# (i.e. ruggedness), but here we can't imagine changing continent so useful to see
muA <- link(m8.4, data=data.frame(cid=1, rugged_std=rugged_seq))
muN <- link(m8.4, data=data.frame(cid=2, rugged_std=rugged_seq))
delta <- muA - muN
pi <- apply(delta, 2, PI)
tibble(rugged_std=rugged_seq, diff = colMeans(delta),
       lower = pi[1, ], upper=pi[2, ]) |>
    ggplot(aes(x=rugged_seq, y=diff)) +
        geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
        geom_hline(yintercept=0, linetype="dashed") +
        annotate("text", x=-Inf, y=0.05, label="Africa higher GDP", hjust=0) +
        annotate("text", x=-Inf, y=-0.05, label="Africa lower GDP", hjust=0) +
        geom_line() +
        theme_minimal() +
        labs(x="Rugged std", y="Expected difference between African and non-African GDP (Africa - NonAfrica)")

# 8.3. Continuous interactions --------------------------------------------
# Aim is to demonstrate how hard it is to interpret interactions between continuous vars
# "triptych" plot for vizualising interactions
data(tulips)
d <- tulips
str(d)
# Want to predict blooms based on amount of water and shade (both ordinal 1-3)
# Bed is a grouping variable
# Makes sense both light and water affect blooms, but is there an interaction here? I.e. with little light the amount of water makes less difference
# then if had full light

# 2 models:
#   water & shade
#   water & shade & interaction
# Causal scenario is W -> B <- S
# But this just means that B = f(W, S). It doesn't specify the nature of the functional form!
# Will center water and shade and scale blooms by max
# NB: does it make sense to center ordinal vars?!
# And we want to keep 0 blooms as a reference point as we have it in the dataset 
#   - Also makes it easier to assign a meaningful prior without knowing much about the blooms values themselves
summary(d)

d$blooms_std <- d$blooms / max(d$blooms)
d$water_cent <- d$water - mean(d$water)
d$shade_cent <- d$shade - mean(d$shade)

# What does 'mean' mean here? Really should be modal value?
table(d$water)
mean(d$water) # Mean is 2
# So centered values are -1, 0, and 1
table(d$water_cent)

# same for shade
table(d$shade)
mean(d$shade) 
table(d$shade_cent)

# Initial priors guess:
#  - interecpt: when have average water and average shade, would expect average blooms. Don't guarantee that 0.5 is average but should be
# close enough.
# In fact 0.36 is avg but 0.5 is ok, so N(0.5, maybe 1?)
mean(d$blooms_std)
#   - What would expect for the slopes?
#   - The maximum range is 2 (-1 to 1), which should encompass 0-1 output, equal to a slope of +/-0.5. So N(0, 0.25) would be sensible imo
# Richard has gone with a ~ N(0.5, 1), b ~ N(0, 1) - can't understand the slopes being this flat

# He then goes through it more and for the intercept, we know that all values must be [0, 1], so he tightens the sd to be 0.25, so only 5% of
# mass is outside this interval. I.e. centered on 0.5, we want max distance of 0.5, which is equal to ~2sd, hence 0.25
# Likewise for slope we want no more than 0.5 difference either side of 0 (as this is the max rise/run in the dataset)

m8.5 <- quap(
    alist(
        blooms_std ~ dnorm(mu, sigma),
        mu <- a + bw*water_cent + bs*shade_cent,
        a ~ dnorm(0.5, 0.25),
        bw ~ dnorm(0, 0.25),
        bs ~ dnorm(0, 0.25),
        sigma ~ dexp(1)
    ),
    data=d
)

# Will look at priors and posteriors later, but first will quickly fit the interaction model to have a comparison point
# Priors for interactions are hard to grok!
# Strongest effect is one where max shade means water has zero effect
# Basically try same prior as the same coefficients
# NB: Interactions are *always* bi-directional! I.e. a*b works in both directions as the increase in a for a 1 unit increase in b and vice-versa
m8.6 <- quap(
    alist(
        blooms_std ~ dnorm(mu, sigma),
        mu <- a + bw*water_cent + bs*shade_cent + bws*water_cent*shade_cent,
        a ~ dnorm(0.5, 0.25),
        bw ~ dnorm(0, 0.25),
        bs ~ dnorm(0, 0.25),
        bws ~ dnorm(0, 0.25),
        sigma ~ dexp(1)
    ),
    data=d
)

# Here going to condition on one variable (say shade) and plot the bivariate relationship between the other predictor (shade) and blooms
# Will do this for different values of shade. Here it's easy as there's just 3 values, but can use the same approach for a fully continuous variable
# Just choose 3 appropriate values (low, mean, high) say

# Will create my own plot rather than using Richard's code here just to ensure I grok this
shades <- c(-1, 0, 1)
waters <- c(-1, 0, 1)
models <- list("No interaction"=m8.5, "Interaction"=m8.6)
stats <- map_dfr(shades, function(s) {
    map_dfr(models, function(mod) {
        post <- link(mod, data=tibble(water_cent=waters, shade_cent=s))
        mu <- colMeans(post)
        PI <- apply(post, 2, PI)
        tibble(water_cent=waters, shade_cent=s, mu=mu, lower=PI[1, ], upper=PI[2, ])
    },
    .id="model"
    )
})

# This looks the same as the book
# The no interaction model believes there's always a strong effect of water on blooms and that shade always hurts
# Note how the blue line decreases from left to right, but is always parallel because this is just a variable intercept
# With the interaction the model believes that as the amount of shade increases, water has less impact (shallower slope)
# So the strongest effect on growth is the left panel, right x-axis where shade is -1 and water is +1.
d |>
    ggplot(aes(x=water_cent, y=blooms_std)) +
        geom_ribbon(aes(ymin=lower, y=mu, ymax=upper), data=stats,
                    alpha=0.2) +
        geom_line(aes(y=mu), data=stats) +
        geom_point() +
        facet_grid(model~shade_cent) +
        theme_bw() +
        scale_colour_discrete("") +
        scale_fill_discrete("") +
        theme(legend.position = "bottom")

# Also going to plot prior predictions!
# Just use the link function but with post=prior
# And also plot the actual lines here rather than shaded regions as we're looking at the priors of regression lines
# Not output shaded regions where individual predicted points might fall
stats <- map_dfr(models, function(mod) {
    prior <- extract.prior(mod)
    map_dfr(shades, function(s) {
        post <- link(mod, post=prior, data=tibble(water_cent=waters, shade_cent=s))
        colnames(post) <- paste0("water_", waters)
        post <- as_tibble(post)
        post |>
            mutate(obs = row_number()) |>
            pivot_longer(-obs, names_pattern="water_([0-9-]+)", names_to="water_cent", values_to="blooms_std") |>
            mutate(water_cent = as.numeric(water_cent),
                   shade_cent=s) |>
            filter(obs <= 20)
    }
    )
}, .id="model")

# Ensuring that use the same prior draws for all 3 links within the same model (i.e. swap the loops over
# shade and models) and it looks the same as the book.
# All priors are fine, but "only weekly realistic" as while most of the lines stay within the valid outcome space, there are still unfeasible trends
# (seen at least 4 here and have only plotted 20)
d |>
    ggplot(aes(x=water_cent, y=blooms_std)) +
        geom_line(aes(group=obs), alpha=0.3, data=stats) +
        geom_line(aes(group=obs), alpha=1, linewidth=1, data=stats |> filter(obs==17)) +
        geom_hline(yintercept = 1, linetype="dashed") +
        geom_hline(yintercept = 0, linetype="dashed") +
        geom_point() +
        facet_grid(model~shade_cent) +
        theme_bw() +
        scale_colour_discrete("") +
        scale_fill_discrete("") +
        theme(legend.position = "bottom")

# Nb: The interaction model is preferred for predictive accuracy with a huge weight, but with low dWAIC and dSE
compare(m8.5, m8.6)

# Practice ----------------------------------------------------------------

# 8E1
# 1) temperature (yeast has a stronger effect at higher temps)
# 2) education and job type, the more educated you are in a higher paying field the more you earn
# 3) 


