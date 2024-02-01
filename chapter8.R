library(rethinking)
library(tidyverse)
library(patchwork)

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
# 3) throttle - the more throttle the stronger the effect of gas?

# 8E2
# 1, don't think any others (more ORs rather than ANDs)

# 8E3
# caramelising = heat + moisture + heat * moisture

# 8M1
# Need cool temperature AND water AND shade

# 8M2
# blooms = (1 - Hot) * (alpha + b_1 * light + b_2 * moisture + b_12 * light * moisture)

# 8M3
# Not really an interaction as don't have a third variable, just a bivariate relationship
# Could simulate a dataset as wolves <- rnorm(10, 50, 10), ravens <- rnorm(10, wolves*10, 20)
# Probably not linear as many biological phenomena are multiplicative. Can think concretely
# about this as saying that the addition of each wolf increases implicatively the number
# of ravens

# 8M4
data(tulips)
d <- tulips
d$blooms_std <- d$blooms / max(d$blooms)
d$water_cent <- d$water - mean(d$water)
d$shade_cent <- d$shade - mean(d$shade)

# 2 models:
#   water & shade
#   water & shade & interaction
# Original models
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
# Had to increase bws variance to fit
set.seed(17)
m8.6_constrained <- quap(
    alist(
        blooms_std ~ dnorm(mu, sigma),
        mu <- a + bw*water_cent - bs*shade_cent + bws*water_cent*shade_cent,
        a ~ dnorm(0.5, 0.25),
        bw ~ dexp(1),
        bs ~ dexp(1),
        bws ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ),
    data=d
)
precis(m8.6_constrained)

# Prior predictive
shades <- c(-1, 0, 1)
waters <- c(-1, 0, 1)
models <- list("No interaction"=m8.5, "Interaction"=m8.6, "Interaction constrained"=m8.6_constrained)
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

# Has much stronger effect now!
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

# Previously the interaction was negative only, while now it is forced to be positive
# I think this is because I've forced a negative b_shade by forcing the coefficient itself positive
# but sticking a negative sign in front of it in the regression
# And indeed bs is -0.11 in the free model and 0.11 in the constrained model (with negative sign)
# The shade and water slopes are pretty much the exact same value in the constrained
# So why is this happening?

plot(coeftab(m8.6, m8.6_constrained))
# a + bw*water_cent - bs*shade_cent + bws*water_cent*shade_cent
# If water is +ve but shade is -ve then for an average water & shade we get average blooms (a)
# I.e. b0 = a
# If water increases by 1 and shade decreases by 1 then our output is:
# b1 = a + bw*1 - bs*-1 + bws*1*-1  
# b1 = a + bw + bs - bws  
# b1 > b0
# a + bw + bs - bws > a
# bw + bs - bws > 0
# bw + bs > bws
post <- extract.samples(m8.6_constrained)
# And yep this is always met!
mean((post$bw + post$bs) > post$bws)

# Can we determine that the interaction prior will always be negative from this?
# Or is this just enough for the prior, that it will be less than bw + bs?

# 8H1
summary(d$bed)
# Include bed as categorical intercept
d$bid <- as.integer(as.factor(d$bed))
m8h1 <- quap(
    alist(
        blooms_std ~ dnorm(mu, sigma),
        mu <- a[bid] + bw*water_cent + bs*shade_cent + bws*water_cent*shade_cent,
        a[bid] ~ dnorm(0.5, 0.25),
        bw ~ dnorm(0, 0.25),
        bs ~ dnorm(0, 0.25),
        bws ~ dnorm(0, 0.25),
        sigma ~ dexp(1)
    ),
    data=d
)

# 8H2
# Slightly better predictions when using bed, although massive standard error (because so few points?)
compare(m8h1, m8.6)

# The bed coefficients aren't that different, the slopes are the same, sigma is the same
# The original intercept looks the average of the 3 beds
# So the main benefit is just adding a bit more observation specific variance, but not a huge difference
plot(coeftab(m8h1, m8.6))

# 8H3
# a)
# Book says m8.5 but think should be m8.4
data(rugged)
d <- rugged
d$log_gdp <- log(d$rgdppc_2000)
dd <- d[complete.cases(d$rgdppc_2000), ]
dd$log_gdp_std <- dd$log_gdp / mean(dd$log_gdp)
dd$rugged_std <- dd$rugged / max(dd$rugged)
dd$cid <- ifelse(dd$cont_africa == 1, 1, 2)
dd$cid_fact <- factor(dd$cid, levels=c(1, 2), labels=c("Africa", "Not Africa")) 
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

psis <- PSIS(m8.4, pointwise=TRUE)
# Other outliers are:
#  - Switzerland (very high rugged for its log GDP outside of Africa which has 
# negative relationship)
#  - Lesotho (lower GDP than would expect from ruggedness)
#  - Tajikstan (much lower GDP than would expect from ruggedness outside of Africa)
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
    cid_fact = factor(rep(c(2, 1), each=50), levels=c(1, 2), labels=c("Africa", "Not Africa")),
    rugged_std = rep(rugged.seq, 2)
)

p1 <- dd |> 
    select(country, cid_fact, rugged_std, log_gdp_std) |>
    ggplot(aes(x=rugged_std, colour=cid_fact, fill=cid_fact)) +
        geom_ribbon(aes(ymin=lower, ymax=upper), data=stats, alpha=0.3) +
        geom_point(aes(y=log_gdp_std)) +
        geom_line(aes(y=mean), data=stats) +
        geom_text(aes(y=log_gdp_std, label=country), data=dd |> mutate(k = psis$k) |> filter(k > 0.4)) +
        scale_colour_discrete("") +
        scale_fill_discrete("") +
        labs(title="Index variable on intercept") +
        theme_minimal() +
        theme(legend.position = "bottom")
p1

# b) Using robust regression:
m8.4_robust <- quap(
    alist(
        log_gdp_std ~ dstudent(2, mu, sigma),
        mu <- a[cid] + b[cid]*(rugged_std - rugged_bar),
        a[cid] ~ dnorm(1, 0.1),
        b[cid] ~ dnorm(0, 0.3),
        sigma ~ dexp(1)
    ),
    data=dd
)
psis <- PSIS(m8.4_robust, pointwise=TRUE)
mu.NotAfrica <- link(m8.4_robust, data=data.frame(cid=2, rugged_std=rugged.seq))
mu.Africa <- link(m8.4_robust, data=data.frame(cid=1, rugged_std=rugged.seq))

mu.NotAfrica_mu <- colMeans(mu.NotAfrica)
mu.NotAfrica_ci <- apply(mu.NotAfrica, 2, PI, prob=0.97)
mu.Africa_mu <- colMeans(mu.Africa)
mu.Africa_ci <- apply(mu.Africa, 2, PI, prob=0.97)

stats <- tibble(
    mean = c(mu.NotAfrica_mu, mu.Africa_mu),
    lower = c(mu.NotAfrica_ci[1, ], mu.Africa_ci[1, ]),
    upper = c(mu.NotAfrica_ci[2, ], mu.Africa_ci[2, ]),
    cid_fact = factor(rep(c(2, 1), each=50), levels=c(1, 2), labels=c("Africa", "Not Africa")),
    rugged_std = rep(rugged.seq, 2)
)

# Now there are no outliers!
p2 <- dd |> 
    select(country, cid_fact, rugged_std, log_gdp_std) |>
    ggplot(aes(x=rugged_std, colour=cid_fact, fill=cid_fact)) +
        geom_ribbon(aes(ymin=lower, ymax=upper), data=stats, alpha=0.3) +
        geom_point(aes(y=log_gdp_std)) +
        geom_line(aes(y=mean), data=stats) +
        geom_text(aes(y=log_gdp_std, label=country), data=dd |> mutate(k = psis$k) |> filter(k > 0.4)) +
        scale_colour_discrete("") +
        scale_fill_discrete("") +
        labs(title="Index variable on intercept") +
        theme_minimal() +
        theme(legend.position = "bottom")
p2
# Now all k values are < 0.32 
# Where even is Seychelles on this?
dd |>
    select(country) |>
    mutate(k = psis$k) |>
    arrange(desc(k)) |>
    as_tibble()
# Down to 0.118 now!
dd |>
    as_tibble() |>
    select(country) |>
    mutate(k = psis$k) |>
    filter(country == 'Seychelles')

# Notice that the Not Africa slope is more extreme now, while
# the Africa one hasn't changed much
p1 | p2

# Looking on this scale and actually the slopes have barely changed, so
# robust regression hasn't substantially changed the interpretation
plot(coeftab(m8.4, m8.4_robust))

# H84
data("nettle")
d <- as_tibble(nettle)
d
summary(d)
# Hypothesis: language diversity is function of food security
# Languages per capita as outcome
d$lang_per_capita <- d$num.lang / d$k.pop
d$loglang <- log(d$lang_per_capita)
d$loglang_scale <- d$loglang / mean(d$loglang)
# This gives value in [0-1.7] that is roughly normally distributed
dens(d$loglang_scale)
# Mean growing season
# Z-score should be sufficient
dens(scale(d$mean.growing.season))
d$mean_scale <- scale(d$mean.growing.season)
# Growing season is very skewed, even under Z-transform!
dens(scale(d$sd.growing.season))
# See here's the raw value
dens(d$sd.growing.season)
# I think a log transform is better here
# Ah but 2 values have SD of 0 so log becomes infinite
dens(log(d$sd.growing.season))
# So maybe just scaling 0-1 is best here
dens(d$sd.growing.season / max(d$sd.growing.season)) 
d$sd_scale <- d$sd.growing.season / max(d$sd.growing.season)
# Now for Area
# Log Area gives very big values, the units of area must be such that they are very big
# and there are several orders of magnitude difference
dens(log(d$area))
# Yep can see that's the case in the raw plot
dens(d$area)
# standardized log area should be useful
dens(scale(log(d$area)))
d$area_scale <- scale(log(d$area))

# Use log of this as outcome to constrain positive
# mean growing season & sd of this, both main effects & interaction
# Need to decide sensible priors
# Whether to use WAIC
# a) lang ~ mean growing
# For priors then, when mean growing is at 0 (mean), we'd expect outcome
# to be mean too (1).  Outcome has range 0-2 so N(1, 0.5) is quite tight but gives realistic values
# For slope Growing season can value from -2 to 2, so max rise/run is 2/4 or 0.5. 
# So N(0, 0.25) should be fine
m_h84_1 <- quap(
    alist(
        loglang_scale ~ dnorm(mu, sigma),
        mu <- a + bM * mean_scale,
        a ~ dnorm(1, 0.5),
        bM ~ dnorm(0, 0.25),
        sigma ~ dexp(1)
    ), data=d
)
# The intercept is around 1 as expected
# The slope is actually negative, so as the farming season increases, the number
# of languages spoken per capita decreases. This is counter to what the hypothesis was
plot(precis(m_h84_1))

# b) Accounting for area
# prior on area is scale[-3,3] so max rise/run of 2/6 or 0.3,
# Will just use N(0, 0.25) to be safe and consistent
m_h84_2 <- quap(
    alist(
        loglang_scale ~ dnorm(mu, sigma),
        mu <- a + bM * mean_scale + bA * area_scale,
        a ~ dnorm(1, 0.5),
        bM ~ dnorm(0, 0.25),
        bA ~ dnorm(0, 0.25),
        sigma ~ dexp(1)
    ), data=d
)
# This doesn't change the effect of growing season but area is slighly positive
# So when we condition on the size of the country, the length of the growing season
# still has a negative impact on the number of languages
plot(coeftab(m_h84_1, m_h84_2))

# c) lang ~ sd growing
# SD of growing season. Here the hypothesis is that this is negatively associated,
# as the more variable the growing season, the more you need bigger communities for resilience
# and hence fewer languages
# For a prior this just goes from 0-1 so max rise/run is 2/1 =2, so will use prior of N(0, 1)
m_h84_3 <- quap(
    alist(
        loglang_scale ~ dnorm(mu, sigma),
        mu <- a + bS * sd_scale,
        a ~ dnorm(1, 0.5),
        bS ~ dnorm(0, 1),
        sigma ~ dexp(1)
    ), data=d
)
# Lots of uncertainty here, but bS is positive, the opposite of what was theorized!
plot(precis(m_h84_3))

# d) lang ~ sd growing + log(area)
m_h84_4 <- quap(
    alist(
        loglang_scale ~ dnorm(mu, sigma),
        mu <- a + bS * sd_scale + bA * area_scale,
        a ~ dnorm(1, 0.5),
        bS ~ dnorm(0, 1),
        bA ~ dnorm(0, 0.25),
        sigma ~ dexp(1)
    ), data=d
)

# Again area has a slight positive effect, but this hasn't really done much to bS except to
# make it less certain and overlap with 0
plot(coeftab(m_h84_3, m_h84_4))

# e) lang ~ mean growing * sd growing
m_h84_5 <- quap(
    alist(
        loglang_scale ~ dnorm(mu, sigma),
        mu <- a + bM * mean_scale + bS * sd_scale + bMS * mean_scale * sd_scale,
        a ~ dnorm(1, 0.5),
        bM ~ dnorm(0, 0.25),
        bS ~ dnorm(0, 1),
        bMS ~ dnorm(0, 1),
        sigma ~ dexp(1)
    ), data=d
)

# This doesn't change the signs on the Mean and SD params
# The interaction is positive, so if the mean growing season increases, then you get more languages
# if there is more variance in growing seasons (how does this work in practice?!)
plot(coeftab(m_h84_1, m_h84_3, m_h84_5))

# The model with interactions gives best WAIC by some way (look at penalty)
# SD has lowest WAIC, seems that Mean contributes more
compare(m_h84_1, m_h84_2, m_h84_3, m_h84_4, m_h84_5)

# 8H5
data("Wines2012") 
d <- Wines2012 |> as_tibble()
summary(d)
head(d)
nrow(d)
# model score ~ judge + wine
# Score is normal so will just z-transform
dens(d$score)
d$score_scale <- scale(d$score)
# Create index variables for judge and wine
d$judge_ind <- as.integer(as.factor(d$judge))
d$wine_ind <- as.integer(as.factor(d$wine))
# Scores range from -3 - 3
# Priors are intercepts, so want average on 0 (average standardized score)
# And max range -3 - 3 so N(0, 1.5)
dens(d$score_scale)
m_h85 <- quap(
    alist(
        score_scale ~ dnorm(mu, sigma),
        mu <- J[judge_ind] + W[wine_ind],
        J[judge_ind] ~ dnorm(0, 1.5),
        W[wine_ind] ~ dnorm(0, 1.5),
        sigma ~ dexp(1)
    ),
    data=d
)

# Wine 18 looks to be the worst, with wine 4 just about the best
plot(precis(m_h85, depth=2, pars=c("W")))
# There is a lot of variance among judges! 8 is the harshest, nearly giving 1SD below the 
# average judge, with 5 the highest scoring.
plot(precis(m_h85, depth=2, pars=c("J")))
# Looks like more variation between wines than judges!

# 8H6
# Create index variables
d$colour_ind <- as.integer(as.factor(d$flight))  # Red=1, White=2
d$wine_country_ind <- as.integer(as.factor(d$wine.amer))  # French=1, American=2
d$judge_country_ind <- as.integer(as.factor(d$judge.amer))  # French=1, American=2
# Model with same priors
m_h86 <- quap(
    alist(
        score_scale ~ dnorm(mu, sigma),
        mu <- C[colour_ind] + WC[wine_country_ind] + JC[judge_country_ind],
        C[colour_ind] ~ dnorm(0, 1.5),
        WC[wine_country_ind] ~ dnorm(0, 1.5),
        JC[judge_country_ind] ~ dnorm(0, 1.5),
        sigma ~ dexp(1)
    ),
    data=d
)
# No real effect!
# I.e. judges scores are not biased in terms of where they are from, the colour of the wine
# Or the location of the wine. This is generally good and what we'd hope to find!
# This gives us confidence in the previous model that some wines are really better than others
# And the judges just have their own natural 'average' score that is independent of these factors
plot(precis(m_h86, depth=2))

# NB: The model is very sensitive to priors!
# If I use smaller priors, say 0.1, then we get smaller standard errors on this scale
# This is because there is relatively limited data
# But the megnitudes are generally the same, and there isn't any coefficient
# that doesn't have zero in its interval
m_h86_tighter <- quap(
    alist(
        score_scale ~ dnorm(mu, sigma),
        mu <- C[colour_ind] + WC[wine_country_ind] + JC[judge_country_ind],
        C[colour_ind] ~ dnorm(0, 0.1),
        WC[wine_country_ind] ~ dnorm(0, 0.1),
        JC[judge_country_ind] ~ dnorm(0, 0.1),
        sigma ~ dexp(1)
    ),
    data=d
)
plot(coeftab(m_h86, m_h86_tighter))

# If we were to read into these small effect sizes, we'd see that the American judges (JC=2)
# give slightly higher scores than french judges (JC=1), but American wines get slightly lower
# scores (WC=2) than French wines (WC=1)

# 8H7
# All 2-way interactions
# Will revert back to index variables for this
# Think model will struggle without sensible priors!
d$colour_red <- as.integer(d$flight == 'red')
# Index prior all have the same interpretation - the _additional_ impact on score
# of the '1' value. So really this should be between -3, 3 (i.e. N(0, 1.5) as used before)
# but probably want to keep a bit tighter to help estimates so will use 0.5
m_h87 <- quap(
    alist(
        score_scale ~ dnorm(mu, sigma),
        mu <- C * colour_red + WC * wine.amer + JC * judge.amer + CWC * colour_red * wine.amer + CJC * colour_red * judge.amer + WCJC * wine.amer * judge.amer,
        C ~ dnorm(0, 0.5),
        WC ~ dnorm(0, 0.5),
        JC ~ dnorm(0, 0.5),
        CWC ~ dnorm(0, 0.5),
        CJC ~ dnorm(0, 0.5),
        WCJC ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ),
    data=d
)
# Ok interesting...
# Red wine seems to score slightly higher than white (not seen before)
# American wine seems to score slightly lower than French (seen in last model)
# American judges seems to score slightly higher than French (seen in last model)
# As for the interactions, CWC is the most significant by far!
# So the effect of red wine and american is extra negative, so red wine scores well,
# American scores badly, but red AND american do extra worse
plot(precis(m_h87))

# Quickly note that the Prior doesn't seem to impact indicator variables as much!
# Below I show that with N(0, 1.5) there's barely any difference
# I suspect this is due to needing to estimate fewer values
# Ideally if you were using index variables for a lot of categorical data you'd use 
# a multi-level model to pool the priors
m_h87_wider <- quap(
    alist(
        score_scale ~ dnorm(mu, sigma),
        mu <- C * colour_red + WC * wine.amer + JC * judge.amer + CWC * colour_red * wine.amer + CJC * colour_red * judge.amer + WCJC * wine.amer * judge.amer,
        C ~ dnorm(0, 1.5),
        WC ~ dnorm(0, 1.5),
        JC ~ dnorm(0, 1.5),
        CWC ~ dnorm(0, 1.5),
        CJC ~ dnorm(0, 1.5),
        WCJC ~ dnorm(0, 1.5),
        sigma ~ dexp(1)
    ),
    data=d
)
plot(coeftab(m_h87, m_h87_wider))

# Now to look at the model through the lens of predictions, not just the coefficients
# Colour by colour
# x-axis is wine country
# Facet is judge country
in_data <- expand_grid(colour_red=c(0, 1), wine.amer=c(0, 1), judge.amer=c(0, 1))
mu <- link(m_h87, data=in_data)

# Can see this visually:
# Red wines tend to do better than white for France,
# But American whites fare better than their reds
# And American judges are slightly more positive all round
in_data |>
    mutate(mean = colMeans(mu),
           lower = apply(mu, 2, PI)[1, ],
           upper = apply(mu, 2, PI)[2, ],
           colour=factor(colour_red, labels=c("White", "Red")),
           country=factor(wine.amer, labels=c("France", "America")),
           judge=factor(judge.amer, labels=c("France", "America"))
    ) |>
    ggplot(aes(x=country, y=mean, colour=colour)) +
        geom_point() +
        geom_errorbar(aes(ymin=lower, ymax=upper)) +
        facet_wrap(~judge) +
        theme_bw() +
        scale_colour_manual("", values=c("Green", "Red")) +
        guides(colour="none") +
        labs(x="Wine country", y="Predicted score (z)")
