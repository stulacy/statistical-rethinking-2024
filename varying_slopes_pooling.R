library(rethinking)
library(tidyverse)
data("chimpanzees")

dlist <- list( 
    pulled_left = chimpanzees$pulled_left ,
    prosoc_left = chimpanzees$prosoc_left ,
    condition = as.integer( 2 - chimpanzees$condition ) ,
    actor = as.integer( chimpanzees$actor ) 
)


cluster_variance <- ulam(
    alist(
        pulled_left ~ bernoulli(theta),
        logit(theta) <- a + aj[actor] + bp[condition]*prosoc_left,
        aj[actor] ~ normal( 0 , sigma_actor ),
        a ~ normal(0,4),
        bp[condition] ~ normal(0,1),
        sigma_actor ~ exponential(1)
    ) ,
    data=dlist, chains=2 , cores=1 , sample=TRUE, log_lik = TRUE )

cluster_mean_variance <- ulam(
    alist(
        pulled_left ~ bernoulli(theta),
        logit(theta) <- aj[actor] + bp[condition]*prosoc_left,
        aj[actor] ~ normal( abar , sigma_actor ),
        abar ~ normal(0, 1),
        bp[condition] ~ normal(0,1),
        sigma_actor ~ exponential(1)
    ) ,
    data=dlist, chains=2 , cores=1 , sample=TRUE, log_lik=TRUE )

# And also make unpooled model for comparison
m2_unpooled <- ulam(
    alist(
        pulled_left ~ bernoulli(theta),
        logit(theta) <- aj[actor] + bp[condition]*prosoc_left,
        aj[actor] ~ normal(0, 1),
        bp[condition] ~ normal(0,1),
        sigma_actor ~ exponential(1)
    ) ,
    data=dlist, chains=2 , cores=1 , sample=TRUE, log_lik=TRUE )

# No difference in WAIC between the 2 varying slopes model, although the unpooled model fares slightly worse
compare(cluster_variance, cluster_mean_variance, m2_unpooled)

# Compare on posterior predictions
ndata <- expand_grid(
    prosoc_left = unique(dlist$prosoc_left),
    condition=unique(dlist$condition),
    actor=unique(dlist$actor)
)

models <- list(
    "Cluster variance"=cluster_variance,
    "Cluster mean + variance"=cluster_mean_variance,
    "Unpooled"=m2_unpooled
)

preds <- map_dfr(models, function(mod) {
    preds <- link(mod, data=ndata)
    ndata |>
        mutate(
            pred_mean = colMeans(preds),
            pred_lower = apply(preds, 2, PI)[1, ],
            pred_upper = apply(preds, 2, PI)[2, ],
        )
}, .id="model")

# Calculate observed rates to plot
observed <- chimpanzees |>
    mutate(
        condition = as.integer(2 - condition),
        prosoc_left = as.factor(prosoc_left)
    ) |>
    group_by(actor, condition, prosoc_left) |>
    summarise(prop_obs = mean(pulled_left, na.rm=T)) |>
    ungroup()

# And plot!
# There's some real differences between the unpooled model and the 2 partially pooled models,
# particularly in actors 2 & 7, but otherwise it seems to make no difference whether the cluster
# have means & variances or just variances
preds |>
    mutate(prosoc_left = as.factor(prosoc_left)) |>
    ggplot(aes(x=prosoc_left, y=pred_mean, colour=model)) +
        geom_point(position=position_dodge(width=0.8)) +
        geom_point(aes(y=prop_obs), data=observed, colour="black") +
        geom_errorbar(aes(ymin=pred_lower, ymax=pred_upper), position=position_dodge(width=0.8)) +
        facet_grid(actor ~ condition, labeller = label_both, scales="free_y") +
        theme_bw() +
        scale_colour_brewer(palette="Dark2") +
        theme(
            legend.position = "bottom"
        )

# Comparing the coefficients themselves and abar is actually quite close to 0 and its CI covers 0 quite comfortably, so it appears
# to have a relatively weak effect. sigma_actor is largely the same, so is there where most of the pooling is done?
# The only other difference is that the CI on the cluster variance model only is that the CIs are much wider for each actor
# and indeed there is a bit more pooling
plot(coeftab(cluster_variance, cluster_mean_variance))
