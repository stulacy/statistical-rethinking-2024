library(rethinking)
library(tidyverse)
library(ggpubr)

data("Monks")
Monks <- as_tibble(Monks)
# Each row is a dyad (dyad_id), i.e. a Monk-Monk relationship of 2 monks:
#  - A (integer Monk A ID)
#  - B (integer Monk B ID)
#  - A_name (name of Monk A - unneeded)
#  - B_name (name of Monk B - unneeded)
# And the outcomes are whether the monks like each other:
#  - like_AB: Number (/3) times Monk A likes Monk B
#  - like_BA: Number (/3) times Monk B likes Monk A
#  - dislike_AB: Number (/3) times Monk A dislikes Monk B
#  - dislike_BA: Number (/3) times Monk B dislikes Monk A
Monks

# Question 1 --------------------------------------------------------------
# The first question asks us to model the amount of reciprocity from each monk,
# i.e. if Monk A likes Monk B a lot, is that reciprocated?
# Going to ignore dislikes for now
# So can just model the liking slope as univariate
# And fit the dyad multi-variate effect
# Will use the example social network code that has dyads
# non-centered, and will remove the dislike bit
# The final question is whether to model number of likes as ordinal or Poisson
# And strictly we should use ordinal since it's a discrete choice (0-3)
m_q1 <- ulam(
    alist(
        # I think can use the same kappa here since it's symmetric
        like_AB ~ ordered_logistic(phiAB, kappa),
        like_BA ~ ordered_logistic(phiBA, kappa),
        # So the number of likes from A->B is a function of:
        #  - overall mean likes (a)
        #  - how many likes A likes to give away in general (likes_gr[A, 1])
        #  - how many likes B tends to receive in general (likes_gr[B, 2])
        #  - the specific number of likes from A to B (d[dyad_id, 1])
        phiAB <- a + likes_gr[A, 1] + likes_gr[B, 2] + d[dyad_id, 1],
        phiBA <- a + likes_gr[B, 1] + likes_gr[A, 2] + d[dyad_id, 2],
        a ~ normal(0, 1),
        kappa ~ normal(0, 1.5),
        
        # likes/dislikes matrix of varying effects. first col is how many likes monk A
        # gives away, second column is how many likes B receives
        transpars> matrix[18, 2]:likes_gr <-
            compose_noncentered(sigma_l, L_Rho_l, z_l),
        matrix[2, 18]:z_l ~ normal(0, 1),
        # Relatively broad prior
        cholesky_factor_corr[2]:L_Rho_l ~ lkj_corr_cholesky(4),
        vector[2]:sigma_l ~ exponential(1),
        
        # dyad effects - non-centered
        # Have 153 dyads
        transpars> matrix[153, 2]:d <-
            compose_noncentered(rep_vector(sigma_d, 2), L_Rho_d, z_d),
        matrix[2, 153]:z_d ~ normal(0, 1),
        # Relatively broad prior
        cholesky_factor_corr[2]:L_Rho_d ~ lkj_corr_cholesky(4),
        sigma_d ~ exponential(1),
        
        # compute centered correlation matrix for dyad effects
        gq> matrix[2,2]:Rho_d <<- Chol_to_Corr(L_Rho_d),
        gq> matrix[2,2]:Rho_l <<- Chol_to_Corr(L_Rho_l)
    ),
    # NB: for ordered logistic the outcome must be positive, think factor levels
    data=Monks |> mutate(like_AB = like_AB + 1, like_BA = like_BA + 1),
    chains=4, cores=4, iter=2000
)
# Good Rhat and N_eff
precis(m_q1, pars=c("likes_gr", "sigma_l", "sigma_d", "Rho_d", "a"), depth=3)

# Not a huge amount of variability in how many likes each monk gives (the first 18 rows
# and also a low sigma_l[1]), but a lot of variability in how many likes a monk receives!
# This is shown by sigma_l[2] > sigma_l[1], and also more general variability in the second column of
# likes_gr
plot(precis(m_q1, pars=c("likes_gr", "sigma_l"), depth=3))

# Is there any evidence of reciprocity?
# NB: Does this mean the general reciprocity? I.e. the correlation between giving likes out and receiving them?
# Or does it refer to the dyad specific effects?
# I think the former..
# In which case there is very little reciprocity, effectively zero correlation between
# someone who gives likes out and receives them
precis(m_q1, pars=c("Rho_l"), depth=3)
# If it's referring to the dyad-specific effects then there is a much stronger correlation
# This is the correlation between the giving and receiving for each dyad
precis(m_q1, depth=3, pars=c("Rho_d", "sigma_d"))

# Plotting mean dyad effects and can see a fair trend
# Since this is ordinal data can see the bulk of the density is around 0,0, i.e. no likes given
# But when likes are given they tend to be given in both directions
post <- extract.samples(m_q1)
dy1 <- colMeans(post$d[, , 1])
dy2 <- colMeans(post$d[, , 2])
plot(dy1, dy2)

# Question 2 --------------------------------------------------------------
# Now want to model dislikes as well
# So will include this as separate likelihoods
# Can I use any of the existing priors or do I need to create a dislikes_gr matrix?
# Could argue that the dislikes recieving is the opposite of the likes receiving?
# Firstly will try with a completely separate matrix for dislikes
# And then will try a 4-column matrix with giving out likes, giving out dislikes, receiving likes, receiving dislikes
# And then will try a 3-column matrix with giving out likes, giving out dislikes, receiving likes
#   and in this will take receiving dislikes as the opposite of receiving likes
m_q2_a <- ulam(
    alist(
        # I think can use the same kappa here since it's symmetric
        like_AB ~ ordered_logistic(phiAB_l, kappa_l),
        like_BA ~ ordered_logistic(phiBA_l, kappa_l),
        dislike_AB ~ ordered_logistic(phiAB_d, kappa_d),
        dislike_BA ~ ordered_logistic(phiBA_d, kappa_d),
        # So the number of likes from A->B is a function of:
        #  - overall mean likes (a)
        #  - how many likes A likes to give away in general (likes_gr[A, 1])
        #  - how many likes B tends to receive in general (likes_gr[B, 2])
        #  - the specific number of likes from A to B (d[dyad_id, 1])
        phiAB_l <- a + likes_gr[A, 1] + likes_gr[B, 2] + dyl[dyad_id, 1],
        phiBA_l <- a + likes_gr[B, 1] + likes_gr[A, 2] + dyl[dyad_id, 2],
        phiAB_d <- b + dislikes_gr[A, 1] + dislikes_gr[B, 2] + dyd[dyad_id, 1],
        phiBA_d <- b + dislikes_gr[B, 1] + dislikes_gr[A, 2] + dyd[dyad_id, 2],
        a ~ normal(0, 1),
        b ~ normal(0, 1),
        kappa_l ~ normal(0, 1.5),
        kappa_d ~ normal(0, 1.5),
        
        # likes matrix of varying effects. first col is how many likes monk A
        # gives away, second column is how many likes B receives
        transpars> matrix[18, 2]:likes_gr <-
            compose_noncentered(sigma_l, L_Rho_l, z_l),
        matrix[2, 18]:z_l ~ normal(0, 1),
        # Relatively broad prior
        cholesky_factor_corr[2]:L_Rho_l ~ lkj_corr_cholesky(4),
        vector[2]:sigma_l ~ exponential(1),
        
        # Dislikes giving receiving
        transpars> matrix[18, 2]:dislikes_gr <-
            compose_noncentered(sigma_dl, L_Rho_dl, z_dl),
        matrix[2, 18]:z_dl ~ normal(0, 1),
        cholesky_factor_corr[2]:L_Rho_dl ~ lkj_corr_cholesky(4),
        vector[2]:sigma_dl ~ exponential(1),
        
        # dyad effects - non-centered
        # Have 153 dyads
        transpars> matrix[153, 2]:dyl <-
            compose_noncentered(rep_vector(sigma_dyl, 2), L_Rho_dyl, z_dyl),
        matrix[2, 153]:z_dyl ~ normal(0, 1),
        # Relatively broad prior
        cholesky_factor_corr[2]:L_Rho_dyl ~ lkj_corr_cholesky(4),
        sigma_dyl ~ exponential(1),
        
        # Dyads for dislikes
        transpars> matrix[153, 2]:dyd <-
            compose_noncentered(rep_vector(sigma_dyd, 2), L_Rho_dyd, z_dyd),
        matrix[2, 153]:z_dyd ~ normal(0, 1),
        cholesky_factor_corr[2]:L_Rho_dyd ~ lkj_corr_cholesky(4),
        sigma_dyd ~ exponential(1),
        
        # compute centered correlation matrix for dyad effects
        gq> matrix[2,2]:Rho_dyl <<- Chol_to_Corr(L_Rho_dyl),
        gq> matrix[2,2]:Rho_dyd <<- Chol_to_Corr(L_Rho_dyd),
        gq> matrix[2,2]:Rho_l <<- Chol_to_Corr(L_Rho_l),
        gq> matrix[2,2]:Rho_dl <<- Chol_to_Corr(L_Rho_dl)
    ),
    # NB: for ordered logistic the outcome must be positive, think factor levels
    data=Monks |> mutate(
        like_AB = like_AB + 1,
        like_BA = like_BA + 1,
        dislike_AB = dislike_AB + 1,
        dislike_BA = dislike_BA + 1,
    ),
    chains=4, cores=4, iter=2000
)

# Is there any evidence of reciprocity?
# Again I don't understand precisely what this is asking for, whether it's the general
# or dyad-specific effects we're after
# There is again little correlation between someone who gives likes and someone who receives them (Rho_l)
# Although there is much stronger correlation between someone who gives out dislikes and receives them!
precis(m_q2_a, pars=c("Rho_l", "Rho_dl"), depth=3)
# Can get the posterior contrast
post <- extract.samples(m_q2_a)
dislike_minus_like <- post$Rho_dl[, 2, 1] - post$Rho_l[, 2, 1]
# Overall there is a slight trend to more reciprocity among dislikes than likes, although there is 
# a wide uncertainty
dens(dislike_minus_like)
precis(dislike_minus_like)

# If the question is asking about reciprocity amongst dyads then we get the opposite interpretation!
# There is a weak negative slant, i.e. more reciprocity amongst likes in dyads
precis(m_q2_a, pars=c("Rho_dyl", "Rho_dyd"), depth=3)
dislike_minus_like_dyad <- post$Rho_dyd[, 2, 1] - post$Rho_dyl[, 2, 1]
dens(dislike_minus_like_dyad)
precis(dislike_minus_like_dyad)

# Question 3 --------------------------------------------------------------
# I've already added general like/dislike varying effects
# I mentioned I was going to try 2 more models:
#   - A 4-column matrix with giving out likes, giving out dislikes, receiving likes, receiving dislikes
#   - A 3-column matrix with giving out likes, giving out dislikes, receiving likes
# I'll fit these extra models now and see if they agree on the most liked/disliked individuals
m_q2_b <- ulam(
    alist(
        # I think can use the same kappa here since it's symmetric
        like_AB ~ ordered_logistic(phiAB_l, kappa_l),
        like_BA ~ ordered_logistic(phiBA_l, kappa_l),
        dislike_AB ~ ordered_logistic(phiAB_d, kappa_d),
        dislike_BA ~ ordered_logistic(phiBA_d, kappa_d),
        # So the number of likes from A->B is a function of:
        #  - overall mean likes (a)
        #  - how many likes A likes to give away in general (likes_gr[A, 1])
        #  - how many likes B tends to receive in general (likes_gr[B, 2])
        #  - the specific number of likes from A to B (d[dyad_id, 1])
        phiAB_l <- a + likes_gr[A, 1] + likes_gr[B, 2] + dyl[dyad_id, 1],
        phiBA_l <- a + likes_gr[B, 1] + likes_gr[A, 2] + dyl[dyad_id, 2],
        phiAB_d <- b + likes_gr[A, 3] + likes_gr[B, 4] + dyd[dyad_id, 1],
        phiBA_d <- b + likes_gr[B, 3] + likes_gr[A, 4] + dyd[dyad_id, 2],
        a ~ normal(0, 1),
        b ~ normal(0, 1),
        kappa_l ~ normal(0, 1.5),
        kappa_d ~ normal(0, 1.5),
        
        # likes matrix of varying effects. 
        # first col is how many likes monk A gives away
        # second column is how many likes B receives
        # third col is how many dislikes A gives away
        # third col is how many dislikes A receives
        transpars> matrix[18, 4]:likes_gr <-
            compose_noncentered(sigma_l, L_Rho_l, z_l),
        matrix[4, 18]:z_l ~ normal(0, 1),
        # Relatively broad prior
        cholesky_factor_corr[4]:L_Rho_l ~ lkj_corr_cholesky(4),
        vector[4]:sigma_l ~ exponential(1),
        
        
        # dyad effects - non-centered
        # Have 153 dyads
        transpars> matrix[153, 2]:dyl <-
            compose_noncentered(rep_vector(sigma_dyl, 2), L_Rho_dyl, z_dyl),
        matrix[2, 153]:z_dyl ~ normal(0, 1),
        # Relatively broad prior
        cholesky_factor_corr[2]:L_Rho_dyl ~ lkj_corr_cholesky(4),
        sigma_dyl ~ exponential(1),
        
        # Dyads for dislikes
        transpars> matrix[153, 2]:dyd <-
            compose_noncentered(rep_vector(sigma_dyd, 2), L_Rho_dyd, z_dyd),
        matrix[2, 153]:z_dyd ~ normal(0, 1),
        cholesky_factor_corr[2]:L_Rho_dyd ~ lkj_corr_cholesky(4),
        sigma_dyd ~ exponential(1),
        
        # compute centered correlation matrix for dyad effects
        gq> matrix[2,2]:Rho_dyl <<- Chol_to_Corr(L_Rho_dyl),
        gq> matrix[2,2]:Rho_dyd <<- Chol_to_Corr(L_Rho_dyd),
        gq> matrix[4,4]:Rho_l <<- Chol_to_Corr(L_Rho_l)
    ),
    # NB: for ordered logistic the outcome must be positive, think factor levels
    data=Monks |> mutate(
        like_AB = like_AB + 1,
        like_BA = like_BA + 1,
        dislike_AB = dislike_AB + 1,
        dislike_BA = dislike_BA + 1,
    ),
    chains=4, cores=4, iter=2000
)
m_q2_c <- ulam(
    alist(
        # I think can use the same kappa here since it's symmetric
        like_AB ~ ordered_logistic(phiAB_l, kappa_l),
        like_BA ~ ordered_logistic(phiBA_l, kappa_l),
        dislike_AB ~ ordered_logistic(phiAB_d, kappa_d),
        dislike_BA ~ ordered_logistic(phiBA_d, kappa_d),
        # So the number of likes from A->B is a function of:
        #  - overall mean likes (a)
        #  - how many likes A likes to give away in general (likes_gr[A, 1])
        #  - how many likes B tends to receive in general (likes_gr[B, 2])
        #  - the specific number of likes from A to B (d[dyad_id, 1])
        phiAB_l <- a + likes_gr[A, 1] + likes_gr[B, 2] + dyl[dyad_id, 1],
        phiBA_l <- a + likes_gr[B, 1] + likes_gr[A, 2] + dyl[dyad_id, 2],
        phiAB_d <- b + likes_gr[A, 3] - likes_gr[B, 2] + dyd[dyad_id, 1],
        phiBA_d <- b + likes_gr[B, 3] - likes_gr[A, 2] + dyd[dyad_id, 2],
        a ~ normal(0, 1),
        b ~ normal(0, 1),
        kappa_l ~ normal(0, 1.5),
        kappa_d ~ normal(0, 1.5),
        
        # likes matrix of varying effects. 
        # first col is how many likes monk A gives away
        # second column is how many likes B receives
        # third col is how many dislikes A gives away
        transpars> matrix[18, 3]:likes_gr <-
            compose_noncentered(sigma_l, L_Rho_l, z_l),
        matrix[3, 18]:z_l ~ normal(0, 1),
        # Relatively broad prior
        cholesky_factor_corr[3]:L_Rho_l ~ lkj_corr_cholesky(4),
        vector[3]:sigma_l ~ exponential(1),
        
        
        # dyad effects - non-centered
        # Have 153 dyads
        transpars> matrix[153, 2]:dyl <-
            compose_noncentered(rep_vector(sigma_dyl, 2), L_Rho_dyl, z_dyl),
        matrix[2, 153]:z_dyl ~ normal(0, 1),
        # Relatively broad prior
        cholesky_factor_corr[2]:L_Rho_dyl ~ lkj_corr_cholesky(4),
        sigma_dyl ~ exponential(1),
        
        # Dyads for dislikes
        transpars> matrix[153, 2]:dyd <-
            compose_noncentered(rep_vector(sigma_dyd, 2), L_Rho_dyd, z_dyd),
        matrix[2, 153]:z_dyd ~ normal(0, 1),
        cholesky_factor_corr[2]:L_Rho_dyd ~ lkj_corr_cholesky(4),
        sigma_dyd ~ exponential(1),
        
        # compute centered correlation matrix for dyad effects
        gq> matrix[2,2]:Rho_dyl <<- Chol_to_Corr(L_Rho_dyl),
        gq> matrix[2,2]:Rho_dyd <<- Chol_to_Corr(L_Rho_dyd),
        gq> matrix[3,3]:Rho_l <<- Chol_to_Corr(L_Rho_l)
    ),
    # NB: for ordered logistic the outcome must be positive, think factor levels
    data=Monks |> mutate(
        like_AB = like_AB + 1,
        like_BA = like_BA + 1,
        dislike_AB = dislike_AB + 1,
        dislike_BA = dislike_BA + 1,
    ),
    chains=4, cores=4, iter=2000
)
# Look at Rho for each
# So in the original (separate) model, there is no correlation between giving out
# likes and receiving likes, but there is a weak positive correlation between giving out dislikes and
# receiving them
precis(m_q2_a, pars=c("Rho_l", "Rho_dl"), depth=3)
# In the second model, there is still little correlation with giving out likes and receiving likes (1,2)
# There is little correlation between giving out likes and giving out dislikes (1,3)
# There is still a weak positive correlation between giving out dislikes and receiving dislikes (3, 4)
# Little correlation between giving out likes and receiving dislikes (1,4)
# Extremely weak negative correlation between receiving likes and giving out dislikes (2, 3)
# Also extremely weak negative correlation receiving likes and receiving dislikes (2, 4) (would have thought would be stronger)
precis(m_q2_b, pars=c("Rho_l"), depth=3)

# In the 3-way model the only correlation is a negative one between
# receiving likes and giving out dislikes (2,3), which is stronger than it was in the previous model
precis(m_q2_c, pars=c("Rho_l"), depth=3)

# Let's rank the monks by their likeability/dislikability using each of these 3 models
models <- list(
    "Separate"=m_q2_a,
    "Four-way"=m_q2_b,
    "Three-way"=m_q2_c
)

# Extract from posterior
res <- map_dfr(models, function(mod) {
    post <- extract.samples(mod)
    likeability <- post$likes_gr[, , 2]
    if (dim(post$likes_gr)[3] == 2) {
        # The 2 separate models
        dislikeability <- post$dislikes_gr[, , 2]
    } else if (dim(post$likes_gr)[3] == 3) {
        # Model with 3 cols, dislikeability is negative likeability
        dislikeability <- -post$likes_gr[, , 2]
    } else if (dim(post$likes_gr)[3] == 4) {
        # Model with 4 cols, have separate dislikeability
        dislikeability <- post$likes_gr[, , 4]
    }
    tibble(
        monk=1:18,
        likeability_mean=colMeans(likeability),
        likeability_lower=apply(likeability, 2, PI)[1, ],
        likeability_upper=apply(likeability, 2, PI)[2, ],
        dislikeability_mean=colMeans(dislikeability),
        dislikeability_lower=apply(dislikeability, 2, PI)[1, ],
        dislikeability_upper=apply(dislikeability, 2, PI)[2, ]
    )
}, .id="model")

# The models tend to agree a bit...
res |>
    pivot_longer(-c(model, monk), names_pattern="(.+)_(.+)", names_to=c("type", "statistic")) |>
    pivot_wider(names_from=statistic, values_from=value) |>
    mutate(
        monk=as.factor(monk),
        type=factor(type, levels=c("likeability", "dislikeability"))
    ) |>
    ggplot(aes(x=monk, y=mean, colour=model)) +
        geom_point() +
        geom_errorbar(aes(ymin=lower, ymax=upper)) +
        facet_grid(type~model, scales="free") +
        theme_classic() +
        theme(
            legend.position = "bottom"
        ) +
        scale_colour_brewer(palette="Dark2")

# On Likeability they all agree that Monk 2 is the most likeable
res |>
    group_by(model) |>
    top_n(1, likeability_mean)
# And 10 has the lowest 'likeability', which means least likely to receive likes, rather than least
# likely to receive dislikes
res |>
    group_by(model) |>
    top_n(1, -likeability_mean)
# The separate models agree on 5 being the most dislikeable, while the three-way model has monk 10
# (by design, as it's saying that least likely to receive likes is the same as most likely to receive dislike)
# Although 5 is also a highly disliked monk by the three-way model
res |>
    group_by(model) |>
    top_n(1, dislikeability_mean)

# The separate model is probably better here (IMO as a four-way with correlated features)
# since just because you don't receive a lot of likes, doesn't also mean you're likely to receive dislikes!


# Question 4 --------------------------------------------------------------
# Adding factions!
# So this will be another varying-effects term with a giving and receiving column for each of the
# 4 factions
# I'll use my four-way model for this
# Firstly add factions to the dataset
Monk_factions <- tibble(
    MonkID=1:18,
    Faction=c(1, 1, 3, 2, 2, 2, 1, 4, 2, 4, 2, 1, 4, 1, 1, 1, 3, 3)
)
Monks <- Monks |>
    inner_join(Monk_factions |> rename(faction_A=Faction), by=c("A"="MonkID")) |>
    inner_join(Monk_factions |> rename(faction_B=Faction), by=c("B"="MonkID"))

# Firstly combine dyads into a single matrix for individuals
# This will allow correlations between if a dyad are more likely to give likes then they are
# less likely to give dislikes, which would expect
m_q4_a <- ulam(
    alist(
        # I think can use the same kappa here since it's symmetric
        like_AB ~ ordered_logistic(phiAB_l, kappa_l),
        like_BA ~ ordered_logistic(phiBA_l, kappa_l),
        dislike_AB ~ ordered_logistic(phiAB_d, kappa_d),
        dislike_BA ~ ordered_logistic(phiBA_d, kappa_d),
        # So the number of likes from A->B is a function of:
        #  - overall mean likes (a)
        #  - how many likes A likes to give away in general (likes_gr[A, 1])
        #  - how many likes B tends to receive in general (likes_gr[B, 2])
        #  - the specific number of likes from A to B (d[dyad_id, 1])
        phiAB_l <- a + likes_gr[A, 1] + likes_gr[B, 2] +  dyi[dyad_id, 1],
        phiBA_l <- a + likes_gr[B, 1] + likes_gr[A, 2] +  dyi[dyad_id, 2],
        phiAB_d <- b + likes_gr[A, 3] + likes_gr[B, 4] +  dyi[dyad_id, 3],
        phiBA_d <- b + likes_gr[B, 3] + likes_gr[A, 4] +  dyi[dyad_id, 4],
        a ~ normal(0, 1),
        b ~ normal(0, 1),
        kappa_l ~ normal(0, 1.5),
        kappa_d ~ normal(0, 1.5),
        
        # likes matrix of varying effects. 
        # first col is how many likes monk A gives away
        # second column is how many likes B receives
        # third col is how many dislikes A gives away
        # third col is how many dislikes A receives
        transpars> matrix[18, 4]:likes_gr <-
            compose_noncentered(sigma_l, L_Rho_l, z_l),
        matrix[4, 18]:z_l ~ normal(0, 1),
        # Relatively broad prior
        cholesky_factor_corr[4]:L_Rho_l ~ lkj_corr_cholesky(4),
        vector[4]:sigma_l ~ exponential(1),
        
        # dyad effects - non-centered
        # Have 153 dyads
        transpars> matrix[153, 4]:dyi <-
            compose_noncentered(rep_vector(sigma_dyi, 4), L_Rho_dyi, z_dyi),
        matrix[4, 153]:z_dyi ~ normal(0, 1),
        # Relatively broad prior
        cholesky_factor_corr[4]:L_Rho_dyi ~ lkj_corr_cholesky(4),
        sigma_dyi ~ exponential(1),
        
        # compute centered correlation matrix for dyad effects
        gq> matrix[4,4]:Rho_dyi <<- Chol_to_Corr(L_Rho_dyi),
        gq> matrix[4,4]:Rho_l <<- Chol_to_Corr(L_Rho_l)
    ),
    # NB: for ordered logistic the outcome must be positive, think factor levels
    data=Monks |> mutate(
        like_AB = like_AB + 1,
        like_BA = like_BA + 1,
        dislike_AB = dislike_AB + 1,
        dislike_BA = dislike_BA + 1,
    ),
    chains=4, cores=4, iter=2000
)
# Strongest correlation is still the reciprocity, i.e. if A gives B likes then B gives likes to A
# If A gives likes to B (col 1), then this somehow has a positive correlation with A giving dislikes to B (col 3), although the CI covers 0
# The only correlations that don't cross 0 are the reciprocity of likes (1-2) and dislikes (3-4)
plot(precis(m_q4_a, pars=c("Rho_dyi"), depth=3))

# Now add in the faction dyad
# Since need to have at least a 2D matrix (one dimension for each individual's faction), I could either add
# 4, 4x4 matrices, or 1 4x4 matrix. If I had 4 it would be the usual effect of A giving B likes, receiving likes, B giving A likes & receiving likes
# If I use a single 4x4 matrix it's a simplified way of saying "how good is the relationship between B & A" and assumes that the effects are symmetrical
# (i.e. A -> B has the same effect of B -> A), and also mirrored (dislikes are equivalent to negative likes)
# I'll write it as 2 2Ds, so fac_likes[fac_a, fac_b] is effect of A giving B likes, while fac_likes[fac_b, fac_a] is B giving A likes
# And then I'll have a separate matrix for dislikes
# This means there's no mechanism for receiving effects, but that's ok as that's basically the same
m_q4_b <- ulam(
    alist(
        # I think can use the same kappa here since it's symmetric
        like_AB ~ ordered_logistic(phiAB_l, kappa_l),
        like_BA ~ ordered_logistic(phiBA_l, kappa_l),
        dislike_AB ~ ordered_logistic(phiAB_d, kappa_d),
        dislike_BA ~ ordered_logistic(phiBA_d, kappa_d),
        # So the number of likes from A->B is a function of:
        #  - overall mean likes (a)
        #  - how many likes A likes to give away in general (likes_gr[A, 1])
        #  - how many likes B tends to receive in general (likes_gr[B, 2])
        #  - the specific number of likes from A to B (d[dyad_id, 1])
        phiAB_l <- a + likes_gr[A, 1] + likes_gr[B, 2] +  dyi[dyad_id, 1] + facs_likes[faction_A, faction_B],
        phiBA_l <- a + likes_gr[B, 1] + likes_gr[A, 2] +  dyi[dyad_id, 2] + facs_likes[faction_B, faction_A],
        phiAB_d <- b + likes_gr[A, 3] + likes_gr[B, 4] +  dyi[dyad_id, 3] + facs_dislikes[faction_A, faction_B],
        phiBA_d <- b + likes_gr[B, 3] + likes_gr[A, 4] +  dyi[dyad_id, 4] + facs_dislikes[faction_B, faction_A],
        a ~ normal(0, 1),
        b ~ normal(0, 1),
        kappa_l ~ normal(0, 1.5),
        kappa_d ~ normal(0, 1.5),
        
        # likes matrix of varying effects. 
        # first col is how many likes monk A gives away
        # second column is how many likes B receives
        # third col is how many dislikes A gives away
        # third col is how many dislikes A receives
        transpars> matrix[18, 4]:likes_gr <-
            compose_noncentered(sigma_l, L_Rho_l, z_l),
        matrix[4, 18]:z_l ~ normal(0, 1),
        # Relatively broad prior
        cholesky_factor_corr[4]:L_Rho_l ~ lkj_corr_cholesky(4),
        vector[4]:sigma_l ~ exponential(1),
        
        # dyad effects - non-centered
        # Have 153 dyads
        transpars> matrix[153, 4]:dyi <-
            compose_noncentered(rep_vector(sigma_dyi, 4), L_Rho_dyi, z_dyi),
        matrix[4, 153]:z_dyi ~ normal(0, 1),
        # Relatively broad prior
        cholesky_factor_corr[4]:L_Rho_dyi ~ lkj_corr_cholesky(4),
        sigma_dyi ~ exponential(1),
        
        # faction matrices - non-centered
        transpars> matrix[4, 4]:facs_likes <-
            compose_noncentered(rep_vector(sigma_fl, 4), L_Rho_fl, z_fl),
        matrix[4, 4]:z_fl ~ normal(0, 1),
        cholesky_factor_corr[4]:L_Rho_fl ~ lkj_corr_cholesky(4),
        sigma_fl ~ exponential(1),
        # faction dislikes
        transpars> matrix[4, 4]:facs_dislikes <-
            compose_noncentered(rep_vector(sigma_fd, 4), L_Rho_fd, z_fd),
        matrix[4, 4]:z_fd ~ normal(0, 1),
        cholesky_factor_corr[4]:L_Rho_fd ~ lkj_corr_cholesky(4),
        sigma_fd ~ exponential(1),
        
        # compute centered correlation matrix for dyad effects
        gq> matrix[4,4]:Rho_dyi <<- Chol_to_Corr(L_Rho_dyi),
        gq> matrix[4,4]:Rho_fl <<- Chol_to_Corr(L_Rho_fl),
        gq> matrix[4,4]:Rho_fd <<- Chol_to_Corr(L_Rho_fd),
        gq> matrix[4,4]:Rho_l <<- Chol_to_Corr(L_Rho_l)
    ),
    # NB: for ordered logistic the outcome must be positive, think factor levels
    data=Monks |> mutate(
        like_AB = like_AB + 1,
        like_BA = like_BA + 1,
        dislike_AB = dislike_AB + 1,
        dislike_BA = dislike_BA + 1,
    ),
    chains=4, cores=4, iter=2000
)

# Look at the correlation between dyads (the like reciprocity in Rho_dyi[2, 1])
# This shrinks from 0.64 to 0.14 when adding factions, suggesting factions explain a lot of the variance in likes
precis(m_q4_a, pars=c("Rho_dyi"), depth=4)
precis(m_q4_b, pars=c("Rho_dyi"), depth=4)

# Here's the faction effects
# As expected, factions tend to give each other more likes ([1,1], [2,2], [3,3]), except for Faction 4 whose effect is zero.
# Every other faction is less likely to give Faction 1 likes ([1,1] down to [4,1])
# But this isn't necessarily the case for every faction, i.e. Faction 4 give a lot of likes to Faction 2, moreso than it gives itself!
# And while Faction 4 doesn't give likes to Faction 1, faction 1 is less likely to give Faction 4 likes and is pretty average for
# everyone else, so maybe there's a bit of a frosty relationship here
plot(precis(m_q4_b, pars=c("facs_likes"), depth=4))

plot(precis(m_q4_b, pars=c("facs_dislikes"), depth=4))

# Does that change the monk's likeability?
models <- list(
    "Separate"=m_q2_a,
    "Four-way"=m_q2_b,
    "Three-way"=m_q2_c,
    "Factions"=m_q4_b
)

# Extract from posterior
res <- map_dfr(models, function(mod) {
    post <- extract.samples(mod)
    likeability <- post$likes_gr[, , 2]
    if (dim(post$likes_gr)[3] == 2) {
        # The 2 separate models
        dislikeability <- post$dislikes_gr[, , 2]
    } else if (dim(post$likes_gr)[3] == 3) {
        # Model with 3 cols, dislikeability is negative likeability
        dislikeability <- -post$likes_gr[, , 2]
    } else if (dim(post$likes_gr)[3] == 4) {
        # Model with 4 cols, have separate dislikeability
        dislikeability <- post$likes_gr[, , 4]
    }
    tibble(
        monk=1:18,
        likeability_mean=colMeans(likeability),
        likeability_lower=apply(likeability, 2, PI)[1, ],
        likeability_upper=apply(likeability, 2, PI)[2, ],
        dislikeability_mean=colMeans(dislikeability),
        dislikeability_lower=apply(dislikeability, 2, PI)[1, ],
        dislikeability_upper=apply(dislikeability, 2, PI)[2, ]
    )
}, .id="model")

# Adding factions in doesn't change the order of the results
# Although it does decrease the likeability a bit as some of this likeability is now due
# to the faction and not the individual
res |>
    group_by(model) |>
    top_n(1, likeability_mean)
res |>
    group_by(model) |>
    top_n(1, -likeability_mean)
# Ditto the dislikeability
res |>
    group_by(model) |>
    top_n(1, dislikeability_mean)