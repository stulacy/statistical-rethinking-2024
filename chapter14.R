library(rethinking)
library(tidyverse)
library(patchwork)
library(ellipse)
library(MASS)
library(dagitty)
library(ggrepel)
library(ape)

# 14.1 Varying slopes by construction ------------------------------------------
a <- 3.5       # average morning wait time (intercept)
b <- -1        # average difference with afternoon wait time (slope)
sigma_a <- 1   # std dev of intercepts
sigma_b <- 0.5 # std dev of slopes
rho <- -0.7    # correlation between intercepts and slopes

# form multivariate Gaussian parameterisation
Mu <- c(a, b)
cov_ab <- sigma_a * sigma_b*rho
Sigma <- matrix(c(sigma_a**2, cov_ab, cov_ab, sigma_b**2), nrow=2)
Mu
Sigma
# Can also encode Sigma as correlation matrix & variance matrix
sigmas <- c(sigma_a, sigma_b)
Rho <- matrix(c(1, rho, rho, 1), nrow=2)
diag(sigmas)  # Converts the 2-length vector into 2x2 matrix with 0s on off-diagonals
Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)

# MASS provides multivariate normal
# Draw 20 cafes from this multivariate normal
N_cafes <- 20
set.seed(5)
vary_effects <- mvrnorm(N_cafes, Mu, Sigma)
a_cafe <- vary_effects[, 1]
b_cafe <- vary_effects[, 2]

plot(a_cafe, b_cafe, col=rangi2, xlab="intercepts", ylab="slopes")
# Overlay population distribution
for (l in c(0.1, 0.3, 0.5, 0.8, 0.99)) {
    # Ellipse takes covariance of multivariate norm as main arg, can also pass
    # centre, and which level for the CI
    lines(ellipse(Sigma, centre=Mu, level=l), col=col.alpha("black", 0.2))
}

# NB: correlation is just scaled covariance to be between [-1, 1]

# Visit each of these 20 cafes 10 times, 5x morning & 5x afternoon
set.seed(22)
N_visits <- 10
afternoon <- rep(0:1, N_visits*N_cafes/2)
cafe_id <- rep(1:N_cafes, each=N_visits)
mu <- a_cafe[cafe_id] + b_cafe[cafe_id] * afternoon
sigma <- 0.5  # stdev within cafes
wait <- rnorm(N_visits*N_cafes, mu, sigma)
d <- tibble(cafe=cafe_id, afternoon, wait)
d

# Have multiple clusters in data (cafes), and each cluster is viewed under different conditions
# (morning & afternoon) so can more accurately estimate slope (does this make a difference from
# the conditions that allow us to fit a varying intercept model? I doubt it, although slope is
# harder to fit)

# Now will fit a varying slopes model using covariance
# First, LKJ prior on the R matrix (correlations, so 1 on diags, and rho on off-diags)
R <- rlkjcorr(1e4, K=2, eta=2)
dim(R)
# NB: the off diags are the same
dens(R[, 1, 2], xlab="Correlation")
dens(R[, 2, 1], xlab="Correlation")
# While the on-diags are always 1
mean(R[, 1, 1] >0.99)
mean(R[, 2, 2] >0.99)

# Eta=1 is basically uniform
dens(rlkjcorr(1e4, K=2, eta=1)[, 1, 2])
# Greater than 1 becomes more peaked around 0 (no correlation), so higher eta are more regularizing
dens(rlkjcorr(1e4, K=2, eta=4)[, 1, 2])

set.seed(867530)
m14.1 <- ulam(
    alist(
        wait ~ normal(mu, sigma),
        mu <- a_cafe[cafe] + b_cafe[cafe] * afternoon,
        c(a_cafe, b_cafe)[cafe] ~ multi_normal(c(a, b), Rho, sigma_cafe),
        a ~ normal(5, 2),
        b ~ normal(-1, 0.5),
        sigma_cafe ~ exponential(1),
        sigma ~ exponential(1),
        Rho ~ lkj_corr(2)
    ),
    data=d, chains=4, cores=4
)

# sigma_cafe is vector of 2, Rho is a 'corr_matrix' of 2 dims
stancode(m14.1)

# Posterior of varying effects
# Model has learned the negative correlation, despite the relatively uniform prior
post <- extract.samples(m14.1)
dens(post$Rho[, 1, 2], xlim=c(-1, 1))
dens(rlkjcorr(1e4, K=2, eta=2)[, 1, 2], add=TRUE, lty=2)

# Compare the partially pooled modelled estimates to unpooled
# Unpooled first
a1 <- d |> group_by(cafe) |> filter(afternoon == 0) |> summarise(foo = mean(wait)) |> pull(foo)
b1 <- d |> group_by(cafe) |> filter(afternoon == 1) |> summarise(foo = mean(wait)) |> pull(foo)
b1 <- b1 - a1

# Partially pooled
post <- extract.samples(m14.1)
a2 <- apply(post$a_cafe, 2, mean)
b2 <- apply(post$b_cafe, 2, mean)

# Posterior mean bivariate Gaussian
Mu_est <- c(mean(post$a), mean(post$b))
rho_est <- mean(post$Rho[, 1, 2])
sa_est <- mean(post$sigma_cafe[, 1])
sb_est <- mean(post$sigma_cafe[, 2])
cov_ab <- sa_est * sb_est * rho_est
Sigma_est <- matrix(c(sa_est**2, cov_ab, cov_ab, sb_est**2), ncol=2)
# Is this the best way to show the variance of cafes?
# Rather than using the ellipse function at specified levels, could instead draw from the posterior of cafes
# (mvnorm using the full posterior for a, b, rho, sigma_cafe etc...) and underlay as heatmap
contours <- map_dfr(c(0.1, 0.3, 0.5, 0.8, 0.99), function(l) as_tibble(ellipse(Sigma_est, centre=Mu_est, level=l)), .id="level")

# More pooling the further away from the centre!
p_1 <- tibble(
    cafe=1:N_cafes,
    a_unpooled=a1,
    a_pooled=a2,
    b_unpooled=b1,
    b_pooled=b2
) |>
    pivot_longer(-cafe, names_pattern="(a|b)_([a-z]+)", names_to=c("parameter", "model")) |>
    pivot_wider(names_from=parameter, values_from=value) |>
    ggplot(aes(x=a, y=b)) +
        geom_path(aes(x=x, y=y, group=level), data=contours, alpha=0.5) +
        geom_point(aes(colour=model)) +
        geom_line(aes(group=cafe), alpha=0.5) +
        theme_classic() +
        theme(legend.position = "bottom") +
        scale_colour_brewer("", palette="Dark2")
p_1

# Plot same information on the outcome scale, i.e. waiting times
# To get the Sigma for waiting times it:
#  - simulate intercepts & slopes using the mean posterior values
#  - convert to afternoon wait
#  - take the cov to get sigma (why not apply this directly on our waiting time distributions?)
# Mu is a direct conversion from the original mu
# Again, couldn't just take the posterior of values from these distributions and underlay contours or so?
v <- mvrnorm(1e4, Mu_est, Sigma_est)
# Calculate afternoon wait
v[, 2] <- v[, 1] + v[, 2]
Sigma_est2 <- cov(v)
Mu_est2 <- Mu_est
Mu_est2[2] <- Mu_est[1] + Mu_est[2]
contours2 <- map_dfr(c(0.1, 0.3, 0.5, 0.8, 0.99), function(l) as_tibble(ellipse(Sigma_est2, centre=Mu_est2, level=l)), .id="level")

p_2 <- tibble(
    cafe=1:N_cafes,
    morning_unpooled = a1,
    afternoon_unpooled = a1 + b1,
    morning_pooled = a2,
    afternoon_pooled=a2 + b2
) |>
    pivot_longer(-cafe, names_pattern="(morning|afternoon)_([a-z]+)", names_to=c("parameter", "model")) |>
    pivot_wider(names_from=parameter, values_from=value) |>
    ggplot(aes(x=morning, y=afternoon)) +
        geom_path(aes(x=x, y=y, group=level), data=contours2, alpha=0.5) +
        geom_point(aes(colour=model)) +
        geom_line(aes(group=cafe), alpha=0.5) +
        theme_classic() +
        theme(legend.position = "bottom") +
        scale_colour_brewer("", palette="Dark2")
p_2

# Just trying my approach to contours
contour_raw <- map_dfr(1:nrow(post$a), function(s) {
    sigmas <- diag(post$sigma_cafe[s, ])
    covar <- sigmas %*% post$Rho[s, , ] %*% sigmas
    res <- rmvnorm(1, c(post$a[s, ], post$b[s, ]))
    tibble(a=res[1, 1], b=res[1, 2])
}, .id="sample")


p_3 <- tibble(
    cafe=1:N_cafes,
    a_unpooled=a1,
    a_pooled=a2,
    b_unpooled=b1,
    b_pooled=b2
) |>
    pivot_longer(-cafe, names_pattern="(a|b)_([a-z]+)", names_to=c("parameter", "model")) |>
    pivot_wider(names_from=parameter, values_from=value) |>
    ggplot(aes(x=a, y=b)) +
        geom_density_2d(aes(x=a, y=b), data=contour_raw) +
        geom_point(aes(colour=model)) +
        geom_line(aes(group=cafe), alpha=0.5) +
        theme_classic() +
        theme(legend.position = "bottom") +
        scale_colour_brewer("", palette="Dark2")
# Certainly not as clean as using the ellipse function directly, but that loses some of
# the uncertainty in the cafes by only using the mean cafe parameters
p_1 + p_3

# Can do the same for getting results on the original scale
contour_raw <- map_dfr(1:nrow(post$a), function(s) {
    sigmas <- diag(post$sigma_cafe[s, ])
    covar <- sigmas %*% post$Rho[s, , ] %*% sigmas
    res <- rmvnorm(1, c(post$a[s, ], post$b[s, ]))
    tibble(morning=res[1, 1], afternoon=morning + res[1, 2])
}, .id="sample")

p_4 <- tibble(
    cafe=1:N_cafes,
    morning_unpooled = a1,
    afternoon_unpooled = a1 + b1,
    morning_pooled = a2,
    afternoon_pooled=a2 + b2
) |>
    pivot_longer(-cafe, names_pattern="(morning|afternoon)_([a-z]+)", names_to=c("parameter", "model")) |>
    pivot_wider(names_from=parameter, values_from=value) |>
    ggplot(aes(x=morning, y=afternoon)) +
        geom_density_2d(data=contour_raw) +
        geom_point(aes(colour=model)) +
        geom_line(aes(group=cafe), alpha=0.5) +
        theme_classic() +
        theme(legend.position = "bottom") +
        scale_colour_brewer("", palette="Dark2")
# Again I think I prefer my version
p_2 + p_4

# There's also the filled version which is quite nice
tibble(
    cafe=1:N_cafes,
    morning_unpooled = a1,
    afternoon_unpooled = a1 + b1,
    morning_pooled = a2,
    afternoon_pooled=a2 + b2
) |>
    pivot_longer(-cafe, names_pattern="(morning|afternoon)_([a-z]+)", names_to=c("parameter", "model")) |>
    pivot_wider(names_from=parameter, values_from=value) |>
    ggplot(aes(x=morning, y=afternoon)) +
        geom_density2d_filled(data=contour_raw, alpha=0.7) +
        geom_point(aes(colour=model)) +
        geom_line(aes(group=cafe), alpha=0.5) +
        theme_classic() +
        theme(legend.position = "bottom") +
        scale_colour_brewer("", palette="Dark2")

# Back to the show!

# 14.2 Advanced varying slopes --------------------------------------------
data("chimpanzees")
chimpanzees <- as_tibble(chimpanzees) |>
    mutate(
        block_id = as.integer(block),
        treatment = 1L + prosoc_left + 2L*condition
    ) |>
    rename(L=pulled_left, tid=treatment)

set.seed(4387510)
m14.2 <- ulam(
    alist(
        L ~ dbinom(1, p),
        logit(p) <- g[tid] + alpha[actor,tid] + beta[block_id, tid],
        
        # adaptive priors (aka hyperparameter priors)
        vector[4]:alpha[actor] ~ multi_normal(0, Rho_actor, sigma_actor),
        vector[4]:beta[block_id] ~ multi_normal(0, Rho_block, sigma_block),
        
        # Fixed priors
        g[tid] ~ dnorm(0, 1),
        sigma_actor ~ dexp(1),
        Rho_actor ~ dlkjcorr(4),
        sigma_block ~ dexp(1),
        Rho_block ~ dlkjcorr(4)
    ),data=chimpanzees |> dplyr::select(L, tid, actor, block_id),
    cores=4, chains=4
)

# Divergences!

# So alpha is 7,4, but not as a matrix, but a 2D array where each actor 
# is distributed around the same Rho and Sigma prior
# I.e. there is a hyper-parameter prior of treatment effects across actors
# and each actor has its 4 treatment effects drawn from this. This is a multi-variate
# equivalent of the hyper-parameter priors we saw in last chapter.
# NB: I don't see how this is varying slopes yet, surely this is still an intercept?
# We also only got this multi-variate nature because we allowed an interaction turn
stancode(m14.2)

# Non-centered. FAR quicker to fit and no divergent warnings
set.seed(4387510)
m14.3 <- ulam(
    alist(
        L ~ dbinom(1, p),
        logit(p) <- g[tid] + alpha[actor,tid] + beta[block_id, tid],
        
        # adaptive priors (aka hyperparameter priors) - non-centered
        transpars> matrix[actor, 4]:alpha <-
            compose_noncentered(sigma_actor, L_Rho_actor, z_actor),
        transpars> matrix[block_id, 4]:beta <-
            compose_noncentered(sigma_block, L_Rho_block, z_block),
        matrix[4, actor]:z_actor ~ normal(0, 1),
        matrix[4, block_id]:z_block ~ normal(0, 1),
        
        # Fixed priors
        g[tid] ~ dnorm(0, 1),
        # TODO why is the vector declaration needed here now? sigma_actor & block
        # were length 4 vectors previously as well
        vector[4]:sigma_actor ~ dexp(1),
        cholesky_factor_corr[4]:L_Rho_actor ~ lkj_corr_cholesky(2),
        vector[4]:sigma_block ~ dexp(1),
        cholesky_factor_corr[4]:L_Rho_block ~ lkj_corr_cholesky(2),
        
         # generated quantities: ordinary correlation matrices
        gq> matrix[4, 4]:Rho_actor <<- Chol_to_Corr(L_Rho_actor),
        gq> matrix[4, 4]:Rho_block <<- Chol_to_Corr(L_Rho_block)
    ),data=chimpanzees |> dplyr::select(L, tid, actor, block_id),
    cores=4, chains=4, log_lik=TRUE
)

stancode(m14.3)

# Pretty bad neff!
# Why also doesn't it sample Rho diagonals?
precis(m14.2, depth=3)
# Far better Rhat, always at 1.0 for non-centered version
precis(m14.3, depth=3)

# Thoughts: is this varying intercepts or varying slopes?
# is it a fair assumption that intercepts for each actor have a multi-variate nature?
# When to use pooling and when to use multivariate? I.e. couldn't we equally model
# the actor-treatment interactions as each actor-treatment is drawn from a distribution of
# treatments unique to that actor? I guess that doesn't make sense though, as you'd expect
# some (all?) treatments to have a different distribution entirely, rather than from the same
# underlying one
# what are both DLKJ priors, and Cholesky factor priors?

# The 76 parameters are 27 effective parameters, pretty heavy regularization! 
WAIC(m14.3)

# These standard deviations are pretty small, especially for the block cluster
# This illustrates how much shrinkage there is
plot(precis(m14.3, depth=2, pars=c("sigma_actor", "sigma_block")))

# Plot posterior predictions
# This is the same code as model m11.4 from 3 chapters ago
# That model had separate intercept for each actor and separate slope for each treatment
# BUT NOT MULTI-LEVEL PRIORS! That would have helped with shrinkage
# Now have separate intercept for each treatment, and each actor~treatment pair has its
# own slope, BUT we have both multi-level priors on these, as well as correlated treatment
# effects within actors
pl <- chimpanzees |> 
        filter(block_id == 5) |>
        group_by(actor, tid) |>
        summarise(L=mean(L)) |>
        ungroup()

ndata <- expand_grid(actor=unique(chimpanzees$actor), 
                     tid=unique(chimpanzees$tid), 
                     block_id=5)
post <- link(m14.3, data=ndata)
ndata$propSim <- colMeans(post)

d_plt <- chimpanzees |>
    group_by(actor, tid) |>
    summarise(propActual = mean(L)) |>
    ungroup() |>
    inner_join(ndata, by=c("actor", "tid"))

d_plt |>
    pivot_longer(c(propActual, propSim)) |>
    mutate(hand = ifelse(tid %% 2 == 1, 'Right', 'Left'),
           partner = ifelse(tid > 2, "Present", "Absent")) |>
    ggplot(aes(x=partner, y=value, colour=hand)) +
        geom_point() +
        geom_line(aes(group=hand)) +
        facet_grid(name ~ actor) +
        geom_hline(yintercept=0.5, linetype="dashed") +
        theme_bw() +
        theme(legend.position="bottom") +
        labs(x="Partner Present?", y="Proportion pulled left lever")


# 14.3 Instruments and Causal Designs -------------------------------------

set.seed(73)
N <- 500
U_sim <- rnorm(N)
Q_sim <- sample(1:4, size=N, replace=TRUE)
E_sim <- rnorm(N, U_sim + Q_sim)
W_sim <- rnorm(N, U_sim + 0 * E_sim)
dat_sim <- list(
    W=standardize(W_sim),
    E=standardize(E_sim),
    Q=standardize(Q_sim)
)

m14.4 <- ulam(
    alist(
        W ~ dnorm(mu, sigma),
        mu <- aW + bEW*E,
        aW ~ dnorm(0, 0.2),
        bEW ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ),
    data=dat_sim, chains=4, cores=4
)
# Because of confounding unobserved U, we get an association between E and W
# that shouldn't exist according to our data generating model
precis(m14.4)

# Add Q as predictor
m14.5 <- ulam(
    alist(
        W ~ dnorm(mu, sigma),
        mu <- aW + bEW*E + bQW*Q,
        aW ~ dnorm(0, 0.2),
        bEW ~ dnorm(0, 0.5),
        bQW ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ),
    data=dat_sim, chains=4, cores=4
)
# Even worse! Now Q is associated with W too when we know from the DGP that it shouldn't be
# E's association has got stronger too
precis(m14.5)

# Instead model the DAG more fully, i.e. W ~ E + U and E ~ Q + U
# But remove U since it's unobserved (could model as latent?)
# Model as multi-variate
# NB: just because modelling something as multi-variate DOESN'T mean that we think there IS
# significant correlation between them
m14.6 <- ulam(
    alist(
        c(W,E) ~ multi_normal(c(muW, muE), Rho, Sigma),
        muW <- aW + bEW*E,
        muE <- aE + bQE*E,
        # These are the 'standard' linear regression priors for
        # standardized data. intercept allows for [-0.5, 0.5] for average predictors
        # Slope allows for slope to be in [-1, 1], which makes sense as max rise/run = 4/4 = 1
        c(aW, aE) ~ dnorm(0, 0.2),
        c(bEW, bQE) ~ dnorm(0, 0.5),  
        Rho ~ lkj_corr(2),
        Sigma ~ exponential(1)  # NB: Sigma is 2-dimensional vector
    ),
    data=dat_sim, cores=4, chains=4
)
# Now bEQ is around 0 as expected from the DGP
# Rho[1,2] (and Rho[2, 1] for that matter) show a correlation of 0.54, this is due to U
# NB: this correlation isn't in the raw data, BUT conditional on E (for W) and Q (for E)
precis(m14.6, depth=3)

# Dagitty can find IVs if they are present in a dag
dagIV <- dagitty("dag{ Q -> E <- U -> W <- E }")
drawdag(dagIV)
instrumentalVariables(dagIV, exposure="E", outcome="W")

# 14.4 Social relations as correlated varying effects ---------------------
data(KosterLeckie)
kl_data <- list(
    N=nrow(kl_dyads),
    N_households=max(kl_dyads$hidB),
    did=kl_dyads$did,
    hidA=kl_dyads$hidA,
    hidB=kl_dyads$hidB,
    giftsAB=kl_dyads$giftsAB,
    giftsBA=kl_dyads$giftsBA
)
# Could imagine doing this as a football / sport model! gr is effectively the attack/defense
# parameters per team, which are allowed to correlate, and the dyad effects allow for specific
# matchups to go certain ways. I wonder if this can be done in a SSM? I doubt it it because of the
# indexing required
m14.7 <- ulam(
    alist(
        giftsAB ~ poisson(lambdaAB),
        giftsBA ~ poisson(lambdaBA),
        log(lambdaAB) <- a + gr[hidA, 1] + gr[hidB, 2] + d[did, 1],
        log(lambdaBA) <- a + gr[hidB, 1] + gr[hidA, 2] + d[did, 2],
        a ~ normal(0, 1),
        
        # gr matrix of varying effects. first col is giving, second is receving
        vector[2]:gr[N_households] ~ multi_normal(0, Rho_gr, sigma_gr),
        Rho_gr ~ lkj_corr(4),  # Stronger prior than seen before
        sigma_gr ~ exponential(1),
        
        # dyad effects - non-centered
        transpars> matrix[N, 2]:d <-
            compose_noncentered(rep_vector(sigma_d, 2), L_Rho_d, z),
        matrix[2, N]:z ~ normal(0, 1),
        # Even stronger prior! Don't expect much correlation between dyads
        cholesky_factor_corr[2]:L_Rho_d ~ lkj_corr_cholesky(8),
        sigma_d ~ exponential(1),
        
        # compute centered correlation matrix for dyad effects
        gq> matrix[2,2]:Rho_d <<- Chol_to_Corr(L_Rho_d)
    ),
    data=kl_data,
    chains=4, cores=4, iter=2000
)
# Surprisingly quick sampling from 300 samples! I guess the non-centered dyad matrix makes
# this possible
precis(m14.7, depth=3, pars=c("Rho_gr", "sigma_gr"))

# Want to plot posterior predictive, i.e for each (ok a small subset of) household,
# how much does it give and receive?
post <- extract.samples(m14.7)
giving <- sapply(1:25, function(i) post$a + post$gr[, i, 1])
receiving <- sapply(1:25, function(i) post$a + post$gr[, i, 2])
# Get the means on outcome scale too
Eg_mu <- apply(exp(giving), 2, mean)
Er_mu <- apply(exp(receiving), 2, mean)
# Want to show uncertainty in both of these. Could display horizontal and vertical error
# bars, but that is messy and besides it would hide the correlated nature

# This shows one household, could do this for all 25?
plot(exp(giving[, 1]), exp(receiving[, 1]))
# Or can approximate using Gaussian and draw contour lines using the ellipse function
# NB: ideally would draw individual contours or heatmaps PER house using some kind
# of density estimation
plot(NULL, xlim=c(0, 8.6), ylim=c(0, 8.6), xlab="generalized giving", 
     ylab="generalized receiving", lwd=1.5)
abline(a=0, b=1, lty=2)
for (i in 1:25) {
    Sigma <- cov(cbind(giving[, i], receiving[, i]))
    Mu <- c(mean(giving[, i]), mean(receiving[, i]))
    for (l in 0.5) {
        el <- ellipse(Sigma, centre=Mu, level=l)
        lines(exp(el), col=col.alpha("black", 0.5))
    }
}
points(Eg_mu, Er_mu, pch=21, bg="white", lwd=1.5)

# Can _kind_ of emulate this by using KDE to draw contours, although would need to play
# with geom_contours settings to only show the median contour
# Doesn't look like it's readily possible, instead would be more explicit to determine
# contours ahead of time using MASS::kde2d, as here: https://stackoverflow.com/a/46010566/1020006
r_df <- as_tibble(receiving)
colnames(r_df) <- paste0('H', 1:ncol(r_df))
r_df$sample <- 1:nrow(r_df)
r_df <- r_df |> pivot_longer(-sample, names_pattern="H([0-9]+)", names_to="house",
                             values_to="receiving")

g_df <- as_tibble(giving)
colnames(g_df) <- paste0('H', 1:ncol(g_df))
g_df$sample <- 1:nrow(g_df)
g_df <- g_df |> pivot_longer(-sample, names_pattern="H([0-9]+)", names_to="house",
                             values_to="giving")
comb_df <- r_df |> inner_join(g_df, by=c("sample", "house"))
comb_df |>
    ggplot(aes(x=exp(receiving), y=exp(giving), group=house)) +
        geom_density_2d(bins=10)

# Visualising dyad effects
# Fair amount of variance amongst dyad effects
# Very correlated!
# Sigma is higher than sigma from g/r, so more variation amongst dyads than
# between households gr rates
# The positive correlation means that if one household in a dyad, so will the other
# So in the football analogy would expect a negative here, so that if one team scores lots
# would expect the other in the dyad to score fewer
precis(m14.7, depth=3, pars=c("Rho_d", "sigma_d"))

# Plotting means, so each point is a dyad of household A & household B
# Very balanced!
# Again in football would be interested to see any outlier matchups!
dy1 <- colMeans(post$d[, , 1])
dy2 <- colMeans(post$d[, , 2])
plot(dy1, dy2)

# 14.5 Continuous categories and the Gaussian process ---------------------
data(islandsDistMatrix)
Dmat <- islandsDistMatrix
colnames(Dmat) <- c("Ml", "Ti", "SC", "Ya", "Fi", "Tr", "Ch", "Mn", "To", "Ha")
round(Dmat, 1)

# Kernel function
# Linear is standard exponential decay
# Quadratic penalizes islands further apart
tibble(
    x = seq(0, 4, length.out=100)
) |>
    mutate(
        linear = exp(-1 * x),
        squared = exp(-1*x**2)
   ) |>
    pivot_longer(-x) |>
    ggplot(aes(x=x, y=value, colour=name)) +
        geom_line()

# Now fit model
data(Kline2)
d <- Kline2
d$society <- 1:10  # Index 

dat_list <- list(
    T = d$total_tools,
    P=d$population,
    society=d$society,
    Dmat=islandsDistMatrix
)
m14.8 <- ulam(
    alist(
        T ~ dpois(lambda),
        lambda <- (a*P^b/g)*exp(k[society]),  # Varying intercept k
        vector[10]:k ~ multi_normal(0, SIGMA), # But now K is multi-norm with covarmat from kernel function
        matrix[10, 10]:SIGMA <- cov_GPL2(Dmat, etasq, rhosq, 0.01),
        c(a,b,g) ~ dexp(1),
        etasq ~ dexp(2),
        rhosq ~ dexp(0.5)
    ), data=dat_list, chains=4, cores=4, iter=2000
)
# Good rhat and neff
precis(m14.8, depth=3)

# Now plotting posterior samples
post <- extract.samples(m14.8)
x_seq <- seq(0, 10, length.out=100)
pmcov <- sapply(x_seq, function(x) post$etasq*exp(-post$rhosq*x^2))
as_tibble(pmcov) |>
    mutate(sample=1:nrow(pmcov)) |>
    pivot_longer(-sample) |>
    inner_join(tibble(x=x_seq, vname=paste0("V", seq_along(x_seq))), by=c("name"="vname")) |>
    filter(sample <= 100) |>
    ggplot(aes(x=x, y=value)) +
        geom_line(aes(group=sample), alpha=0.3) +
        theme_classic() +
        labs(y="covariance", x="Distance between islands (km)") +
        geom_smooth()

# Want to look at the covariances themselves but hard to interpret as their effects
# are on the log-count scale. Instead will look at the corelations, although for this
# need to use the posterior median (or some other scalar)
K <- matrix(0, nrow=10, ncol=10)
for (i in 1:10) {
    for (j in 1:10) {
        K[i,j] <- median(post$etasq) * exp(-median(post$rhosq) * islandsDistMatrix[i, j]^2)
    }
}

# Now convert to correlation
# NB: cov2cor is in base stats package
Rho <- round(cov2cor(K), 2)
colnames(Rho) <- d$culture
rownames(Rho) <- d$culture
# Hawaii has 0 correlation with anywhere else as it's so far away
Rho_df <- Rho |>
    as_tibble() |>
    mutate(
        island1 = rownames(Rho)
    ) |>
    pivot_longer(-island1, names_to="island2", values_to="correlation") |>
    inner_join(d |> dplyr::select(culture, lat1=lat, lon1=lon2), by=c("island1"="culture")) |>
    inner_join(d |> dplyr::select(culture, lat2=lat, lon2=lon2), by=c("island2"="culture")) 

# Plot map of islands with their correlations
# scale point size to logpop
d |>
    mutate(
        psize = logpop / max(logpop),
        psize = exp(psize*1.5)-2
    ) |>
    ggplot(aes(x=lon2, y=lat)) +
        geom_point(aes(size=psize), colour="blue") +
        geom_text_repel(aes(label=culture)) +
        geom_segment(aes(x=lon1, xend=lon2, y=lat1, yend=lat2, linewidth=correlation),
                     data=Rho_df) +
        scale_linewidth_continuous(range=c(0, 0.4)) +
        guides(size="none", linewidth="none") +
        labs(x="lon", y="lat") +
        theme_classic()

# Plot alongside the population & tools relationship
# Calculate the relationship between population and tools, ignoring distance
# I.e. for an 'average' island intercept
logpop.seq <- seq(6, 14, length.out=30)
lambda <- sapply(logpop.seq, function(lp) post$a * (exp(lp)^post$b)/post$g)
vals_sim <- tibble(
    logpop = logpop.seq,
    lambda_mean = colMeans(lambda),
    lambda_lower = apply(lambda, 2, PI)[1, ],
    lambda_upper = apply(lambda, 2, PI)[2, ]
) 

Rho_df <- Rho |>
    as_tibble() |>
    mutate(
        island1 = rownames(Rho)
    ) |>
    pivot_longer(-island1, names_to="island2", values_to="correlation") |>
    inner_join(d |> dplyr::select(culture, pop1=logpop, tools1=total_tools), by=c("island1"="culture")) |>
    inner_join(d |> dplyr::select(culture, pop2=logpop, tools2=total_tools), by=c("island2"="culture")) 

# According to the book, whenever there are islands above or below the expected line,
# then the difference is likely accounted for by their covariance
# I.e. Tonga has more tools than expected for an island of its size, but it is connected to
# Fiji which is connected to the 3 Tikopa, Malekula, Santa Cruz so "Tonga's proximity to Fiji
# counteracts some of the tug Fiji's smaller neighbours exert on it. So the model seems to think Fiji
# would have fewer tools, if it weren't for Tonga"
# So maybe it's saying that Tonga is a bigger island with lots of tools, and Fiji being closer
# means it benefits from this, and if it werent close to Tonga, it would have fewer tools
# because it would be "sharing" them more with the 3 smaller islands?
# (I don't really follow this rationale).
# This makes it sound like it's a zero-sum game, rather than being closer to another island means
# have more tools in general
# The book also states that Tikopa, Malekula, and Santa Cruz have lower tools than expected,
# but only Malekula does - the other 2 are right on the expected line
d |>
    mutate(
        psize = logpop / max(logpop),
        psize = exp(psize*1.5)-2
    ) |>
    ggplot(aes(x=logpop, y=total_tools)) +
        geom_ribbon(aes(ymin=lambda_lower, y=lambda_upper, ymax=lambda_upper), data=vals_sim, alpha=0.3) +
        geom_line(aes(x=logpop, y=lambda_mean), data=vals_sim, colour="steelblue", linewidth=1) +
        geom_point(aes(size=psize), colour="blue") +
        geom_text_repel(aes(label=culture)) +
        geom_segment(aes(x=pop1, xend=pop2, y=tools1, yend=tools2, linewidth=correlation),
                     data=Rho_df) +
        ylim(0, 100) +
        scale_linewidth_continuous(range=c(0, 0.4)) +
        guides(size="none", linewidth="none") +
        labs(x="logpop", y="tools") +
        theme_classic()

# Non-centering Gaussian process
# NB: We were using a CENTERED parameterisation of the varying slopes k, which
# were defined as a function of sigma
# vector[10]:k ~ multi_normal(0, SIGMA)
# So to optimise fitting can use a Cholesky decomposition
m14.8nc <- ulam(
    alist(
        T ~ dpois(lambda),
        lambda <- (a*P^b/g)*exp(k[society]),  # Varying intercept k
        
        transpars> vector[10]: k <<- L_SIGMA*z,
        vector[10]:z ~ normal(0, 1),
        transpars> matrix[10,10]: L_SIGMA <<- cholesky_decompose(SIGMA),
        transpars> matrix[10, 10]: SIGMA <- cov_GPL2(Dmat, etasq, rhosq, 0.01),
        
        c(a,b,g) ~ dexp(1),
        etasq ~ dexp(2),
        rhosq ~ dexp(0.5)
    ), data=dat_list, chains=4, cores=4, iter=2000
)
precis(m14.8, depth=2)
# Much better n_eff! have other 2,000!
precis(m14.8nc, depth=2)

# 14.5.2 Phylogenetic distance --------------------------------------------
data("Primates301")
data("Primates301_nex")
# Plot phylogenetic distance between apes
plot(ladderize(Primates301_nex), type="fan", font=1, no.margin=TRUE, 
     label.offset=1, cex=0.5)
# Primates has 16 variables from 301 primate species
Primates30 |> dim()
# the nex version is a phylogenetic tree
Primates301_nex
# Class 'phylo'
class(Primates301_nex)

# Now standard linear model of brain size ~ group size + body mass
d <- as_tibble(Primates301)
dstan <- d |>
        mutate(name = as.character(name)) |>
        filter(!is.na(group_size), !is.na(body), !is.na(brain)) |>
        mutate(
            N_spp = n(),
            M = standardize(log(body)),
            B = standardize(log(brain)),
            G = standardize(log(group_size)),
            Imat = diag(n())
)

dat_list <- dstan |> dplyr::select(B, G, M, Imat) |> as.list()
dat_list$N_spp <- nrow(dstan)

# Why log all these variables?
# Ah yes, all quite skewed
dstan |>
    dplyr::select(body, brain, group_size, species=name) |>
    pivot_longer(-species) |>
    ggplot(aes(x=value)) +
        facet_wrap(~name, scales = "free") +
        geom_density() +
        theme_classic()

# NB: Won't let me calculate log_lik!
m14.9 <- ulam(
    alist(
        B ~ multi_normal(mu, SIGMA),
        mu <- a + bM*M + bG*G,
        matrix[N_spp, N_spp]: SIGMA <- Imat * sigma_sq,
        a ~ normal(0, 1),
        c(bM,bG) ~ normal(0, 0.5),
        sigma_sq ~ exponential(1)
    ), data=dat_list, chains=4, cores=4
)
# Just to check that this would give the same results as a non-multivariate model
precis(m14.9)

m14.9_uni <- ulam(
    alist(
        B ~ normal(mu, sigma),
        mu <- a + bM*M + bG*G,
        a ~ normal(0, 1),
        c(bM,bG) ~ normal(0, 0.5),
        sigma ~ exponential(1)
    ), data=dat_list, chains=4, cores=4
)
# Yep, the only difference is in sigma vs sigma_sq!
# But I imagine this is some subtley of how multivariates are parameterised in terms
# of sigma**2 rather than sigma
precis(m14.9_uni)

# Now want to replace identity matrix with something that provides covariance as a function
# of phylogenetic distance, which we can do using the ape package
# Firstly trim the tree to the species we're modelling (i.e. the 151 that had complete data)
tree_trimmed <- keep.tip(Primates301_nex, dstan$name)
plot(tree_trimmed)
# Now form a correlation matrix from this using Brownian motion, i.e. the covariance depends
# only on the phylogenetic distance
Rbm <- corBrownian(phy=tree_trimmed) 
# Can't seem to do anything with this
Rbm
# Until initialise it, at which point we get a covariance matrix
# Ignore the warning for now...
V <- vcv(Rbm)
dim(V)
# Can plot this covariance vs distance and see that it's linear
Dmat <- cophenetic(tree_trimmed)
plot(Dmat, V, xlab="phylognetic distance", ylab="covariance")

# Order the covariance matrix. As that warning said, it was using the order in the data frame
# but I think since we didn't provide one it went with a random order
dat_list$V <- V[dstan$name, dstan$name]
# Convert from covariance to correlation matrix (Is this right?, this just goes from 
# 0 to 1 rather than -1 to 1)
dat_list$R <- cov2cor(dat_list$V)
# Now replace the identity matrix with this and we're good to go
m14.10 <- ulam(
    alist(
        B ~ multi_normal(mu, SIGMA),
        mu <- a + bM*M + bG*G,
        matrix[N_spp, N_spp]: SIGMA <- R * sigma_sq,
        a ~ normal(0, 1),
        c(bM,bG) ~ normal(0, 0.5),
        sigma_sq ~ exponential(1)
    ), data=dat_list, chains=4, cores=4
)
# The group size parameter that was previously important has now shrunk to 0
# So brain size was co-clustering with group-size in the phylogenetic tree
precis(m14.10)

# Previously our kernel was basically just a linear function of phylogenetic distance
# But a scientifically more reasonable model is the exponential distance kernel (rather
# than the quadratic kernel we used for islands). This is also called the L1 kernel, where
# L2 is the quadratic.
# Just need a distance matrix!
# I don't understand this, but we DON'T use R here
# Instead want to use Dmat, which arose from the cophenetic pairwise distances
# But since the Brownian motion is basically linear with cophenetic distances
# (c.f. that scatter plot above), would it really make a difference if use R or Dmat here?

# I'll use Dmat as the book does at first. NB: they scale Dmat so that it's between 0-1,
# but I think this is just for ease of computation as usual, as the Dmat we used in the 
# islands example wasn't scaled
dat_list$Dmat <- Dmat[dstan$name, dstan$name] / max(Dmat)
m14.11 <- ulam(
    alist(
        B ~ multi_normal(mu, SIGMA),
        mu <- a + bM*M + bG*G,
        matrix[N_spp, N_spp]: SIGMA <- cov_GPL1(Dmat, etasq, rhosq, 0.01),
        a ~ normal(0, 1),
        c(bM,bG) ~ normal(0, 0.5),
        etasq ~ half_normal(1, 0.25),
        rhosq ~ half_normal(3, 0.25)
    ), data=dat_list, chains=4, cores=4
)
# Group size is weakly associated with brain size again!
# How can we explain this just based on a different kernel?
# Maybe because we're losing the effect that distance species have on each other
precis(m14.11)

# Also just want to see what happens if we use R as the distance matrix
# Rather than Dmat
# Huh, that fails completely..
# Ah of course! R is correlation, D is distance.
# We want distance!
# So could try doing 1 - R
#ulam(
#    alist(
#        B ~ multi_normal(mu, SIGMA),
#        mu <- a + bM*M + bG*G,
#        matrix[N_spp, N_spp]: SIGMA <- cov_GPL1(R, etasq, rhosq, 0.01),
#        a ~ normal(0, 1),
#        c(bM,bG) ~ normal(0, 0.5),
#        etasq ~ half_normal(1, 0.25),
#        rhosq ~ half_normal(3, 0.25)
#    ), data=dat_list, chains=4, cores=4
#) |> precis()

# Yep that basically gives the same results so didn't need to calculate Dmat
# Lesson learned: first argument to cov_GPLx is the DISTANCE matrix not similarity
dat_list$R_dist <- 1 - dat_list$R
ulam(
    alist(
        B ~ multi_normal(mu, SIGMA),
        mu <- a + bM*M + bG*G,
        matrix[N_spp, N_spp]: SIGMA <- cov_GPL1(R_dist, etasq, rhosq, 0.01),
        a ~ normal(0, 1),
        c(bM,bG) ~ normal(0, 0.5),
        etasq ~ half_normal(1, 0.25),
        rhosq ~ half_normal(3, 0.25)
    ), data=dat_list, chains=4, cores=4
) |> precis()

# Plot predictions from this kernel
# Basically the prior allowed for quite a high effect from close neighbours
# But the posterior has drastically shrunk this effect
# The Brownian motion didn't do this because we weren't using a Gaussian process
# i.e. we just gave it a correlation matrix and allowed it to find sigma
# But here we have given it a prior of a kernel function and allowed it to find
# one that actually fits the data, so this _should_ be a better fit...
# And therefore it does seem like group size is associated with brain size
post <- extract.samples(m14.11)
x_seq <- seq(0, 1, length.out=100)
post_draws <- sapply(x_seq, function(x) post$etasq*exp(-post$rhosq*x))
eta_draws <- abs(rnorm(1000, 1, 0.25))
rho_draws <- abs(rnorm(1000, 3, 0.25))
prior_draws <- sapply(x_seq, function(x) eta_draws * exp(-rho_draws * x))
tibble(
    x=x_seq,
    post_mean = colMeans(post_draws),
    post_lower = apply(post_draws, 2, PI)[1, ],
    post_upper = apply(post_draws, 2, PI)[2, ],
    prior_mean = colMeans(prior_draws),
    prior_lower = apply(prior_draws, 2, PI)[1, ],
    prior_upper = apply(prior_draws, 2, PI)[2, ]
) |>
    pivot_longer(-x, names_pattern="(.+)_(.+)",
                                     names_to=c("dist", "stat")) |>
    pivot_wider(names_from=stat, values_from=value) |>
    ggplot(aes(x=x, y=mean, colour=dist)) +
        geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
        geom_line(linewidth=1) +
        theme_classic() +
        theme(legend.position="bottom") +
        labs(y="covariance", x="Phylogenetic distance") +
        geom_smooth()

# Ask Richard -------------------------------------------------------------

# When to use multi-level prior and when to use multivariate?
# I.e. how to decide to use a multivariate prior on treatments, and THAT has a 
# multi-level prior across actors?
# Instead could do either each treatment has its own multi-level prior but without correlation?
# Or even a double-level prior, so that get treatment effects have their own dist per actor
# and actors have their own dist

