library(rethinking)
library(dagitty)
library(tidyverse)
library(ggpubr)

# Question 1 --------------------------------------------------------------

dag_1 <- dagitty("dag {
    A -> F -> G -> W <- F
                   }")
coordinates(dag_1) <- list(x=c(F=0, A=1, W=1, G=2),
                               y=c(F=1, A=0, W=2, G=1))
drawdag(dag_1)

# Total causal effect of A on F, only path is A->F
# No need to adjust for anything
adjustmentSets(dag_1, exposure = "A", outcome="F")

# Just a bivariate model is needed
data(foxes)
d <- foxes
d$F <- scale(d$avgfood)
d$A <- scale(d$area)

m1 <- quap(
    alist(
        F ~ dnorm(mu, sigma),
        mu <- a + bA * A,
        a ~ dnorm(0, 0.2),
        c(bA) ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
    ), 
    data=d
)
# A 1 std change in area is associated with a 0.88 std increase in average food
precis(m1)
# Or in measurement units, 0.928 increase in area is associated with a 0.174 increase in average food
sd(d$area)
0.88 * sd(d$avgfood)

# Question 2 --------------------------------------------------------------
# Total causal effect of F -> W
# This acts through 2 paths: F -> W and F->G->W
# Fit a joint model with both paths that can simulate from
d$G <- scale(d$groupsize)
d$W <- scale(d$weight)
m2 <- quap(
    alist(
        ## F -> W <- G
        W ~ dnorm(mu, sigma),
        mu <- a + bF*F + bG*G,
        a ~ dnorm(0, 0.2),
        bF ~ dnorm(0, 0.5),
        bG ~ dnorm(0, 0.5),
        sigma ~ dexp(1),
        ## F -> G
        G ~ dnorm(mu_G, sigma_G),
        mu_G <- aG + bFG*F,
        aG ~ dnorm(0, 0.2),
        bFG ~ dnorm(0, 0.5),
        sigma_G ~ dexp(1)
    ),
    data=d
)

# Can now simulate from this model, by sweeping over F, simulating G, then simulating W
F_seq <- seq(from=min(d$F), to=max(d$F), length.out=50)
sim_dat <- data.frame(F = F_seq)
s <- sim(m2, data=sim_dat, vars=c("G", "W"))
pi <- apply(s$W, 2, PI)
pi_tbl <- tibble(ymin=pi[1, ], ymax=pi[2, ], F=sim_dat$F)

# Plot counter-factual:
# So 1 sd increase in F, corresponds to a 0.033 sd decrease in W
tibble(F=sim_dat$F, W=colMeans(s$W)) |>
    ggplot(aes(x=F)) +
        geom_ribbon(aes(ymin=ymin, ymax=ymax), data=pi_tbl,
                    alpha=0.4) +
        geom_line(aes(y=W)) +
        ylim(c(-2, 2)) +
        labs(x="Manipulated F", y="Counterfactual W") +
        stat_regline_equation(aes(y=W)) +
        theme_minimal()

# So a 0.198 increase in average food is associated with a 0.039 decrease in weight
# Even without knowing the measurement units, this seems a very weak association, and is in the
# oppposite direction to what you'd expect
sd(d$avgfood)
0.033 * sd(d$weight)

# Question 3 --------------------------------------------------------------
# Total causal effect of F on W
# Need to adjust for G as have a masked relationship
adjustmentSets(dag_1, exposure = "F", outcome="W", effect="direct")
# This is just bF from the first model
# So 0.48 +/- 0.18 units
precis(m2)

# The finding from these questions is that food does directly increase wolves weights
# after accounting for the change in group size, but the food size also directly
# causes an increase in wolf group size which then has a dentrimental effect
# on weight, stronger than the direct increase in weight from adjusting for group

# Question 4 --------------------------------------------------------------
dag_2 <- dagitty("dag {
    U[u]
    G <- U -> F
    A -> F -> G -> W <- F
                   }")
coordinates(dag_2) <- list(x=c(F=0, A=0, W=1, G=2, U=2),
                               y=c(F=1, A=0, W=2, G=1, U=0))
drawdag(dag_2)
adjustmentSets(dag_2, exposure = "F", outcome = "W", effect="direct")

# How can we simulate an unobserved variable?
# Could either use the original data and add in a random variable that isn't
# directly causing these
# Or could simulate an entire dataset using the DAG
# Or can use the DAG and account for it somehow?

# Now need to simulate a random variable U affecting F and G, and
# then can model the usual relationships
d$U <- rnorm(nrow(d))
m4 <- quap(
    alist(
        # U
        U ~ dnorm(0, 2),
        ## F -> W <- G
        W ~ dnorm(mu_W, sigma_W),
        mu_W <- aW + bWF*F + bWG*G,
        aW ~ dnorm(0, 0.2),
        bWF ~ dnorm(0, 0.5),
        bWG ~ dnorm(0, 0.5),
        sigma_W ~ dexp(1),
        ## F -> G <- U
        G ~ dnorm(mu_G, sigma_G),
        mu_G <- aG + bGF*F + bGU * U,
        aG ~ dnorm(0, 0.2),
        bGF ~ dnorm(0, 0.5),
        bGU ~ dnorm(0, 0.5),
        sigma_G ~ dexp(1),
        # U -> F
        F ~ dnorm(mu_F, sigma_F),
        mu_F <- U,
        #mu_F <- aF,
        #aF ~ dnorm(0, 0.2),
        sigma_F ~ dexp(1)
    ),
    data=d,
    start=list(U=rep(0.1, nrow(d)))
)

# Now simulate just providing U
# (Should I have simulated the entiredataset so that )