library(rethinking)
library(tidyverse)
library(dagitty)


# Question 1 --------------------------------------------------------------
data("NWOGrants")
NWOGrants <- NWOGrants |> as_tibble()
NWOGrants 
# Have dept/discipline, gender, applications, and number of awards
# "What are the total and indirect causal effects of gender on grant awards?"
# Consider mediation path (aka pipe) through discipline
dag <- dagitty("dag { 
  awards <- gender -> discipline -> awards;
}")
coordinates(dag) <- list(x=c(gender=0, discipline=1, awards=2),
                         y=c(gender=1, discipline=0, awards=1))
drawdag(dag)
# For the total effect we don't need to condition on anything
adjustmentSets(dag, exposure="gender", outcome="awards", effect="total")

# Create index variables
NWOGrants <- NWOGrants |>
    mutate(
        gid = as.numeric(factor(gender, levels=c('m', 'f'))),
        did = as.numeric(discipline)
    )

# For prior, want to centre on 0.5, and +/- 2SD means 12%-88% prob, which sounds about right
# as I doubt the actual award probability will be outside of these
m1 <- ulam(
    alist(
        awards ~ dbinom(applications, p),
        logit(p) <- b[gid],
        b[gid] ~ dnorm(0, 1)
    ),
    data=NWOGrants, cores=4, chains=4
)
# Rhat=1 and decent neff
# Has given men a slightly higher log-odds than women overall
# But we need to evaluate the _constrasts_ to identify the difference
precis(m1, depth=2)

post <- extract.samples(m1)
# A 3% difference in all (in favour of men), which although the 89% CI is outside of 0, is not a huge
# effect size
total_effect <- inv_logit(post$b[, 1]) - inv_logit(post$b[, 2])
precis(total_effect)
dens(total_effect)

# Question 2 --------------------------------------------------------------
# For the direct effect we have a straight forward mediation analysis so we need
# to condition on discipline. Dagitty agrees
adjustmentSets(dag, exposure="gender", outcome="awards", effect="direct")

m2 <- ulam(
    alist(
        awards ~ dbinom(applications, p),
        logit(p) <- b[gid] + d[did],
        b[gid] ~ dnorm(0, 1),
        d[did] ~ dnorm(0, 1)
    ),
    data=NWOGrants, cores=4, chains=4
)
# neff a little low at ~230, rhat at 1.01 is ok
precis(m2, depth=2)
post <- extract.samples(m2)
direct_effect <- inv_logit(post$b[, 1]) - inv_logit(post$b[, 2])
precis(direct_effect)
# Little difference between total and direct effect, slightly reduced from 3% to 2% in favour of men
# But still showing a male advantage
dens(direct_effect)
# 2.5% male DIRECT advantage
mean(direct_effect)
# 2.8% male TOTAL advantage
mean(total_effect)

# HOWEVER! The question now asks for "average direct causal effect of gender,
# weighting each discipline in proportion to the number of applications in each sample".
# This _isn't_ covered in the book, only the video lecture 9.
# Basically rather than assume that the gender effect is the same in every dept, we should
# assume it can be different
m3 <- ulam(
    alist(
        awards ~ dbinom(applications, p),
        logit(p) <- a[gid, did],
        matrix[gid,did]:a ~ normal(0, 1)
    ),
    data=NWOGrants, cores=4, chains=4
)
# Much better rhat and neff than the far simpler m2! How is that the case?
precis(m3, depth=3)

# To calculate the average direct causal effect over gender we would simply
# take the average effect
# But the question asks for this to be weighted by each depts # applicants

# So we'll simulate that number of both male and female and see what the 
# overall difference in acceptance probability is
ndata_men <- NWOGrants |> 
    group_by(did) |> 
    summarise(applications=sum(applications)) |>
    mutate(gid=1)
ndata_women <- NWOGrants |> 
    group_by(did) |> 
    summarise(applications=sum(applications)) |>
    mutate(gid=2)
probs_men <- link(m3, data=ndata_men)
# NB: getting 1 column per row of dataset, i.e. 1 column per dept
dim(probs_men)
probs_women <- link(m3, data=ndata_women)
# Can directly subtract them, this will do elementwise subtraction so it will get
# the direct causal effect per dept
direct_effect_2 <- probs_men - probs_women
# And then we can just flatten that into a single distribution to marginalize it
dens(as.numeric(direct_effect_2))
# Giving us an effect of -0.01 on average, with 89% CI from -0.16 - 0.12, i.e. nothing
precis(as.numeric(direct_effect_2))
mean(as.numeric(direct_effect_2))
median(as.numeric(direct_effect_2))

# Marking -----------------------------------------------------------------
# NB: What we did here was called POST-STRATIFYING
# Where reweight counter-factuals according to the distribution
# of samples in the 'real' world. Think how polls work (fit model of voter
# preference by certain characteristics from small sample, then simulate 
# how the entire voting population would vote) using their known characteristics

# The difference between this model (A ~ G * D) and the one in the book (A ~ G + D) is based
# on how you model the functional relationship between Gender and Discipline. Are they independent (book)
# or interacting (homework)? They both answer the same question of whether there is a direct relationship, but
# they differ in their assumptions (does each dept have a different gender imbalance?). Would be interesting to
# see what would happen if had a partially pooled version of this!
# The reason for needing to do the post-stratifying, is unlike the book model where you get a single coefficient
# the direct causal effect estimate, here you get one estimate per discipline and so need to marginalize over
# discipline. The way this is done is using a WEIGHTED average, which makes sense as we have data that the applications
# per discipline aren't uniform.
# Using interactions more explicitly like this also means you can pull out more interesting observations, and provided
# you have sufficient data, can understand your problem much more in depth. Week 6's homework shows a use-case where
# an interesting relationship would be hidden if a full interaction wasn't used. This approach also scales well to
# partial pooling, i.e. could have a male effect across disciplines that each discipline is drawn from and ditto for
# female, and then can use THESE hyper-parameters as the direct effect estimates (I think... sounds right in my head!)

# Question 3 --------------------------------------------------------------
data("UFClefties")
# What is the average advantage of lefties over righties?
UFClefties <- as_tibble(UFClefties)
# Columns:
#   - fight: match id
#   - episode: episode id
#   - fight.in.episode: order of fight in episode
#   - fighter1.win: outcome
#   - fighter1: fighter id
#   - fighter2: fighter id
#   - fighter1: if fighter 1 was lefty
#   - fighter2: if fighter 2 was lefty
# I can't see that any of these factors could influence the outcome beyond the handedness of both fighters
UFClefties

# In which case we could have 4 inputs:
# L/L
# L/R
# R/L
# R/R
UFClefties <- UFClefties |>
    mutate(matchup = 2*fighter1.lefty + fighter2.lefty + 1)

m4 <- ulam(
    alist(
        winner ~ dbern(p),
        logit(p) <- a[matchup],
        a[matchup] ~ dnorm(0, 1)
    ), data=UFClefties |> select(winner=fighter1.win, matchup), cores=4, chains=4
)
# Decent rhat and neff
precis(m4, depth=2)

# None of these results look particularly strong, except for matchup 1, where no one
# is a lefty and the fighter 1 has a slightly higher chance of winning

# Let's average the causal effect over all proportions of fighters (as not conditioning on 
# anything like we were in the last model)
post <- extract.samples(m4)
# Interesting, when both are righties, fighter 1 has a slightly higher chance of winning
# than fighter 2 (V1), but when both are lefties fighter 2 has a slightly higher chance (V4)
# But when there is one leftie and one not( V2, V3), the prob is ~0.48 in both directions
# which doesn't seem like a leftie advantage!
precis(as_tibble(inv_logit(post$a)))

# Based on this, why do I think lefties are overrepresented in the UFC?
# Potentially just old wives' tales that they are better fighters, without looking at the actual data?

# Marking -----------------------------------------------------------------
# The solutions says that we "need a model that is invariant to which
# figher is listed as number 1 or 2", which isn't the case for my model as
# I explicitly differentiate between fighter 1 being a leftie and fighter 2 not
# and vice-versa
# TODO: Why is this required?
# So the straight forward solution is rather than having a categorical predictor
# with all 4 matchups, to just encode this as continuous so get a single 
# coefficient, e.g. matchup = fighter1_lefty - fighter2_lefty
# So if both are lefties or righties this is 0
# If F1 is leftie and F2 is rightie then this gives the log-odds of a F1 win,
# and if F1 is rightie and F2 is leftie then it's the -ve log-odds of a F1 win,
# i.e. the log-odds of an F2 win. This means the model is exchangeable to fighter
# labels.
# I.e. say matchup = F1_leftie - F2_leftie = M
# Then prob of F1 win is inv_logit(M)
# and prob of F2 win is inv_logit(-M) = 1 - inv_logit(M)
# Demonstrate this holds:
# 80% chance of F1 win if F1 leftie and F2 rightie
M <- 1.4
inv_logit(M)

# 20% chance if the opposite
inv_logit(-M)

inv_logit(-M) == (1 - inv_logit(M))
# Basically the same except for float precision
inv_logit(-M) - (1 - inv_logit(M))

m5 <- ulam(
    alist(
        winner ~ dbern(p),
        logit(p) <- b*matchup,
        b ~ dnorm(0, 1)
    ), data=UFClefties |> mutate(matchup=fighter1.lefty - fighter2.lefty) |>
        select(matchup, winner=fighter1.win), cores=4, chains=4
)
# Decent rhat and neff
# Effect is centered on 0!
# Ok rhat and neff, so looks like no leftie effect!
# The solutions discusses possible explanations for this common belief
precis(m5)

