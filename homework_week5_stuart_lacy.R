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
