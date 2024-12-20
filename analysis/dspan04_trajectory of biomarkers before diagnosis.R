

# 1. Restrict to individuals with newly diagnosed T2D or no T2D. Set max_age as dmagediag if new T2D or max(age) if no T2D

# 2. Take their data retrospectively up to 15 years (i.e. max_age - age <= 15)
# Create a 't' variable: max_age - age

# 3. Include individuals who have at least one wave before T2D diagnosis (if newly diagnosed) or have at least two waves (no T2D)
# 3a. Calculate inverse probability weights for those individuals you had to drop based on their socio-demographic characteristics
# ipw = inverse of logit (Pr[included = 1 | X]) ~ sex + age + study + race 
# 3b. Create an exposure variable called 'subtype' with 5 categories: SIDD, SIRD, MOD, MARD and NoT2D



# 4. Fit a model that accounts for clustering

library(geepack)
m1 = geeglm(homa2b ~ subtype*ns(t) + age + sex + study + race, family = gaussian(),weights = ipw,data = df,id = study_id,corstr = "exchangeable")

# 5. Estimate predicted trajectories for each subtype over time along with confidence intervals
# https://stats.stackexchange.com/questions/109369/how-can-i-estimate-model-predicted-means-a-k-a-marginal-means-lsmeans-or-em
