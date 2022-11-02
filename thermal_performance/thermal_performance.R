
# Great tutorial here:
# https://ourcodingclub.github.io/tutorials/mixed-models/#what

library(tidyverse)
library(here)

library(lmerTest)
library(report)
library(performance)
# library(stargazer) # only works on lme4::lmer(), not lmerTest::lmer()
library(sjPlot)
library(emmeans)

options(scipen=999) # switch off scientific notation

# Import and prepare data ----------------------------------------

DATA <- read_csv("./DATA_MS.csv")


## Colony as factor --------------
# Since (old) colony names are numbers, let's make sure they're treated
# as categorical
DATA <- DATA %>% mutate(Colony = factor(Colony))
DATA$Colony %>% class()


# Check replication -----------------
DATA %>% 
  group_by(State, Day, Incub_type) %>% 
  summarize(n_polyps = length(Polyp_ID)) # %>% view()

# MAYBE OUTPUT TABLE?


# Plot: Overview of polyps by incubation -----------------------------------
ggplot(data = DATA, aes(x = factor(Day), y = Polyp_ID, color = Colony_ms)) + 
  geom_point(position = position_dodge(width = 0.75), 
             aes(shape = forcats::fct_rev(factor(Incub_type)))
  ) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    axis.title = element_blank()
  ) +
  facet_wrap(~forcats::fct_rev(factor(State)), 
             strip.position = "right", dir = "v", scale = "free" )


## Subsets --------------------------------------

# Remove controls and separate (subset) light and dark incubations

corals <- DATA %>%
  filter(State != "Control")

dark <- corals %>% 
  filter(Incub_type == "Dark")

light <- corals %>% 
  filter(Incub_type == "Light")



# Q1. Effect of **symbiotic state** on R rates (R~symbiotic~ vs R~bleached~) ------

# Consider only **dark** incubations because in light it is already clear 
# that symbiotic and bleached are different (positive and negative 
# rates resp.) and therefore no need to apply a statistical test.  

## Mod state1 ------------------
state1 <- lmerTest::lmer(Pn_ug ~ State + (1|Polyp_ID), 
                         data = dark, 
                         REML = F)

summary(state1)
report::report(state1)

### Checks ------------------
performance::check_normality(state1, type = "qq")
performance::check_normality(state1, type = "qq", effects = "fixed")
performance::check_normality(state1, type = "qq", effects = "random") # 
# all residuals appear normally distributed (p > 0.05) excpet for 
# the random effect (p = 0.042)

# Residuals
plot(state1)

qqnorm(resid(state1)); qqline(resid(state1))
 # not exactly the best ...


# P value
coefficients(summary(state1))[2,5] %>% round(., 4) # 0.0047 * significant


## Mod state2 ------------------
state2 <- lmerTest::lmer(Pn_ug ~ State + Colony_ms + (1|Polyp_ID), 
                         data = dark, 
                         REML = F)

# state2 <- lme4::lmer(Pn_ug ~ State + Colony_ms + (1|Polyp_ID), 
#                          data = dark, 
#                          REML = F)

summary(state2)
summary(state2) %>% coefficients()
report::report(state2)

performance::check_normality(state2, type = "qq")
performance::check_normality(state2, type = "qq", effects = "fixed")
performance::check_normality(state2, type = "qq", effects = "random")
# all residuals appear normally distributed (p > 0.05)

# Residuals
plot(state2)

qqnorm(resid(state2)); qqline(resid(state2))
# not exactly the best, but ok ...


# P value
coefficients(summary(state2))[2,5] %>% round(., 4) # 0.0002 *** very significant


# mod summary table
# stargazer::stargazer(state2, type = "text",
#                      digits = 3,
#                      star.cutoffs = c(0.05, 0.01, 0.001),
#                      digit.separator = "")



## Compare the two models --------------------

AIC(state1) # 1560.608
AIC(state2) # 1547.212 


### plot(performance::compare_performance(
###   state1,
###   state2,
###   rank = F))
### 
### 
### performance::compare_performance(
###   state1,
###   state2,
###   rank = T)

# In both models, 'State' is significant. The second model (state2) is better 
#  in term of distribution of residuals and AIC score.





# Q2. Effect of **colony** identity on Pn and R (RS1 vs RS2 vs RS3) -----------

# Q2: is there an *overall* significant difference between colonies in 
# their physiological performance?      
# Where *overall* means across the whole data set, hence considering both 
# light (Pn) and dark (R) incubations together (in the same model).  
# 
# Note that I have a unique variable `Pn_ug` (= Pn as ug of DO per ...) for 
# both light and dark incubations, which is positive when there is 
# a net production of oxygen (symbiotic in light), while is negative when 
# there is a net consumption of oxygen (symb in dark, and bleached always).


## Colony effect across all incubations (light & dark together) -----

## Mod colony 1 ---------------
colony1 <- lmer(Pn_ug ~ Colony + (1|Polyp_ID), 
                data = corals, REML = F)

# Residuals
plot(colony1)
qqnorm(resid(colony1)); qqline(resid(colony1))

plot(performance::check_normality(colony1, type = "qq"))

performance::check_normality(colony1, type = "qq")
performance::check_normality(colony1, type = "qq", effects = "fixed")
performance::check_normality(colony1, type = "qq", effects = "random")
# Non-normality detected for all effects!



## Mod colony 2 ---------------
colony2 <- lmer(Pn_ug ~ Colony + Incub_type + State + (1|Polyp_ID),
                data = corals, REML = F)

# Residuals
plot(colony2)
qqnorm(resid(colony2)); qqline(resid(colony2))

performance::check_normality(colony2, type = "qq")
performance::check_normality(colony2, type = "qq", effects = "fixed")
performance::check_normality(colony2, type = "qq", effects = "random")
# Non-normality detected for all effects!


## Mod colony 3 ---------------
colony3 <- lmer(Pn_ug ~ Colony + Incub_type * State + (1|Polyp_ID),
                data = corals, REML = F)

summary(colony3)
# Fixed effects:
# Estimate Std. Error       df t value             Pr(>|t|)    
# (Intercept)                    -11.2923     0.7881  38.0252 -14.328 < 0.0000000000000002 ***
# Colony53                        -2.5925     0.8647  24.3632  -2.998              0.00617 ** 
# Colony60                         0.5211     0.8606  24.5553   0.606              0.55040    
# Incub_typeLight                  3.1647     0.7175 442.4623   4.411           0.00001296 ***
# StateSymbiotic                  -4.6462     0.8654  54.6220  -5.369           0.00000167 ***
# Incub_typeLight:StateSymbiotic  22.3724     0.9999 442.6877  22.374 < 0.0000000000000002 ***

# Residuals
plot(colony3)
qqnorm(resid(colony3)); qqline(resid(colony3))

plot(performance::check_normality(colony3, type = "qq")) # normally distr.!

performance::check_normality(colony3, type = "qq")
performance::check_normality(colony3, type = "qq", effects = "fixed")
performance::check_normality(colony3, type = "qq", effects = "random")
# Warning: Non-normality of random effects detected (p = 0.002)

report::report(colony3) 
# We fitted a linear mixed model (estimated using ML and nloptwrap optimizer) to predict Pn_ug with Colony, Incub_type and State (formula: Pn_ug ~ Colony + Incub_type * State). The model included Polyp_ID as random effect (formula: ~1 | Polyp_ID). The model's total explanatory power is substantial (conditional R2 = 0.77) and the part related to the fixed effects alone (marginal R2) is of 0.76. The model's intercept, corresponding to Colony = 6, Incub_type = Dark and State = Bleached, is at -11.29 (95% CI [-12.84, -9.74], t(458) = -14.33, p < .001). Within this model:
#   
# - The effect of Colony [53] is statistically significant and negative 
# (beta = -2.59, 95% CI [-4.29, -0.89], t(458) = -3.00, p = 0.003; Std. beta = -0.23, 95% CI [-0.38, -0.08])
# - The effect of Colony [60] is statistically non-significant and positive 
# (beta = 0.52, 95% CI [-1.17, 2.21], t(458) = 0.61, p = 0.545; Std. beta = 0.05, 95% CI [-0.10, 0.20])
# - The effect of Incub type [Light] is statistically significant and positive 
# (beta = 3.16, 95% CI [1.75, 4.57], t(458) = 4.41, p < .001; Std. beta = 0.28, 95% CI [0.15, 0.40])
# - The effect of State [Symbiotic] is statistically significant and negative 
# (beta = -4.65, 95% CI [-6.35, -2.95], t(458) = -5.37, p < .001; Std. beta = -0.41, 95% CI [-0.56, -0.26])
# - The interaction effect of State [Symbiotic] on Incub type [Light] is statistically 
#   significant and positive (beta = 22.37, 95% CI [20.41, 24.34], t(458) = 22.37, p < .001; Std. beta = 1.97, 95% CI [1.80, 2.15])



## Mod colony 4 ---------------
colony4 <- lmer(Pn_ug ~ Colony + (1|Incub_type) + (1|State) + (1|Polyp_ID),
                data = corals, REML = F)

# Residuals
plot(colony4)
qqnorm(resid(colony4)); qqline(resid(colony4))

performance::check_normality(colony4, type = "qq")
performance::check_normality(colony4, type = "qq", effects = "fixed")
performance::check_normality(colony4, type = "qq", effects = "random")
# Non-normality detected for all effects!


## Conclusions --------------------

# In all models the residuals show clear patterns (bad) and 
#  are non-normally distributed!

# EXCEPT FOR MOD 3 THAT HAS NORMALLY DISTRIBUTED RESIDUALS!

# So, need to consider light and dark incubations separately.



# Q2.1. Colony effect in LIGHT incubations only -------------

## Mod colony light 1 ----------------
colony_light1 <- lmer(Pn_ug ~ Colony + (1|Polyp_ID),
                      data = light, REML = F)

summary(colony_light1)

# Residuals
plot(colony_light1) # clear pattern ...
qqnorm(resid(colony_light1)); qqline(resid(colony_light1))

performance::check_normality(colony_light1, type = "qq")
performance::check_normality(colony_light1, type = "qq", effects = "fixed")
performance::check_normality(colony_light1, type = "qq", effects = "random")
# Warning: Non-normality of random effects detected (p < .001).


## Mod colony light 2 ---------------
colony_light2 <- lmer(Pn_ug ~ Colony + State + (1|Polyp_ID),
                      data = light, REML = F)

summary(colony_light2)

# Residuals
plot(colony_light2) # clear pattern ...
qqnorm(resid(colony_light2)); qqline(resid(colony_light2))

performance::check_normality(colony_light2, type = "qq")
performance::check_normality(colony_light2, type = "qq", effects = "fixed")
performance::check_normality(colony_light2, type = "qq", effects = "random")
# All appear as normally distributed (p > 0.05)




## Mod colony light 3 ---------------
colony_light3 <- lmer(Pn_ug ~ Colony + State + Colony * State + (1|Polyp_ID),
                      data = light, REML = F)

summary(colony_light3)

# Residuals
plot(colony_light3) # clear pattern ...
qqnorm(resid(colony_light3)); qqline(resid(colony_light3))



## Mod colony light 4 ---------------
colony_light4 <- lmer(Pn_ug ~ Colony + (1|State) + (1|Polyp_ID),
                      data = light, REML = F)

summary(colony_light4)

# Residuals
plot(colony_light4) # clear pattern ...
qqnorm(resid(colony_light4)); qqline(resid(colony_light4))




# Q2.2. Colony effect in DARK incubations only -------------





# Q3. Pairwise comparison btw colonies (emmeans) -----------------------

# Test differences between colonies as pairwise comparison 
# of *estimated marginal means*, with *Bonferroni* correction.
ems <- emmeans::emmeans(colony_dark1, list(pairwise ~ Colony_ms), adjust = "bonferroni")

ems$`pairwise differences of Colony_ms` %>% as.data.frame()
