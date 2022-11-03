# Script description --------------------

# What's happening here:

# Q1: is there a significant difference between SYMBIOTIC STATES 
#  (symb vs bleached) in DARK incubations (so looking at R - respiration)?
#   (note: no need to test this for light incubation as that is obvious:
#   symbiotic have net O2 production while bleached have net O2 consumption)

# Q2: is there a significant difference between COLONIES (RS1, RS2, RS3)?
#  For this we look at light incubations (PN, Q2.1) and 
#  dark incubations (R; Q2.2) separately.


# Load packages -------------------------

# R version 4.1.0 (2021-05-18)

library(tidyverse)
library(here)

library(lmerTest)
library(report)
library(performance)
# library(sjPlot)
library(emmeans)


options(scipen = 999) # switch off scientific notation


# Import and prepare data ----------------------------------------

DATA <- read_csv("./DATA_MS.csv")

## Colony as factor 
# Since (old) colony names are numbers, let's make sure they're treated
# as categorical, however I will use 'Colony_ms' to match names with the MS
DATA <- DATA %>% mutate(Colony = factor(Colony))
DATA$Colony %>% class()


# Check replication 
DATA %>% 
  group_by(State, Day, Incub_type) %>% 
  summarize(n_polyps = length(Polyp_ID)) # %>% view()



# Plot: Overview of polyps by incubation
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



## Subsets
# Remove controls and separate (subset) light and dark incubations
corals <- DATA %>%
  filter(State != "Control")

dark <- corals %>% 
  filter(Incub_type == "Dark")

light <- corals %>% 
  filter(Incub_type == "Light")


# Notes on statistical approach --------------------

# 1. In all models: set Polyp_id as random factor (1|Polyp_ID) because the same
#    polyps were repeatedly measured at every incubation (every day)
# 2. Data is not transformed, as normality of the DV is not a requirement, rather
#    residuals should show no clear patterns. See: https://www.biorxiv.org/content/10.1101/305946v1
# 3. Variable such as 'State' (= symbiotic state) and 'Colony_ms' (colony id)
#    are always treated as fixed factors (even when I would consider them a 
#    random factor) because the have < 5 levels.
#    -> as a result, Q1 and Q2.2 can be answered by the same model
#  4. Temperature ('temp_byWB_mean' = mean temperature of the waterbath) is 
#     a continuos vaiable and as such is always treated as fixed factor.



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

# Checks
performance::check_normality(state1, type = "qq", effects = "fixed")
performance::check_normality(state1, type = "qq", effects = "random")

plot(state1)
qqnorm(resid(state1)); qqline(resid(state1))



## Mod state2 ------------------
state2 <- lmerTest::lmer(Pn_ug ~ State + Colony_ms + (1|Polyp_ID), 
                         data = dark, 
                         REML = F)

summary(state2)
report::report(state2)

performance::check_normality(state2, type = "qq", effects = "fixed")
performance::check_normality(state2, type = "qq", effects = "random")

plot(state2)
qqnorm(resid(state2)); qqline(resid(state2))



## Mod state3 ------------------
state3 <- lmerTest::lmer(Pn_ug ~ State + Colony_ms + Temp_byWB_mean + (1|Polyp_ID), 
                         data = dark, 
                         REML = F)

summary(state3)
report::report(state3)

performance::check_normality(state3, type = "qq", effects = "fixed")
performance::check_normality(state3, type = "qq", effects = "random")

plot(state3)
qqnorm(resid(state3)); qqline(resid(state3))



## Compare the models --------------------

AIC(state1) 
AIC(state2) 
AIC(state3) # <- best


 plot(performance::compare_performance(
   state1,
   state2,
   state3, # <- best
   rank = F))
 
 
 performance::compare_performance(
   state1,
   state2,
   state3, # <- best
   rank = T) %>% view()
 

## Q2. Conclusions -----------------------

# In all models, 'State' is highly significant. 
# The model 'state3' (includes temperature) has the best explanatory power and 
#  and performance, so we consider the significance levels of this model.
# Namely, we see that 'State', 'Colony_ms', and 'Temp_byWB_mean' are all 
# significant:
#  Fixed effects:
#                    Estimate Std. Error        df t value             Pr(>|t|)    
#  (Intercept)     31.06005    2.15223 185.29230  14.432 < 0.0000000000000002 ***
#  StateSymbiotic  -4.17309    1.08594  23.95737  -3.843             0.000785 ***
#  Colony_msRS2    -5.83244    1.32940  23.91336  -4.387             0.000199 ***
#  Colony_msRS3    -0.28202    1.32799  23.82782  -0.212             0.833626    
#  Temp_byWB_mean  -1.61009    0.07248 210.75853 -22.216 < 0.0000000000000002 ***

# >>> Therefore this already answers the following question of whether colony 
# identity has a significant effect on O2 evolution (R) in dark incubations
# (Q2.2)
 

# >>> Next, we want to look at **light** (PN) incubations  
#  to see if there are significant differences between colonies 




# Q2.1. Colony effect in LIGHT incubations only -------------

## Mod colony light 1 ----------------
colony_light1 <- lmer(Pn_ug ~ Colony_ms + (1|Polyp_ID),
                      data = light, REML = F)

summary(colony_light1)

performance::check_normality(colony_light1, type = "qq", effects = "fixed")
performance::check_normality(colony_light1, type = "qq", effects = "random")

plot(colony_light1) # clear pattern ...
qqnorm(resid(colony_light1)); qqline(resid(colony_light1))



## Mod colony light 2 ---------------
colony_light2 <- lmer(Pn_ug ~ Colony_ms + State + (1|Polyp_ID),
                      data = light, REML = F)

summary(colony_light2) # again, State is signif, but not Colony_ms

performance::check_normality(colony_light2, type = "qq", effects = "fixed")
performance::check_normality(colony_light2, type = "qq", effects = "random")

plot(colony_light2)
qqnorm(resid(colony_light2)); qqline(resid(colony_light2))



## Mod colony light 3 ---------------

colony_light3 <- lmer(Pn_ug ~ Colony_ms + State + Colony_ms * State + (1|Polyp_ID),
                      data = light, REML = F)

summary(colony_light3)

performance::check_normality(colony_light3, type = "qq", effects = "fixed")
performance::check_normality(colony_light3, type = "qq", effects = "random")

plot(colony_light3) 
qqnorm(resid(colony_light3)); qqline(resid(colony_light3))



## Mod colony light 4 ---------------
colony_light4 <- lmer(Pn_ug ~ Colony_ms + State + Temp_byWB_mean + (1|Polyp_ID),
                      data = light, REML = F)

summary(colony_light4)

performance::check_normality(colony_light4, type = "qq", effects = "fixed")
performance::check_normality(colony_light4, type = "qq", effects = "random")

plot(colony_light4) 
qqnorm(resid(colony_light4)); qqline(resid(colony_light4))




## Compare the models --------------------

AIC(colony_light1) 
AIC(colony_light2)  
AIC(colony_light3) 
AIC(colony_light4) # <- the best


plot(performance::compare_performance(
  colony_light1,
  colony_light2,
  colony_light3,
  colony_light4, # <- the best
  rank = F))


performance::compare_performance(
  colony_light1,
  colony_light2,
  colony_light3,
  colony_light4, # <- the best
  rank = T) %>% view()



## Q2.1. Conclusions ----------------------------------

# In none of the models 'Colony_ms' (= colony id) is significant. 
# The model 'colony_light4' (includes temperature) has the best explanatory 
# power and performance, so we consider the significance levels of this model.
# Namely, we see that 'Colony_ms' is not significant in light incubations 
#  (this answers our question), while, as before, 'State', and 'Temp_byWB_mean' are significant:
# 
#                  Estimate Std. Error        df t value             Pr(>|t|)    
#  (Intercept)     11.30993    1.56905 189.51179   7.208      0.0000000000131 ***
#  Colony_msRS2     0.13288    0.95211  24.58620   0.140                0.890    
#  Colony_msRS3     0.56802    0.94038  25.22523   0.604                0.551    
#  StateSymbiotic  17.83048    0.77183  25.11451  23.101 < 0.0000000000000002 ***
#  Temp_byWB_mean  -0.78535    0.05236 212.04537 -14.999 < 0.0000000000000002 ***







# Q2.2.3. Pairwise comparison btw colonies (emmeans) -----------------------

# Now that we have established that colony id has a significant effect in 
#  **dark** incubations (R), we want to see which colonies are significantly 
# different from each other.

# >>> Test differences between colonies as pairwise comparison 
# of *estimated marginal means*, with *Bonferroni* correction.
ems <- emmeans::emmeans(state3, list(pairwise ~ Colony_ms), adjust = "bonferroni")

ems$`pairwise differences of Colony_ms` %>% as.data.frame()
#             estimate       SE       df    t.ratio     p.value
# RS1 - RS2  5.8324385 1.456848 28.78975  4.0034649 0.001201065
# RS1 - RS3  0.2820232 1.455486 28.71707  0.1937656 1.000000000
# RS2 - RS3 -5.5504153 1.459409 28.94023 -3.8031929 0.002047508


## Q2.2.3. Conclusions ---------------------
# Colony RS2 (old name: '#53') is significantly different from the other two
#  colonies (stats confirm pattern visible in Fig. S4).
