library(plotrix)
library(microbiome)
library(gridExtra)
library(nlme)
library(emmeans)
library(dplyr)
library(ggplot2)
library(Rmisc)
library(multcomp)

###Stats###

### Category ###

#Observed ASVs_ANOVA
Obs <- lm(Observed ~ Category, data = divMeta, na.action = na.omit)
summary(Obs)
anova(Obs) #no significant difference

ObsASVs <- summarySE(divMeta, groupvars = "Category", measurevar = "Observed")

 # Obtain least-squares means and confidence intervals
lsm <- emmeans(Obs, "Category", adjust = "sidak")

# Get Tukey's letters
lsm_tukey <- cld(lsm, alpha = 0.05)

pairs(emmeans(Obs, "Category", adjust="sidak"))

#contrast               estimate   SE df t.ratio p.value
#Bloodworms - Earthworm   -264.7 37.4 29  -7.068  <.0001
#Bloodworms - Mealworm     -13.3 37.4 29  -0.356  0.9964
#Bloodworms - Platypus     -63.2 28.3 29  -2.232  0.1965
#Bloodworms - Water        -94.8 35.0 29  -2.708  0.0772
#Earthworm - Mealworm      251.3 37.4 29   6.712  <.0001
#Earthworm - Platypus      201.5 28.3 29   7.118  <.0001
#Earthworm - Water         169.8 35.0 29   4.849  0.0003
#Mealworm - Platypus       -49.9 28.3 29  -1.761  0.4144
#Mealworm - Water          -81.5 35.0 29  -2.327  0.1654
#Platypus - Water          -31.6 25.0 29  -1.265  0.7141


#Simpsons_ANOVA
Sim <- lm(Simpson ~ Category,
          data = divMeta, na.action = na.omit)
summary(Sim)
anova(Sim) #no significant difference


#Shannon_ANOVA
Shan <- lm(Shannon ~ Category,
           data = divMeta, na.action = na.omit)
summary(Shan)
anova(Shan)

#Df  Sum Sq Mean Sq F value    Pr(>F)    
#Category   4 11.0487 2.76218  8.5987 0.0001054 ***
#Residuals 29  9.3158 0.32123   

# Obtain least-squares means and confidence intervals
lsm <- emmeans(Shan, "Category", adjust = "sidak")

# Get Tukey's letters
lsm_tukey <- cld(lsm, alpha = 0.05)

pairs(emmeans(Shan, "Category", adjust="sidak"))

### Captive ###

#Observed ASVs_ANOVA
Obs <- lm(Observed ~ Captive,
          data = divplatypus, na.action = na.omit)
summary(Obs)
anova(Obs)#no significant difference


#Simpsons_ANOVA
Sim <- lm(Simpson ~ Captive,
          data = divplatypus, na.action = na.omit)
summary(Sim)
anova(Sim)#no significant difference


#Shannon_ANOVA
Shan <- lm(Shannon ~ Captive,
           data = divplatypus, na.action = na.omit)
summary(Shan)
anova(Shan)#no significant difference

rm(Obs)
rm(Shan)
rm(Sim)
rm(divMeans)
rm(divplatypus)