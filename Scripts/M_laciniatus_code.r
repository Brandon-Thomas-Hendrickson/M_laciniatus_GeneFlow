# Download packages
install.packages(c("plyr", "lme4", "lmerTest", "DescTools", "psc1"))
lapply(c("plyr", "lme4", "lmerTest", "DescTools", "psc1"), require, character.only = TRUE)
list.of.packages <- c("plyr", "lme4", "lmerTest", "DescTools", "psc1")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(new.packages)

# Upload files
HILGARD <- read.csv("~/Documents/Github_Projects/M_laciniatus_GeneFlow/data/csv/HILGARD_ALL_ANALYSIS_READY.csv")

###############################
### Generalized Linear Models###
###############################

# Survival trait with binomial family distribution

survivalGLM <- glm(Survived_To_Flower ~ cross_type + as.factor(dam.transect) + cross_type * as.factor(dam.transect), family = binomial(link = "logit"), data = HILGARD)

survivalGLM_2009 <- glm(Survived_To_Flower ~ cross_type + as.factor(dam.transect) + cross_type * as.factor(dam.transect), family = binomial(link = "logit"), data = HILGARD[HILGARD$Year == 2009, ])

survivalGLM_2010 <- glm(Survived_To_Flower ~ cross_type + as.factor(dam.transect) + cross_type * as.factor(dam.transect), family = binomial(link = "logit"), data = HILGARD[HILGARD$Year == 2010,])

# Phenology
phenologyGLM <- glm(june_phenology ~ cross_type + as.factor(dam.transect) + cross_type * as.factor(dam.transect), family = poisson(link = "log"), data = HILGARD)

phenologyGLM_2009 <- glm(june_phenology ~ cross_type + as.factor(dam.transect) + cross_type * as.factor(dam.transect), family = poisson(link = "log"), data = HILGARD[HILGARD$Year == 2009, ])

phenologyGLM_2010 <- glm(june_phenology ~ cross_type + as.factor(dam.transect) + cross_type * as.factor(dam.transect), family = poisson(link = "log"), data = HILGARD[HILGARD$Year == 2010,])


# Lifetime Fitness (Composite Fitness)
fitnessGLM <- glm(Composite_Fitness ~ cross_type + as.factor(dam.transect) + cross_type * as.factor(dam.transect), inverse.gaussian(link = "1/mu^2"), data = HILGARD[HILGARD$Composite_Fitness>0,])

fitnessGLM_2009 <- glm(Composite_Fitness ~ cross_type + as.factor(dam.transect) + cross_type * as.factor(dam.transect), inverse.gaussian(link = "1/mu^2"), data = HILGARD[HILGARD$Year == 2009 & HILGARD$Composite_Fitness>0, ])

fitnessGLM_2010 <- glm(Composite_Fitness ~ cross_type + as.factor(dam.transect) + cross_type * as.factor(dam.transect), inverse.gaussian(link = "1/mu^2"), data = HILGARD[HILGARD$Year == 2010 & HILGARD$Composite_Fitness>0,])

# Fitness measured as seed mass (Only for 2010)

fitnessGLM_2010 <- glm(fruitmass.mg ~ cross_type + as.factor(dam.transect) + cross_type * as.factor(dam.transect), inverse.gaussian(link = "1/mu^2"), data = HILGARD[HILGARD$Year == 2010 & HILGARD$fruitmass.mg>0,])

# Final Height

heightGLM <- glm(final_height ~ cross_type + as.factor(dam.transect) + cross_type * as.factor(dam.transect), gaussian(link = "identity"),data=HILGARD)

heightGLM_2009 <- glm(final_height ~ cross_type + as.factor(dam.transect) + cross_type * as.factor(dam.transect), gaussian(link = "identity"), data = HILGARD[HILGARD$Year == 2009, ])

heightGLM_2010 <- glm(final_height ~ cross_type + as.factor(dam.transect) + cross_type * as.factor(dam.transect), gaussian(link = "identity"),data = HILGARD[HILGARD$Year == 2010,])

#############################################################################
#### Phenotypic Selection Analysis Using Reproductive Output and Survival####
#############################################################################

# Phenology selection
phenologySelection <- lm(Composite_Fitness ~ std_phenology, data = HILGARD)
phenologySelection_BLOCK <- lm(Composite_Fitness ~ std_block_phenology, data = HILGARD)

phenologySelection_2009 <- lm(Composite_Fitness ~ std_phenology, data = HILGARD)
phenologySelection_BLOCK <- lm(Composite_Fitness ~ std_block_phenology, data = HILGARD)

phenologySelection_2010 <- lm(Composite_Fitness ~ std_phenology, data = HILGARD[HILGARD$Year == 2009, ])
phenologySelection_BLOCK <- lm(Composite_Fitness ~ std_block_phenology, data = HILGARD[HILGARD$Year == 2010, ])


#Final Height

final_heightSelection <- lm(Composite_Fitness ~ std__final_height, data = HILGARD)
final_heightSelection_BLOCK <- lm(Composite_Fitness ~ std_block_final_height, data = HILGARD)

final_heightSelection_2009 <- lm(Composite_Fitness ~ std__final_height, data = HILGARD)
final_heightSelection_BLOCK <- lm(Composite_Fitness ~ std_block_final_height, data = HILGARD)

final_heightSelection_2010 <- lm(Composite_Fitness ~ std__final_height, data = HILGARD[HILGARD$Year == 2009, ])
final_heightSelection_BLOCK <- lm(Composite_Fitness ~ std_block_final_height, data = HILGARD[HILGARD$Year == 2010, ])

# Node

nodesSelection <- lm(Composite_Fitness ~ std__nodes, data = HILGARD)
nodesSelection_BLOCK <- lm(Composite_Fitness ~ std_block_nodes, data = HILGARD)

nodesSelection_2009 <- lm(Composite_Fitness ~ std__nodes, data = HILGARD)
nodesSelection_BLOCK <- lm(Composite_Fitness ~ std_block_nodes, data = HILGARD)

nodesSelection_2010 <- lm(Composite_Fitness ~ std__nodes, data = HILGARD[HILGARD$Year == 2009, ])
nodesSelection_BLOCK <- lm(Composite_Fitness ~ std_block_nodes, data = HILGARD[HILGARD$Year == 2010, ])

# # of Nodes at Flower (Only 2010)
numbernodesSelection <- lm(Composite_Fitness ~ std_..nodes, data = HILGARD)
numbernodesSelection_BLOCK <- lm(Composite_Fitness ~ std_block_..nodes, data = HILGARD)

# Vegetative Mass (Only 2010)
vegMassSelection <- lm(Composite_Fitness ~ std_veg.mass..mg., data = HILGARD)
vegMassSelection_BLOCK <- lm(Composite_Fitness ~ std_block_veg.mass..mg., data = HILGARD)

# Height @ Flower (Only 2010)
heightSelection <- lm(Composite_Fitness ~ std_height, data = HILGARD)
heightSelection_BLOCK <- lm(Composite_Fitness ~ std_block_height, data = HILGARD)

# Biomass (Only 2010)
heightSelection <- lm(Composite_Fitness ~  std_total.mass.2, data = HILGARD)
heightSelection_BLOCK <- lm(Composite_Fitness ~ std_block_total.mass.2, data = HILGARD)

# June Phenology Height (Only 2009)
heightSelection <- lm(Composite_Fitness ~  std_june_height, data = HILGARD)
heightSelection_BLOCK <- lm(Composite_Fitness ~ std_block_june_height, data = HILGARD)

############################################################################
#### Phenotypic Selection Analysis using Fruit mass (Only 2010) ############
############################################################################

# Phenology
phenologySelection_mass <- lm(fruitmass.mg ~ std_phenology, data = HILGARD)
phenologySelection_BLOCK_mass <- lm(fruitmass.mg ~ std_block_phenology, data = HILGARD)

# Final Height
finalheightSelection_mass <- lm(fruitmass.mg ~ std__final_height, data = HILGARD)
finalheightSelection_BLOCK_mass <- lm(fruitmass.mg ~ std_block_final_height, data = HILGARD)

# Height @ Flower
heightSelection_mass <- lm(fruitmass.mg ~ std_height, data = HILGARD)
heightSelection_BLOCK_mass <- lm(fruitmass.mg ~ std_block_height, data = HILGARD)

# Nodes 
nodesSelection_mass <- lm(fruitmass.mg ~ std__nodes, data = HILGARD)
nodesSelection_BLOCK_mass <- lm(fruitmass.mg ~ std_block_nodes, data = HILGARD)

# Nodes @ Flower
nodesFlowerSelection_mass <- lm(fruitmass.mg ~ std_..nodes, data = HILGARD)
nodesFlowerSelection_BLOCK_mass <- lm(fruitmass.mg ~ std_block_..nodes, data = HILGARD)

# Vegetative Mass
vegMassSelection_mass <- lm(fruitmass.mg ~ std_veg.mass..mg., data = HILGARD)
vegMassSelection_BLOCK_mass <- lm(fruitmass.mg ~ std_block_veg.mass..mg., data = HILGARD)

# Biomass 
biomassSelection_mass <- lm(fruitmass.mg ~ std_total.mass.2, data = HILGARD)
biomassSelection_BLOCK_mass <- lm(fruitmass.mg ~ std_block_total.mass.2, data = HILGARD)

################################################
#### REML Rank Transformed Traits and Fitness####
################################################

# Fitness
fitnessREML <- lmer(rank_Composite_Fitness ~ cross_type + as.factor(dam.transect) + cross_type * as.factor(dam.transect) + (1 | block), data = HILGARD, REML = TRUE)
anova(fitnessREML)


fitnessREML_2009 <- lmer(rank_Composite_Fitness ~ cross_type + as.factor(dam.transect) + cross_type * as.factor(dam.transect) + (1 | block), data = HILGARD[HILGARD$Year == 2009,], REML = TRUE)
anova(fitnessREML_2009)


fitnessREML_2010 <- lmer(rank_Composite_Fitness ~ cross_type + as.factor(dam.transect) + cross_type * as.factor(dam.transect) + (1 | block), data = HILGARD[HILGARD$Year == 2010,], REML = TRUE)
anova(fitnessREML_2010)


# Final Height

finalheightREML <- lmer(rank_final_height ~ cross_type + as.factor(dam.transect) + cross_type * as.factor(dam.transect) + (1 | block), data = HILGARD, REML = TRUE)
anova(finalheightREML)


finalheightREML_2009 <- lmer(rank_final_height ~ cross_type + as.factor(dam.transect) + cross_type * as.factor(dam.transect) + (1 | block), data = HILGARD[HILGARD$Year == 2009,], REML = TRUE)
anova(finalheightREML_2009)


finalheightREML_2010 <- lmer(rank_final_height ~ cross_type + as.factor(dam.transect) + cross_type * as.factor(dam.transect) + (1 | block), data = HILGARD[HILGARD$Year == 2010,], REML = TRUE)
anova(finalheightREML_2010)


# Phenology
phenologyREML <- lmer(rank_phenology ~ cross_type + as.factor(dam.transect) + cross_type * as.factor(dam.transect) + (1 | block), data = HILGARD, REML = TRUE)
anova(phenologyREML)


phenologyREML_2009 <- lmer(rank_phenology ~ cross_type + as.factor(dam.transect) + cross_type * as.factor(dam.transect) + (1 | block), data = HILGARD[HILGARD$Year == 2009,], REML = TRUE)
anova(phenologyREML_2009)


phenologyREML_2010 <- lmer(rank_phenology ~ cross_type + as.factor(dam.transect) + cross_type * as.factor(dam.transect) + (1 | block), data = HILGARD[HILGARD$Year == 2010,], REML = TRUE)
anova(phenologyREML_2010)


# Height @ Flower

height_flowerREML <- lmer(rank_height ~ cross_type + as.factor(dam.transect) + cross_type * as.factor(dam.transect) + (1 | block), data = HILGARD, REML = TRUE)
anova(height_flowerREML)

################################################################
### Dunnett's multiple comparisons test with Selfed as control###
################################################################

# Phenology
phenologyDunnett <- DescTools::DunnettTest(x = HILGARD$rank_phenology, g = HILGARD$cross_type, control = "SELF", data = HILGARD)

phenologyDunnett_2009 <- DescTools::DunnettTest(x = HILGARD$rank_phenology, g = HILGARD$cross_type, control = "SELF", data = HILGARD[HILGARD$Year == 2009, ])

phenologyDunnett_2010 <- DescTools::DunnettTest(x = HILGARD$rank_phenology, g = HILGARD$cross_type, control = "SELF", data = HILGARD[HILGARD$Year == 2010,])

# Lifetime Fitness
fitnessDunnett <- DescTools::DunnettTest(x = HILGARD$rank_Composite_Fitness, g = HILGARD$cross_type, control = "SELF", data = HILGARD)

fitnessDunnett_2009 <- DescTools::DunnettTest(x = HILGARD$rank_Composite_Fitness, g = HILGARD$cross_type, control = "SELF", data = HILGARD[HILGARD$Year == 2009, ])

fitnessDunnett_2010 <- DescTools::DunnettTest(x = HILGARD$rank_Composite_Fitness, g = HILGARD$cross_type, control = "SELF", data = HILGARD[HILGARD$Year == 2010,])

# Final Height
finalheightDunnett <- DescTools::DunnettTest(x = HILGARD$rank_final_height, g = HILGARD$cross_type, control = "SELF", data = HILGARD)

finalheightDunnett_2009 <- DescTools::DunnettTest(x = HILGARD$rank_final_height, g = HILGARD$cross_type, control = "SELF", data = HILGARD[HILGARD$Year == 2009, ])

finalheightDunnett_2010 <- DescTools::DunnettTest(x = HILGARD$rank_final_height, g = HILGARD$cross_type, control = "SELF", data = HILGARD[HILGARD$Year == 2010,])

# Height @ Flower (Only for 2010)
heightDunnett_2010 <- DescTools::DunnettTest(x = HILGARD$rank_height, g = HILGARD$cross_type, control = "SELF", data = HILGARD[HILGARD$Year == 2010, ])

# Vegetative Mass (Only for 2010)
vegmassDunnett_2010 <- DescTools::DunnettTest(x = HILGARD$rank_veg.mass..mg., g = HILGARD$cross_type, control = "SELF", data = HILGARD[HILGARD$Year == 2010, ])

# Total Biomass (Only for 2010)
biomassDunnett_2010 <- DescTools::DunnettTest(x =  HILGARD$rank_total.mass.2, g = HILGARD$cross_type, control = "SELF", data = HILGARD[HILGARD$Year == 2010, ])

##################################################
### Zero Inflated and Count submodel Regressions###
##################################################

# Zero inflated regression with Poisson family (Using untransformed Composite_Fitness data)

zeroPoisson <- zeroinfl(Composite_Fitness ~ cross_type + as.factor(dam.transect) + block| 1, dist = "poisson", link = "logit", data = HILGARD)

# Zero inflated regression with Poisson family (Using untransformed fruit mass data)

zeroPoisson_mass <- zeroinfl(fruitmass.microgram ~ cross_type + as.factor(dam.transect) + block| 1, dist = "poisson", link = "logit", data = HILGARD)

# Count regression with Poiisson family (Using untransformed Composite_Fitness data)

countPoisson <- glm(Composite_Fitness ~ cross_type + as.factor(dam.transect) + block, family=poisson, data = HILGARD)

# Count regression with Poisson family (Using untransformed fruit mass data)

massPoisson <- glm(fruitmass.microgram ~ cross_type + as.factor(dam.transect) + block, family=poisson, data = HILGARD)

# Model Comparison of Zero-Inflated and Count submodels (Composite Fitness)

vTest <- vuong(zeroPoisson, countPoisson)

# Model Comparison of Zero-Inflated and Count submodels (fruit mass)

vTest_mass <- vuong(zeroPoisson_mass, massPoisson)





