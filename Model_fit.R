library(dplyr)
library(ade4)
library(lme4)
library(glmmTMB)
library(MASS)
library(car)
library(emmeans)

final_dataset$`Time period` <- as.numeric(final_dataset$`Time period`)
final_dataset$Elevation <- as.numeric(final_dataset$Elevation)

group_stats <- final_dataset %>%
  summarise(
    Mean_Number = mean(Number),
    Variance_Number = var(Number),
    Overdispersion = Variance_Number / Mean_Number
  )

group_stats

# Fit a Poisson GLM
poisson_model <- glmer(Number ~ Elevation + Mountain + Wingspan*Dietary+(1|`Time period`)+(1|Species), data = final_dataset,family = poisson(link="log"))

summary(poisson_model)

# Calculate the Pearson chi-square statistic
resid_dev <- sum(residuals(poisson_model, type = "pearson")^2)
resid_df <- df.residual(poisson_model)

overdispersion_ratio <- resid_dev / resid_df
overdispersion_ratio

nb_model <- glmmTMB(
  Number ~ Elevation + Temperature+ Wind +Mountain+ (1 | `Time period`) + (1 | Species),
  data = final_dataset, family = nbinom2(link = "log"))
summary(nb_model)

Anova(nb_model, type = "III")

library(knitr)
tidy_mod <- tidy(nb_model, effects = "fixed", conf.int = TRUE)
kable(tidy_mod, digits = 3, caption = "Model Estimates and 95% Confidence Intervals")

################################################################################################################
#MODEL pro stanovení interakce SpeciesRichness~Movement*Treatment rozdělených do Season!! final!!
species_richness_data <- final_dataset %>% 
  group_by(`Time period`, Elevation, Mountain) %>% 
  summarize(species_richness = n_distinct(Species))

final_dataset2 <- final_dataset %>% 
  left_join(species_richness_data, by = c("Time period", "Elevation", "Mountain"))

mod <- glmmTMB(species_richness ~ Elevation + Mountain + Temperature + Wind + 
                 (1 | `Time period`) + (1 | Species),
               data = final_dataset2, 
               family = nbinom2(link = "log"))

summary(mod)
Anova(mod,type = "III")

library(broom.mixed)
library(ggplot2)

# Tidy the model to get fixed effects with confidence intervals
tidy_mod <- tidy(mod, effects = "fixed", conf.int = TRUE)

# Plot the estimates and their confidence intervals
ggplot(tidy_mod, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(title = "Fixed Effects Estimates with 95% Confidence Intervals",
       x = "", y = "Estimate (log scale)")

library(knitr)
tidy_mod <- tidy(mod, effects = "fixed", conf.int = TRUE)
kable(tidy_mod, digits = 3, caption = "Model Estimates and 95% Confidence Intervals")

################################################################################################################
# Load required libraries
library(dplyr)
library(lme4)
library(car)
library(DHARMa)
library(ade4)
library(reshape2)
library(dplyr)
library(tibble)  # Load tibble for row name conversion

final_dataset <- final_dataset %>%
  mutate(
    Dietary = as.numeric(as.factor(Dietary)),
    Distribution = as.numeric(as.factor(Distribution)),
    `Host species` = as.numeric(as.factor(`Host species`)),
    Overwintering = as.numeric(as.factor(Overwintering)),
    `Leaf action` = as.numeric(as.factor(`Leaf action`))
  )

colnames(final_dataset)

cwm_results <- final_dataset %>%
  group_by(Elevation, Mountain) %>%
  summarize(
    Dietary_cwm = weighted.mean(Dietary, Number, na.rm = TRUE),
    Red_list_cwm = weighted.mean(`Red list species`, Number, na.rm = TRUE),
    Wingspan_cwm = weighted.mean(Wingspan, Number, na.rm = TRUE),
    Distribution_cwm = weighted.mean(Distribution, Number, na.rm = TRUE),
    Host_species_cwm = weighted.mean(`Host species`, Number, na.rm = TRUE),
    Overwintering_cwm = weighted.mean(Overwintering, Number, na.rm = TRUE),
    Leaf_action_cwm = weighted.mean(`Leaf action`, Number, na.rm = TRUE),Abundance=sum(Number)
  )

# Display results
print(cwm_results)

# Save results to a CSV file
write.csv(cwm_results, "cwm_results.csv", row.names = FALSE)

# Convert categorical variables to factors
cwm_results <- cwm_results %>%
  mutate(Mountain = as.factor(Mountain))

# Convert categorical variables to factors
cwm_results <- cwm_results %>%
  mutate(Mountain = as.factor(Mountain),
         `Time period` = as.factor(`Time period`))

# Fit Generalized Linear Mixed Model (GLMM)
mod1 <- lm(Dietary_cwm ~ Elevation + Mountain,
               data = cwm_results)
Anova(mod1,type = "III")

confint(mod1, method = "boot", nsim = 999)

library(vegan)
d <- dist(cwm_results$Distribution_cwm)
mod2 <- adonis2(d ~ Elevation, data = cwm_results,permutations = 999)
print(mod2)


#####
library(pbkrtest)

# Full model (with all predictors)
mod1 <- lm(Dietary_cwm ~ Elevation + Mountain, data = cwm_results)

# Reduced model (without 'Mountain')
mod0 <- lm(Dietary_cwm ~ Mountain, data = cwm_results)

# Run parametric bootstrap test
pb <- PBmodcomp(mod1, mod0, nsim = 999)

# Print results
summary(pb)

isSingular(mod1) #if random are properly estimated and does not collapse to zero variance
VarCorr(mod1) #what variance Time period has

#####
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)

# Summarize total abundance per site (optional but informative)
abundance_summary <- final_dataset %>%
  group_by(Elevation, Mountain) %>%
  summarise(total_abundance = sum(Number), .groups = "drop")

# Create site identifier
final_dataset <- final_dataset %>%
  mutate(Site = paste(Elevation, Mountain, sep = "_"))

# Create community matrix (species x site)
comm_matrix <- final_dataset %>%
  group_by(Site, Species) %>%
  summarise(Abundance = sum(Number), .groups = "drop") %>%
  pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0) %>%
  as.data.frame()

# Set rownames and drop the 'Site' column
rownames(comm_matrix) <- comm_matrix$Site
comm_matrix <- comm_matrix[, -1]

# Compute diversity indices
shannon <- diversity(comm_matrix, index = "shannon")      # Shannon entropy (H')
richness <- specnumber(comm_matrix)                       # Species richness (S)
evenness <- shannon / log(richness)                       # Pielou's evenness

# Combine into one dataframe
evenness_df <- data.frame(
  Site = rownames(comm_matrix),
  Shannon = shannon,
  Richness = richness,
  Evenness = evenness
)

#####
# Calculate eveness for each site 
# Create a site ID
final_dataset <- final_dataset %>%
  mutate(Site = paste(Elevation, Mountain, sep = "_"))

# Build the site-by-species abundance matrix
comm_matrix <- final_dataset %>%
  group_by(Site, Species) %>%
  summarise(Abundance = sum(Number), .groups = "drop") %>%
  pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0) %>%
  as.data.frame()

# Set rownames
rownames(comm_matrix) <- comm_matrix$Site
comm_matrix <- comm_matrix[, -1]  # Remove 'Site' column now in rownames

# Shannon entropy (H')
shannon <- diversity(comm_matrix, index = "shannon")

# Species richness (S)
richness <- specnumber(comm_matrix)

# Pielou's evenness (J′ = H' / log(S))
evenness <- shannon / log(richness)

# Combine into a results table
evenness_df <- data.frame(
  Site = rownames(comm_matrix),
  Shannon = shannon,
  Richness = richness,
  Evenness = evenness
)

evenness_df <- evenness_df %>%
  separate(Site, into = c("Elevation", "Mountain"), sep = "_") %>%
  mutate(Elevation = as.numeric(Elevation))

ggplot(evenness_df, aes(x = Elevation, y = Evenness)) +
  geom_point() +
  geom_smooth(method = "loess") +
  labs(title = "Pielou’s Evenness (J′) across Elevation",
       x = "Elevation", y = "Evenness (J′)")
