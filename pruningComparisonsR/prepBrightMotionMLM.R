#### LOAD PACKAGES ####
# Load necessary libraries
library(lme4)
library(ggplot2)
library(MASS)  # For Box-Cox 
library(car)   # For powerTransform 
# for bootstrapping:
library(foreach)
library(doParallel)
library(dplyr)
library(effectsize)
library(MuMIn)
library(viridis) #for colourblind friendly palette


#### WORKSPACE AND PARAMETERS SETUP ####
#set working directory (for saving outputs - bootstrapping is time-consuming!)
setwd("[path-to...]/pruningComparisonsR")

# Load the data
dataExp <- read.csv("[path-to...]/pruningComparisonsMatlab/stats/bright/overall/pruneChannelInfoTableQT.csv")


# Calculate total, median, mean and standard deviation of Percentage_Motion_Windows
result <- dataExp %>%
  group_by(Age, Cohort) %>%
  summarize(
    total_percentage_motion = sum(Percentage_Motion_Windows, na.rm = TRUE)*100,  
    median_percentage_motion = median(Percentage_Motion_Windows, na.rm = TRUE)*100,
    mean_percentage_motion = mean(Percentage_Motion_Windows, na.rm = TRUE)*100,   
    sd_percentage_motion = sd(Percentage_Motion_Windows, na.rm = TRUE)*100       
  )
print(result)

# rescale variables
dataExp <- dataExp %>%
  mutate(
    Age = scale(Age),
    Channel = scale(Channel),
    Homer_Pruned_Total = scale(Homer_Pruned_Total),
    Average_SCI = scale(Average_SCI),
    Average_PSP = scale(Average_PSP),
    Percentage_Motion_Windows = scale(Percentage_Motion_Windows),
    AP_Orientation = scale(AP_Orientation)
  )

# Bootstrapping params
# Set number of bootstrap iterations
nBoot <- 10000
#nBoot <- 2000 #for testing code - comment out!
# number of cores to use
numCores <- detectCores() - 1
# set seed for reproducibility
set.seed(42)

#### HELPER FUNCTIONS ####
# Function to create model formulas given a set of predictors and outcome variable
create_model <- function(outcome, predictors) {
  formula <- as.formula(paste(outcome, "~", paste(predictors, collapse = " + "), "+ (1 | ID)"))
  model <- lmer(formula, data = dataExp)
  return(model)
}

# Helper function for bootstrapping
bootstrap_iteration <- function(data, model_eqn, i) {
  # Resample data with replacement
  indices <- sample(nrow(data), replace = TRUE)
  bootData <- data[indices, ]
  
  # Fit model
  bootModel <- try(lmer(model_eqn, data = bootData), silent = TRUE)
  
  # Extract fixed effects
  if (!inherits(bootModel, "try-error")) {
    fixefs <- fixef(bootModel)
    seFixefs <- sqrt(diag(vcov(bootModel)))
    
    return(fixefs)
  } else {
    return(rep(NA, length(fixef(model_Motion))))  
  }
}

# Function to calculate standardised coefficients and their CIs
standardiseEffects <- function(model, summary_df) {
  # Get model frame
  mf <- model.frame(model)
  
  # Get response variable
  response <- model.response(mf)
  response_sd <- sd(response)
  
  # Get predictor standard deviations (for numeric variables)
  predictors <- model.matrix(model)
  numeric_vars <- sapply(mf, is.numeric)
  predictor_sds <- sapply(mf[, numeric_vars, drop = FALSE], sd)
  
  # Get original coefficients and SEs
  coefs <- summary_df$Estimate
  ses <- summary_df$SE
  
  # Initialise results df
  results <- data.frame(
    Variable = rownames(summary_df),
    Raw_Estimate = coefs,
    Raw_CI_lower = summary_df$CI_lower,
    Raw_CI_upper = summary_df$CI_upper,
    Raw_SE = ses,
    stringsAsFactors = FALSE
  )
  
  # Calculate standardised effects
  std_effects <- coefs
  std_ses <- ses
  
  for(i in 1:length(coefs)) {
    var_name <- rownames(summary_df)[i]
    
    if(var_name == "(Intercept)") {
      std_effects[i] <- coefs[i] / response_sd
      std_ses[i] <- ses[i] / response_sd
    } else {
      # Check if main effect or interaction
      terms <- strsplit(var_name, ":")[[1]]
      
      if(length(terms) == 1) {
        # Main effect
        if(terms %in% names(predictor_sds)) {
          # Numeric predictor
          std_effects[i] <- coefs[i] * predictor_sds[terms] / response_sd
          std_ses[i] <- ses[i] * predictor_sds[terms] / response_sd
        } else {
          # Categorical predictor
          std_effects[i] <- coefs[i] / response_sd
          std_ses[i] <- ses[i] / response_sd
        }
      } else {
        # Interaction
        multiplier <- 1
        for(term in terms) {
          if(term %in% names(predictor_sds)) {
            multiplier <- multiplier * predictor_sds[term]
          }
        }
        std_effects[i] <- coefs[i] * multiplier / response_sd
        std_ses[i] <- ses[i] * multiplier / response_sd
      }
    }
  }
  
  # Add standardised effects to results
  results$Std_Estimate <- std_effects
  results$Std_SE <- std_ses
  results$Std_CI_lower <- std_effects - 1.96 * std_ses
  results$Std_CI_upper <- std_effects + 1.96 * std_ses
  
  # Calculate Cohen's f^2 for each variable
  f2 <- numeric(length(coefs))
  names(f2) <- rownames(summary_df)
  
  for(i in 1:length(coefs)) {
    var_name <- rownames(summary_df)[i]
    if(var_name != "(Intercept)") {
      # Fit model without this term
      terms <- strsplit(var_name, ":")[[1]]
      if(length(terms) == 1) {
        # For main effects
        formula_without <- update(formula(model), paste(". ~ . -", var_name))
        model_without <- lmer(formula_without, data = mf)
        r2_with <- r.squaredGLMM(model)[1]
        r2_without <- r.squaredGLMM(model_without)[1]
        f2[i] <- (r2_with - r2_without) / (1 - r2_with)
      }
    }
  }
  
  results$Cohens_f2 <- f2
  
  # Add effect size interpretation
  results$Effect_Size <- case_when(
    abs(results$Std_Estimate) < 0.2 ~ "Small",
    abs(results$Std_Estimate) < 0.5 ~ "Medium",
    TRUE ~ "Large"
  )
  
  return(results)
}


# #### DETERMINE MODEL ####
# # Define the list of factors for Percentage_Motion_Windows
# factors <- c("Age", "Task", "Cohort", "AP_Orientation")
# 
# # Initialise an empty list to store models
# models <- list()
# 
# # Define letters for naming models (a, b, c, ...)
# letters_list <- letters
# 
# # Set the outcome variable
# outcome <- "Percentage_Motion_Windows"
# 
# # Loop over the number of factors in each model (1 to 4)
# for (n_factors in 1:4) {
#   # Generate combinations of factors
#   factor_combinations <- combn(factors, n_factors, simplify = FALSE)
#   
#   # Create models for each combination of factors
#   for (i in seq_along(factor_combinations)) {
#     predictors <- factor_combinations[[i]]
#     model_name <- paste0("model", n_factors, letters_list[i], "_Motion") # e.g., model1a_Motion
#     models[[model_name]] <- create_model(outcome, predictors)
#   }
# }
# 
# # Compare all models
# # Collect model summaries
# modelSummaries <- models
# 
# # Extract AIC and BIC for comparison
# aicBic <- sapply(modelSummaries, function(model) c(AIC(model), BIC(model)))
# 
# # Combine results into a data frame for easy viewing
# results <- data.frame(
#   Model = names(modelSummaries),
#   AIC = aicBic[1, ],
#   BIC = aicBic[2, ]
# )
# 
# # Sort results by AIC and BIC in ascending order
# # AIC
# results_sorted <- results[order(results$AIC), ]
# # Display sorted results
# print("Results for Percentage_Motion_Windows Models sorted by AIC:")
# print(results_sorted)
# 
# # BIC
# results_sorted <- results[order(results$BIC), ]
# # Display sorted results
# print("Results for Percentage_Motion_Windows Models sorted by BIC:")
# print(results_sorted)
# 
# # Find models with the lowest AIC and BIC
# lowestAIC <- results_sorted[which.min(results_sorted$AIC), ]
# lowestBIC <- results_sorted[which.min(results_sorted$BIC), ]
# 
# # Display the best models
# print("Lowest AIC Model for Percentage_Motion_Windows:")
# print(lowestAIC)
# print("Lowest BIC Model for Percentage_Motion_Windows:")
# print(lowestBIC)
# 
# # Show formulas for the top 3 models by AIC
# print("Top 3 Models for Percentage_Motion_Windows by AIC:")
# for (i in 1:min(3, nrow(results_sorted))) {
#   model_name <- results_sorted$Model[i]
#   cat(paste0("Model: ", model_name, "\n"))
#   print(formula(models[[model_name]]))
#   cat("\n")
# }
# 
# # Show formulas for the top 3 models by BIC
# print("Top 3 Models for Percentage_Motion_Windows by BIC:")
# for (i in 1:min(3, nrow(results_sorted))) {
#   model_name <- results_sorted$Model[i]
#   cat(paste0("Model: ", model_name, "\n"))
#   print(formula(models[[model_name]]))
#   cat("\n")
# }
# 
# 
# # Add interaction models to the models list
# # Effect on motion by age might differ by task if infants bored due to stimulus 
# # and/or fact that HaND task is run second (infants have been capped and 
# # involved in data recording for a longer amount of time => bored!)
# # Effect on motion may differ by cohort due to e.g. UK infants becoming more
# # familiar with screens and therefore less interested
# # Possible for effect of optode position to change with age if optodes located
# # in positions with more curvature in skull => more susceptible to being displaced
# # due to movement
# modelIntTask_Motion <- lmer(Percentage_Motion_Windows ~ Age * Task + (1 | ID), data = dataExp)
# modelIntCohort_Motion <- lmer(Percentage_Motion_Windows ~ Age * Cohort + (1 | ID), data = dataExp)
# modelIntOrientation_Motion <- lmer(Percentage_Motion_Windows ~ Age * AP_Orientation + (1 | ID), data = dataExp)
# modelIntTaskCohort_Motion <- lmer(Percentage_Motion_Windows ~ Age * Task + Age * Cohort + AP_Orientation + (1 | ID), data = dataExp)
# modelIntCohortOrientation_Motion <- lmer(Percentage_Motion_Windows ~ Age * Cohort + Age * AP_Orientation + Task + (1 | ID), data = dataExp)
# modelIntTaskOrientation_Motion <- lmer(Percentage_Motion_Windows ~ Age * Task + Age * AP_Orientation + Cohort + (1 | ID), data = dataExp)
# modelIntAll_Motion <- lmer(Percentage_Motion_Windows ~ Age * Task + Age * Cohort + Age * AP_Orientation + (1 | ID), data = dataExp)
# 
# models[["modelIntTask_Motion"]] <- modelIntTask_Motion
# models[["modelIntCohort_Motion"]] <- modelIntCohort_Motion
# models[["modelIntOrientation_Motion"]] <- modelIntOrientation_Motion 
# models[["modelIntTaskCohort_Motion"]] <- modelIntTaskCohort_Motion
# models[["modelIntCohortOrientation_Motion"]] <- modelIntCohortOrientation_Motion
# models[["modelIntTaskOrientation_Motion"]] <- modelIntTaskOrientation_Motion
# models[["modelIntAll_Motion"]] <- modelIntAll_Motion
# 
# # Compare all models again
# # Collect model summaries
# modelSummaries <- models
# 
# # Extract AIC and BIC for comparison
# aicBic <- sapply(modelSummaries, function(model) c(AIC(model), BIC(model)))
# 
# # Combine results into a data frame for easy viewing
# results <- data.frame(
#   Model = names(modelSummaries),
#   AIC = aicBic[1, ],
#   BIC = aicBic[2, ]
# )
# 
# # Sort results by AIC and BIC in ascending order
# # AIC
# results_sorted <- results[order(results$AIC), ]
# # Display sorted results
# print("Results for Percentage_Motion_Windows Models sorted by AIC:")
# print(results_sorted)
# 
# # BIC
# results_sorted <- results[order(results$BIC), ]
# # Display sorted results
# print("Results for Percentage_Motion_Windows Models sorted by BIC:")
# print(results_sorted)
# 
# # Find models with the lowest AIC and BIC
# lowestAIC <- results_sorted[which.min(results_sorted$AIC), ]
# lowestBIC <- results_sorted[which.min(results_sorted$BIC), ]
# 
# # Display the best models
# print("Lowest AIC Model for Percentage_Motion_Windows:")
# print(lowestAIC)
# print("Lowest BIC Model for Percentage_Motion_Windows:")
# print(lowestBIC)
# 
# # Show formulas for the top 5 models by AIC
# print("Top 3 Models for Percentage_Motion_Windows by AIC:")
# for (i in 1:min(5, nrow(results_sorted))) {
#   model_name <- results_sorted$Model[i]
#   cat(paste0("Model: ", model_name, "\n"))
#   print(formula(models[[model_name]]))
#   cat("\n")
# }
# 
# # Show formulas for the top 5 models by BIC
# print("Top 5 Models for Percentage_Motion_Windows by BIC:")
# for (i in 1:min(5, nrow(results_sorted))) {
#   model_name <- results_sorted$Model[i]
#   cat(paste0("Model: ", model_name, "\n"))
#   print(formula(models[[model_name]]))
#   cat("\n")
# }
# 
# # Rename the model with best outcome.
# model_MotionEqn = Percentage_Motion_Windows ~ Age * Task + Age * Cohort + AP_Orientation + (1 | ID)
# model_Motion <- lmer(model_MotionEqn, data = dataExp)
# 
# 
# 
# 
# 
# 
# #### Examine residual normality and fit ####
# # interaction model best choice overall, so this used to look @ residuals
# # extract residuals for overall best fitting model(s)
# residuals_Motion <- resid(model_Motion)
# 
# 
# # plot QQ-plot and histogram modelIntTaskCohort_SCI
# # QQ-plot for modelIntTaskCohort_SCI
# par(mfrow = c(1, 1))  # single plot layout
# qqnorm(residuals_Motion, main = "QQ-Plot for Residuals (model_Motion)")
# qqline(residuals_Motion, col = "red", lwd = 2)
# 
# # Histogram for modelIntTaskCohort_SCI
# par(mfrow = c(1, 1))  
# hist(residuals_Motion, breaks = 100, main = "Histogram of Residuals (model_Motion)", xlab = "Residuals", col = "lightblue")
# 
# # Model extremely non-normal. Examine outliers and influential points.
# 
# #### Examine outlier and influential points. ####
# 
# # Extract residuals for overall best fitting model(s)
# hatValues_Motion <- hatvalues(model_Motion)  # Leverage values
# 
# # Cook's Distance (Influential Points)
# cooksD_Motion <- cooks.distance(model_Motion)
# influentialPoints_Motion <- which(cooksD_Motion > (4 / length(cooksD_Motion)))  # Influential points threshold
# influentialPointsPercentage_Motion <- length(influentialPoints_Motion) / length(cooksD_Motion) * 100
# 
# # Identify outliers: Residuals greater than 3 standard deviations
# outliers_Motion <- which(abs(residuals_Motion) > 3 * sd(residuals_Motion))
# outliersPercentage_Motion <- length(outliers_Motion) / length(residuals_Motion) * 100
# 
# # Identify high leverage points: Leverage greater than 2 * mean leverage
# highLeverage_Motion <- which(hatValues_Motion > 2 * mean(hatValues_Motion))
# highLeveragePercentage_Motion <- length(highLeverage_Motion) / length(hatValues_Motion) * 100
# 
# # Output percentage of influential points, outliers, and high leverage points
# cat("Influential Points (Cook's Distance) Percentage:", round(influentialPointsPercentage_Motion, 2), "%\n")
# cat("Outliers Percentage:", round(outliersPercentage_Motion, 2), "%\n")
# cat("High Leverage Points Percentage:", round(highLeveragePercentage_Motion, 2), "%\n\n")
# 

#### DEFINE FINAL MODEL AND GET RESULTS ####
# Model repeated here for ease in case you don't want to run model comparisons 
model_MotionEqn = Percentage_Motion_Windows ~ Age * Task + Age * Cohort + AP_Orientation + (1 | ID)
model_Motion <- lmer(model_MotionEqn, data = dataExp)

# Extract fixed effects and standard errors
fixed_effects_Motion <- fixef(model_Motion)  # Fixed effects coefficients
se_Motion <- sqrt(diag(vcov(model_Motion)))  # Standard errors for fixed effects

# Calculate t-values for each effect
t_values_Motion <- fixed_effects_Motion / se_Motion
p_values_Motion <- 2 * (1 - pt(abs(t_values_Motion), df = 1000))  # Two-tailed p-value

# Combine fixed effects, t-values, and p-values into data frames
fixed_effects_Motion_df <- data.frame(
  Effect = names(fixed_effects_Motion),
  Estimate = fixed_effects_Motion,
  Std.Error = se_Motion,
  t_value = t_values_Motion,
  p_value = p_values_Motion
)

# View the results
print(fixed_effects_Motion_df)
#### BOOTSTRAPPING ####
# Set up parallel processing
registerDoParallel(numCores)

if (!file.exists("stats/bootResults_Motion.rds")) {
  
  # Time the bootstrapping process
  startTime <- Sys.time()
  
  # Initialise variables
  successfulBootstraps <- 0
  seValues <- numeric(nBoot)  # Storage for SE values
  
  # Perform bootstrap using foreach
  bootResults_Motion <- foreach(i = 1:nBoot, .combine = rbind, .packages = c("lme4")) %dopar% {
    bootstrap_iteration(dataExp, model_MotionEqn, i)
  }
  
  # Calculate the time taken
  endTime <- Sys.time()
  timeTaken <- endTime - startTime
  cat("Time taken for Motion bootstrapping:", timeTaken, "hours\n")
  
  # Save results
  saveRDS(bootResults_Motion, "stats/bootResults_Motion.rds")
  
} else {
  cat("Motion bootstrapping results already exist. Loading instead.\n")
  bootResults_Motion <- readRDS("stats/bootResults_Motion.rds")
}

# Stop parallel processing
stopImplicitCluster()

# Diagnostics for failed bootstraps
cat("\nNumber of failed bootstraps:", sum(rowSums(is.na(bootResults_Motion)) > 0), "\n")


#### ---- Display model summary ####
# Calculate confidence intervals for model_Motion
nCoefMotion <- length(fixef(model_Motion))
motionCoefNames <- names(fixef(model_Motion))
motionCIs <- t(sapply(1:nCoefMotion, function(i) {
  quantile(bootResults_Motion[, i], probs = c(0.025, 0.975), na.rm = TRUE)
}))
rownames(motionCIs) <- motionCoefNames

# Calculate standard errors from bootstrapped estimates
motionSEs <- apply(bootResults_Motion[, 1:nCoefMotion], 2, sd, na.rm = TRUE)

# Create summary table for model_Motion
motionSummary <- data.frame(
  Estimate = colMeans(bootResults_Motion),
  CI_lower = motionCIs[, "2.5%"],
  CI_upper = motionCIs[, "97.5%"],
  SE = motionSEs
)

# Print results
cat("\nModel Motion Results:\n")
print(motionSummary)

saveRDS(motionSummary, "stats/Motion_Summary.rds")


# Calculate predicted values for the interaction effects (Age * Task, Age * cohort, Age * AP_Orientation)
# Get unique values for all 4 factors
predictedData_Motion <- expand.grid(
  Age = unique(dataExp$Age),  
  Task = unique(dataExp$Task), 
  Cohort = unique(dataExp$Cohort), 
  AP_Orientation = unique(dataExp$AP_Orientation) 
)

# Get the predicted values and the standard errors
predictedData_Motion$predictedValues <- predict(model_Motion, newdata = predictedData_Motion, re.form = NA)
predictedData_Motion$se <- motionSEs[1]  # Standard error (adjust if needed for interaction terms)

# Calculate the confidence intervals for the predictions
predictedData_Motion$lowerCi <- predictedData_Motion$predictedValues - 1.96 * predictedData_Motion$se
predictedData_Motion$upperCi <- predictedData_Motion$predictedValues + 1.96 * predictedData_Motion$se


#### ---- Calculate effect sizes ####

Motion_Effects = standardiseEffects(model_Motion, motionSummary)

saveRDS(Motion_Effects, "stats/Motion_Effects.rds")
write.csv(Motion_Effects, "stats/Motion_Effects.csv")

#caluclate and save fit
# calculatemodel fit for SNR
Motion_Fit <- r.squaredGLMM(model_Motion)
# save
saveRDS(Motion_Fit, "stats/Motion_Fit.rds")


#### PLOT TRENDS ####
#### ---- Re-define model and equation if necessary ####
# Model repeated here for ease in case you don't want to run model comparisons 
model_MotionEqn = Percentage_Motion_Windows ~ Age * Task + Age * Cohort + AP_Orientation + (1 | ID)
model_Motion <- lmer(model_MotionEqn, data = dataExp)
#### ---- Load variables if necessary ####
# Check and load outputs if they are not already in the environment
if (!exists("bootResults_Motion")) {
  bootResults_Motion <- readRDS("stats/bootResults_Motion.rds")
}

if (!exists("Motion_Summary")) {
  Motion_Summary <- readRDS("stats/Motion_Summary.rds")
}

if (!exists("Motion_Effects")) {
  Motion_Effects <- readRDS("stats/Motion_Effects.rds")
}

if (!exists("Motion_Fit")) {
  Motion_Fit <- readRDS("stats/Motion_Fit.rds")
}

#### ---- Forest plots ####

# New column names for bootResults_Motion for plots
new_colnames <- c("Intercept", "Age", "Task (Social)", "Cohort (UK)", "Channel Location", 
                  "Age * Task (Social)", "Age * Cohort (UK)")

# Assign the new column names to bootResults_Motion
colnames(bootResults_Motion) <- new_colnames

# Prepare data for forest plot of Model_Motion
Motion_ForestPlot_Data <- data.frame(
  Coefficient = colnames(bootResults_Motion),
  Estimate = Motion_Summary$Estimate,
  CI_lower = Motion_Summary$CI_lower,
  CI_upper = Motion_Summary$CI_upper,
  SE = Motion_Summary$SE,
  Effect_Size = Motion_Effects$Effect_Size
)

# Create a new column for effect type based on size and direction
Motion_ForestPlot_Data <- Motion_ForestPlot_Data %>%
  mutate(Effect_Type = case_when(
    CI_lower < 0 & CI_upper > 0 ~ "Not Significant",
    CI_lower > 0 & CI_upper < 0 ~ "Not Significant",
    Effect_Size == "Large" & Estimate > 0 ~ "Large_Positive",
    Effect_Size == "Medium" & Estimate > 0 ~ "Medium_Positive",
    Effect_Size == "Small" & Estimate > 0 ~ "Small_Positive",
    Effect_Size == "Large" & Estimate < 0 ~ "Large_Negative",
    Effect_Size == "Medium" & Estimate < 0 ~ "Medium_Negative",
    Effect_Size == "Small" & Estimate < 0 ~ "Small_Negative"
  ))

# Define color and shape mappings
color_mapping <- c(
  "Not Significant" = "red",
  "Large_Positive" = "darkgreen",
  "Medium_Positive" = "dodgerblue",
  "Small_Positive" = "orange",
  "Large_Negative" = "chartreuse2",
  "Medium_Negative" = "blue3",
  "Small_Negative" = "darkorange3"
)

shape_mapping <- c(
  "Not Significant" = 4, # cross
  "Large_Positive" = 15,  # Square
  "Medium_Positive" = 19, # Circle
  "Small_Positive" = 17,  # Triangle
  "Large_Negative" = 0,   # 
  "Medium_Negative" = 1,  # 
  "Small_Negative" = 2    #
)

# Generate forest plots for Model_Motion

# With intercept
ggplot(Motion_ForestPlot_Data, aes(x = Estimate, y = Coefficient)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray", size = 1) + 
  geom_point(aes(color = Effect_Type, shape = Effect_Type), size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.4) +
  scale_color_manual(values = color_mapping) +
  scale_shape_manual(values = shape_mapping) +
  labs(title = "Bootstrapped Fixed Effects for Percentage of Motion",
       x = "Estimate with 95% CI", y = "Coefficient",
       color = "Standardised Effect Type", shape = "Standardised Effect Type") +
  theme_minimal()

#### ---- ---- Without Intercept ####
# Remove intercept from the SCI data
Motion_ForestPlot_Data_no_intercept <- Motion_ForestPlot_Data %>%
  filter(Coefficient != "Intercept")

ggplot(Motion_ForestPlot_Data_no_intercept, aes(x = Estimate, y = Coefficient)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray", size = 1) + 
  geom_point(aes(color = Effect_Type, shape = Effect_Type), size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.4) +
  scale_color_manual(values = color_mapping) +
  scale_shape_manual(values = shape_mapping) +
  labs(title = "Bootstrapped Fixed Effects for Percentage of Motion",
       x = "Estimate with 95% CI", y = "Coefficient",
       color = "Standardised Effect Type", shape = "Standardised Effect Type") +
  theme_minimal()

#### ---- Trends, effects and interaction effects ####
#### ---- ---- Generate new DF ####
# Create a new data frame for predictions
avgFixef_Motion <- colMeans(bootResults_Motion, na.rm = TRUE)

# Create a grid of combinations for Age, Task, Cohort, etc.
predictedData_Motion <- expand.grid(
  Age = unique(dataExp$Age),
  Task = unique(dataExp$Task),
  Cohort = unique(dataExp$Cohort),
  AP_Orientation = unique(dataExp$AP_Orientation)
)

# Add intercept to design matrix
X <- model.matrix(~ Age * Task + Age * Cohort + AP_Orientation,
                  data = predictedData_Motion)

# Compute predictions using the average fixed effects
predictedData_Motion$predictedValues <- X %*% avgFixef_Motion

# Calculate predictions for all bootstrap iterations to assess uncertainty
bootstrapPredictions_Motion <- apply(bootResults_Motion, 1, function(coefs) X %*% coefs)
predictedData_Motion$LowerCI <- apply(bootstrapPredictions_Motion, 1, quantile, probs = 0.025, na.rm = TRUE)
predictedData_Motion$UpperCI <- apply(bootstrapPredictions_Motion, 1, quantile, probs = 0.975, na.rm = TRUE)


#### ---- ---- Effect of Age on Motion ####
# Aggregate data by Age
predictedData_Motion_summary <- predictedData_Motion %>%
  group_by(Age) %>%
  summarize(predictedValues = mean(predictedValues), 
            LowerCI = mean(LowerCI), 
            UpperCI = mean(UpperCI), .groups = 'drop')

# Plot the effect of Age on the predicted values
ggplot(predictedData_Motion, aes(x = Age, y = predictedValues)) +
  geom_point(color = viridis(3)[1], size = 3, alpha = 0.8) +  
  geom_smooth(method = "lm", color = viridis(3)[2], fill = viridis(3)[2], alpha = 0.2) +
  labs(title = "Predicted Motion Trends by Age (Scatter and Regression)", 
       x = "Age", 
       y = "Predicted Motion %") +
  theme_minimal() +  
  theme(
    panel.grid.major = element_line(color = "grey80"),  
    panel.grid.minor = element_line(color = "grey90"), 
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA)
  )

#### ---- ---- Effect of Cohort on Motion ####
# Violin plot to visualise the distribution of predicted values by Cohort
ggplot(predictedData_Motion, aes(x = Cohort, y = predictedValues, fill = Cohort)) +
  geom_violin(trim = FALSE, alpha = 0.7) +  
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +  
  scale_fill_viridis(discrete = TRUE) +  
  scale_x_discrete(labels = c("gm" = "Gambia", "uk" = "UK")) +  
  labs(title = "Effect of Cohort on Predicted Motion", 
       x = "Cohort", 
       y = "Predicted Motion %") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "grey80"),  
    panel.grid.minor = element_line(color = "grey90"), 
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA)
  )


#### ---- ---- Effect of Age * Task on Motion ####
# Plot the effect of Age and Task on the predicted values
ggplot(predictedData_Motion, aes(x = Age, y = predictedValues, color = Task, shape = Task)) +
  geom_point(size = 3, alpha = 0.8, position = position_dodge(width = 0.5)) +  
  geom_smooth(method = "lm", aes(fill = Task), alpha = 0.2, linetype = "solid", show.legend = FALSE) +  
  scale_color_viridis(discrete = TRUE, labels = c("uk" = "UK", "gm" = "Gambia", "hand" = "HaND", "social" = "SNS")) +  
  scale_fill_viridis(discrete = TRUE) +  # Use viridis for CI fills
  scale_shape_manual(values = c("uk" = 16, "gm" = 17, "hand" = 18, "social" = 19), 
                     labels = c("uk" = "UK", "gm" = "Gambia", "hand" = "HaND", "social" = "SNS")) + 
  labs(title = "Predicted Motion Trends by Age and Task",
       x = "Age", 
       y = "Predicted Motion %") +
  theme_minimal() +  
  theme(
    panel.grid.major = element_line(color = "grey80"),  
    panel.grid.minor = element_line(color = "grey90"), 
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank(),  
    legend.position = "right"  
  )

#### ---- ---- Effect of Age * Cohort on Motion ####
# Plot the effect of Age and Cohort on the predicted values
ggplot(predictedData_Motion, aes(x = Age, y = predictedValues, color = Cohort, shape = Cohort)) +
  geom_point(size = 3, alpha = 0.8, position = position_dodge(width = 0.5)) +  
  geom_smooth(method = "lm", aes(fill = Cohort), alpha = 0.2, linetype = "solid", show.legend = FALSE) +  
  scale_color_viridis(discrete = TRUE, labels = c("uk" = "UK", "gm" = "Gambia", "hand" = "HaND", "social" = "SNS")) +  
  scale_fill_viridis(discrete = TRUE) +  # Use viridis
  scale_shape_manual(values = c("uk" = 16, "gm" = 17, "hand" = 18, "social" = 19), 
                     labels = c("uk" = "UK", "gm" = "Gambia", "hand" = "HaND", "social" = "SNS")) +  
  labs(title = "Predicted Motion Trends by Age and Cohort",
       x = "Age (months)", 
       y = "Predicted Motion %") +
  theme_minimal() +  
  theme(
    panel.grid.major = element_line(color = "grey80"),  
    panel.grid.minor = element_line(color = "grey90"), 
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank(), 
    legend.position = "right"  
  )

#### ---- ---- Effect of Age, Task, and Cohort ####
# Plot for Percentage_Motion_Windows with faceting by Cohort
# ggplot(predictedData_Motion, aes(x = Age, y = predictedValues, color = Cohort, shape = Task, alpha = Task)) +
#   geom_point(size = 3, position = position_dodge(width = 0.2), show.legend = c(shape = TRUE, color = TRUE)) +  
#   geom_smooth(method = "lm", aes(fill = Cohort), linetype = "solid", show.legend = FALSE, alpha = 0.3) +  
#   scale_color_viridis(discrete = TRUE, option = "D", labels = c("uk" = "UK", "gm" = "Gambia")) +  
#   scale_fill_viridis(discrete = TRUE, option = "D", alpha = 0.3) +  
#   scale_alpha_manual(values = c("hand" = 1, "social" = 0.4), guide = "none") +  
#   scale_shape_manual(values = c("hand" = 16, "social" = 17), labels = c("hand" = "HaND", "social" = "SNS")) +  
#   labs(title = "Predicted Percentage of Motion Windows by Age, Task, and Cohort", 
#        x = "Age", 
#        y = "Predicted Motion %") +
#   theme_minimal() +  
#   theme(
#     panel.grid.major = element_line(color = "grey80"),  
#     panel.grid.minor = element_line(color = "grey90"), 
#     panel.background = element_rect(fill = "white", color = NA),  
#     plot.background = element_rect(fill = "white", color = NA),
#     strip.text = element_text(size = 14, face = "bold")  
#   ) +
#   scale_y_continuous(limits = c(-0.75, 0.75)) +
#   guides(shape = guide_legend(title = "Task", override.aes = list(alpha = 1)), 
#          color = "none") +  
#   facet_wrap(~ Cohort, labeller = labeller(Cohort = c("uk" = "UK", "gm" = "Gambia")))  


ggplot(predictedData_Motion, aes(x = Age, y = predictedValues, shape = Task)) +
  geom_point(aes(fill = Cohort, color = Cohort), 
             size = 3, stroke = 0.4, position = position_dodge(width = 0.2), 
             show.legend = c(shape = TRUE, fill = TRUE, color = FALSE)) +  
  
  geom_smooth(method = "lm", aes(fill = Cohort), linetype = "solid", 
              color = NA, show.legend = FALSE, alpha = 0.6) +  # Increased opacity of trendline
  
  # Manually set colors for fill and outline
  scale_fill_manual(values = c("uk" = "#E6AC00", "gm" = "#440154")) +  # UK = lighter yellow, GM = dark blue
  scale_color_manual(values = c("uk" = "#B8860B", "gm" = "white")) +  # Maintains outline color for UK points
  
  scale_alpha_manual(values = c("hand" = 1, "social" = 1), guide = "none") +  # Ensures SNS points are fully opaque
  scale_shape_manual(values = c("hand" = 21, "social" = 24),  # 21 = circle, 24 = triangle
                     labels = c("hand" = "HaND", "social" = "SNS")) +  
  
  labs(title = "Predicted Percentage of Motion Windows by Age, Task, and Cohort", 
       x = "Age", 
       y = "Predicted Motion %") +
  
  theme_minimal() +  
  theme(
    panel.grid.major = element_line(color = "grey80"),  
    panel.grid.minor = element_line(color = "grey90"), 
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 14, face = "bold")  
  ) +
  
  scale_y_continuous(limits = c(-0.75, 0.75)) +
  guides(shape = guide_legend(title = "Task", override.aes = list(alpha = 1)), 
         fill = guide_legend(title = "Cohort")) +  
  
  facet_wrap(~ Cohort, labeller = labeller(Cohort = c("uk" = "UK", "gm" = "Gambia")))  


