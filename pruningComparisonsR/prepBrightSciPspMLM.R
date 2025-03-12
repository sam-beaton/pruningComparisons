#### LOAD PACKAGES ####
# Load necessary libraries
library(lme4)
library(ggplot2)
library(e1071)
# for bootstrapping:
library(foreach)
library(doParallel)
library(dplyr)
library(effectsize)
library(MuMIn)
library(viridis)
library(Matrix)
library(tidyr)
library(gridExtra) #to arrange plots

#### SET PARAMS ####
#set working directory (for saving outputs - bootstrapping is time-consuming!)
setwd("[path-to...]/pruningComparisonsR")

# load the data
dataExp <- read.csv("[path-to...]/pruningComparisonsMatlab/stats/bright/overall/pruneChannelInfoTableQT.csv")

#clone for plotting
dataExpOrig <- dataExp

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

#check first few rows
head(dataExp)
tail(dataExp)
summary(dataExp)

# Define the list of factors
factors <- c("Age", "Task", "Cohort", "Homer_Pruned_Total", "Percentage_Motion_Windows", "AP_Orientation")

# Colourblind-friendly palette 
cbPalette <- viridis(length(unique(dataExp$Task)))

# parallel processing param(s)
numCores <- detectCores() - 1  

# Bootstrapping inference params
# for reproducibility
set.seed(42)
# Set number of bootstrap iterations
nBoot <- 5000
# nBoot <- 250 #for testing code - comment out!
# # divide and conquer:
# chunkSize <- 100  # Number of iterations per chunk
# chunkSize <- 10 #for testing code - comment out!
# #early stopping:
# stabilityThreshold <- 0.00001   # Stability threshold (max SE change for fixed effects)
 
# Bootstrapping DF generation param(s)
chunkSizeDF <- 40 # Number of rows per chunk


#### DEFINE HELPER FUNCTIONS ####
# Function to create model formulas given a set of predictors and outcome variable
create_model <- function(outcome, predictors) {
  formula <- as.formula(paste(outcome, "~", paste(predictors, collapse = " + "), "+ (1 | ID)"))
  model <- lmer(formula, data = dataExp)
  return(model)
}

# Function to compute predicted values, standard errors, and confidence intervals from fitted models
computePredictions <- function(model, data, fixedEffectFormula) {
  vcovMatrix <- vcov(model)  # Variance-covariance matrix
  modelFormula <- as.formula(paste("~", fixedEffectFormula))  
  modelMatrix <- model.matrix(modelFormula, data)  
  
  # Get the rows used in the model (data used for fitting)
  modelData <- model.frame(model)  
  filteredData <- data[rownames(modelData), ]  
  
  # Predictions and standard errors
  predicted <- predict(model, re.form = NA)
  se <- sqrt(rowSums((modelMatrix %*% vcovMatrix)^2))
  
  # Add predictions and CI to the filtered data
  filteredData <- filteredData %>%
    mutate(
      predicted = predicted,
      lowerCi = predicted - 1.96 * se,
      upperCi = predicted + 1.96 * se
    )
  
  return(filteredData)
}

# Helper function for bootstrapping
bootstrap_iteration <- function(data, model_eqn, i) {
  # Resample data with replacement
  indices <- sample(nrow(data), replace = TRUE)
  bootData <- data[indices, ]
  
  # Fit the model
  bootModel <- try(lmer(model_eqn, data = bootData), silent = TRUE)
  
  # Extract fixed effects
  if (!inherits(bootModel, "try-error")) {
    fixefs <- fixef(bootModel)
    seFixefs <- sqrt(diag(vcov(bootModel)))

    return(fixefs)
  } else {
    return(rep(NA, length(fixef(model_Motion))))  # Return NA for errors
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
  
  # Initialise results dataframe
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
      # Check if it's a main effect or interaction
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


#### DETERMINE MODEL ####
#### ---- Initial models (no interaction, non-linearity etc.) ####
# Initialise an empty list to store models
models <- list()

# Define letters for naming models (a, b, c, ...)
letters_list <- letters

# Loop over the outcomes
outcomes <- c("Average_SCI", "Average_PSP")

# Loop over the number of factors in each model
for (outcome in outcomes) {
  for (n_factors in 1:6) {
    # Generate combinations of factors
    factor_combinations <- combn(factors, n_factors, simplify = FALSE)

    # Create models for each combination of factors
    for (i in seq_along(factor_combinations)) {
      predictors <- factor_combinations[[i]]
      model_name <- paste0("model", n_factors, letters_list[i], "_", sub("Average_", "", outcome)) # e.g., model3a_SCI, model4b_PSP
      models[[model_name]] <- create_model(outcome, predictors)
    }
  }
}


# Compare all models
# Collect model summaries by separating SCI and PSP models
modelSummariesSCI <- models[grep("_SCI$", names(models))]
modelSummariesPSP <- models[grep("_PSP$", names(models))]

# Extract AIC and BIC for comparison for both sets
aicBicSCI <- sapply(modelSummariesSCI, function(model) c(AIC(model), BIC(model)))
aicBicPSP <- sapply(modelSummariesPSP, function(model) c(AIC(model), BIC(model)))

# Combine results into data frame
resultsSCI <- data.frame(
  Model = names(modelSummariesSCI),
  AIC = aicBicSCI[1, ],
  BIC = aicBicSCI[2, ]
)

resultsPSP <- data.frame(
  Model = names(modelSummariesPSP),
  AIC = aicBicPSP[1, ],
  BIC = aicBicPSP[2, ]
)

# Sort and display results from worst to best model by AIC
# Sort SCI results by AIC (descending order for worst to best)
resultsSCI <- resultsSCI[order(-resultsSCI$AIC), ]
# Sort PSP results by AIC (descending order for worst to best)
resultsPSP <- resultsPSP[order(-resultsPSP$AIC), ]
# Display (AIC) sorted results for Average_SCI and Average_PSP models
print("\nResults for Average_SCI Models (Worst to Best by AIC):")
print(resultsSCI)
print("\nResults for Average_PSP Models (Worst to Best by AIC):")
print(resultsPSP)

# Now by BIC
# Sort SCI results by BIC (descending order)
resultsSCI <- resultsSCI[order(-resultsSCI$BIC), ]
# Sort PSP results by BIC (descending)
resultsPSP <- resultsPSP[order(-resultsPSP$BIC), ]
# Display (AIC) sorted results for Average_SCI and Average_PSP models
print("\nResults for Average_SCI Models (Worst to Best by BIC):")
print(resultsSCI)
print("\nResults for Average_PSP Models (Worst to Best by BIC):")
print(resultsPSP)


# Function to extract and print model equations
print_model_equations <- function(model_names, models) {
  for (model_name in model_names) {
    model <- models[[model_name]]
    cat("Equation for", model_name, ":\n")
    print(formula(model))  
    cat("\n")
  }
}

# Select the 5 best models based on AIC for each outcome
top5ModelsSCI <- resultsSCI[order(resultsSCI$AIC), ][1:5, ]
top5ModelsPSP <- resultsPSP[order(resultsPSP$AIC), ][1:5, ]
# Print the 5 best model equations for Average_SCI
cat("\nTop 5 Model Equations for Average_SCI according to AIC:\n")
print_model_equations(top5ModelsSCI$Model, models)
# Print the 5 best model equations for Average_PSP
cat("\nTop 5 Model Equations for Average_PSP according to AIC:\n")
print_model_equations(top5ModelsPSP$Model, models)

# Select the 5 best models based on BIC for each outcome
top5ModelsSCI <- resultsSCI[order(resultsSCI$BIC), ][1:5, ]
top5ModelsPSP <- resultsPSP[order(resultsPSP$BIC), ][1:5, ]
# Print the 5 best model equations for Average_SCI
cat("\nTop 5 Model Equations for Average_SCI according to BIC:\n")
print_model_equations(top5ModelsSCI$Model, models)
# Print the 5 best model equations for Average_PSP
cat("\nTop 5 Model Equations for Average_PSP according to BIC:\n")
print_model_equations(top5ModelsPSP$Model, models)


# Model5e seems to be best performing overall, however based on the examination
# of models leading to the eventual final model of
# modelMotion <- lmer(Percentage_Motion_Windows ~ Age * Task + (1 | ID), data = dataExp)
# in
# 'prepBrightMotionMLM.R'
# Task is included and model6a chosen for both outcomes

# Now compare model6a_SCI and model6a_PSP to similar models which utilise
# interaction terms of Age * Task, Age * Cohort, and Age * AP_Orientation again
# based on modelMotion, above, and for similar reasons.
# Also add Age * Percentage_Motion_Windows, based on outcome of modelMotion,
# which showed a change in motion with age.

# Add the new models  to the 'models' list
modelIntTask_SCI <- lmer(Average_SCI ~ Age * Task + Cohort + Homer_Pruned_Total + Percentage_Motion_Windows + AP_Orientation + (1 | ID), data = dataExp)
modelIntTask_PSP <- lmer(Average_PSP ~ Age * Task + Cohort + Homer_Pruned_Total + Percentage_Motion_Windows + AP_Orientation + (1 | ID), data = dataExp)
modelIntCohort_SCI <- lmer(Average_SCI ~ Age * Cohort + Task + Homer_Pruned_Total + Percentage_Motion_Windows + AP_Orientation + (1 | ID), data = dataExp)
modelIntCohort_PSP <- lmer(Average_PSP ~ Age * Cohort + Task + Homer_Pruned_Total + Percentage_Motion_Windows + AP_Orientation + (1 | ID), data = dataExp)
modelIntOrientation_SCI <- lmer(Average_SCI ~ Age * AP_Orientation + Task + Homer_Pruned_Total + Percentage_Motion_Windows + Cohort + (1 | ID), data = dataExp)
modelIntOrientation_PSP <- lmer(Average_PSP ~ Age * AP_Orientation + Task + Homer_Pruned_Total + Percentage_Motion_Windows + Cohort + (1 | ID), data = dataExp)
modelIntMotion_SCI <- lmer(Average_SCI ~ Age * Percentage_Motion_Windows + Task + Cohort + Homer_Pruned_Total + AP_Orientation + (1 | ID), data = dataExp)
modelIntMotion_PSP <- lmer(Average_PSP ~ Age * Percentage_Motion_Windows + Task + Cohort + Homer_Pruned_Total + AP_Orientation + (1 | ID), data = dataExp)


modelIntTaskCohort_SCI <- lmer(Average_SCI ~ Age * Cohort + Age * Task + Homer_Pruned_Total + Percentage_Motion_Windows + AP_Orientation + (1 | ID), data = dataExp)
modelIntTaskCohort_PSP <- lmer(Average_PSP ~ Age * Cohort + Age * Task + Homer_Pruned_Total + Percentage_Motion_Windows + AP_Orientation + (1 | ID), data = dataExp)
modelIntTaskOrientation_SCI <- lmer(Average_SCI ~ Age * AP_Orientation + Age * Task + Homer_Pruned_Total + Percentage_Motion_Windows + Cohort + (1 | ID), data = dataExp)
modelIntTaskOrientation_PSP <- lmer(Average_PSP ~ Age * AP_Orientation + Age * Task + Homer_Pruned_Total + Percentage_Motion_Windows + Cohort + (1 | ID), data = dataExp)
modelIntTaskMotion_SCI <- lmer(Average_SCI ~ Age * Task + Age * Percentage_Motion_Windows + Cohort + Homer_Pruned_Total + AP_Orientation + (1 | ID), data = dataExp)
modelIntTaskMotion_PSP <- lmer(Average_PSP ~ Age * Task + Age * Percentage_Motion_Windows + Cohort + Homer_Pruned_Total + AP_Orientation + (1 | ID), data = dataExp)
modelIntCohortOrientation_SCI <- lmer(Average_SCI ~ Age * Cohort + Age * AP_Orientation + Homer_Pruned_Total + Percentage_Motion_Windows + Task + (1 | ID), data = dataExp)
modelIntCohortOrientation_PSP <- lmer(Average_PSP ~ Age * Cohort + Age * AP_Orientation + Homer_Pruned_Total + Percentage_Motion_Windows + Task + (1 | ID), data = dataExp)
modelIntCohortMotion_SCI <- lmer(Average_SCI ~ Age * Cohort + Age * Percentage_Motion_Windows + Task + Homer_Pruned_Total + AP_Orientation + (1 | ID), data = dataExp)
modelIntCohortMotion_PSP <- lmer(Average_PSP ~ Age * Cohort + Age * Percentage_Motion_Windows + Task + Homer_Pruned_Total + AP_Orientation + (1 | ID), data = dataExp)
modelIntOrientationMotion_SCI <- lmer(Average_SCI ~ Age * AP_Orientation + Age * Percentage_Motion_Windows + Task + Cohort + Homer_Pruned_Total + (1 | ID), data = dataExp)
modelIntOrientationMotion_PSP <- lmer(Average_PSP ~ Age * AP_Orientation + Age * Percentage_Motion_Windows + Task + Cohort + Homer_Pruned_Total + (1 | ID), data = dataExp)

modelIntTaskCohortOrientation_SCI <- lmer(Average_SCI ~ Age * Cohort + Age * Task + Age * AP_Orientation + Homer_Pruned_Total + Percentage_Motion_Windows + AP_Orientation + (1 | ID), data = dataExp)
modelIntTaskCohortOrientation_PSP <- lmer(Average_PSP ~ Age * Cohort + Age * Task + Age * AP_Orientation + Homer_Pruned_Total + Percentage_Motion_Windows + AP_Orientation + (1 | ID), data = dataExp)
modelIntTaskCohortMotion_SCI <- lmer(Average_SCI ~ Age * Cohort + Age * Task + Age * Percentage_Motion_Windows + Homer_Pruned_Total + AP_Orientation + (1 | ID), data = dataExp)
modelIntTaskCohortMotion_PSP <- lmer(Average_PSP ~ Age * Cohort + Age * Task + Age * Percentage_Motion_Windows + Homer_Pruned_Total + AP_Orientation + (1 | ID), data = dataExp)
modelIntTaskOrientationMotion_SCI <- lmer(Average_SCI ~ Age * AP_Orientation + Age * Task + Age * Percentage_Motion_Windows + Homer_Pruned_Total + Cohort + (1 | ID), data = dataExp)
modelIntTaskOrientationMotion_PSP <- lmer(Average_PSP ~ Age * AP_Orientation + Age * Task + Age * Percentage_Motion_Windows + Homer_Pruned_Total + Cohort + (1 | ID), data = dataExp)
modelIntCohortOrientationMotion_SCI <- lmer(Average_SCI ~ Age * Cohort + Age * AP_Orientation + Age * Percentage_Motion_Windows + Homer_Pruned_Total + Task + (1 | ID), data = dataExp)
modelIntCohortOrientationMotion_PSP <- lmer(Average_PSP ~ Age * Cohort + Age * AP_Orientation + Age * Percentage_Motion_Windows + Homer_Pruned_Total + Task + (1 | ID), data = dataExp)

modelIntAll_SCI <- lmer(Average_SCI ~ Age * Cohort + Age * Task + Age * AP_Orientation + Age * Percentage_Motion_Windows + Homer_Pruned_Total + (1 | ID), data = dataExp)
modelIntAll_PSP <- lmer(Average_PSP ~ Age * Cohort + Age * Task + Age * AP_Orientation + Age * Percentage_Motion_Windows + Homer_Pruned_Total + (1 | ID), data = dataExp)

#### ---- Consider interaction models ####
# Add these models to the existing 'models' list (which already contains the other models)
models[["modelIntTask_SCI"]] <- modelIntTask_SCI
models[["modelIntTask_PSP"]] <- modelIntTask_PSP
models[["modelIntCohort_SCI"]] <- modelIntCohort_SCI
models[["modelIntCohort_PSP"]] <- modelIntCohort_PSP
models[["modelIntOrientation_SCI"]] <- modelIntOrientation_SCI
models[["modelIntOrientation_PSP"]] <- modelIntOrientation_PSP
models[["modelIntMotion_SCI"]] <- modelIntMotion_SCI
models[["modelIntMotion_PSP"]] <- modelIntMotion_PSP

models[["modelIntTaskCohort_SCI"]] <- modelIntTaskCohort_SCI
models[["modelIntTaskCohort_PSP"]] <- modelIntTaskCohort_PSP
models[["modelIntTaskOrientation_SCI"]] <- modelIntTaskOrientation_SCI
models[["modelIntTaskOrientation_PSP"]] <- modelIntTaskOrientation_PSP
models[["modelIntTaskMotion_SCI"]] <- modelIntTaskMotion_SCI
models[["modelIntTaskMotion_PSP"]] <- modelIntTaskMotion_PSP
models[["modelIntCohortOrientation_SCI"]] <- modelIntCohortOrientation_SCI
models[["modelIntCohortOrientation_PSP"]] <- modelIntCohortOrientation_PSP
models[["modelIntCohortMotion_SCI"]] <- modelIntCohortMotion_SCI
models[["modelIntCohortMotion_PSP"]] <- modelIntCohortMotion_PSP
models[["modelIntOrientationMotion_SCI"]] <- modelIntOrientationMotion_SCI
models[["modelIntOrientationMotion_PSP"]] <- modelIntOrientationMotion_PSP

models[["modelIntTaskCohortOrientation_SCI"]] <- modelIntTaskCohortOrientation_SCI
models[["modelIntTaskCohortOrientation_PSP"]] <- modelIntTaskCohortOrientation_PSP
models[["modelIntTaskCohortMotion_SCI"]] <- modelIntTaskCohortMotion_SCI
models[["modelIntTaskCohortMotion_PSP"]] <- modelIntTaskCohortMotion_PSP
models[["modelIntTaskOrientationMotion_SCI"]] <- modelIntTaskOrientationMotion_SCI
models[["modelIntTaskOrientationMotion_PSP"]] <- modelIntTaskOrientationMotion_PSP
models[["modelIntCohortOrientationMotion_SCI"]] <- modelIntCohortOrientationMotion_SCI
models[["modelIntCohortOrientationMotion_PSP"]] <- modelIntCohortOrientationMotion_PSP

models[["modelIntAll_SCI"]] <- modelIntAll_SCI
models[["modelIntAll_PSP"]] <- modelIntAll_PSP

# Compare all models again, including the new ones
# Collect model summaries for both SCI and PSP models
modelSummariesSCI <- models[grep("_SCI$", names(models))]
modelSummariesPSP <- models[grep("_PSP$", names(models))]

# Extract AIC and BIC for comparison for both sets
aicBicSCI <- sapply(modelSummariesSCI, function(model) c(AIC(model), BIC(model)))
aicBicPSP <- sapply(modelSummariesPSP, function(model) c(AIC(model), BIC(model)))

# Combine results into data frames for easy viewing
resultsSCI <- data.frame(
  Model = names(modelSummariesSCI),
  AIC = aicBicSCI[1, ],
  BIC = aicBicSCI[2, ]
)

resultsPSP <- data.frame(
  Model = names(modelSummariesPSP),
  AIC = aicBicPSP[1, ],
  BIC = aicBicPSP[2, ]
)

# Sort and display results from worst to best model by AIC
resultsSCI <- resultsSCI[order(resultsSCI$AIC), ]
resultsPSP <- resultsPSP[order(resultsPSP$AIC), ]
cat("\nResults for Average_SCI Models (Worst to Best by AIC):\n")
print(resultsSCI)
cat("\nResults for Average_PSP Models (Worst to Best by AIC):\n")
print(resultsPSP)

# Now by BIC
resultsSCI <- resultsSCI[order(resultsSCI$BIC), ]
resultsPSP <- resultsPSP[order(resultsPSP$BIC), ]
cat("\nResults for Average_SCI Models (Worst to Best by BIC):\n")
print(resultsSCI)
cat("\nResults for Average_PSP Models (Worst to Best by BIC):\n")
print(resultsPSP)

# Function to extract and print model equations
print_model_equations <- function(model_names, models) {
  for (model_name in model_names) {
    model <- models[[model_name]]
    cat("Equation for", model_name, ":\n")
    print(formula(model))  
    cat("\n")
  }
}

# Select the 5 best models based on AIC for each outcome
top5ModelsSCI <- resultsSCI[order(resultsSCI$AIC), ][1:5, ]
top5ModelsPSP <- resultsPSP[order(resultsPSP$AIC), ][1:5, ]
cat("\nTop 5 Model Equations for Average_SCI according to AIC:\n")
print_model_equations(top5ModelsSCI$Model, models)
cat("\nTop 5 Model Equations for Average_PSP according to AIC:\n")
print_model_equations(top5ModelsPSP$Model, models)

# Select the 5 best models based on BIC for each outcome
top5ModelsSCI <- resultsSCI[order(resultsSCI$BIC), ][1:5, ]
top5ModelsPSP <- resultsPSP[order(resultsPSP$BIC), ][1:5, ]
cat("\nTop 5 Model Equations for Average_SCI according to BIC:\n")
print_model_equations(top5ModelsSCI$Model, models)
cat("\nTop 5 Model Equations for Average_PSP according to BIC:\n")
print_model_equations(top5ModelsPSP$Model, models)




# interaction model best choice overall, so this used to look @ residuals
# extract residuals for overall best fitting model(s)
residuals_SCI <- resid(models[["modelIntTaskCohort_SCI"]])
residuals_PSP <- resid(models[["modelIntTaskCohort_PSP"]])

# plot QQ-plot and histogram for modelIntTaskCohort_SCI
# QQ-plot for modelIntTaskCohort_SCI
par(mfrow = c(1, 1))  
qqnorm(residuals_SCI, main = "QQ-Plot for Residuals (modelIntTaskCohort_SCI)")
qqline(residuals_SCI, col = "red", lwd = 2)

# Histogram for modelIntTaskCohort_SCI
par(mfrow = c(1, 1))
hist(residuals_SCI, breaks = 20, main = "Histogram of Residuals (modelIntTaskCohort_SCI)", xlab = "Residuals", col = "lightblue")

# plot QQ-plot and histogram for modelIntTaskCohort_PSP
# QQ-plot for modelIntTaskCohort_PSP
par(mfrow = c(1, 1))
qqnorm(residuals_PSP, main = "QQ-Plot for Residuals (modelIntTaskCohort_PSP)")
qqline(residuals_PSP, col = "red", lwd = 2)

# Histogram for modelIntTaskCohort_PSP
par(mfrow = c(1, 1))
hist(residuals_PSP, breaks = 20, main = "Histogram of Residuals (modelIntTaskCohort_PSP)", xlab = "Residuals", col = "lightblue")


# Residuals are left-skewed in a manner very similar to the SCI and PSP parameter models
# => run bootstrapping, as likely that best case scenario is a Box-Cox transformation
# which results in *slightly* less non-normal residuals

#### FINAL MODELS ####
#### ---- Redefine final models and fit them####
# Define the model equations
model_SCI_Eqn <- Average_SCI ~ Age * Cohort + Age * Task + Age * AP_Orientation + Age * Percentage_Motion_Windows + Homer_Pruned_Total + (1 | ID)
model_PSP_Eqn <- Average_PSP ~ Age * Cohort + Age * Task + Age * AP_Orientation + Age * Percentage_Motion_Windows + Homer_Pruned_Total + (1 | ID)
# Fit the models
model_SCI <- lmer(model_SCI_Eqn, data = dataExp)
model_PSP <- lmer(model_PSP_Eqn, data = dataExp)

#### ---- Fit models and generate estimates for comparison ####
# Get the data used in the models
model_SCI_Data <- model.frame(model_SCI)
model_PSP_Data <- model.frame(model_PSP)
# Filter the original data to match the model rows
data_SCI <- dataExp[rownames(model_SCI_Data), ]
data_PSP <- dataExp[rownames(model_PSP_Data), ]
# Compute predictions for model_SCI and model_PSP with the filtered data
data_SCI <- computePredictions(model_SCI, data_SCI, "Age * Cohort + Age * Task + Age * AP_Orientation + Age * Percentage_Motion_Windows + Homer_Pruned_Total")
data_PSP <- computePredictions(model_PSP, data_PSP, "Age * Cohort + Age * Task + Age * AP_Orientation + Age * Percentage_Motion_Windows + Homer_Pruned_Total")

# Extract fixed effects and standard errors
fixedEffects_SCI <- fixef(model_SCI)  # Fixed effects coefficients
se_SCI <- sqrt(diag(vcov(model_SCI)))  # Standard errors for fixed effects

fixedEffects_PSP <- fixef(model_PSP)  # Fixed effects coefficients
se_PSP <- sqrt(diag(vcov(model_PSP)))  # Standard errors for fixed effects

# Calculate t-values for each effect
tValues_SCI <- fixedEffects_SCI / se_SCI
tValues_PSP <- fixedEffects_PSP / se_PSP
pValues_SCI <- 2 * (1 - pt(abs(tValues_SCI), df = 1000))  # two-tailed p-value
pValues_PSP <- 2 * (1 - pt(abs(tValues_PSP), df = 1000))  # two-tailed p-value

# Combine fixed effects, t-values, and p-values into data frames
fixedEffects_SCI_df <- data.frame(
  Effect = names(fixedEffects_SCI),
  Estimate = fixedEffects_SCI,
  Std.Error = se_SCI,
  t_value = tValues_SCI,
  p_value = pValues_SCI
)

fixedEffects_PSP_df <- data.frame(
  Effect = names(fixedEffects_PSP),
  Estimate = fixedEffects_PSP,
  Std.Error = se_PSP,
  t_value = tValues_PSP,
  p_value = pValues_PSP
)

# View the results for both models
print(fixedEffects_SCI_df)
print(fixedEffects_PSP_df)


#### BOOTSTRAPPING ####
#### ---- Run bootstrapping ####
## SCI ##
# Set up parallel processing
registerDoParallel(numCores)

if (!file.exists("stats/bootResults_SCI.rds")) {

  # Time the bootstrapping process
  startTime <- Sys.time()


  # Perform bootstrap 
  bootResults_SCI <- foreach(i = 1:nBoot, .combine = rbind, .packages = c("lme4")) %dopar% {
    bootstrap_iteration(dataExp, model_SCI_Eqn, i)
  }

  # Calculate time taken
  endTime <- Sys.time()
  timeTaken <- endTime - startTime
  cat("Time taken for SCI bootstrapping:", timeTaken, "hours\n")

  # Save results
  saveRDS(bootResults_SCI, "stats/bootResults_SCI.rds")

} else {
  cat("SCI bootstrapping results already exist. Loading instead.\n")
  bootResults_SCI <- readRDS("stats/bootResults_SCI.rds")
}

# Stop parallel processing
stopImplicitCluster()


## PSP ##
# Set up parallel processing
registerDoParallel(numCores)

if (!file.exists("stats/bootResults_PSP.rds")) {

  # Time the bootstrapping process
  startTime <- Sys.time()


  # Perform bootstrap 
  bootResults_PSP <- foreach(i = 1:nBoot, .combine = rbind, .packages = c("lme4")) %dopar% {
    bootstrap_iteration(dataExp, model_PSP_Eqn, i)
  }

  # Calculate time taken
  endTime <- Sys.time()
  timeTaken <- endTime - startTime
  cat("Time taken for PSP bootstrapping:", timeTaken, "hours\n")

  # Save results
  saveRDS(bootResults_PSP, "stats/bootResults_PSP.rds")

} else {
  cat("PSP bootstrapping results already exist. Loading instead.\n")
  bootResults_PSP <- readRDS("stats/bootResults_PSP.rds")
}

# Stop parallel processing
stopImplicitCluster()

# Diagnostics: count the number of failed bootstraps (where NA is returned)
cat("\nNumber of failed bootstraps (SCI):", sum(rowSums(is.na(bootResults_SCI)) > 0), "\n")
cat("\nNumber of failed bootstraps (PSP):", sum(rowSums(is.na(bootResults_PSP)) > 0), "\n")


#### ---- Examine residuals ####
## SCI
# Compute average fixed effects from bootstrapping
avgFixef_Sci <- colMeans(bootResults_SCI, na.rm = TRUE)
# Create a grid of combinations of Age, Task, Cohort, etc.
unique_Age <- unique(dataExp$Age)
unique_Task <- unique(dataExp$Task)
unique_Cohort <- unique(dataExp$Cohort)
unique_Motion <- unique(dataExp$Percentage_Motion_Windows)
unique_AP <- unique(dataExp$AP_Orientation)
# Limit the number of unique values to a smaller subset (in this case, 10 random levels)
unique_Age <- sample(unique_Age, min(10, length(unique_Age)), replace = FALSE)
unique_Task <- sample(unique_Task, min(10, length(unique_Task)), replace = FALSE)
unique_Cohort <- sample(unique_Cohort, min(10, length(unique_Cohort)), replace = FALSE)
unique_Motion <- sample(unique_Motion, min(10, length(unique_Motion)), replace = FALSE)
unique_AP <- sample(unique_AP, min(10, length(unique_AP)), replace = FALSE)
# Create the grid with a reduced number of combinations
predictedData_SCI <- expand.grid(
  Age = unique_Age,
  Task = unique_Task,
  Cohort = unique_Cohort,
  Homer_Pruned_Total = mean(dataExp$Homer_Pruned_Total, na.rm = TRUE),  # Use mean to avoid large matrix
  Percentage_Motion_Windows = unique_Motion,
  AP_Orientation = unique_AP
)
# Add intercept to design matrix
X <- model.matrix(~ Age * Cohort + Age * Task + Age * AP_Orientation +
                    Age * Percentage_Motion_Windows + Homer_Pruned_Total,
                  data = predictedData_SCI)
# Compute predictions using the average fixed effects
predictedData_SCI$predictedValues <- X %*% avgFixef_Sci
# Merge predictedData_SCI with the true values from dataExp based on common columns (so rows match)
predictedData_SCI <- merge(predictedData_SCI, dataExp[, c("Age", "Task", "Cohort", "Percentage_Motion_Windows", "AP_Orientation", "Average_SCI")],
                           by = c("Age", "Task", "Cohort", "Percentage_Motion_Windows", "AP_Orientation"),
                           all.x = TRUE)
# Compute residuals (absolute differences between predicted and true values)
predictedData_SCI$residuals <- abs(predictedData_SCI$predictedValues - predictedData_SCI$Average_SCI)
# Summary of residuals
summary(predictedData_SCI$residuals)
# Histogram of residuals for SCI
ggplot(predictedData_SCI, aes(x = residuals)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Residuals for SCI Predictions",
       x = "Residuals",
       y = "Frequency") +
  theme_minimal()
# QQ plot of residuals for SCI
ggplot(predictedData_SCI, aes(sample = residuals)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "QQ Plot of Residuals for SCI Predictions") +
  theme_minimal()

## PSP
# Compute average fixed effects from bootstrapping
avgFixef_Psp <- colMeans(bootResults_PSP, na.rm = TRUE)
# Create a grid of combinations of Age, Task, Cohort, etc.
# Subset the data for selected variables if necessary
# Reduce the number of unique values for each factor variable to avoid excessive memory usage
unique_Age <- unique(dataExp$Age)
unique_Task <- unique(dataExp$Task)
unique_Cohort <- unique(dataExp$Cohort)
unique_Motion <- unique(dataExp$Percentage_Motion_Windows)
unique_AP <- unique(dataExp$AP_Orientation)
# Limit the number of unique values to a smaller subset (in this case, 10 random levels)
set.seed(42)  # Set seed for reproducibility
unique_Age <- sample(unique_Age, min(10, length(unique_Age)), replace = FALSE)
unique_Task <- sample(unique_Task, min(10, length(unique_Task)), replace = FALSE)
unique_Cohort <- sample(unique_Cohort, min(10, length(unique_Cohort)), replace = FALSE)
unique_Motion <- sample(unique_Motion, min(10, length(unique_Motion)), replace = FALSE)
unique_AP <- sample(unique_AP, min(10, length(unique_AP)), replace = FALSE)
# Create the grid with a reduced number of combinations
predictedData_PSP <- expand.grid(
  Age = unique_Age,
  Task = unique_Task,
  Cohort = unique_Cohort,
  Homer_Pruned_Total = mean(dataExp$Homer_Pruned_Total, na.rm = TRUE),  # Use mean to avoid large matrix
  Percentage_Motion_Windows = unique_Motion,
  AP_Orientation = unique_AP
)
# Add intercept to design matrix
X <- model.matrix(~ Age * Cohort + Age * Task + Age * Percentage_Motion_Windows +
                    Homer_Pruned_Total + AP_Orientation,
                  data = predictedData_PSP)
# Compute predictions using the average fixed effects
predictedData_PSP$predictedValues <- X %*% avgFixef_Psp
# Merge predictedData_PSP with the true values from dataExp based on common columns (so rows match)
predictedData_PSP <- merge(predictedData_PSP, dataExp[, c("Age", "Task", "Cohort", "Percentage_Motion_Windows", "AP_Orientation", "Average_PSP")],
                           by = c("Age", "Task", "Cohort", "Percentage_Motion_Windows", "AP_Orientation"),
                           all.x = TRUE)
# Compute residuals (absolute differences between predicted and true values)
predictedData_PSP$residuals <- abs(predictedData_PSP$predictedValues - predictedData_PSP$Average_PSP)
# Summary of residuals
summary(predictedData_PSP$residuals)
# Histogram of residuals for PSP
ggplot(predictedData_PSP, aes(x = residuals)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Residuals for PSP Predictions",
       x = "Residuals",
       y = "Frequency") +
  theme_minimal()
# QQ plot of residuals for PSP
ggplot(predictedData_PSP, aes(sample = residuals)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "QQ Plot of Residuals for PSP Predictions") +
  theme_minimal()

## Calculate skewness values
# between -0.5 and 0.5 OK; abs(value > 1 indicates severe skewness)
# Compute skewness of residuals while ignoring NA values
skewValue_SCI <- skewness(predictedData_SCI$residuals, na.rm = TRUE)
print(skewValue_SCI)

skewValue_PSP <- skewness(predictedData_PSP$residuals, na.rm = TRUE)
print(skewValue_PSP)

# Skewness of PSP model is fine; skewness of SCI model more extreme so will
# examine stability of the model

#### ---- Examine stability of SCI model ####

## Visualise the distribution of bootstrapped estimates
#  for 'Age'
ggplot(data = as.data.frame(bootResults_SCI), aes(x = Age)) +
  geom_density(fill = "lightblue", color = "darkblue", alpha = 0.5) +
  labs(title = "Distribution of Bootstrapped Estimates for Age",
       x = "Age Parameter Estimate", y = "Density")
# for Task
ggplot(data = as.data.frame(bootResults_SCI), aes(x = Tasksocial)) +
  geom_density(fill = "lightblue", color = "darkblue", alpha = 0.5) +
  labs(title = "Distribution of Bootstrapped Estimates for Task",
       x = "Age Parameter Estimate", y = "Density")
# for Cohort
ggplot(data = as.data.frame(bootResults_SCI), aes(x = Cohortuk)) +
  geom_density(fill = "lightblue", color = "darkblue", alpha = 0.5) +
  labs(title = "Distribution of Bootstrapped Estimates for Cohort",
       x = "Age Parameter Estimate", y = "Density")
#for Motion
ggplot(data = as.data.frame(bootResults_SCI), aes(x = Percentage_Motion_Windows)) +
  geom_density(fill = "lightblue", color = "darkblue", alpha = 0.5) +
  labs(title = "Distribution of Bootstrapped Estimates for Motion",
       x = "Age Parameter Estimate", y = "Density")
#for channel location
ggplot(data = as.data.frame(bootResults_SCI), aes(x = AP_Orientation)) +
  geom_density(fill = "lightblue", color = "darkblue", alpha = 0.5) +
  labs(title = "Distribution of Bootstrapped Estimates for Channel Location",
       x = "Age Parameter Estimate", y = "Density")
# for over/undersaturation
ggplot(data = as.data.frame(bootResults_SCI), aes(x = Homer_Pruned_Total)) +
  geom_density(fill = "lightblue", color = "darkblue", alpha = 0.5) +
  labs(title = "Distribution of Bootstrapped Estimates for Prior Channels Removed",
       x = "Age Parameter Estimate", y = "Density")

## Check skewness of bootstrapped estimates
# Check if the structure is correct
str(bootResults_SCI)
#  for 'Age'
if("Age" %in% colnames(bootResults_SCI)) {
  skew_age <- skewness(bootResults_SCI[, "Age"], na.rm = TRUE)
  print(paste("Skewness of Age bootstrapped estimates: ", skew_age))
} else {
  print("Column 'Age' not found in bootResults_SCI.")
}
#  for 'Task'
if("Age" %in% colnames(bootResults_SCI)) {
  skew_age <- skewness(bootResults_SCI[, "Tasksocial"], na.rm = TRUE)
  print(paste("Skewness of Task bootstrapped estimates: ", skew_age))
} else {
  print("Column 'Task' not found in bootResults_SCI.")
}
# for Cohort
if("Age" %in% colnames(bootResults_SCI)) {
  skew_age <- skewness(bootResults_SCI[, "Cohortuk"], na.rm = TRUE)
  print(paste("Skewness of Cohort bootstrapped estimates: ", skew_age))
} else {
  print("Column 'Cohort' not found in bootResults_SCI.")
}
# for Motion
if("Age" %in% colnames(bootResults_SCI)) {
  skew_age <- skewness(bootResults_SCI[, "Percentage_Motion_Windows"], na.rm = TRUE)
  print(paste("Skewness of Motion bootstrapped estimates: ", skew_age))
} else {
  print("Column 'Cohort' not found in bootResults_SCI.")
}
# for channel location
if("Age" %in% colnames(bootResults_SCI)) {
  skew_age <- skewness(bootResults_SCI[, "AP_Orientation"], na.rm = TRUE)
  print(paste("Skewness of Channel location bootstrapped estimates: ", skew_age))
} else {
  print("Column 'Age' not found in bootResults_SCI.")
}
# for channel signal extrema
if("Age" %in% colnames(bootResults_SCI)) {
  skew_age <- skewness(bootResults_SCI[, "Homer_Pruned_Total"], na.rm = TRUE)
  print(paste("Skewness of Homer_Pruned_Channel bootstrapped estimates: ", skew_age))
} else {
  print("Column 'Age' not found in bootResults_SCI.")
}

## compare original (fitted) results wiht bootstrapped results
# 'FALSE' shows that original estimates sit outside CI for bootstrapped results, possibly indicating instability
# assign variable to store original results
original_results <- fixef(Model_SCI)
# Check if the original estimates fall within the CI range
within_ci <- sapply(names(original_results), function(param) {
  original_val <- original_results[param]
  lower <- lower_ci_bootstrapped_params[param]
  upper <- upper_ci_bootstrapped_params[param]
  return(original_val >= lower && original_val <= upper)
})
# Print which parameters are within the CI range
print(within_ci)

#### ---- Display bootstrapping results ####
# Calculate confidence intervals for model_SCI
nCoef_SCI <- length(fixef(model_SCI))
SCI_CoefNames <- names(fixef(model_SCI))
SCI_CIs <- t(sapply(1:nCoef_SCI, function(i) {
  quantile(bootResults_SCI[, i], probs = c(0.025, 0.975), na.rm = TRUE)
}))
rownames(SCI_CIs) <- SCI_CoefNames

# Calculate standard errors from bootstrapped estimates for model_SCI
SCI_SEs <- apply(bootResults_SCI[, 1:nCoef_SCI], 2, sd, na.rm = TRUE)

# Create summary table for model_SCI
SCI_Summary <- data.frame(
  Estimate = colMeans(bootResults_SCI),
  CI_lower = SCI_CIs[, "2.5%"],
  CI_upper = SCI_CIs[, "97.5%"],
  SE = SCI_SEs
)

# Print results for model_SCI
cat("\nModel Interaction SCI Results:\n")
print(SCI_Summary)

# Calculate confidence intervals for model_PSP
nCoef_PSP <- length(fixef(model_PSP))
PSP_CoefNames <- names(fixef(model_PSP))
PSP_CIs <- t(sapply(1:nCoef_PSP, function(i) {
  quantile(bootResults_PSP[, i], probs = c(0.025, 0.975), na.rm = TRUE)
}))
rownames(PSP_CIs) <- PSP_CoefNames

# Calculate standard errors from bootstrapped estimates for model_PSP
PSP_SEs <- apply(bootResults_PSP[, 1:nCoef_PSP], 2, sd, na.rm = TRUE)

# Create summary table for model_PSP
PSP_Summary <- data.frame(
  Estimate = colMeans(bootResults_PSP),
  CI_lower = PSP_CIs[, "2.5%"],
  CI_upper = PSP_CIs[, "97.5%"],
  SE = PSP_SEs
)


# Print results for model_PSP
cat("\nModel Interaction PSP Results:\n")
print(PSP_Summary)

# save summaries
saveRDS(PSP_Summary, "stats/PSP_Summary.rds")
saveRDS(SCI_Summary, "stats/SCI_Summary.rds")


#### ---- Calculate and display effect sizes ####
# Calculate effect sizes for both models
SCI_Effects <- standardiseEffects(model_SCI, SCI_Summary)
PSP_Effects <- standardiseEffects(model_PSP, PSP_Summary)

# calculate models fits
SCI_Fit <- r.squaredGLMM(model_SCI)
PSP_Fit <- r.squaredGLMM(model_PSP)

# Print overall model fit
cat("\nSCI Model Overall Fit:\n")
print(SCI_Fit)
cat("\nChannels Retained Model Overall Fit:\n")
print(PSP_Fit)

# Print detailed results
cat("\nSCI Model Effect Sizes:\n")
print(SCI_Effects)
cat("\nPSP Model Effect Sizes:\n")
print(PSP_Effects)

#save effect sizes
saveRDS(SCI_Effects, "stats/SCI_Effects.rds")
saveRDS(PSP_Effects, "stats/PSP_Effects.rds")
write.csv(SCI_Effects, "stats/SCI_Effects.csv")
write.csv(PSP_Effects, "stats/PSP_Effects.csv")

# save fits
saveRDS(SCI_Fit, "stats/SCI_Fit.rds")
saveRDS(PSP_Fit, "stats/PSP_Fit.rds")

#### PLOT TRENDS ####
#### ---- Load variables if necessary ####
# Check and load outputs if they are not already in the environment
if (!exists("bootResults_SCI")) {
  bootResults_SCI <- readRDS("stats/bootResults_SCI.rds")
}
if (!exists("bootResults_PSP")) {
  bootResults_PSP <- readRDS("stats/bootResults_PSP.rds")
}

if (!exists("SCI_Summary")) {
  SCI_Summary <- readRDS("stats/SCI_Summary.rds")
}
if (!exists("PSP_Summary")) {
  PSP_Summary <- readRDS("stats/PSP_Summary.rds")
}

if (!exists("SCI_Effects")) {
  SCI_Effects <- readRDS("stats/SCI_Effects.rds")
}
if (!exists("PSP_Effects")) {
  PSP_Effects <- readRDS("stats/PSP_Effects.rds")
}

if (!exists("SCI_Fit")) {
  SCI_Fit <- readRDS("stats/SCI_Fit.rds")
}
if (!exists("PSP_Fit")) {
  PSP_Fit <- readRDS("stats/PSP_Fit.rds")
}
#### ---- Forest plots ----

# New column names for bootResults_SCI for plots
newColnames <- c("Intercept", "Age", "Cohort (UK)", "Task (Social)", "Channel Location", "Percentage of Motion", "CSE",
                 "Age * Cohort (UK)", "Age * Task (Social)", "Age * Channel Location", "Age * Percentage of Motion")
# Assign the new column names to bootResults_SCI
colnames(bootResults_SCI) <- newColnames

# New column names for bootResults_PSP for plots
newColnames <- c("Intercept", "Age", "Cohort (UK)", "Task (Social)", "Percentage of Motion", "CSE", "Channel Location",
                 "Age * Cohort (UK)", "Age * Task (Social)","Age * Channel Location", "Age * Percentage of Motion")
# Assign the new column names to bootResults_PSP
colnames(bootResults_PSP) <- newColnames

# Prepare data for forest plot of Model_SCI
SCI_ForestPlot_Data <- data.frame(
  Coefficient = colnames(bootResults_SCI),
  Estimate = SCI_Summary$Estimate,
  CI_lower = SCI_Summary$CI_lower,
  CI_upper = SCI_Summary$CI_upper,
  SE = SCI_Summary$SE,
  Effect_Size = SCI_Effects$Effect_Size
)

# Prepare data for forest plot of Model_PSP
PSP_ForestPlot_Data <- data.frame(
  Coefficient = colnames(bootResults_PSP),
  Estimate = PSP_Summary$Estimate,
  CI_lower = PSP_Summary$CI_lower,
  CI_upper = PSP_Summary$CI_upper,
  SE = PSP_Summary$SE,
  Effect_Size = PSP_Effects$Effect_Size
)

# Create a new column for effect type based on size and direction
# SNR
SCI_ForestPlot_Data <- SCI_ForestPlot_Data %>%
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

# CR
PSP_ForestPlot_Data <- PSP_ForestPlot_Data %>%
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
  "Not Significant" = 4, 
  "Large_Positive" = 15,  
  "Medium_Positive" = 19, 
  "Small_Positive" = 17,  
  "Large_Negative" = 0,   
  "Medium_Negative" = 1,  
  "Small_Negative" = 2    
)

#### ---- ---- with intercept ----

# Forest plot with the intercept (default)
ggplot(SCI_ForestPlot_Data, aes(x = Estimate, y = Coefficient)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray", size = 1) + 
  geom_point(aes(color = Effect_Type, shape = Effect_Type), size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.4) +
  scale_color_manual(values = color_mapping) +
  scale_shape_manual(values = shape_mapping) +
  labs(title = "Bootstrapped Fixed Effects for Average SCI",
       x = "Estimate with 95% CI", y = "Coefficient",
       color = "Standardised Effect Type", shape = "Standardised Effect Type") +
  theme_minimal()

# Forest plot with the intercept for PSP model
ggplot(PSP_ForestPlot_Data, aes(x = Estimate, y = Coefficient)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray", size = 1) + 
  geom_point(aes(color = Effect_Type, shape = Effect_Type), size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.4) +
  scale_color_manual(values = color_mapping) +
  scale_shape_manual(values = shape_mapping) +
  labs(title = "Bootstrapped Fixed Effects for Average PSP",
       x = "Estimate with 95% CI", y = "Coefficient",
       color = "Standardised Effect Type", shape = "Standardised Effect Type") +
  theme_minimal()

#### ---- ---- without intercept ----

# Remove intercept from the SCI data
SCI_ForestPlot_Data_no_intercept <- SCI_ForestPlot_Data %>%
  filter(Coefficient != "Intercept")

# Forest plot without the intercept for SCI
ggplot(SCI_ForestPlot_Data_no_intercept, aes(x = Estimate, y = Coefficient)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray", size = 1) + 
  geom_point(aes(color = Effect_Type, shape = Effect_Type), size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.4) +
  scale_color_manual(values = color_mapping) +
  scale_shape_manual(values = shape_mapping) +
  labs(title = "Bootstrapped Fixed Effects for Average SCI",
       x = "Estimate with 95% CI", y = "Coefficient",
       color = "Standardised Effect Type", shape = "Standardised Effect Type") +
  theme_minimal()

# Remove intercept from the PSP data
PSP_ForestPlot_Data_no_intercept <- PSP_ForestPlot_Data %>%
  filter(Coefficient != "Intercept")

# Forest plot without the intercept for PSP model
ggplot(PSP_ForestPlot_Data_no_intercept, aes(x = Estimate, y = Coefficient)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray", size = 1) + 
  geom_point(aes(color = Effect_Type, shape = Effect_Type), size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.4) +
  scale_color_manual(values = color_mapping) +
  scale_shape_manual(values = shape_mapping) +
  labs(title = "Bootstrapped Fixed Effects for Average PSP",
       x = "Estimate with 95% CI", y = "Coefficient",
       color = "Standardised Effect Type", shape = "Standardised Effect Type") +
  theme_minimal()


#### ---- Effects and interaction effects  ####
#### ---- ---- Create DF with effects (N.B. ALWAYS RUN AGAIN FIRST!) ####

outputDir <- "stats"

# Check and load predicted SCI results if they exist; create them if not
if (file.exists(paste0(outputDir, "/predictedData_SCI.rds"))) {
  cat("\nSCI predicted data already exists. Loading...\n")
  predictedData_SCI <- readRDS(paste0(outputDir, "/predictedData_SCI.rds"))
} else {
  
  # Create a new data frame with all necessary variables
  # Compute average fixed effects from bootstrapping
  avgFixef_Sci <- colMeans(bootResults_SCI, na.rm = TRUE)
  
  # Create a grid of combinations of Age, Task, Cohort, etc.
  predictedData_SCI <- expand.grid(
    Age = unique(dataExp$Age),
    Task = unique(dataExp$Task),
    Cohort = unique(dataExp$Cohort),
    Homer_Pruned_Total = mean(dataExp$Homer_Pruned_Total, na.rm = TRUE), #uses too much RAM otherwise
    Percentage_Motion_Windows = unique(dataExp$Percentage_Motion_Windows),
    AP_Orientation = mean(dataExp$AP_Orientation, na.rm = TRUE) #uses too much RAM otherwise
  )
  
  # Initialise output directory for saving chunks
  if (!dir.exists(outputDir)) {
    dir.create(outputDir)
  }
  
  # Set up parallel backend for SCI
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  
  # Calculate chunks and split predictedData_SCI into chunks
  nChunks_SCI <- ceiling(nrow(predictedData_SCI) / chunkSizeDF)
  chunks_SCI <- split(predictedData_SCI, rep(1:nChunks_SCI, each = chunkSizeDF, length.out = nrow(predictedData_SCI)))
  
  foreach(chunkIdx = seq_along(chunks_SCI), .packages = c("stats")) %dopar% {
    cat("\nProcessing chunk", chunkIdx, "of", nChunks_SCI, "\n")
    
    # Extract current chunk
    chunk <- chunks_SCI[[chunkIdx]]
    
    # Create design matrix for SCI
    X_SCI <- model.matrix(~ Age * Cohort + Age * Task + Age * AP_Orientation + 
                            Age * Percentage_Motion_Windows + Homer_Pruned_Total,
                          data = chunk)
    
    # Calculate predictions for all bootstrap iterations
    chunk$predictedValues <- apply(bootResults_SCI, 1, function(coefs) X_SCI %*% coefs)
    
    # Add confidence intervals to the chunk
    chunk$LowerCI <- apply(chunk$predictedValues, 1, quantile, probs = 0.025, na.rm = TRUE)
    chunk$UpperCI <- apply(chunk$predictedValues, 1, quantile, probs = 0.975, na.rm = TRUE)
    
    # Save chunk results to file to save memory
    saveRDS(chunk, file = paste0(outputDir, "/predictedData_SCI_chunk_", chunkIdx, ".rds"))
  }
  
  # Combine all saved chunks into a single data frame
  predictedData_SCI <- do.call(rbind, lapply(
    list.files(outputDir, pattern = "predictedData_SCI_chunk_.*\\.rds", full.names = TRUE),
    readRDS
  ))
  
  # Save combined results as "predictedData_PSP.rds"
  saveRDS(predictedData_SCI, file = paste0(outputDir, "/predictedData_SCI.rds"))
  
  # Clean up (delete) all chunk files once combined file is saved
  file.remove(list.files(outputDir, pattern = "predictedData_SCI_chunk_.*\\.rds", full.names = TRUE))
  
  # stop par. processing
  stopCluster(cl)
}

# Check and load predicted PSP results if they exist; create them if not
if (file.exists(paste0(outputDir, "/predictedData_PSP.rds"))) {
  cat("\nPSP predicted data already exists. Loading...\n")
  predictedData_PSP <- readRDS(paste0(outputDir, "/predictedData_PSP.rds"))
} else {
  # Create a new data frame with all necessary variables
  # Compute average fixed effects from bootstrapping
  avgFixef_Psp <- colMeans(bootResults_PSP, na.rm = TRUE)
  
  # Create a grid of combinations of Age, Task, Cohort, etc.
  predictedData_PSP <- expand.grid(
    Age = unique(dataExpOrig$Age),
    Task = unique(dataExpOrig$Task),
    Cohort = unique(dataExpOrig$Cohort),
    Homer_Pruned_Total = mean(dataExpOrig$Homer_Pruned_Total, na.rm = TRUE), #uses too much RAM otherwise
    Percentage_Motion_Windows = unique(dataExpOrig$Percentage_Motion_Windows),
    AP_Orientation = mean(dataExpOrig$AP_Orientation, na.rm = TRUE) #uses too much RAM otherwise
  )
  
  # Initialise output directory for saving chunks
  if (!dir.exists(outputDir)) {
    dir.create(outputDir)
  }
  
  # Set up parallel backend for PSP
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  
  # Calculate chunks and split predictedData_PSP into chunks
  nChunks <- ceiling(nrow(predictedData_PSP) / chunkSizeDF)
  chunks <- split(predictedData_PSP, rep(1:nChunks, each = chunkSizeDF, length.out = nrow(predictedData_PSP)))
  
  foreach(chunkIdx = seq_along(chunks), .packages = c("stats")) %dopar% {
    cat("\nProcessing chunk", chunkIdx, "of", nChunks, "\n")
    
    # Extract current chunk
    chunk <- chunks[[chunkIdx]]
    
    # Create design matrix for PSP
    X_PSP <- model.matrix(~ ~ Age * Cohort + Age * Task + Age * AP_Orientation + 
                            Age * Percentage_Motion_Windows + Homer_Pruned_Total,
                          data = chunk)
    
    # Calculate predictions for all bootstrap iterations
    chunk$predictedValues <- apply(bootResults_PSP, 1, function(coefs) X_PSP %*% coefs)
    
    # Add confidence intervals to the chunk
    chunk$LowerCI <- apply(chunk$predictedValues, 1, quantile, probs = 0.025, na.rm = TRUE)
    chunk$UpperCI <- apply(chunk$predictedValues, 1, quantile, probs = 0.975, na.rm = TRUE)
    
    # Save chunk results to file to save memory
    saveRDS(chunk, file = paste0(outputDir, "/predictedData_PSP_chunk_", chunkIdx, ".rds"))
  }
  
  # Combine all saved chunks into a single data frame
  predictedData_PSP <- do.call(rbind, lapply(
    list.files(outputDir, pattern = "predictedData_PSP_chunk_.*\\.rds", full.names = TRUE),
    readRDS
  ))
  
  # Save combined results as "predictedData_PSP.rds"
  saveRDS(predictedData_PSP, file = paste0(outputDir, "/predictedData_PSP.rds"))
  
  # Clean up (delete) all chunk files once combined file is saved
  file.remove(list.files(outputDir, pattern = "predictedData_PSP_chunk_.*\\.rds", full.names = TRUE))
  
  # Clean up
  stopCluster(cl)
}

#### ---- ---- Effect of Age on SCI and PSP ####
# SCI
# Aggregate data by Age
predictedData_SCI_summary <- predictedData_SCI %>%
  group_by(Age) %>%
  summarise(predictedValues = mean(predictedValues), 
            LowerCI = mean(LowerCI), 
            UpperCI = mean(UpperCI), .groups = 'drop')

# Plot
ggplot(predictedData_SCI_summary, aes(x = Age, y = predictedValues)) +
  geom_line(color = viridis(20)[5]) +  
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI), fill = viridis(20)[5], alpha = 0.2) +  
  labs(title = "Predicted SCI Trend by Age", 
       x = "Age", 
       y = "Average SCI") +
  theme_minimal() +  
  theme(
    panel.grid.major = element_line(color = "grey80"),  
    panel.grid.minor = element_line(color = "grey90"), 
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA)
  )

# PSP
# Aggregate data by Age
predictedData_PSP_summary <- predictedData_PSP %>%
  group_by(Age) %>%
  summarise(predictedValues = mean(predictedValues), 
            LowerCI = mean(LowerCI), 
            UpperCI = mean(UpperCI), .groups = 'drop')

# Plot
ggplot(predictedData_PSP_summary, aes(x = Age, y = predictedValues)) +
  geom_line(color = viridis(20)[17]) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI), fill = viridis(20)[17], alpha = 0.2) + 
  labs(title = "Predicted PSP Trend by Age", 
       x = "Age", 
       y = "Average PSP") +
  theme_minimal() +  
  theme(
    panel.grid.major = element_line(color = "grey80"),  
    panel.grid.minor = element_line(color = "grey90"), 
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA)
  )

#### ---- ---- Effect of Motion on SCI and PSP ####
# SCI
# Aggregate data by Age
predictedData_SCI_summary <- predictedData_SCI %>%
  group_by(Percentage_Motion_Windows) %>%
  summarise(predictedValues = mean(predictedValues), 
            LowerCI = mean(LowerCI), 
            UpperCI = mean(UpperCI), .groups = 'drop')

# Plot
ggplot(predictedData_SCI_summary, aes(x = Percentage_Motion_Windows, y = predictedValues)) +
  geom_line(color = viridis(20)[5]) +  
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI), fill = viridis(20)[5], alpha = 0.2) + 
  labs(title = "Predicted SCI Trend by motion", 
       x = "Percentage of windows with motion", 
       y = "Average SCI") +
  theme_minimal() +  
  theme(
    panel.grid.major = element_line(color = "grey80"),  
    panel.grid.minor = element_line(color = "grey90"), 
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA)
  )

# PSP
# Aggregate data by Age
predictedData_PSP_summary <- predictedData_PSP %>%
  group_by(Percentage_Motion_Windows) %>%
  summarise(predictedValues = mean(predictedValues), 
            LowerCI = mean(LowerCI), 
            UpperCI = mean(UpperCI), .groups = 'drop')

# Plot
ggplot(predictedData_PSP_summary, aes(x = Percentage_Motion_Windows, y = predictedValues)) +
  geom_line(color = viridis(20)[15]) +  
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI), fill = viridis(20)[15], alpha = 0.2) + 
  labs(title = "Predicted PSP Trend by Motion", 
       x = "Percentage of windows containing motion", 
       y = "Average PSP") +
  theme_minimal() +  
  theme(
    panel.grid.major = element_line(color = "grey80"),  
    panel.grid.minor = element_line(color = "grey90"), 
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA)
  )

#### ---- ---- Effect of Age * Task ####
# SCI
# Aggregate data by Age and Task
predictedData_SCI_summary <- predictedData_SCI %>%
  group_by(Age, Task) %>%
  summarise(predictedValues = mean(predictedValues), 
            LowerCI = mean(LowerCI), 
            UpperCI = mean(UpperCI), .groups = 'drop')

# Plot
ggplot(predictedData_SCI_summary, aes(x = Age, y = predictedValues, color = Task)) +
  geom_line(size = 2) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = Task), alpha = 0.2) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  labs(title = "Predicted SCI Trends by Age and Task", 
       x = "Age", 
       y = "Average SCI") +
  theme_minimal() +  
  theme(
    panel.grid.major = element_line(color = "grey80"),  
    panel.grid.minor = element_line(color = "grey90"), 
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA)
  )

# PSP
# Aggregate data by Age and Task
predictedData_PSP_summary <- predictedData_PSP %>%
  group_by(Age, Task) %>%
  summarise(predictedValues = mean(predictedValues), 
            LowerCI = mean(LowerCI), 
            UpperCI = mean(UpperCI), .groups = 'drop')

# Plot
ggplot(predictedData_PSP_summary, aes(x = Age, y = predictedValues, color = Task)) +
  geom_line(size = 2) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = Task), alpha = 0.2) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  labs(title = "Predicted PSP Trends by Age and Task", 
       x = "Age", 
       y = "Average PSP") +
  theme_minimal() +  
  theme(
    panel.grid.major = element_line(color = "grey80"),  
    panel.grid.minor = element_line(color = "grey90"), 
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA)
  )
#### ---- ---- Effect of Age * Cohort ####
# SCI
# Aggregate data by Age and Cohort
predictedData_SCI_summary <- predictedData_SCI %>%
  group_by(Age, Cohort) %>%
  summarise(predictedValues = mean(predictedValues), 
            LowerCI = mean(LowerCI), 
            UpperCI = mean(UpperCI), .groups = 'drop')

# Plot
ggplot(predictedData_SCI_summary, aes(x = Age, y = predictedValues, color = Cohort)) +
  geom_line(size = 2) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = Cohort), alpha = 0.2) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  labs(title = "Predicted SCI Trends by Age and Cohort", 
       x = "Age", 
       y = "Average SCI") +
  theme_minimal() +  
  theme(
    panel.grid.major = element_line(color = "grey80"),  
    panel.grid.minor = element_line(color = "grey90"), 
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA)
  )

# PSP
# Aggregate data
predictedData_PSP_summary <- predictedData_PSP %>%
  group_by(Age, Cohort) %>%
  summarise(predictedValues = mean(predictedValues), 
            LowerCI = mean(LowerCI), 
            UpperCI = mean(UpperCI), .groups = 'drop')

# Plot
ggplot(predictedData_PSP_summary, aes(x = Age, y = predictedValues, color = Cohort)) +
  geom_line(size = 2) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = Cohort), alpha = 0.2) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  labs(title = "Predicted PSP Trends by Age and Cohort", 
       x = "Age", 
       y = "Average PSP") +
  theme_minimal() +  
  theme(
    panel.grid.major = element_line(color = "grey80"),  
    panel.grid.minor = element_line(color = "grey90"), 
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA)
  )
#### ---- ---- Effect of Age * Motion ####
# SCI
# Aggregate data by Age and Motion
predictedData_SCI_summary <- predictedData_SCI %>%
  group_by(Age, Percentage_Motion_Windows) %>%
  summarise(predictedValues = mean(predictedValues), 
            LowerCI = mean(LowerCI), 
            UpperCI = mean(UpperCI), .groups = 'drop')

# Plot (treats Age as a factor for colour-coding purposes)
ggplot(predictedData_SCI_summary, aes(x = Percentage_Motion_Windows, y = predictedValues, color = factor(Age))) +
  geom_line(size = 2) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = factor(Age)), alpha = 0.2) +
  scale_color_viridis(discrete = TRUE) + 
  scale_fill_viridis(discrete = TRUE) +
  labs(title = "Predicted SCI Trends by Percentage of Motion and Age", 
       x = "Percentage of windows containing motion", 
       y = "Average SCI") +
  theme_minimal() +  
  theme(
    panel.grid.major = element_line(color = "grey80"),  
    panel.grid.minor = element_line(color = "grey90"), 
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA)
  )


# PSP
# Aggregate data by Age and Motion
predictedData_PSP_summary <- predictedData_PSP %>%
  group_by(Age, Percentage_Motion_Windows) %>%
  summarise(predictedValues = mean(predictedValues), 
            LowerCI = mean(LowerCI), 
            UpperCI = mean(UpperCI), .groups = 'drop')

# Plot (treats Age as a factor for colour-coding purposes)
ggplot(predictedData_PSP_summary, aes(x = Percentage_Motion_Windows, y = predictedValues, color = factor(Age))) +
  geom_line(size = 2) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = factor(Age)), alpha = 0.1) +
  scale_color_viridis(discrete = TRUE) + 
  scale_fill_viridis(discrete = TRUE) +
  labs(title = "Predicted PSP Trends by Percentage of Motion and Age", 
       x = "Percentage of windows containing motion", 
       y = "Average PSP") +
  theme_minimal() +  
  theme(
    panel.grid.major = element_line(color = "grey80"),  
    panel.grid.minor = element_line(color = "grey90"), 
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA)
  )
#### ---- ---- Effect of Age, Task and Cohort ####
# Aggregate data by Age, Task, and Cohort
predictedData_SCI_summary <- predictedData_SCI %>%
  group_by(Age, Task, Cohort) %>%
  summarise(predictedValues = mean(predictedValues), 
            LowerCI = mean(LowerCI), 
            UpperCI = mean(UpperCI), .groups = 'drop')

# Plot effect of Age by Task on SCI with faceting by Cohort
ggplot(predictedData_SCI_summary, aes(x = Age, y = predictedValues, color = Task)) +
  geom_line(size = 2) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = Task), alpha = 0.2) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  labs(title = "Predicted SCI Trends by Age and Task, per Cohort", 
       x = "Age", 
       y = "Average SCI") +
  theme_minimal() +  
  theme(
    panel.grid.major = element_line(color = "grey80"),  
    panel.grid.minor = element_line(color = "grey90"), 
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  facet_wrap(~ Cohort)  # (creates two plots in same figure)

# Plot effect of Age by Cohort on SCI with faceting by Task
ggplot(predictedData_SCI_summary, aes(x = Age, y = predictedValues, color = Cohort)) +
  geom_line(size = 2) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = Cohort), alpha = 0.2) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  labs(title = "Predicted SCI Trends by Age and Cohort, per Task", 
       x = "Age", 
       y = "Average SCI") +
  theme_minimal() +  
  theme(
    panel.grid.major = element_line(color = "grey80"),  
    panel.grid.minor = element_line(color = "grey90"), 
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  facet_wrap(~ Task)  # (creates two plots in same figure)

# Aggregate data by Age, Task, and Cohort
predictedData_PSP_summary <- predictedData_PSP %>%
  group_by(Age, Task, Cohort) %>%
  summarise(predictedValues = mean(predictedValues), 
            LowerCI = mean(LowerCI), 
            UpperCI = mean(UpperCI), .groups = 'drop')

# Plot for PSP with faceting by Cohort
ggplot(predictedData_PSP_summary, aes(x = Age, y = predictedValues, color = Task)) +
  geom_line(size = 2) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = Task), alpha = 0.2) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  labs(title = "Predicted PSP Trends by Age and Task, per Cohort", 
       x = "Age", 
       y = "Average PSP") +
  theme_minimal() +  
  theme(
    panel.grid.major = element_line(color = "grey80"),  
    panel.grid.minor = element_line(color = "grey90"), 
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  facet_wrap(~ Cohort)  # (creates two plots in same figure)

# Plot for PSP with faceting by Task
ggplot(predictedData_PSP_summary, aes(x = Age, y = predictedValues, color = Cohort)) +
  geom_line(size = 2) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = Cohort), alpha = 0.2) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  labs(title = "Predicted PSP Trends by Age and Cohort, per Task", 
       x = "Age", 
       y = "Average PSP") +
  theme_minimal() +  
  theme(
    panel.grid.major = element_line(color = "grey80"),  
    panel.grid.minor = element_line(color = "grey90"), 
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  facet_wrap(~ Task)  # (creates two plots in same figure)













#### EXAMINE EFFECT OF AGE PER TASK & COHORT ####
# define new models wthout Task or Cohort factors
Model_SCIEqn_NoTaskCohort <- Average_SCI ~ Age * AP_Orientation + 
  Age * Percentage_Motion_Windows + Homer_Pruned_Total + (1 | ID)
Model_PSPEqn_NoTaskCohort <- Average_PSP ~ Age * Percentage_Motion_Windows + 
  Homer_Pruned_Total + AP_Orientation + (1 | ID)

# Subset data by Task and Cohort
gm_hand <- dataExpOrig %>% filter(Cohort == "gm", Task == "hand")
gm_social <- dataExpOrig %>% filter(Cohort == "gm", Task == "social")
uk_hand <- dataExpOrig %>% filter(Cohort == "uk", Task == "hand")
uk_social <- dataExpOrig %>% filter(Cohort == "uk", Task == "social")

# Fit models
# SCI
Model_SCI_gm_hand <- lmer(Model_SCIEqn_NoTaskCohort, data = gm_hand)
Model_SCI_gm_social <- lmer(Model_SCIEqn_NoTaskCohort, data = gm_social)
Model_SCI_uk_hand <- lmer(Model_SCIEqn_NoTaskCohort, data = uk_hand)
Model_SCI_uk_social <- lmer(Model_SCIEqn_NoTaskCohort, data = uk_social)
#PSP
Model_PSP_gm_hand <- lmer(Model_PSPEqn_NoTaskCohort, data = gm_hand)
Model_PSP_gm_social <- lmer(Model_PSPEqn_NoTaskCohort, data = gm_social)
Model_PSP_uk_hand <- lmer(Model_PSPEqn_NoTaskCohort, data = uk_hand)
Model_PSP_uk_social <- lmer(Model_PSPEqn_NoTaskCohort, data = uk_social)

# Summarise results
#SCI
summary(Model_SCI_gm_hand)
summary(Model_SCI_gm_social)
summary(Model_SCI_uk_hand)
summary(Model_SCI_uk_social)
#PSP
summary(Model_PSP_gm_hand)
summary(Model_PSP_gm_social)
summary(Model_PSP_uk_hand)
summary(Model_PSP_uk_social)
#### PLOT AVG SCI AND PSP BY AGE ####
#### ---- Filter negative values ####
dataExp_filtered <- dataExpOrig %>%
  filter(Average_SCI > 0, Average_PSP > 0)
#### ---- Jittered scatter plots with density ####
# Scatter plot for Average SCI with density-based transparency and black crosses
scatterplot_SCI_density <- ggplot(dataExp_filtered, aes(x = Age, y = Average_SCI)) +
  stat_density_2d(aes(fill = ..density..), geom = "tile", contour = FALSE, alpha = 0.6) +
  geom_point(shape = 4, size = 1, color = "black", position = position_jitter(width = 0.3), alpha = 0.05) +
  geom_smooth(method = "lm", color = "darkblue", se = FALSE) +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Scatter Plot of Average SCI by Age",
       x = "Age (Months)",
       y = "Average SCI") +
  theme_minimal() +
  ylim(0, NA) +
  theme(plot.title = element_text(size = 10))

# Scatter plot for Average PSP with density-based transparency and black crosses
scatterplot_PSP_density <- ggplot(dataExp_filtered, aes(x = Age, y = Average_PSP)) +
  stat_density_2d(aes(fill = ..density..), geom = "tile", contour = FALSE, alpha = 0.6) +
  geom_point(shape = 4, size = 1, color = "black", position = position_jitter(width = 0.3), alpha = 0.05) +
  geom_smooth(method = "lm", color = "darkred", se = FALSE) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Scatter Plot of Average PSP by Age",
       x = "Age (Months)",
       y = "Average PSP") +
  theme_minimal() +
  ylim(0, NA) +
  theme(plot.title = element_text(size = 10))  


print(scatterplot_SCI_density)
print(scatterplot_PSP_density)

#### ---- Density plots ####
# Transposed Density plot for Average SCI by Age 
density_SCI <- ggplot(dataExpOrig, aes(y = Average_SCI, x = ..density.., fill = factor(Age))) +
  geom_density(alpha = 0.6) +
  labs(title = "Density Plot of Average SCI by Age",
       x = "Density",
       y = "Average SCI") +
  theme_minimal() +
  ylim(0.98, NA) +
  theme(plot.title = element_text(size = 10))  

# Transposed Density plot for Average PSP by Age
density_PSP <- ggplot(dataExp, aes(y = Average_PSP, x = ..density.., fill = factor(Age))) +
  geom_density(alpha = 0.6) +
  labs(title = "Density Plot of Average PSP by Age",
       x = "Density",
       y = "Average PSP") +
  theme_minimal() +
  ylim(0, NA) +  
  theme(plot.title = element_text(size = 10))  

# print(density_SCI)
# print(density_PSP)


#### ---- ---- Predicted, scatter and density plots on same figure ####
# Aggregate data by Age, Task, and Cohort
predictedData_SCI_summary <- predictedData_SCI %>%
  group_by(Age, Task, Cohort) %>%
  summarise(predictedValues = mean(predictedValues), 
            LowerCI = mean(LowerCI), 
            UpperCI = mean(UpperCI), .groups = 'drop')
# PSP
# Aggregate data by Age, Task, and Cohort
predictedData_PSP_summary <- predictedData_PSP %>%
  group_by(Age, Task, Cohort) %>%
  summarise(predictedValues = mean(predictedValues), 
            LowerCI = mean(LowerCI), 
            UpperCI = mean(UpperCI), .groups = 'drop')

# Arrange the plots in a 3x2 grid (3 rows, 2 columns)
grid.arrange(
  # SCI Plots (Top row)
  ggplot(predictedData_SCI_summary, aes(x = Age, y = predictedValues, color = Cohort)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = Cohort), alpha = 0.2) +
    scale_color_viridis(discrete = TRUE) +
    scale_fill_viridis(discrete = TRUE) +
    labs(title = "Predicted SCI Trends by Age and Cohort, per Task", 
         x = "Age", 
         y = "Average SCI") +
    theme_minimal() +  
    theme(
      panel.grid.major = element_line(color = "grey80"),  
      panel.grid.minor = element_line(color = "grey90"), 
      panel.background = element_rect(fill = "white", color = NA),  
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    theme(plot.title = element_text(size = 10)) +
    facet_wrap(~ Task),  
  scatterplot_SCI_density,
  density_SCI,
  # PSP plots (bottow row)
  ggplot(predictedData_PSP_summary, aes(x = Age, y = predictedValues, color = Cohort)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = Cohort), alpha = 0.2) +
    scale_color_viridis(discrete = TRUE) +
    scale_fill_viridis(discrete = TRUE) +
    labs(title = "Predicted PSP Trends by Age and Cohort, per Task", 
         x = "Age", 
         y = "Average PSP") +
    theme_minimal() +  
    theme(
      panel.grid.major = element_line(color = "grey80"),  
      panel.grid.minor = element_line(color = "grey90"), 
      panel.background = element_rect(fill = "white", color = NA),  
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    theme(plot.title = element_text(size = 10)) + 
    facet_wrap(~ Task),  
  scatterplot_PSP_density,   
  density_PSP,               
  
  ncol = 3  # layout should have 3 columns so arrangement is 2 x 3
)


#### PROPORTIONS OF LOW RAW VALUES ####
#### ---- Proportion of values below adult thresholds ####
#### ---- ---- Overall ####
# Calculate proportions for each Age group
dataProportions <- dataExpOrig %>%
  group_by(Age) %>%
  summarise(
    SCI_below_0_8 = mean(Average_SCI < 0.8),
    PSP_below_0_1 = mean(Average_PSP < 0.1)
  ) %>%
  gather(key = "Condition", value = "Proportion", SCI_below_0_8, PSP_below_0_1)

# Plot the bar chart with the viridis color map and custom labels
ggplot(dataProportions, aes(x = factor(Age), y = Proportion, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d(
    labels = c("Avg SCI below 0.8", "Avg PSP below 0.1") 
  ) +
  labs(
    x = "Age",
    y = "Proportion",
    title = "Proportion of Average SCI and PSP values below adult thresholds",
    fill = "Condition"
  ) +
  theme_minimal()

#### ---- ---- Split by Cohort ####
# Calculate proportions for each Age and Cohort group
dataProportions <- dataExpOrig %>%
  group_by(Age, Cohort) %>%
  summarise(
    SCI_below_0_8 = mean(Average_SCI < 0.8),
    PSP_below_0_1 = mean(Average_PSP < 0.1)
  ) %>%
  gather(key = "Condition", value = "Proportion", SCI_below_0_8, PSP_below_0_1)

# Create a new variable for grouping
dataProportions$Condition_Cohort <- paste(dataProportions$Condition, dataProportions$Cohort, sep = " - ")

# Plot the bar chart with four bars per Age and Cohort and custom key labels
ggplot(dataProportions, aes(x = factor(Age), y = Proportion, fill = Condition_Cohort)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d(
    labels = c("Avg SCI below 0.8 - Gambia", "Avg PSP below 0.1 - Gambia", "Avg SCI below 0.8 - UK", "Avg PSP below 0.1 - UK")  # Custom labels
  ) +
  labs(
    x = "Age",
    y = "Proportion",
    title = "Adult thresholds",
    fill = "Condition - Cohort"
  ) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 0.35)) +
  scale_x_discrete(labels = c("5", "8", "12", "18", "24"))


#### ---- Proportion of values below HALF OF adult thresholds ####
#### ---- ---- Overall ####
# Calculate proportions for each Age group
dataProportions <- dataExpOrig %>%
  group_by(Age) %>%
  summarise(
    SCI_below_0_4 = mean(Average_SCI < 0.4),
    PSP_below_0_05 = mean(Average_PSP < 0.05)
  ) %>%
  gather(key = "Condition", value = "Proportion", SCI_below_0_4, PSP_below_0_05)

# Plot the bar chart with the viridis color map and custom labels
ggplot(dataProportions, aes(x = factor(Age), y = Proportion, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d(
    labels = c("Avg SCI below 0.4", "Avg PSP below 0.05")
  ) +
  labs(
    x = "Age",
    y = "Proportion",
    title = "Proportion of Average SCI and PSP values below half of adult thresholds",
    fill = "Condition"
  ) +
  theme_minimal()


#### ---- ---- Split by Cohort ####
# Calculate proportions for each Age and Cohort group
dataProportions <- dataExpOrig %>%
  group_by(Age, Cohort) %>%
  summarise(
    SCI_below_0_4 = mean(Average_SCI < 0.4),
    PSP_below_0_05 = mean(Average_PSP < 0.05)
  ) %>%
  gather(key = "Condition", value = "Proportion", SCI_below_0_4, PSP_below_0_05)

# Create a new variable for grouping
dataProportions$Condition_Cohort <- paste(dataProportions$Condition, dataProportions$Cohort, sep = " - ")

# Plot the bar chart with four bars per Age and Cohort and custom key labels
ggplot(dataProportions, aes(x = factor(Age), y = Proportion, fill = Condition_Cohort)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d(
    labels = c("Avg SCI below 0.4 - Gambia", "Avg PSP below 0.05 - Gambia", "Avg SCI below 0.4 - UK", "Avg PSP below 0.05 - UK")  # Custom labels
  ) +
  labs(
    x = "Age",
    y = "Proportion",
    title = "Proportion of Average SCI and PSP values below half of adult thresholds, by Cohort",
    fill = "Condition - Cohort"
  ) +
  theme_minimal() +
  scale_y_continuous(limits = c(0.0, 0.35)) + 
  scale_x_discrete(labels = c("5", "8", "12", "18", "24"))



#### PLOT AVG SCI AND PSP BY Channel Location ####
# Boxplot for SCI 
ggplot(dataExp_filtered %>% filter(AP_Orientation != 0), aes(x = factor(AP_Orientation), y = Average_SCI, fill = factor(AP_Orientation))) +
  geom_boxplot(outlier.alpha = 0.3, outlier.size = 0.5) + 
  scale_fill_viridis_d(option = "D", name = "Channel Location") + 
  labs(title = "Box Plot of Average SCI by Channel Location",
       x = "Channel Location",
       y = "Average SCI") +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none", 
    plot.title = element_text(size = 12) 
  )

# Boxplot for PSP 
ggplot(dataExp_filtered %>% filter(AP_Orientation != 0), aes(x = factor(AP_Orientation), y = Average_PSP, fill = factor(AP_Orientation))) +
  geom_boxplot(outlier.alpha = 0.3, outlier.size = 0.5) +
  scale_fill_viridis_d(option = "D", name = "Channel Location") +
  labs(title = "Box Plot of Average PSP by Channel Location",
       x = "Channel Location",
       y = "Average PSP") +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none", 
    plot.title = element_text(size = 12) 
  )





