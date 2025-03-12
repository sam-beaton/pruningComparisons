#### LOAD PACKAGES ####
library(lme4)
library(foreach)
library(doParallel)
library(dplyr)
library(effectsize)
library(MuMIn)
library(ggplot2)
library(viridis)
library(arules)
library(data.table)
library(effects)
library(tidyr)
library(lmerTest)
library(scales)

#### SET UP PARAMS AND LOAD DATA ####

#set working directory (for saving outputs - bootstrapping is time-consuming!)
setwd("[path-to...]/pruningComparisonsR")

# Create a folder for saving outputs
if (!dir.exists("stats")) {
  dir.create("stats")
}

# load the data
data <- read.csv("[path-to...]/pruningComparisonsMatlab/stats/bright/overall/pruneMLMInputTable.csv")

# Scale the continuous variables
data <- data %>%
  mutate(
    SCI = scale(SCI),
    PSP = scale(PSP),
    Age = scale(Age),
    Percentage_Motion_Windows = scale(Percentage_Motion_Windows),
    ROI_SNR = scale(ROI_SNR),
    Channels_Retained = scale(Channels_Retained),
    Homer_Pruned_Total = scale(Homer_Pruned_Total)
  )


#check first few rows
head(data)
tail(data)
summary(data)

#model fitting
# Define the list of factors (excluding SCI and PSP as they are always included)
factors <- c("Age", "Cohort", "Task", "Homer_Pruned_Total", "Percentage_Motion_Windows")

# Bootstrapping params
# for reproducibility
set.seed(42)
# number of cores for parallel backend
numCores <- detectCores() - 1
# Set number of bootstrap iterations
nBoot <- 5000
# nBoot <- 500 #for testing code - comment out!
# divide and conquer:
chunkSize <- 250  # Number of iterations per chunk
# chunkSize  <- 25 #for testing code - comment out!
#early stopping:
stabilityThreshold <- 0.0000001   # Stability threshold (max SE change for fixed effects)

# Bootstrapping DF generation param(s)
outputDir <- "stats" # directory to save DF files
chunkSizeDF <- 100 # Number of rows per chunk
#can't use higher than 2250 with the models utilised here as exceeds RAM:
subsetSizeSNR <- 2000 #number of files to use to generate predictedData_SNR for plotting
subsetSizeCR <- 2000 #number of files to use to generate predictedData_CR for plotting
subsetSizeSNR <- 1500 #for testing
subsetSizeCR <- 1500 #for testing


#### DEFINE HELPER FUNCTIONS ####
# Helper function to generate interaction models
generateInteractionModels <- function(baseModel, interactionTerms, data, outcomeName) {
  models <- list()
  interactionCombinations <- unlist(
    lapply(1:length(interactionTerms), function(n) combn(interactionTerms, n, simplify = FALSE)),
    recursive = FALSE
  )
  for (i in seq_along(interactionCombinations)) {
    interactions <- paste(interactionCombinations[[i]], collapse = " + ")
    modelFormula <- as.formula(paste(baseModel, "+", interactions))
    modelName <- paste0(outcomeName, "InteractionModel", i)
    models[[modelName]] <- lmer(modelFormula, data = data)
  }
  return(models)
}

# Helper function to evaluate and rank models
evaluateModels <- function(models) {
  # Extract AIC and BIC
  aicBic <- sapply(models, function(model) c(AIC(model), BIC(model)))
  
  # Create results dataframe
  results <- data.frame(
    Model = names(models),
    AIC = aicBic[1, ],
    BIC = aicBic[2, ]
  )
  
  # Sort results
  resultsSortedAIC <- results[order(results$AIC), ]
  resultsSortedBIC <- results[order(results$BIC), ]
  
  return(list(
    resultsSortedAIC = resultsSortedAIC,
    resultsSortedBIC = resultsSortedBIC
  ))
}

# Function to compute predicted values, standard errors, and confidence intervals from fitted models
compute_predictions <- function(model, data, fixed_effect_formula) {
  vcov_matrix <- vcov(model)  
  model_formula <- as.formula(paste("~", fixed_effect_formula))  
  model_matrix <- model.matrix(model_formula, data)  
  
  # Get the rows used in the model (data used for fitting)
  model_data <- model.frame(model)  
  filtered_data <- data[rownames(model_data), ]  
  
  # Predictions and standard errors
  predicted <- predict(model, re.form = NA)
  se <- sqrt(rowSums((model_matrix %*% vcov_matrix)^2))
  
  # Add predictions and CI to the filtered data
  filtered_data <- filtered_data %>%
    mutate(
      predicted = predicted,
      lower_ci = predicted - 1.96 * se,
      upper_ci = predicted + 1.96 * se
    )
  
  return(filtered_data)
}


# Function to run bootstrapping with chunking and early stopping
runBootstrap <- function(modelEqn, modelBase, modelName, nBoot, chunkSize, stabilityThreshold) {
  nChunks <- ceiling(nBoot / chunkSize)  
  prevSe <- NULL  
  stable <- FALSE  
  
  for (chunk in 1:nChunks) {
    chunkStart <- (chunk - 1) * chunkSize + 1
    chunkEnd <- min(chunk * chunkSize, nBoot)
    cat("\nRunning chunk", chunk, "of", nChunks, "for", modelName, "\n")
    
    # Perform bootstrapping for the current chunk
    bootResultsChunk <- foreach(i = chunkStart:chunkEnd, .combine = rbind) %dopar% {
      indices <- sample(nrow(data), replace = TRUE)
      bootData <- data[indices, ]
      
      # Fit the model
      bootModel <- try(lmer(modelEqn, data = bootData), silent = TRUE)
      
      # Extract fixed effects;return NA if there's an error
      if (!inherits(bootModel, "try-error")) {
        fixef(bootModel)
      } else {
        rep(NA, length(fixef(modelBase)))
      }
    }
    
    # Save chunk results
    saveRDS(bootResultsChunk, paste0(outputDir, "/bootResults_", modelName, "_chunk_", chunk, ".rds"))
    
    # Combine results so far
    bootResults <- do.call(rbind, lapply(
      list.files(outputDir, pattern = paste0("bootResults_", modelName, "_chunk_.*\\.rds"), full.names = TRUE),
      readRDS
    ))
    
    # Calculate mean and SE of fixed effects
    fixefSummary <- apply(bootResults, 2, function(x) c(mean = mean(x, na.rm = TRUE), se = sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))))
    
    # Check stability based on SE change
    if (!is.null(prevSe)) {
      maxSeChange <- max(abs(fixefSummary["se", ] - prevSe))
      cat("Max SE change for", modelName, ":", maxSeChange, "\n")
      if (maxSeChange < stabilityThreshold) {
        stable <- TRUE
        cat("Estimates for", modelName, "have stabilised. Stopping early at chunk", chunk, "\n")
        break
      }
    }
    
    prevSe <- fixefSummary["se", ]  # Save SE for the next iteration
  }
  
  # Combine all results after bootstrapping
  finalResults <- do.call(rbind, lapply(
    list.files(outputDir, pattern = paste0("bootResults_", modelName, "_chunk_.*\\.rds"), full.names = TRUE),
    readRDS
  ))
  
  saveRDS(finalResults, paste0(outputDir, "/bootResults_", modelName, ".rds"))
  return(finalResults)
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
  results$Std_CI_lower <- std_effects - 1.96 * std_ses
  results$Std_CI_upper <- std_effects + 1.96 * std_ses
  results$Std_SE <- std_ses
  
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


#### DETERMINE ADDITIVE MODELS ####
#### ---- Fit simple models for SCI, PSP, and both combined #####
model1a_SNR <- lmer(ROI_SNR ~ SCI + (1 | ID), data = data) 
model1b_SNR <- lmer(ROI_SNR ~ PSP + (1 | ID), data = data) 
model2a_SNR <- lmer(ROI_SNR ~ SCI + PSP + (1 | ID), data = data)

model1a_CR <- lmer(Channels_Retained ~ SCI + (1 | ID), data = data) 
model1b_CR <- lmer(Channels_Retained ~ PSP + (1 | ID), data = data) 
model2a_CR <- lmer(Channels_Retained ~ SCI + PSP + (1 | ID), data = data) 

# Compare models based on AIC (ROI_SNR outcome)
AIC(model1a_SNR, model1b_SNR, model2a_SNR)
# Compare models based on AIC (Channels_Retained outcome)
AIC(model1a_CR, model1b_CR, model2a_CR)

# According to AIC and BIC, Model2a (SCI + PSP) is the best for both outcomes
# So use that as basis moving forwards

#### ---- Fit models with combinations of the other factors ####
# Initialise an empty list to store models
models_SNR <- list()
models_CR <- list()

# Set the outcome variables
outcome_SNR <- "ROI_SNR"
outcome_CR <- "Channels_Retained"

# Loop over the number of additional factors in each model (0 to length of factors)
for (n_factors in 0:length(factors)) {
  # Generate combinations of factors
  factor_combinations <- combn(factors, n_factors, simplify = FALSE)

  # Create models for each combination of factors
  for (i in seq_along(factor_combinations)) {
    # Extract predictors for this combination
    additional_predictors <- factor_combinations[[i]]
    all_predictors <- c("SCI", "PSP", additional_predictors)

    # Generate model names 
    model_name_SNR <- paste0("model", n_factors + 2, letters[i], "_SNR")
    model_name_CR <- paste0("model", n_factors + 2, letters[i], "_CR")  

    # Define the formula 
    formula_SNR <- as.formula(paste(outcome_SNR, "~", paste(all_predictors, collapse = " + "), "+ (1 | ID)"))
    formula_CR <- as.formula(paste(outcome_CR, "~", paste(all_predictors, collapse = " + "), "+ (1 | ID)"))

    # Fit the models and store them
    models_SNR[[model_name_SNR]] <- lmer(formula_SNR, data = data)
    models_CR[[model_name_CR]] <- lmer(formula_CR, data = data)
  }
}

#### ---- Compare SNR additive models ####
# Extract AIC and BIC for SNR models
aicBic_SNR <- sapply(models_SNR, function(model) c(AIC(model), BIC(model)))

# Combine results into a data frame 
results_SNR <- data.frame(
  Model = names(models_SNR),
  AIC = aicBic_SNR[1, ],
  BIC = aicBic_SNR[2, ]
)

# Sort results by AIC and BIC
results_sorted_AIC_SNR <- results_SNR[order(results_SNR$AIC), ]
results_sorted_BIC_SNR <- results_SNR[order(results_SNR$BIC), ]

# Print sorted results
cat("Results for ROI_SNR Models sorted by AIC:\n")
print(results_sorted_AIC_SNR)

cat("\nResults for ROI_SNR Models sorted by BIC:\n")
print(results_sorted_BIC_SNR)

# find models with lowest AIC and BIC
lowestAIC_SNR <- results_sorted_AIC_SNR[1, ]
lowestBIC_SNR <- results_sorted_BIC_SNR[1, ]

# Display the best models
cat("\nLowest AIC Model for ROI_SNR:\n")
print(lowestAIC_SNR)

cat("\nLowest BIC Model for ROI_SNR:\n")
print(lowestBIC_SNR)

# Show formulas for the top 5 models by AIC
cat("\nTop 5 Models for ROI_SNR by AIC:\n")
for (i in 1:min(5, nrow(results_sorted_AIC_SNR))) {
  model_name <- results_sorted_AIC_SNR$Model[i]
  cat(paste0("Model: ", model_name, "\n"))
  print(formula(models_SNR[[model_name]]))
  cat("\n")
}

# top 5 models by BIC
cat("\nTop 5 Models for ROI_SNR by BIC:\n")
for (i in 1:min(5, nrow(results_sorted_BIC_SNR))) {
  model_name <- results_sorted_BIC_SNR$Model[i]
  cat(paste0("Model: ", model_name, "\n"))
  print(formula(models_SNR[[model_name]]))
  cat("\n")
}

#### ---- Compare Channels retained additive models ####
# Extract AIC and BIC for CR models
aicBic_CR <- sapply(models_CR, function(model) c(AIC(model), BIC(model)))

# Combine results into a data frame 
results_CR <- data.frame(
  Model = names(models_CR),
  AIC = aicBic_CR[1, ],
  BIC = aicBic_CR[2, ]
)

# Sort results by AIC and BIC
results_sorted_AIC_CR <- results_CR[order(results_CR$AIC), ]
results_sorted_BIC_CR <- results_CR[order(results_CR$BIC), ]

# Print sorted results
cat("Results for Channels_Retained Models sorted by AIC:\n")
print(results_sorted_AIC_CR)

cat("\nResults for Channels_Retained Models sorted by BIC:\n")
print(results_sorted_BIC_CR)

# Find models with the lowest AIC and BIC
lowestAIC_CR <- results_sorted_AIC_CR[1, ]
lowestBIC_CR <- results_sorted_BIC_CR[1, ]

# Display the best models
cat("\nLowest AIC Model for Channels_Retained:\n")
print(lowestAIC_CR)

cat("\nLowest BIC Model for Channels_Retained:\n")
print(lowestBIC_CR)

# Show formulas for the top 5 models by AIC
cat("\nTop 5 Models for Channels_Retained by AIC:\n")
for (i in 1:min(5, nrow(results_sorted_AIC_CR))) {
  model_name <- results_sorted_AIC_CR$Model[i]
  cat(paste0("Model: ", model_name, "\n"))
  print(formula(models_CR[[model_name]]))
  cat("\n")
}

# formulas for the top 5 models by BIC
cat("\nTop 5 Models for Channels_Retained by BIC:\n")
for (i in 1:min(5, nrow(results_sorted_BIC_CR))) {
  model_name <- results_sorted_BIC_CR$Model[i]
  cat(paste0("Model: ", model_name, "\n"))
  print(formula(models_CR[[model_name]]))
  cat("\n")
}







#### ---- Compare models with interaction terms added ######
# Add appropriate interaction effects with Age to the best additive models and
# compare.
# Interaction terms taken from work modelling effects on Motion, Average SCI
# and Average PSP
# According to both the AIC and BIC above, the best models were:
# SNR Model: model6d_SNR
# ROI_SNR ~ SCI + PSP +
#          Age + Task + Homer_Pruned_Total + Percentage_Motion_Windows + (1 | ID)
# Model: model7a_CR
#Channels_Retained ~ SCI + PSP +
#                   Age + Cohort + Task + Homer_Pruned_Total + Percentage_Motion_Windows + (1 | ID)
# so use them as a basis

# Clear memory (except 'data' dataframe)
rm(list = setdiff(ls(), c("data", "evaluateModels", "generateInteractionModels", "runBootstrap")))

# Define the base models
additivemodel_SNR <- "ROI_SNR ~ SCI + PSP + Age + Task + Homer_Pruned_Total + Percentage_Motion_Windows + (1 | ID)"
additivemodel_CR <- "Channels_Retained ~ SCI + PSP + Age + Cohort + Task + Homer_Pruned_Total + Percentage_Motion_Windows + (1 | ID)"

# Interaction terms based on results from intermediary analyses into effects on
# Motion, Average SCI and Average PSP, plus consideration of appropriate
# hypothetical interactions (SCI * PSP)
# Cohort not included for SNR because not in additive model; not included for CR
# because no effect on SCI



# Define the interaction terms
interactionTermsSNR <- c("SCI * PSP", "Age * SCI", "Age * PSP", "Age * Percentage_Motion_Windows", "Age * Task")
interactionTermsCR <- c("SCI * PSP", "Age * SCI", "Age * PSP", "Age * Percentage_Motion_Windows", "Age * Task", "Age * Cohort")

# Generate models for ROI_SNR
modelsSNR <- generateInteractionModels(additivemodel_SNR, interactionTermsSNR, data, "SNR")

# Generate models for Channels_Retained 
modelsCR <- generateInteractionModels(additivemodel_CR, interactionTermsCR, data, "CR")


# Evaluate models for ROI_SNR
evaluationSNR <- evaluateModels(modelsSNR)

# Evaluate models for Channels_Retained
evaluationCR <- evaluateModels(modelsCR)

# Display results for ROI_SNR
cat("\n### Results for ROI_SNR Models ###\n")
cat("Models sorted by AIC:\n")
print(evaluationSNR$resultsSortedAIC)
cat("\nModels sorted by BIC:\n")
print(evaluationSNR$resultsSortedBIC)

# Display results for Channels_Retained
cat("\n### Results for Channels_Retained Models ###\n")
cat("Models sorted by AIC:\n")
print(evaluationCR$resultsSortedAIC)
cat("\nModels sorted by BIC:\n")
print(evaluationCR$resultsSortedBIC)

# Display best models and top 5 for each outcome
displayBestModels <- function(evaluation, models, outcomeName) {
  # Best by AIC
  bestAIC <- evaluation$resultsSortedAIC[1, ]
  cat(paste("\nLowest AIC Model for", outcomeName, ":\n"))
  print(bestAIC)
  print(formula(models[[bestAIC$Model]]))

  # Best by BIC
  bestBIC <- evaluation$resultsSortedBIC[1, ]
  cat(paste("\nLowest BIC Model for", outcomeName, ":\n"))
  print(bestBIC)
  print(formula(models[[bestBIC$Model]]))

  # Top 5 models by AIC
  cat(paste("\nTop 5 Models for", outcomeName, "by AIC:\n"))
  for (i in 1:min(5, nrow(evaluation$resultsSortedAIC))) {
    modelName <- evaluation$resultsSortedAIC$Model[i]
    cat(paste0("Model: ", modelName, "\n"))
    print(formula(models[[modelName]]))
    cat("\n")
  }

  # Top 5 models by BIC
  cat(paste("\nTop 5 Models for", outcomeName, "by BIC:\n"))
  for (i in 1:min(5, nrow(evaluation$resultsSortedBIC))) {
    modelName <- evaluation$resultsSortedBIC$Model[i]
    cat(paste0("Model: ", modelName, "\n"))
    print(formula(models[[modelName]]))
    cat("\n")
  }
}

# Display best models for ROI_SNR
displayBestModels(evaluationSNR, modelsSNR, "ROI_SNR")

# Display best models for Channels_Retained
displayBestModels(evaluationCR, modelsCR, "Channels_Retained")

# ALTHOUGH THESE MODELS HAVE LOWER BIC THAN THOSE WITHOUT INTERACTIONS, THE
# COMPLEXITY INTRODUCED AND LOSS OF BREVITY IS NOT WORTH THE SLIGHTLY BETTER MODEL
# FITS. THEREFORE BUILD MODELS BASED ONLY ON PRIMARY AIMS, THEORETICAL ASSUMPTIONS,
# AND RESULTS FROM SUBSIDIARY INVESTIGATIONS INTO MOTION, AVERAGE SCI AND AVERAGE PSP

#### ---- Add terms based on theoretical justifications from previous models ####

# Priority is placed on fixed effects of interest, their interaction terms,
# random intercepts and random slopes. Model complexity allows for the
# inclusion of covariates as fixed effects but increases number of effects to an
# unmanageable degree and risks overfitting, so are not included.

# Approach, in Order of importance:
#   (1)	Add all of Age*PSP, Age*SCI, SCI*PSP, as this is what we are interested
#       in primarily
# Then, working backwards back thru investigations into Avg SCI & PSP, then
# Motion, in order of importance/effect size of predictors:
#   (2)	Add Motion * PSP, as Motion has medium effect on Average PSP and is a
#       predictor which can be extrapolated to wider populations/research
#   (3)	Add Motion * Age, as Age has small effect on Motion BUT is a primary
#       predictor which can be extrapolated to wider populations/research
#   (4)	Add Motion * SCI, as Motion and SCI are primary
#       predictors which can be extrapolated to wider populations/research
#   (5) Add random effects for SCI and PSP with Cohort as Cohort has medium
#       effect on SCI and large effect on SCI, altering slopes
#   (6) Try random slope for PSP:Task, as Average PSP varies in slope for each task
#   (7)	Try random slope: (SCI | Percentage_Motion_Windows) due to wide
#       confidence interval in Avg SCI investigation and possibility that latent
#       motion falsely inflating correlation as well as acting as proxy for worse coupling

#   (1) Add all of Age*PSP, Age*SCI, SCI*PSP, as this is what we are interested
#       in primarily
model_SNR_Eqn_TEST <- "ROI_SNR ~ SCI * PSP + Age * SCI + Age * PSP +
                      Task + Cohort + Homer_Pruned_Total +
                      (1 | ID)"
model_CR_Eqn_TEST <- "Channels_Retained ~ SCI * PSP + Age * SCI + Age * PSP +
                      Task + Cohort + Homer_Pruned_Total +
                      (1 | ID)"
# Fit the models
model_SNR_TEST <- lmer(model_SNR_Eqn_TEST, data = data)
model_CR_TEST <- lmer(model_CR_Eqn_TEST, data = data)
# check model fit
BIC(model_SNR_TEST) # 1605156
BIC(model_CR_TEST) # 1220298

#   (2)	Add Motion * PSP, as Motion has medium effect on Average PSP and is a
#       predictor which can be extrapolated to wider populations/research
model_SNR_Eqn_TEST <- "ROI_SNR ~ SCI * PSP + Age * SCI + Age * PSP +
                      Percentage_Motion_Windows * PSP +
                      Task + Cohort + Homer_Pruned_Total +(1 | ID)"
model_CR_Eqn_TEST <- "Channels_Retained ~ SCI * PSP + Age * SCI + Age * PSP +
                      Percentage_Motion_Windows * PSP +
                      Task + Cohort + Homer_Pruned_Total + (1 | ID)"
# Fit the models
model_SNR_TEST <- lmer(model_SNR_Eqn_TEST, data = data)
model_CR_TEST <- lmer(model_CR_Eqn_TEST, data = data)
# check model fit
BIC(model_SNR_TEST) # 1452076. Better - keep.
BIC(model_CR_TEST) # 1104364. Better- keep.

#   (3) Add Motion * Age, as Motion has small effect on Motion BUT is a primary
#       predictor which can be extrapolated to wider populations/research
model_SNR_Eqn_TEST <- "ROI_SNR ~ SCI * PSP + Age * SCI + Age * PSP +
                      Percentage_Motion_Windows * PSP + Percentage_Motion_Windows * Age +
                      Task + Cohort + Homer_Pruned_Total + (1 | ID)"
model_CR_Eqn_TEST <- "Channels_Retained ~ SCI * PSP + Age * SCI + Age * PSP +
                      Percentage_Motion_Windows * PSP + Percentage_Motion_Windows * Age +
                      Task + Cohort + Homer_Pruned_Total + (1 | ID)"
# Fit the models
model_SNR_TEST <- lmer(model_SNR_Eqn_TEST, data = data)
model_CR_TEST <- lmer(model_CR_Eqn_TEST, data = data)
# check model fit
BIC(model_SNR_TEST) # 1450169. Better - keep.
BIC(model_CR_TEST) # 1099041. Better - keep.

#   (4)	Add Motion * SCI, as Motion and SCI are primary
#       predictors which can be extrapolated to wider populations/research
model_SNR_Eqn_TEST <- "ROI_SNR ~ SCI * PSP + Age * SCI + Age * PSP +
                      Percentage_Motion_Windows * SCI + Percentage_Motion_Windows * PSP + Percentage_Motion_Windows * Age +
                      Task + Cohort + Homer_Pruned_Total + (1 | ID)"
model_CR_Eqn_TEST <- "Channels_Retained ~ SCI * PSP + Age * SCI + Age * PSP +
                      Percentage_Motion_Windows * SCI + Percentage_Motion_Windows * PSP + Percentage_Motion_Windows * Age +
                      Task + Cohort + Homer_Pruned_Total + (1 | ID)"
# Fit the models
model_SNR_TEST <- lmer(model_SNR_Eqn_TEST, data = data)
model_CR_TEST <- lmer(model_CR_Eqn_TEST, data = data)
# check model fit
BIC(model_SNR_TEST) # 1450186. Worse but keep due to small change and advantages for looking into primary effects
BIC(model_CR_TEST) # 1098959. Better - keep.

#   (5a) Add random slopes for SCI and PSP with Cohort as Cohort has medium
#       effect on SCI and large effect on SCI, altering slopes
model_SNR_Eqn_TEST <- "ROI_SNR ~ SCI * PSP + Age * SCI + Age * PSP +
                      Percentage_Motion_Windows * SCI + Percentage_Motion_Windows * PSP + Percentage_Motion_Windows * Age +
                      Task + Cohort + Homer_Pruned_Total +
                      (1 + SCI + PSP | Cohort) +
                      (1 | ID)"
model_CR_Eqn_TEST <- "Channels_Retained ~ SCI * PSP + Age * SCI + Age * PSP +
                      Percentage_Motion_Windows * SCI + Percentage_Motion_Windows * PSP + Percentage_Motion_Windows * Age +
                      Task + Cohort + Homer_Pruned_Total +
                      (1 + SCI + PSP | Cohort) +
                      (1 | ID)"
# Fit the models
model_SNR_TEST <- lmer(model_SNR_Eqn_TEST, data = data) # didn't converge
model_CR_TEST <- lmer(model_CR_Eqn_TEST, data = data) # didn't converge

#   (5b) Add random intercepts for SCI and PSP with Cohort
model_SNR_Eqn_TEST <- "ROI_SNR ~ SCI * PSP + Age * SCI + Age * PSP +
                      Percentage_Motion_Windows * SCI + Percentage_Motion_Windows * PSP + Percentage_Motion_Windows * Age +
                      (1 | Cohort) + (1 | SCI:Cohort) + (1 | PSP:Cohort) +
                      Task + Cohort + Homer_Pruned_Total + (1 | ID)"
model_CR_Eqn_TEST <- "Channels_Retained ~ SCI * PSP + Age * SCI + Age * PSP +
                      Percentage_Motion_Windows * SCI + Percentage_Motion_Windows * PSP + Percentage_Motion_Windows * Age +
                      (1 | Cohort) + (1 | SCI:Cohort) + (1 | PSP:Cohort) +
                      Task + Cohort + Homer_Pruned_Total + (1 | ID)"
# Fit the models
model_SNR_TEST <- lmer(model_SNR_Eqn_TEST, data = data) # didn't converge
model_CR_TEST <- lmer(model_CR_Eqn_TEST, data = data) # didn't converge

#   (5c) Add random intercepts for *JUST* SCI and PSP with Cohort
model_SNR_Eqn_TEST <- "ROI_SNR ~ SCI * PSP + Age * SCI + Age * PSP +
                      Percentage_Motion_Windows * SCI + Percentage_Motion_Windows * PSP + Percentage_Motion_Windows * Age +
                      Task + Cohort + Homer_Pruned_Total +
                      (1 | SCI:Cohort) + (1 | PSP:Cohort) + (1 | ID)"
model_CR_Eqn_TEST <- "Channels_Retained ~ SCI * PSP + Age * SCI + Age * PSP +
                      Percentage_Motion_Windows * SCI + Percentage_Motion_Windows * PSP + Percentage_Motion_Windows * Age +
                      Task + Cohort + Homer_Pruned_Total +
                      (1 | SCI:Cohort) + (1 | PSP:Cohort) + (1 | ID)"
# Fit the models
model_SNR_TEST <- lmer(model_SNR_Eqn_TEST, data = data)
model_CR_TEST <- lmer(model_CR_Eqn_TEST, data = data)
# check model fit
BIC(model_SNR_TEST) # 1449667. Better - keep.
BIC(model_CR_TEST) # 1043455. Better - keep.

#   (6) Try random slope for PSP:Task, as Average PSP varies in slope for each task
model_SNR_Eqn_TEST <- "ROI_SNR ~ SCI * PSP + Age * SCI + Age * PSP +
                      Percentage_Motion_Windows * SCI + Percentage_Motion_Windows * PSP + Percentage_Motion_Windows * Age +
                      Task + Cohort + Homer_Pruned_Total +
                      (1 | SCI:Cohort) + (1 | PSP:Cohort) + (1 | ID) +
                      (1 + PSP | Task)"
model_CR_Eqn_TEST <- "Channels_Retained ~ SCI * PSP + Age * SCI + Age * PSP +
                      Percentage_Motion_Windows * SCI + Percentage_Motion_Windows * PSP + Percentage_Motion_Windows * Age +
                      Task + Cohort + Homer_Pruned_Total +
                      (1 | SCI:Cohort) + (1 | PSP:Cohort) + (1 | ID) +
                      (1 + PSP | Task)"
# Fit the models
model_SNR_TEST <- lmer(model_SNR_Eqn_TEST, data = data) # didn't converge
model_CR_TEST <- lmer(model_CR_Eqn_TEST, data = data) # didn't converge

#   (7)	Add random slope: (SCI | Percentage_Motion_Windows) due to wide
#       confidence interval in Avg SCI investigation and possibility that latent
#       motion falsely inflating correlation as well as acting as proxy for worse coupling
model_SNR_Eqn_TEST <- "ROI_SNR ~ SCI * PSP + Age * SCI + Age * PSP +
                      Percentage_Motion_Windows * SCI + Percentage_Motion_Windows * PSP + Percentage_Motion_Windows * Age +
                      Task + Cohort + Homer_Pruned_Total +
                      (1 | SCI:Cohort) + (1 | PSP:Cohort) + (1 | ID) +
                      (SCI | Percentage_Motion_Windows)"
model_CR_Eqn_TEST <- "Channels_Retained ~ SCI * PSP + Age * SCI + Age * PSP +
                      Percentage_Motion_Windows * SCI + Percentage_Motion_Windows * PSP + Percentage_Motion_Windows * Age +
                      Task + Cohort + Homer_Pruned_Total +
                      (1 | SCI:Cohort) + (1 | PSP:Cohort) + (1 | ID) +
                      (SCI | Percentage_Motion_Windows)"
# Fit the models
model_SNR_TEST <- lmer(model_SNR_Eqn_TEST, data = data)
model_CR_TEST <- lmer(model_CR_Eqn_TEST, data = data)
# check model fit
BIC(model_SNR_TEST) # # didn't converge
BIC(model_CR_TEST) # didn't converge


#### BEST MODELS ####
#### ---- Redefine final models and fit them ####
# 1st row: main outcomes of interest
# 2nd row: interaction effects, based on predictors from models for
# Motion, Avg SCI and Avg PSP
# 3rd row: covariates
# 4th row (CR): random intercepts for SCI and PSP per cohort, and participant

# Best fitting model eqns:
model_SNR_Eqn <- "ROI_SNR ~ SCI * PSP + Age * SCI + Age * PSP +
                      Percentage_Motion_Windows * SCI + Percentage_Motion_Windows * PSP + Percentage_Motion_Windows * Age +
                      Task + Cohort + Homer_Pruned_Total +
                      (1 | SCI:Cohort) + (1 | PSP:Cohort) + (1 | ID)"
model_CR_Eqn <- "Channels_Retained ~ SCI * PSP + Age * SCI + Age * PSP +
                      Percentage_Motion_Windows * SCI + Percentage_Motion_Windows * PSP + Percentage_Motion_Windows * Age +
                      Task + Cohort + Homer_Pruned_Total +
                      (1 | SCI:Cohort) + (1 | PSP:Cohort) + (1 | ID)"
# Fit the models
model_SNR <- lmer(model_SNR_Eqn, data = data)
model_CR<- lmer(model_CR_Eqn, data = data)


#### ---- Generate estimates for comparison later ####
# Extract fixed effects and standard errors
fixed_effects_SNR <- fixef(model_SNR)  # Fixed effects coefficients
se_SNR <- sqrt(diag(vcov(model_SNR)))  # Standard errors for fixed effects

fixed_effects_CR <- fixef(model_CR)  
se_CR <- sqrt(diag(vcov(model_CR))) 

# Calculate t-values for each effect
t_values_SNR <- fixed_effects_SNR / se_SNR
t_values_CR <- fixed_effects_CR / se_CR

# calculate approximate p-values from the t-distribution (more precise methods could be used)
# degrees of freedom tricky with mixed models: assume that df = n - number of fixed effects
# df = large sample size assumption
p_values_SNR <- 2 * (1 - pt(abs(t_values_SNR), df = 1000))  # Two-tailed p-value
p_values_CR <- 2 * (1 - pt(abs(t_values_CR), df = 1000))  

# Combine fixed effects, t-values, and p-values into data frames
fixed_effects_SNR_df <- data.frame(
  Effect = names(fixed_effects_SNR),
  Estimate = fixed_effects_SNR,
  Std.Error = se_SNR,
  t_value = t_values_SNR,
  p_value = p_values_SNR
)

fixed_effects_CR_df <- data.frame(
  Effect = names(fixed_effects_CR),
  Estimate = fixed_effects_CR,
  Std.Error = se_CR,
  t_value = t_values_CR,
  p_value = p_values_CR
)

# View the results for both models
print(fixed_effects_SNR_df)
print(fixed_effects_CR_df)


#### EXAMINE RESIDUALS ####
# Extract residuals
residuals_SNR <- residuals(model_SNR)
residuals_CR <- residuals(model_CR)

# Plot QQ-plot and histogram for model_SNR residuals
par(mfrow = c(1, 1))  # single plot layout

# QQ-plot for model_SNR
qqnorm(residuals_SNR, main = "QQ-Plot for Residuals (model_SNR)")
qqline(residuals_SNR, col = "red", lwd = 2)

# Histogram for model_SNR
hist(residuals_SNR, breaks = 100, main = "Histogram of Residuals (model_SNR)", xlab = "Residuals", col = "lightblue")

# Plot QQ-plot and histogram for model_CR residuals
par(mfrow = c(1, 1))  # single plot layout

# QQ-plot for model_CR
qqnorm(residuals_CR, main = "QQ-Plot for Residuals (model_CR)")
qqline(residuals_CR, col = "red", lwd = 2)

# Histogram for model_CR
hist(residuals_CR, breaks = 100, main = "Histogram of Residuals (model_CR)", xlab = "Residuals", col = "lightblue")

# Residual vs. Fitted Plot
plot(fitted(model_SNR), resid(model_SNR), main = "Residuals vs. Fitted values for SNR",
     xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

plot(fitted(model_CR), resid(model_CR), main = "Residuals vs. Fitted values for CR",
     xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

# Find influential points
cooks.distance <- cooks.distance(model_SNR)
plot(cooks.distance, main = "Cook's Distance (SNR model)",
     ylab = "Cook's Distance", xlab = "Index")
abline(h = 4/length(cooks.distance), col = "red")  # Common threshold for Cook's distance used

cooks.distance <- cooks.distance(model_CR)
plot(cooks.distance, main = "Cook's Distance (CR model)",
     ylab = "Cook's Distance", xlab = "Index")
abline(h = 4/length(cooks.distance), col = "red")  # Common threshold as above

# There are outliers as seen in the above plots. Transformations did not work
# (Box Cox the best but still doesn't fully solve issue)
# Try robust LMMs and bootstrapping

#### ROBUST LMMs (convergence and performance issues) ####

# SNR
robustModel_SNR <- rlmer(model_SNR_Eqn, data = data)
summary(robustModel_SNR)

# SNR
robustModel_CR <- rlmer(model_CR_Eqn, data = data)
summary(robustModel_CR)



#### BOOTSTRAPPING ####
#### ---- Run Bootstrapping ####

# --- SNR --- #
# Set up parallel processing
registerDoParallel(numCores)

# Bootstrapping
if (!file.exists("stats/bootResults_SNR.rds")) {
  # Time the SNR bootstrapping process
  startTime <- Sys.time() 

  # Run bootstrapping for SNR model
  bootResults_SNR <- runBootstrap(model_SNR_Eqn, model_SNR_Eqn, "SNR", nBoot, chunkSize, stabilityThreshold)

  endTime <- Sys.time() 
  timeTaken <- endTime - startTime  # Calculate time taken
  timeTaken
} else {
  cat("SNR bootstrapping results already exist. Skipping SNR bootstrapping.\n")
  bootResults_SNR <- readRDS("stats/bootResults_SNR.rds")
}

# Stop parallel processing
stopImplicitCluster()

# Clean up (delete) all chunk files with the specified pattern
file.remove(list.files("stats", pattern = "bootResults_.*_chunk_.*\\.rds", full.names = TRUE))

# --- CR --- #

# Set up parallel processing
registerDoParallel(numCores)

# Bootstrapping
if (!file.exists("stats/bootResults_CR.rds")) {
  # Time the CR bootstrapping process
  startTime <- Sys.time()  

  # Run bootstrapping for CR model
  bootResults_CR <- runBootstrap(model_CR_Eqn, model_CR_Eqn, "CR", nBoot, chunkSize, stabilityThreshold)

  endTime <- Sys.time()  
  timeTaken <- endTime - startTime  # Calculate time taken
  timeTaken
} else {
  cat("CR bootstrapping results already exist. Skipping CR bootstrapping.\n")
  bootResults_CR <- readRDS("stats/bootResults_CR.rds")
}

# Stop parallel processing
stopImplicitCluster()

# Clean up (delete) all chunk files with the specified pattern
file.remove(list.files("stats", pattern = "bootResults_.*_chunk_.*\\.rds", full.names = TRUE))



# Diagnostics: count the number of failed bootstraps (where NA is returned)
if (exists("bootResults_SNR")) {
  cat("\nNumber of failed bootstraps (SNR):", sum(rowSums(is.na(bootResults_SNR)) > 0), "\n")
}
if (exists("bootResults_CR")) {
  cat("\nNumber of failed bootstraps (CR):", sum(rowSums(is.na(bootResults_CR)) > 0), "\n")
}

#### ---- Display bootstrapping results for model_SNR ####
#### ---- Calculate confidence intervals for model_SNR
nCoef_SNR <- length(fixef(model_SNR))
SNR_CoefNames <- names(fixef(model_SNR))
SNR_CIs <- t(sapply(1:nCoef_SNR, function(i) {
  quantile(bootResults_SNR[, i], probs = c(0.025, 0.975), na.rm = TRUE)
}))
rownames(SNR_CIs) <- SNR_CoefNames

# Calculate standard errors from bootstrapped estimates for model_SNR
SNR_SEs <- apply(bootResults_SNR[, 1:nCoef_SNR], 2, sd, na.rm = TRUE)

# Calculate bootstrapped fixed effects (mean of bootstrapped estimates)
SNR_Estimates <- apply(bootResults_SNR, 2, mean, na.rm = TRUE)

# Create summary table for model_SNR
SNR_Summary <- data.frame(
  Estimate = SNR_Estimates,
  CI_lower = SNR_CIs[, "2.5%"],
  CI_upper = SNR_CIs[, "97.5%"],
  SE = SNR_SEs
)


# Print results for model_SNR
cat("\nModel SNR Results:\n")
print(SNR_Summary)

# Save estimates, confidence intervals, standard errors
saveRDS(SNR_CIs, "stats/SNR_CIs.rds")
saveRDS(SNR_SEs, "stats/SNR_SEs.rds")
saveRDS(SNR_Estimates, "stats/SNR_Estimates.rds")

# Save summary table
saveRDS(SNR_Summary, "stats/SNR_Summary.rds")


#### ---- Display bootstrapping results for model_CR ####
# Calculate confidence intervals for model_CR
nCoef_CR <- length(fixef(model_CR))
CR_CoefNames <- names(fixef(model_CR))
CR_CIs <- t(sapply(1:nCoef_CR, function(i) {
  quantile(bootResults_CR[, i], probs = c(0.025, 0.975), na.rm = TRUE)
}))
rownames(CR_CIs) <- CR_CoefNames

# Calculate standard errors from bootstrapped estimates for model_CR
CR_SEs <- apply(bootResults_CR[, 1:nCoef_CR], 2, sd, na.rm = TRUE)

# Calculate bootstrapped fixed effects
CR_Estimates <- apply(bootResults_CR, 2, mean, na.rm = TRUE)

# Create summary table for model_CR
CR_Summary <- data.frame(
  Estimate = CR_Estimates,
  CI_lower = CR_CIs[, "2.5%"],
  CI_upper = CR_CIs[, "97.5%"],
  SE = CR_SEs
)

# Print results for model_CR
cat("\nModel CR Results:\n")
print(CR_Summary)

# Save estimates, confidence intervals and standard errors
saveRDS(CR_CIs, "stats/CR_CIs.rds")
saveRDS(CR_SEs, "stats/CR_SEs.rds")
saveRDS(CR_Estimates, "stats/CR_Estimates.rds")

# Save summary table
saveRDS(CR_Summary, "stats/CR_Summary.rds")

#### ---- Compare to estimates from fitted model ####
# Model estimates and effect names for SNR
SNR_model_estimates <- fixef(model_SNR)
SNR_CoefNames <- names(fixef(model_SNR))
# Model estimates and effect names for CR
CR_model_estimates <- fixef(model_CR)
CR_CoefNames <- names(fixef(model_CR))

# Compare model estimates (fixed effects) to bootstrapped estimates
comparison_SNR <- data.frame(
  Coefficient = SNR_CoefNames,
  Model_Estimate = SNR_model_estimates,
  Bootstrap_Estimate = SNR_Estimates,
  Difference = SNR_model_estimates - SNR_Estimates
)

# Print comparison
print(comparison_SNR)

# plot the comparison
ggplot(comparison_SNR, aes(x = Model_Estimate, y = Bootstrap_Estimate)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  theme_minimal() +
  labs(title = "Model vs Bootstrap Estimates for SNR",
       x = "Model Estimate", y = "Bootstrap Estimate")

# Compare model estimates to bootstrapped estimates for CR
comparison_CR <- data.frame(
  Coefficient = CR_CoefNames,
  Model_Estimate = CR_model_estimates,
  Bootstrap_Estimate = CR_Estimates,
  Difference = CR_model_estimates - CR_Estimates
)

# Print comparison
print(comparison_CR)

# plot the comparison
ggplot(comparison_CR, aes(x = Model_Estimate, y = Bootstrap_Estimate)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  theme_minimal() +
  labs(title = "Model vs Bootstrap Estimates for CR",
       x = "Model Estimate", y = "Bootstrap Estimate")


# Check if model estimates for SNR are within the CI range from bootstrapping
within_CI_SNR <- sapply(1:length(SNR_CoefNames) , function(i) {
  SNR_model_estimates[i] >= SNR_CIs[i, "2.5%"] && SNR_model_estimates[i] <= SNR_CIs[i, "97.5%"]
})

cat("SNR Model estimates within CI: ", sum(within_CI_SNR), "/", length(SNR_CoefNames), " coefficients.\n")

# For CR
within_CI_CR <- sapply(1:length(CR_CoefNames) , function(i) {
  CR_model_estimates[i] >= CR_CIs[i, "2.5%"] && CR_model_estimates[i] <= CR_CIs[i, "97.5%"]
})

cat("CR Model estimates within CI: ", sum(within_CI_CR), "/", length(CR_CoefNames), " coefficients.\n")

#### ---- Calculate and display effect sizes for model_SNR ####
if (!file.exists("stats/SNR_Effects.rds")) {

  # Calculate effect sizes for model_SNR
  SNR_Effects <- standardiseEffects(model_SNR, SNR_Summary)

  # calculatemodel fit for SNR
  SNR_Fit <- r.squaredGLMM(model_SNR)

  # Print overall model fit for model_SNR
  cat("\nSNR Model Overall Fit:\n")
  print(SNR_Fit)

  # Print detailed results for model_SNR
  cat("\nSNR Model Effect Sizes:\n")
  print(SNR_Effects)

  # Save standardised effects
  saveRDS(SNR_Effects, "stats/SNR_Effects.rds")
  write.csv(SNR_Effects, "stats/SNR_Effects.csv")

  # Save model fit metric
  saveRDS(SNR_Fit, "stats/SNR_Fit.rds")

} else {

  cat("SNR Effect Size calculations already exist. Skipping calculation and displaying SNR Effects.\n")
  SNR_Effects <- readRDS("stats/SNR_Effects.rds")
  print(SNR_Effects)
}

#### ---- Calculate and display effect sizes for model_CR ####
if (!file.exists("stats/CR_Effects.rds")) {

  # Calculate effect sizes for model_CR
  CR_Effects <- standardiseEffects(model_CR, CR_Summary)

  # calculatemodel fit for CR
  CR_Fit <- r.squaredGLMM(model_CR)

  # Print overall model fit for model_CR
  cat("\nChannels Retained Model Overall Fit:\n")
  print(CR_Fit)

  # Print detailed results for model_CR
  cat("\nChannels Retained Model Effect Sizes:\n")
  print(CR_Effects)

  # Save standardised effects
  saveRDS(CR_Effects, "stats/CR_Effects.rds")
  write.csv(CR_Effects, "stats/CR_Effects.csv")

  # Save model fit metric
  saveRDS(CR_Fit, "stats/CR_Fit.rds")

} else {

  cat("CR Effect Size calculationss already exist. Skipping calculation and displaying  CR Effects.\n")
  CR_Effects <- readRDS("stats/CR_Effects.rds")
  print(CR_Effects)
}

#### PLOT RESULTS ####
#### ---- Load variables if necessary ####
# Check and load outputs if they are not already in the environment
if (!exists("bootResults_SNR")) {
  bootResults_SNR <- readRDS("stats/bootResults_SNR.rds")
}
if (!exists("bootResults_CR")) {
  bootResults_CR <- readRDS("stats/bootResults_CR.rds")
}
if (!exists("SNR_Summary")) {
  SNR_Summary <- readRDS("stats/SNR_Summary.rds")
}
if (!exists("CR_Summary")) {
  CR_Summary <- readRDS("stats/CR_Summary.rds")
}
if (!exists("SNR_Effects")) {
  SNR_Effects <- readRDS("stats/SNR_Effects.rds")
}
if (!exists("CR_Effects")) {
  CR_Effects <- readRDS("stats/CR_Effects.rds")
}
# 
# if (!exists("SNR_CIs")) {
#   SNR_CIs <- readRDS("stats/SNR_CIs.rds")
# }
# if (!exists("CR_CIs")) {
#   CR_CIs <- readRDS("stats/CR_CIs.rds")
# }
# 
# if (!exists("SNR_SEs")) {
#   SNR_SEs <- readRDS("stats/SNR_SEs.rds")
# }
# if (!exists("CR_SEs")) {
#   CR_SEs <- readRDS("stats/CR_SEs.rds")
# }
# 
# if (!exists("SNR_Estimates")) {
#   SNR_Estimates <- readRDS("stats/SNR_Estimates.rds")
# }
# if (!exists("CR_Estimates")) {
#   CR_Estimates <- readRDS("stats/CR_Estimates.rds")
# }
# 
# if (!exists("SNR_Fit")) {
#   SNR_Fit <- readRDS("stats/SNR_Fit.rds")
# }
# if (!exists("CR_Fit")) {
#   CR_Fit <- readRDS("stats/CR_Fit.rds")
# }
#### ---- FOREST PLOTS ####
#### ---- ---- Prepare Forest plot data ####

# New column names for bootResults_SCI for plots
newColnames <- c("Intercept", "SCI Threshold", "PSP Threshold", "Age", "Percentage of Motion", "Task (Social)", "Cohort (UK)", "CSE", 
                 "SCI Threshold * PSP Threshold", "Age * SCI Threshold", "Age * PSP Threshold", 
                 "SCI Threshold * Percentage of Motion", "PSP Threshold * Percentage of Motion", "Age * Percentage of Motion")
# Assign the new column names to bootResults_SCI
colnames(bootResults_SNR) <- newColnames

# New column names for bootResults_PSP for plots
newColnames <- c("Intercept", "SCI Threshold", "PSP Threshold", "Age", "Percentage of Motion", "Task (Social)", "Cohort (UK)", "CSE", 
                 "SCI Threshold * PSP Threshold", "Age * SCI Threshold", "Age * PSP Threshold", 
                 "SCI Threshold * Percentage of Motion", "PSP Threshold * Percentage of Motion", "Age * Percentage of Motion")

# Assign the new column names to bootResults_PSP
colnames(bootResults_CR) <- newColnames

#SNR
SNR_ForestPlot_Data <- data.frame(
  Coefficient = colnames(bootResults_SNR),
  Estimate = SNR_Summary$Estimate,
  CI_lower = SNR_Summary$CI_lower,
  CI_upper = SNR_Summary$CI_upper,
  SE = SNR_Summary$SE,
  Effect_Size = SNR_Effects$Effect_Size
)

#CR
CR_ForestPlot_Data <- data.frame(
  Coefficient = colnames(bootResults_CR),
  Estimate = CR_Summary$Estimate,
  CI_lower = CR_Summary$CI_lower,
  CI_upper = CR_Summary$CI_upper,
  SE = CR_Summary$SE, 
  Effect_Size = CR_Effects$Effect_Size
)

# Define the order of coefficients
orderPredictors_SNR <- c("Intercept", "SCI Threshold", "PSP Threshold", "SCI Threshold * PSP Threshold", 
                         "Percentage of Motion", "PSP Threshold * Percentage of Motion", "SCI Threshold * Percentage of Motion",
                         "Age", "Age * SCI Threshold", "Age * PSP Threshold", "Age * Percentage of Motion",
                         "CSE", "Cohort (UK)", "Task (Social)")


# Change the order of the factor levels in the SNR forest plot data
SNR_ForestPlot_Data$Coefficient <- factor(SNR_ForestPlot_Data$Coefficient, 
                                          levels = rev(orderPredictors_SNR))

# For CR data, change the order of coefficients
orderPredictors_CR <- c("Intercept", "SCI Threshold", "PSP Threshold", "SCI Threshold * PSP Threshold", 
                        "Percentage of Motion", "PSP Threshold * Percentage of Motion", "SCI Threshold * Percentage of Motion",
                        "Age", "Age * SCI Threshold", "Age * PSP Threshold", "Age * Percentage of Motion",
                        "CSE", "Cohort (UK)", "Task (Social)")

CR_ForestPlot_Data$Coefficient <- factor(CR_ForestPlot_Data$Coefficient, 
                                         levels = rev(orderPredictors_CR))

# Create a new column for effect type based on size and direction
# SNR
SNR_ForestPlot_Data <- SNR_ForestPlot_Data %>%
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
CR_ForestPlot_Data <- CR_ForestPlot_Data %>%
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

gc()
#### ---- ---- With Intercept ####
# SNR
ggplot(SNR_ForestPlot_Data, aes(x = Estimate, y = Coefficient)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray", linewidth = 1) + 
  geom_point(aes(color = Effect_Type, shape = Effect_Type), size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.4) +
  scale_color_manual(values = color_mapping) +
  scale_shape_manual(values = shape_mapping) +
  labs(title = "Bootstrapped Fixed Effects for ROI SNR",
       x = "Estimate with 95% CI", y = "Coefficient",
       color = "Standardised Effect Type", shape = "Standardised Effect Type") +
  theme_minimal()

# CR
ggplot(CR_ForestPlot_Data, aes(x = Estimate, y = Coefficient)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray", size = 1) + 
  geom_point(aes(color = Effect_Type, shape = Effect_Type), size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.4) +
  scale_color_manual(values = color_mapping) +
  scale_shape_manual(values = shape_mapping) +
  labs(title = "Bootstrapped Fixed Effects for Channels Retained",
       x = "Estimate with 95% CI", y = "Coefficient",
       color = "Standardised Effect Type", shape = "Standardised Effect Type") +
  theme_minimal()
#### ---- ---- Without Intercept ####

# Remove intercept from the SNR data
SNR_ForestPlot_Data_no_intercept <- SNR_ForestPlot_Data %>%
  filter(Coefficient != "Intercept")

# Forest plot without the intercept for SNR
ggplot(SNR_ForestPlot_Data_no_intercept, aes(x = Estimate, y = Coefficient)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray", size = 1) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray", size = 1) + 
  geom_point(aes(color = Effect_Type, shape = Effect_Type), size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.4) +
  scale_color_manual(values = color_mapping) +
  scale_shape_manual(values = shape_mapping) +
  labs(title = "Bootstrapped Fixed Effects for ROI SNR",
       x = "Estimate with 95% CI", y = "Coefficient",
       color = "Standardised Effect Type", shape = "Standardised Effect Type") +
  theme_minimal()

# Remove intercept from the PSP data
CR_ForestPlot_Data_no_intercept <- CR_ForestPlot_Data %>%
  filter(Coefficient != "Intercept")

# Forest plot without the intercept for PSP model
ggplot(CR_ForestPlot_Data_no_intercept, aes(x = Estimate, y = Coefficient)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray", size = 1) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray", size = 1) + 
  geom_point(aes(color = Effect_Type, shape = Effect_Type), size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.4) +
  scale_color_manual(values = color_mapping) +
  scale_shape_manual(values = shape_mapping) +
  labs(title = "Bootstrapped Fixed Effects for Channels Retained",
       x = "Estimate with 95% CI", y = "Coefficient",
       color = "Standardised Effect Type", shape = "Standardised Effect Type") +
  theme_minimal()

#### ---- ---- Without Percentage of Motion and PSP * Percentage of Motion (SNR) ####
# Remove PSP*Perc and Perc. from the SNR data
SNR_ForestPlot_Data_no_PSP_Perc <- SNR_ForestPlot_Data_no_intercept %>%
  filter(Coefficient != "PSP Threshold * Percentage of Motion") %>%
  filter(Coefficient != "Percentage of Motion")

# Forest plot without the intercept for SNR
ggplot(SNR_ForestPlot_Data_no_PSP_Perc, aes(x = Estimate, y = Coefficient)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray", size = 1) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray", size = 1) + 
  geom_point(aes(color = Effect_Type, shape = Effect_Type), size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.4) +
  scale_color_manual(values = color_mapping) +
  scale_shape_manual(values = shape_mapping) +
  labs(title = "Bootstrapped Fixed Effects for ROI SNR",
       x = "Estimate with 95% CI", y = "Coefficient",
       color = "Standardised Effect Type", shape = "Standardised Effect Type") +
  theme_minimal()
#### ---- ---- Without PSP * Percentage of Motion (CR) ####
# Remove PSP*Perc from the SNR data
CR_ForestPlot_Data_no_PSP_Perc <- CR_ForestPlot_Data_no_intercept %>%
  filter(Coefficient != "PSP Threshold * Percentage of Motion")

# Forest plot without the intercept for PSP model
ggplot(CR_ForestPlot_Data_no_PSP_Perc, aes(x = Estimate, y = Coefficient)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray", size = 1) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray", size = 1) + 
  geom_point(aes(color = Effect_Type, shape = Effect_Type), size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.4) +
  scale_color_manual(values = color_mapping) +
  scale_shape_manual(values = shape_mapping) +
  labs(title = "Bootstrapped Fixed Effects for Channels Retained",
       x = "Estimate with 95% CI", y = "Coefficient",
       color = "Standardised Effect Type", shape = "Standardised Effect Type") +
  theme_minimal() 

#### ----BOOTSTRAPPED ROI SNR PLOTS ####
#### ---- ---- Remove previous data ####
rm(list = c("CR_ForestPlot_Data", "CR_ForestPlot_Data_no_intercept", "CR_ForestPlot_Data_no_PSP_Perc"))
rm(list = c("SNR_ForestPlot_Data", "SNR_ForestPlot_Data_no_intercept", "SNR_ForestPlot_Data_no_PSP_Perc"))
gc()
#### ---- ---- Prepare data ####
if (file.exists(paste0(outputDir, "/predictedData_SNR_subset.rds"))) {
  cat("\nSNR predicted data already exists")
  if (!exists("predictedData_SNR_subset")) {
    cat(". Loading...\n")
    predictedData_SNR_subset <- readRDS(paste0(outputDir, "/predictedData_SNR_subset.rds"))
  } else {
    cat(" and loaded\n")
  }
} else {
  
  # Create a new data frame with all necessary variables
  # Compute average fixed effects from bootstrapping
  avgFixef_SNR <- colMeans(bootResults_SNR, na.rm = TRUE)
  
  # Create a grid of combinations of Age, Task, Cohort, etc. 
  predictedData_SNR <- expand.grid(
    SCI = unique(data$SCI),
    PSP = unique(data$PSP),
    Age = unique(data$Age),
    Task = unique(data$Task),
    Cohort = unique(data$Cohort),
    Homer_Pruned_Total = mean(data$Homer_Pruned_Total, na.rm = TRUE), #uses too much RAM otherwise
    Percentage_Motion_Windows = unique(data$Percentage_Motion_Windows)
  )
  
  # Initialise output directory for saving chunks
  if (!dir.exists(outputDir)) {
    dir.create(outputDir)
  }
  
  # Set up parallel backend for SNR
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  
  # Calculate chunks and split predictedData_SNR into chunks
  cat("Splitting data into chunks...\n")
  nChunks_SNR <- ceiling(nrow(predictedData_SNR) / chunkSizeDF)
  chunks_SNR <- split(predictedData_SNR, rep(1:nChunks_SNR, each = chunkSizeDF, length.out = nrow(predictedData_SNR)))
  cat("Data successfully split into", length(chunks_SNR), "chunks.\n")
  
  # Randomly select subset of chunks 
  selectedChunkIndices <- sample(seq_along(chunks_SNR), size = min(subsetSizeSNR, length(chunks_SNR)))
  selected_chunks_SNR <- chunks_SNR[selectedChunkIndices]
  
  # Process only selected chunks in parallel, saving each chunk
  foreach(chunkIdx = seq_along(selected_chunks_SNR), .packages = c("stats")) %dopar% {
    cat("\nProcessing chunk", chunkIdx, "of", length(selected_chunks_SNR), "\n")
    
    # Extract current chunk
    chunk <- selected_chunks_SNR[[chunkIdx]]
    
    # Create design matrix for SCI
    X_SNR <- model.matrix(~ SCI * PSP + Age * SCI + Age * PSP +
                            Percentage_Motion_Windows * SCI + Percentage_Motion_Windows * PSP + Percentage_Motion_Windows * Age +
                            Task + Cohort + Homer_Pruned_Total,
                          data = chunk)
    
    # Calculate predictions for all bootstrap iterations
    chunk$predictedValues <- apply(bootResults_SNR, 1, function(coefs) X_SNR %*% coefs)
    
    # Add confidence intervals to the chunk
    chunk$LowerCI <- apply(chunk$predictedValues, 1, quantile, probs = 0.025, na.rm = TRUE)
    chunk$UpperCI <- apply(chunk$predictedValues, 1, quantile, probs = 0.975, na.rm = TRUE)
    
    # Save chunk results to file to save memory
    saveRDS(chunk, file = paste0(outputDir, "/predictedData_SNR_chunk_", chunkIdx, ".rds"))
  }
  
  # Get all chunk files 
  chunkFiles <- list.files(outputDir, pattern = "predictedData_SNR_chunk_.*\\.rds", full.names = TRUE)
  
  # Combine all chunk files into a single data frame
  predictedData_SNR_subset <- do.call(rbind, lapply(chunkFiles, readRDS))
  
  # Save the combined data frame
  saveRDS(predictedData_SNR_subset, file = paste0(outputDir, "/predictedData_SNR_subset.rds"))
  
  # Clean up (delete) all chunk files once combined file is saved
  file.remove(chunkFiles)
  
  # Clean up
  stopCluster(cl)
}
#### ---- ---- (1a) Effect of SCI ####
# Aggregate data for ROI_SNR by SCI
predictedData_SNR_summary <- predictedData_SNR_subset %>%
  group_by(SCI) %>%
  summarise(
    predictedValues = mean(predictedValues[,1], na.rm = TRUE),
    LowerCI = mean(LowerCI, na.rm = TRUE),
    UpperCI = mean(UpperCI, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot ROI_SNR by SCI
ggplot(predictedData_SNR_summary, aes(x = factor(SCI), y = predictedValues)) +
  geom_point(size = 3, color = viridis(3)[1]) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, color = viridis(3)[2]) +
  labs(title = "Predicted ROI SNR by SCI Threshold values",
       x = "SCI",
       y = "ROI SNR") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )


#### ---- ---- (1b) Effect of PSP ####
# Aggregate data for ROI_SNR by PSP
predictedData_SNR_summary <- predictedData_SNR_subset %>%
  group_by(PSP) %>%
  summarise(
    predictedValues = mean(predictedValues[,1], na.rm = TRUE),
    LowerCI = mean(LowerCI, na.rm = TRUE),
    UpperCI = mean(UpperCI, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot ROI_SNR by PSP
ggplot(predictedData_SNR_summary, aes(x = factor(PSP), y = predictedValues)) +
  geom_point(size = 3, color = viridis(3)[1]) +  
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, color = viridis(3)[2]) +
  labs(title = "Predicted ROI SNR by PSP Threshold values", 
       x = "PSP Threshold", 
       y = "ROI SNR") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

#### ---- ---- (1c) Effect of PSP Threshold * SCI Threshold ####
# Aggregate data grouped by Age and SCI
predictedData_SNR_summary <- predictedData_SNR_subset %>%
  group_by(SCI, PSP) %>%
  summarise(
    predictedValues = mean(predictedValues[,1], na.rm = TRUE),
    LowerCI = mean(LowerCI, na.rm = TRUE),
    UpperCI = mean(UpperCI, na.rm = TRUE),
    .groups = 'drop'
  )

# Heatmap
ggplot(predictedData_SNR_subset, aes(x = factor(SCI), y = factor(PSP), fill = predictedValues[,1])) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma") +
  labs(title = "Heatmap of Predicted ROI SNR by SCI and PSP",
       x = "SCI Threshold",
       y = "PSP Threshold",
       fill = "Predicted ROI SNR") +
  theme_minimal()

# Plot with trend lines behind the scatterplot and slightly more faint
ggplot(predictedData_SNR_summary, aes(x = PSP, y = predictedValues, group = SCI, color = SCI)) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = SCI), alpha = 0.1) +  
  geom_point(alpha = 0.6) +  
  scale_color_viridis_c() +  
  scale_fill_viridis_c() +   
  labs(title = "Predicted SNR Trend by PSP Threshold, grouped by SCI Threshold", 
       x = "PSP Threshold", 
       y = "ROI SNR",
       color = "SCI Threshold",
       fill = "SCI Threshold") +
  theme_minimal() +
  theme(legend.position = "right")

# Plot with trend lines behind the scatterplot and slightly more faint
ggplot(predictedData_SNR_summary, aes(x = SCI, y = predictedValues, group = PSP, color = PSP)) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = PSP), alpha = 0.1) +  
  geom_point(alpha = 0.6) +  
  scale_color_viridis_c() + 
  scale_fill_viridis_c() +   
  labs(title = "Predicted SNR Trend by PSP Threshold, grouped by SCI Threshold", 
       x = "SCI Threshold", 
       y = "ROI SNR",
       color = "PSP Threshold",
       fill = "PSP Threshold") +
  theme_minimal() +
  theme(legend.position = "right")



#### ---- ---- (2a) Effect of Motion ####
# Aggregate data
predictedData_SNR_summary <- predictedData_SNR_subset %>%
  group_by(Percentage_Motion_Windows) %>%
  summarise(
    predictedValues = mean(predictedValues[,1], na.rm = TRUE),
    LowerCI = mean(LowerCI, na.rm = TRUE),
    UpperCI = mean(UpperCI, na.rm = TRUE)
  )

# Plot
ggplot(predictedData_SNR_summary, aes(x = Percentage_Motion_Windows, y = predictedValues)) +
  geom_point(alpha = 0.6) +  
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI), alpha = 0.2) +
  labs(title = "Predicted SNR Trend by Motion", 
       x = "Percentage of Windows with Motion", 
       y = "ROI SNR") +
  theme_minimal()

#### ---- ---- (2b) Effect of Motion and SCI Threshold####
# Aggregate data grouped by PSP and Motion Windows
predictedData_SNR_summary <- predictedData_SNR_subset %>%
  group_by(SCI, Percentage_Motion_Windows) %>%
  summarise(
    predictedValues = mean(predictedValues[,1], na.rm = TRUE),
    LowerCI = mean(LowerCI, na.rm = TRUE),
    UpperCI = mean(UpperCI, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot with different colors for each PSP level
ggplot(predictedData_SNR_summary, aes(x = Percentage_Motion_Windows, y = predictedValues, color = factor(SCI))) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = factor(SCI))) +
  scale_color_viridis(discrete = TRUE) +  
  scale_fill_viridis(discrete = TRUE) +   
  labs(title = "Predicted SNR Trend by Motion and SCI Threshold", 
       x = "Percentage of Windows with Motion", 
       y = "ROI SNR",
       color = "SCI Threshold",
       fill = "SCI Threshold") +
  theme_minimal() +
  theme(legend.position = "right")

#### ---- ---- (2c) Effect of Motion and PSP Threshold####
# Aggregate data grouped by PSP and Motion Windows
predictedData_SNR_summary <- predictedData_SNR_subset %>%
  group_by(PSP, Percentage_Motion_Windows) %>%
  summarise(
    predictedValues = mean(predictedValues[,1], na.rm = TRUE),
    LowerCI = mean(LowerCI, na.rm = TRUE),
    UpperCI = mean(UpperCI, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot with different colors for each PSP level
ggplot(predictedData_SNR_summary, aes(x = Percentage_Motion_Windows, y = predictedValues, color = factor(PSP))) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = factor(PSP))) +
  scale_color_viridis(discrete = TRUE) +  
  scale_fill_viridis(discrete = TRUE) +  
  labs(title = "Predicted SNR Trend by Motion and PSP Threshold", 
       x = "Percentage of Windows with Motion", 
       y = "ROI SNR",
       color = "PSP Threshold",
       fill = "PSP Threshold") +
  theme_minimal() +
  theme(legend.position = "right")


#### ---- ---- (3a) Effect of Age ####
# Aggregate data grouped by Age
predictedData_SNR_summary <- predictedData_SNR_subset %>%
  group_by(Age) %>%
  summarise(
    predictedValues = mean(predictedValues[,1], na.rm = TRUE),
    LowerCI = mean(LowerCI, na.rm = TRUE),
    UpperCI = mean(UpperCI, na.rm = TRUE)
  )

# Plot ROI_SNR by Age
ggplot(predictedData_SNR_summary, aes(x = Age, y = predictedValues)) +
  geom_point(alpha = 0.6, color = "blue") +  
  geom_smooth(method = "lm", color = "blue", se = FALSE) +  
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI), alpha = 0.2, fill = "blue") +  
  labs(
    title = "Predicted SNR Trend by Age", 
    x = "Age (months)", 
    y = "ROI SNR"
  ) +
  theme_minimal()

#### ---- ---- (3b) Effect of SCI Threshold * Age ####
# Aggregate data grouped by Age and SCI
predictedData_SNR_summary <- predictedData_SNR_subset %>%
  group_by(Age, SCI) %>%
  summarise(
    predictedValues = mean(predictedValues[,1], na.rm = TRUE),
    LowerCI = mean(LowerCI, na.rm = TRUE),
    UpperCI = mean(UpperCI, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot with different colors for each Age group
ggplot(predictedData_SNR_summary, aes(x = SCI, y = predictedValues, color = factor(Age))) +
  geom_point(alpha = 0.6) +  
  geom_smooth(method = "lm", se = TRUE, aes(fill = factor(Age))) +
  scale_color_viridis(discrete = TRUE) +  
  scale_fill_viridis(discrete = TRUE) +   
  labs(title = "Predicted SNR Trend by SCI Threshold and Age", 
       x = "SCI Threshold", 
       y = "ROI SNR",
       color = "Age",
       fill = "Age") +
  theme_minimal() +
  theme(legend.position = "right")

#### ---- ---- (3c) Effect of PSP Threshold * Age ####
# Aggregate data grouped by Age and PSP
predictedData_SNR_summary <- predictedData_SNR_subset %>%
  group_by(Age, PSP) %>%
  summarise(
    predictedValues = mean(predictedValues[,1], na.rm = TRUE),
    LowerCI = mean(LowerCI, na.rm = TRUE),
    UpperCI = mean(UpperCI, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot with different colors for each Age group
ggplot(predictedData_SNR_summary, aes(x = PSP, y = predictedValues, color = factor(Age))) +
  geom_point(alpha = 0.6) +  
  geom_smooth(method = "lm", se = TRUE, aes(fill = factor(Age))) +
  scale_color_viridis(discrete = TRUE) +  
  scale_fill_viridis(discrete = TRUE) +   
  labs(title = "Predicted SNR Trend by PSP Threshold and Age", 
       x = "PSP Threshold", 
       y = "ROI SNR",
       color = "Age",
       fill = "Age") +
  theme_minimal() +
  theme(legend.position = "right")
#### ---- ---- (3d) Effect of Motion and Age ####
# Aggregate data grouped by PSP and Motion Windows
predictedData_SNR_summary <- predictedData_SNR_subset %>%
  group_by(Age, Percentage_Motion_Windows) %>%
  summarise(
    predictedValues = mean(predictedValues[,1], na.rm = TRUE),
    LowerCI = mean(LowerCI, na.rm = TRUE),
    UpperCI = mean(UpperCI, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot with different colors for each Age
ggplot(predictedData_SNR_summary, aes(x = Percentage_Motion_Windows, y = predictedValues, color = factor(Age))) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = factor(Age)), alpha = 0.2) +
  scale_color_viridis(discrete = TRUE) +  
  scale_fill_viridis(discrete = TRUE) +   
  labs(title = "Predicted SNR Trend by Motion and Age", 
       x = "Percentage of Windows with Motion", 
       y = "ROI SNR",
       color = "Age",
       fill = "Age") +
  theme_minimal() +
  theme(legend.position = "right")


# #### ---- ---- (4) Effect of CSE ####
# # Aggregate data
# predictedData_SNR_summary <- predictedData_SNR_subset %>%
#   group_by(Homer_Pruned_Total) %>%
#   summarise(
#     predictedValues = mean(predictedValues[,1], na.rm = TRUE),
#     LowerCI = mean(LowerCI, na.rm = TRUE),
#     UpperCI = mean(UpperCI, na.rm = TRUE)
#   )
# 
# # Plot
# ggplot(predictedData_SNR_summary, aes(x = Homer_Pruned_Total, y = predictedValues)) +
#   geom_smooth(method = "lm", color = "blue", se = TRUE) +
#   geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI), alpha = 0.2) +
#   labs(title = "Predicted ROI SNR Trend by CSE", 
#        x = "CSE", 
#        y = "Channels Retained") +
#   theme_minimal()
#### ---- CR PLOTS ####
#### ---- ---- Remove previous data ####
rm(list = c("predictedData_SNR_subset", "predictedData_SNR_summary"))
gc()
#### ---- ---- Prepare data ####
# Check and load predicted CR results if they exist; create them if not
if (file.exists(paste0(outputDir, "/predictedData_CR_subset.rds"))) {
  cat("\nCR predicted data already exists")
  if (!exists("predictedData_CR_subset")) {
    cat(". Loading...\n")
    predictedData_CR_subset <- readRDS(paste0(outputDir, "/predictedData_CR_subset.rds"))
  } else {
    cat(" and loaded\n")
  }
} else {
  
  # Create a new data frame with all necessary variables
  # Compute average fixed effects from bootstrapping
  avgFixef_CR <- colMeans(bootResults_CR, na.rm = TRUE)
  
  # Create a grid of combinations of Age, Task, Cohort, etc.
  predictedData_CR <- expand.grid(
    SCI = unique(data$SCI),
    PSP = unique(data$PSP),
    Age = unique(data$Age),
    Task = unique(data$Task),
    Cohort = unique(data$Cohort),
    Homer_Pruned_Total = mean(data$Homer_Pruned_Total, na.rm = TRUE), #uses too much RAM otherwise
    Percentage_Motion_Windows = unique(data$Percentage_Motion_Windows)
  )
  
  # Initialise output directory for saving chunks
  if (!dir.exists(outputDir)) {
    dir.create(outputDir)
  }
  
  # Set up parallel backend for CR
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  
  # Calculate chunks and split predictedData_CR into chunks
  cat("Splitting data into chunks...\n")
  nChunks_CR <- ceiling(nrow(predictedData_CR) / chunkSizeDF)
  chunks_CR <- split(predictedData_CR, rep(1:nChunks_CR, each = chunkSizeDF, length.out = nrow(predictedData_CR)))
  cat("Data successfully split into", length(chunks_CR), "chunks.\n")
  
  # Randomly select a subset of chunks BEFORE processing
  selectedChunkIndices <- sample(seq_along(chunks_CR), size = min(subsetSizeCR, length(chunks_CR)))
  selected_chunks_CR <- chunks_CR[selectedChunkIndices]
  
  # Process only selected chunks in parallel, saving each chunk to a file
  foreach(chunkIdx = seq_along(selected_chunks_CR), .packages = c("stats")) %dopar% {
    cat("\nProcessing chunk", chunkIdx, "of", length(selected_chunks_CR), "\n")
    
    # Extract current chunk
    chunk <- selected_chunks_CR[[chunkIdx]]
    
    # Create design matrix for CR
    X_CR <- model.matrix(~ SCI * PSP + Age * SCI + Age * PSP +
                           Percentage_Motion_Windows * SCI + Percentage_Motion_Windows * PSP + Percentage_Motion_Windows * Age +
                           Task + Cohort + Homer_Pruned_Total,
                         data = chunk)
    
    # Calculate predictions for all bootstrap iterations
    chunk$predictedValues <- apply(bootResults_CR, 1, function(coefs) X_CR %*% coefs)
    
    # Add confidence intervals to the chunk
    chunk$LowerCI <- apply(chunk$predictedValues, 1, quantile, probs = 0.025, na.rm = TRUE)
    chunk$UpperCI <- apply(chunk$predictedValues, 1, quantile, probs = 0.975, na.rm = TRUE)
    
    # Save chunk results to file to save memory
    saveRDS(chunk, file = paste0(outputDir, "/predictedData_CR_chunk_", chunkIdx, ".rds"))
  }
  
  # Get all chunk files (which are now only from the pre-selected chunks)
  chunkFiles <- list.files(outputDir, pattern = "predictedData_CR_chunk_.*\\.rds", full.names = TRUE)
  
  # Combine the selected chunk files into a single data frame
  predictedData_CR_subset <- do.call(rbind, lapply(chunkFiles, readRDS))
  
  # Save the combined data frame
  saveRDS(predictedData_CR_subset, file = paste0(outputDir, "/predictedData_CR_subset.rds"))
  
  # Clean up (delete) all chunk files once combined file is saved
  file.remove(chunkFiles)
  
  # Clean up
  stopCluster(cl)
}

gc()
#### ---- ---- (1a) Effect of SCI ####
# Aggregate data for ROI_SNR by PSP
predictedData_CR_summary <- predictedData_CR_subset %>%
  group_by(SCI) %>%
  summarise(
    predictedValues = mean(predictedValues[,1], na.rm = TRUE),
    LowerCI = mean(LowerCI, na.rm = TRUE),
    UpperCI = mean(UpperCI, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot ROI_SNR by PSP
ggplot(predictedData_CR_summary, aes(x = factor(SCI), y = predictedValues)) +
  geom_point(size = 3, color = viridis(3)[1]) +  
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, color = viridis(3)[2]) +
  labs(title = "Predicted Channels Retained by SCI Threshold values", 
       x = "SCI", 
       y = "Channels Retained") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

#### ---- ---- (1b) Effect of PSP ####
# Aggregate data for ROI_SNR by PSP
predictedData_CR_summary <- predictedData_CR_subset %>%
  group_by(PSP) %>%
  summarise(
    predictedValues = mean(predictedValues[,1], na.rm = TRUE),
    LowerCI = mean(LowerCI, na.rm = TRUE),
    UpperCI = mean(UpperCI, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot ROI_SNR by PSP
ggplot(predictedData_CR_summary, aes(x = factor(PSP), y = predictedValues)) +
  geom_point(size = 3, color = viridis(3)[1]) +  
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, color = viridis(3)[2]) +
  labs(title = "Predicted Channels Retained by PSP Threshold values", 
       x = "PSP", 
       y = "Channels Retained") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

#### ---- ---- (1c) Effect of PSP Threshold * SCI Threshold ####
# Aggregate data grouped by Age and SCI
predictedData_CR_summary <- predictedData_CR_subset %>%
  group_by(SCI, PSP) %>%
  summarise(
    predictedValues = mean(predictedValues[,1], na.rm = TRUE),
    LowerCI = mean(LowerCI, na.rm = TRUE),
    UpperCI = mean(UpperCI, na.rm = TRUE),
    .groups = 'drop'
  )

# Heatmap
ggplot(predictedData_CR_subset, aes(x = factor(SCI), y = factor(PSP), fill = predictedValues[,1])) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma") +
  labs(title = "Heatmap of Predicted Channels Retained by SCI and PSP",
       x = "SCI Threshold",
       y = "PSP Threshold",
       fill = "Predicted Channels Retained") +
  theme_minimal()

# Plot with trend lines behind the scatterplot and slightly more faint
ggplot(predictedData_CR_summary, aes(x = PSP, y = predictedValues, group = SCI, color = SCI)) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = SCI), alpha = 0.1) +  
  geom_point(alpha = 0.6) + 
  scale_color_viridis_c() +  
  scale_fill_viridis_c() +   
  labs(title = "Predicted Channels Retained Trend by PSP Threshold, grouped by SCI Threshold", 
       x = "PSP Threshold", 
       y = "Channels Retained",
       color = "SCI Threshold",
       fill = "SCI Threshold") +
  theme_minimal() +
  theme(legend.position = "right")

# Plot with trend lines behind the scatterplot and slightly more faint
ggplot(predictedData_CR_summary, aes(x = SCI, y = predictedValues, group = PSP, color = PSP)) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = PSP), alpha = 0.1) +  
  geom_point(alpha = 0.6) + 
  scale_color_viridis_c() + 
  scale_fill_viridis_c() +   
  labs(title = "Predicted Channels Retained Trend by SCI Threshold, grouped by PSP Threshold", 
       x = "SCI Threshold", 
       y = "Channels Retained",
       color = "PSP Threshold",
       fill = "PSP Threshold") +
  theme_minimal() +
  theme(legend.position = "right")



#### ---- ---- (2a) Effect of Motion ####
# Aggregate data
predictedData_CR_summary <- predictedData_CR_subset %>%
  group_by(Percentage_Motion_Windows) %>%
  summarise(
    predictedValues = mean(predictedValues[,1], na.rm = TRUE),
    LowerCI = mean(LowerCI, na.rm = TRUE),
    UpperCI = mean(UpperCI, na.rm = TRUE)
  )

# Plot
ggplot(predictedData_CR_summary, aes(x = Percentage_Motion_Windows, y = predictedValues)) +
  geom_point(alpha = 0.6) +  
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI), alpha = 0.2) +
  labs(title = "Predicted Channels Retained Trend by Motion", 
       x = "Percentage of Windows with Motion", 
       y = "Channels Retained") +
  theme_minimal()

#### ---- ---- (2b) Effect of Motion and SCI Threshold####
# Aggregate data grouped by PSP and Motion Windows
predictedData_CR_summary <- predictedData_CR_subset %>%
  group_by(SCI, Percentage_Motion_Windows) %>%
  summarise(
    predictedValues = mean(predictedValues[,1], na.rm = TRUE),
    LowerCI = mean(LowerCI, na.rm = TRUE),
    UpperCI = mean(UpperCI, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot with different colors for each PSP level
ggplot(predictedData_CR_summary, aes(x = Percentage_Motion_Windows, y = predictedValues, color = factor(SCI))) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = factor(SCI))) +
  scale_color_viridis(discrete = TRUE) +  
  scale_fill_viridis(discrete = TRUE) +  
  labs(title = "Predicted Channels Retained Trend by Motion and SCI Threshold", 
       x = "Percentage of Windows with Motion", 
       y = "Channels Retained",
       color = "SCI Threshold",
       fill = "SCI Threshold") +
  theme_minimal() +
  theme(legend.position = "right")

#### ---- ---- (2c) Effect of Motion and PSP Threshold####
# Aggregate data grouped by PSP and Motion Windows
predictedData_CR_summary <- predictedData_CR_subset %>%
  group_by(PSP, Percentage_Motion_Windows) %>%
  summarise(
    predictedValues = mean(predictedValues[,1], na.rm = TRUE),
    LowerCI = mean(LowerCI, na.rm = TRUE),
    UpperCI = mean(UpperCI, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot with different colors for each PSP level
ggplot(predictedData_CR_summary, aes(x = Percentage_Motion_Windows, y = predictedValues, color = factor(PSP))) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = factor(PSP))) +
  scale_color_viridis(discrete = TRUE) +  
  scale_fill_viridis(discrete = TRUE) +   
  labs(title = "Predicted Channels Retained Trend by Motion and PSP Threshold", 
       x = "Percentage of Windows with Motion", 
       y = "Channels Retained",
       color = "PSP Threshold",
       fill = "PSP Threshold") +
  theme_minimal() +
  theme(legend.position = "right")

#### ---- ---- (2d) Effect of Age ####
# Aggregate data grouped by Age
predictedData_CR_summary <- predictedData_CR_subset %>%
  group_by(Age) %>%
  summarise(
    predictedValues = mean(predictedValues[,1], na.rm = TRUE),
    LowerCI = mean(LowerCI, na.rm = TRUE),
    UpperCI = mean(UpperCI, na.rm = TRUE)
  )

# Plot ROI_SNR by Age
ggplot(predictedData_CR_summary, aes(x = Age, y = predictedValues)) +
  geom_point(alpha = 0.6, color = "blue") +  
  geom_smooth(method = "lm", color = "blue", se = FALSE) + 
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI), alpha = 0.2, fill = "blue") + 
  labs(
    title = "Predicted CR Trend by Age", 
    x = "Age (months)", 
    y = "ROI SNR"
  ) +
  theme_minimal()
#### ---- ---- (3a) Effect of SCI Threshold * Age ####
# Aggregate data grouped by Age and SCI
predictedData_CR_summary <- predictedData_CR_subset %>%
  group_by(Age, SCI) %>%
  summarise(
    predictedValues = mean(predictedValues[,1], na.rm = TRUE),
    LowerCI = mean(LowerCI, na.rm = TRUE),
    UpperCI = mean(UpperCI, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot with different colors for each Age group
ggplot(predictedData_CR_summary, aes(x = SCI, y = predictedValues, color = factor(Age))) +
  geom_point(alpha = 0.6) +  
  geom_smooth(method = "lm", se = TRUE, aes(fill = factor(Age))) +
  scale_color_viridis(discrete = TRUE) +  
  scale_fill_viridis(discrete = TRUE) +  
  labs(title = "Predicted Channels Retained Trend by SCI Threshold and Age", 
       x = "SCI Threshold", 
       y = "Channels Retained",
       color = "Age",
       fill = "Age") +
  theme_minimal() +
  theme(legend.position = "right")
#### ---- ---- (3b) Effect of PSP Threshold * Age ####
# Aggregate data grouped by Age and PSP
predictedData_CR_summary <- predictedData_CR_subset %>%
  group_by(Age, PSP) %>%
  summarise(
    predictedValues = mean(predictedValues[,1], na.rm = TRUE),
    LowerCI = mean(LowerCI, na.rm = TRUE),
    UpperCI = mean(UpperCI, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot with different colors for each Age group
ggplot(predictedData_CR_summary, aes(x = PSP, y = predictedValues, color = factor(Age))) +
  geom_point(alpha = 0.6) +  
  geom_smooth(method = "lm", se = TRUE, aes(fill = factor(Age))) +
  scale_color_viridis(discrete = TRUE) +  
  scale_fill_viridis(discrete = TRUE) +  
  labs(title = "Predicted Channels Retained Trend by PSP Threshold and Age", 
       x = "PSP Threshold", 
       y = "ROI SNR",
       color = "Age",
       fill = "Age") +
  theme_minimal() +
  theme(legend.position = "right")
#### ---- ---- (3c) Effect of Motion and Age ####
# Aggregate data grouped by Age and Motion Windows
predictedData_CR_summary <- predictedData_CR_subset %>%
  group_by(Age, Percentage_Motion_Windows) %>%
  summarise(
    predictedValues = mean(predictedValues[,1], na.rm = TRUE),
    LowerCI = mean(LowerCI, na.rm = TRUE),
    UpperCI = mean(UpperCI, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot with different colors for each Age
ggplot(predictedData_CR_summary, aes(x = Percentage_Motion_Windows, y = predictedValues, color = factor(Age))) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = factor(Age))) +
  scale_color_viridis(discrete = TRUE) +  
  scale_fill_viridis(discrete = TRUE) +  
  labs(title = "Predicted Channels Retained Trend by Motion and Age", 
       x = "Percentage of Windows with Motion", 
       y = "Channels Retained",
       color = "Age",
       fill = "Age") +
  theme_minimal() +
  theme(legend.position = "right")
