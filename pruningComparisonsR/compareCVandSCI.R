
# install.packages("R.matlab")
library(R.matlab)
library(dplyr)
library(tidyr)
library(broom)
library(lmerTest)
library(ggplot2)
library(ggpubr)
library(car) # for Yeo-Johnson transformation



##### SET PARAMS #######
#Set to 1 to generate new tables:
genTables = 0
# Specify ages, tasks, and cohorts
ages <- c("05", "08", "12", "18", "24")
tasks <- c("hand", "social")
cohorts <- c("gm", "uk")

#CV directory
fullPath <- "[path-to...]/pruningComparisonsMatlab/stats/bright"
# SCI directory
sciPath <- "[path-to...]/pruningComparisonsMatlab/stats/brightSCIOnly"

#total channel number
totalChans = 68
#min % of channels required for infant inclusion (as decimal)
infInclThresh = 0.6
# alpha for significance
pThreshold = 0.05

numCombinations = length(cohorts)*length(ages)*length(tasks)


##### HELPER FuNCTIONS #####
# Function to assign significance stars in polots
assignStars <- function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("")  # No stars for p >= 0.05
  }
}

#### GENERATE TABLES ###################
#### ---- Generate tables, if necessary ####
if (genTables == 0) {
  # Check if the tables already exist
  loadCV <- function(tableDir) {
    tryCatch(
      {
      cvTable = read.table(paste(tableDir, 'cvTableCompareCvSciQt', sep = "/"))
      return(cvTable)
      },
      error = function(e) {
        message('Unable to load all tables')
      }
    )
  }
  loadSCI <- function(tableDir) {
    tryCatch(
      {
        sciTable = read.table(paste(tableDir, 'sciTableCompareCvSciQt', sep = "/"))
        return(sciTable)
      },
      error = function(e) {
        message('Unable to load all tables')
      }
    )
  }
  loadQT <- function(tableDir) {
    tryCatch(
      {
        qtTable = read.table(paste(tableDir, 'qtTableCompareCvSciQt', sep = "/"))
        return(qtTable)
      },
      error = function(e) {
        message('Unable to load all tables')
      }
    )
  }
  
  tableDir = '/Users/sambe/Documents/R_Projects/PreP/data'
  cvTable <- loadCV(tableDir)
  sciTable <- loadSCI(tableDir)
  qtTable <- loadQT(tableDir)
}

# Generate tables if specified by user or not able to load them
if (!exists("cvTable") && !exists("sciTable") && !exists("qtTable")) {

  #create CV table
  cvTable <- data.frame(
    Participant = numeric(),
    Cohort = character(),
    Task = character(),
    Age = numeric(),
    Channels = numeric(),
    ROI_SNR = numeric(),
    Infant_Excluded = numeric(),
    stringsAsFactors = FALSE # Prevent automatic conversion of strings to factors
  )
  
  #create SCI files table
  sciTable <- data.frame(
    Participant = numeric(),
    Cohort = character(),
    Task = character(),
    Age = numeric(),
    SCI_Param = numeric(),
    Channels = numeric(),
    ROI_SNR = numeric(),
    Infant_Excluded = numeric(),
    stringsAsFactors = FALSE # Prevent automatic conversion of strings to factors
  )
  
  #create QT files table
  qtTable <- data.frame(
    Participant = numeric(),
    Cohort = character(),
    Task = character(),
    Age = numeric(),
    SCI_Param = numeric(),
    PSP_Param = numeric(),
    Channels = numeric(),
    ROI_SNR = numeric(),
    Infant_Excluded = numeric(),
    stringsAsFactors = FALSE # Prevent automatic conversion of strings to factors
  )
  
  #start cohort loop
  for (iCohort in 1:length(cohorts)) {
    #start task loop
    for (iTask in 1:length(tasks)) {
      #start age loop
      for (iAge in 1:length(ages)) {
  
        ### CV data ###
        # load CV files
        cvSnrFilePath <- paste(fullPath, cohorts[iCohort], paste(ages[iAge], "mo", sep = ""), tasks[iTask], "prune", paste(tasks[iTask], "FullPruneSNR_CV_20.mat", sep=""), sep = "/")
        cvChanFilePath <- paste(fullPath, cohorts[iCohort], paste(ages[iAge], "mo", sep = ""), tasks[iTask], "prune", paste(tasks[iTask], "FullPruneChan_CV_20.mat", sep=""), sep = "/")
        cvSnrData <- readMat(cvSnrFilePath)
        cvChanData <- readMat(cvChanFilePath)
  
        # find number of participants
        numParts = length(cvChanData$channelsRetained)

        #create CV table
        cvCurrTable <- data.frame(
          Participant = rep(NA, numParts),
          Cohort = rep(NA, numParts),
          Task = rep(NA, numParts),
          Age = rep(NA, numParts),
          Channels = rep(NA, numParts),
          ROI_SNR = rep(NA, numParts),
          Infant_Excluded = rep(NA, numParts),
          stringsAsFactors = FALSE # Prevent automatic conversion of strings to factors
        )

        #pre-calculate row means so that they can be inserted into the table
        # inside the participant loop
        rowMeans <- apply(cvSnrData$snrMat, 1, function(x) mean(as.numeric(x), na.rm = TRUE))

        for (iPart in 1:numParts) {
          cvCurrTable$Participant[iPart] = iPart
          cvCurrTable$Cohort[iPart] = cohorts[iCohort]
          cvCurrTable$Task[iPart] = tasks[iTask]
          cvCurrTable$Age[iPart] = as.numeric(ages[iAge])
          cvCurrTable$Channels[iPart] = cvChanData$channelsRetained[iPart]
          if (cvCurrTable$Channels[iPart] <= floor(infInclThresh*totalChans)) {
            cvCurrTable$Infant_Excluded[iPart] = 1
          } else {
            cvCurrTable$Infant_Excluded[iPart] = 0
          }
          cvCurrTable$ROI_SNR[iPart] = rowMeans[iPart]
        }

        #append this version of the cvTable to the overall table
        cvTable <- rbind(cvTable, cvCurrTable)

        ### SCI Data ###
        #define directory
        sciDirectory <- paste(sciPath, cohorts[iCohort], paste(ages[iAge], "mo", sep = ""), tasks[iTask], "prune", sep = "/")
        # List all files with the specified suffix in the directory
        allMatFiles <- list.files(sciDirectory, pattern = paste0("*", ".mat"), full.names = TRUE)
        # Filter files that contain SNR values
        qtSnrFiles <- allMatFiles[grepl("SNR", basename(allMatFiles))]
        # Filter files that contain chan numbers
        qtChanFiles <- allMatFiles[grepl("Chan", basename(allMatFiles))]
        # Count the number of matching files
        numFiles <- length(qtSnrFiles)

        sciCurrTable <- data.frame(
          Participant = rep(NA, numParts*numFiles),
          Cohort = rep(NA, numParts*numFiles),
          Task = rep(NA, numParts*numFiles),
          Age = rep(NA, numParts*numFiles),
          SCI_Param = rep(NA, numParts*numFiles),
          Channels = rep(NA, numParts*numFiles),
          ROI_SNR = rep(NA, numParts*numFiles),
          Infant_Excluded = rep(NA, numParts*numFiles),
          stringsAsFactors = FALSE # Prevent automatic conversion of strings to factors
        )

        # initialise count variables for rows in tables
        iSCI = 1

        # start SCI files loop
        for (iFile in 1:numFiles) {
          qtSnrData <- readMat(qtSnrFiles[iFile])
          qtChanData <- readMat(qtChanFiles[iFile])

          #pre-calculate row means so that they can be inserted into the table
          # inside the participant loop
          rowMeans <- apply(qtSnrData$snrMat, 1, function(x) mean(as.numeric(x), na.rm = TRUE))

          for (iPart in 1:numParts) {

            sciCurrTable$Participant[iSCI] = iPart
            sciCurrTable$Cohort[iSCI] = cohorts[iCohort]
            sciCurrTable$Task[iSCI] = tasks[iTask]
            sciCurrTable$Age[iSCI] = as.numeric(ages[iAge])
            sciCurrTable$Channels[iSCI] = qtChanData$channelsRetained[iPart]
            if (sciCurrTable$Channels[iSCI] <= floor(infInclThresh*totalChans)) {
              sciCurrTable$Infant_Excluded[iSCI] = 1
            } else {
              sciCurrTable$Infant_Excluded[iSCI] = 0
            }
            sciCurrTable$ROI_SNR[iSCI] = rowMeans[iPart]
            #insert SCI param
            # Get the start and end positions for the SCI parameter
            startIndex <- gregexpr("_SCI", qtSnrFiles[iFile])[[1]][1] + nchar("_SCI")
            endIndex <- gregexpr("_PSP", qtSnrFiles[iFile])[[1]][1] - 1
            # Extract the SCI parameter
            sciParam <- substr(qtSnrFiles[iFile], startIndex, endIndex)
            # Replace the character at position 2 with a decimal point
            dotPos <- 2
            part1 <- substr(sciParam, 1, dotPos - 1)
            part2 <- substr(sciParam, dotPos, nchar(sciParam))
            sciParam <- paste0(part1, ".", part2)
            # Convert to numeric
            sciCurrTable$SCI_Param[iSCI] <- sciParam

            #increment counter
            iSCI = iSCI + 1
          }
        }
        #append current table to overall table
        sciTable <- rbind(sciTable, sciCurrTable)


        ### Both Params Data ###
        #define directory
        qtDirectory <- paste(fullPath, cohorts[iCohort], paste(ages[iAge], "mo", sep = ""), tasks[iTask], "prune", sep = "/")
        # List all files with the specified suffix in the directory
        allMatFiles <- list.files(qtDirectory, pattern = paste0("*", ".mat"), full.names = TRUE)
        # Filter files that contain SNR values
        qtSnrFiles <- allMatFiles[grepl("SNR_Prune", basename(allMatFiles))]
        # Filter files that contain chan numbers
        qtChanFiles <- allMatFiles[grepl("Chan_Prune", basename(allMatFiles))]
        # Count the number of matching files
        numFiles <- length(qtSnrFiles)

        qtCurrTable <- data.frame(
          Participant = rep(NA, numParts*numFiles),
          Cohort = rep(NA, numParts*numFiles),
          Task = rep(NA, numParts*numFiles),
          Age = rep(NA, numParts*numFiles),
          SCI_Param = rep(NA, numParts*numFiles),
          PSP_Param = rep(NA, numParts*numFiles),
          Channels = rep(NA, numParts*numFiles),
          ROI_SNR = rep(NA, numParts*numFiles),
          Infant_Excluded = rep(NA, numParts*numFiles),
          stringsAsFactors = FALSE # Prevent automatic conversion of strings to factors
        )

        # initialise count variables for rows in tables
        iQT = 1

        # start SCI files loop
        for (iFile in 1:numFiles) {
          qtSnrData <- readMat(qtSnrFiles[iFile])
          qtChanData <- readMat(qtChanFiles[iFile])

          #pre-calculate row means so that they can be inserted into the table
          # inside the participant loop
          rowMeans <- apply(qtSnrData$snrMat, 1, function(x) mean(as.numeric(x), na.rm = TRUE))

          for (iPart in 1:numParts) {

            qtCurrTable$Participant[iQT] = iPart
            qtCurrTable$Cohort[iQT] = cohorts[iCohort]
            qtCurrTable$Task[iQT] = tasks[iTask]
            qtCurrTable$Age[iQT] = as.numeric(ages[iAge])
            qtCurrTable$Channels[iQT] = qtChanData$channelsRetained[iPart]
            qtCurrTable$ROI_SNR[iQT] = rowMeans[iPart]
            if (qtCurrTable$Channels[iQT] <= floor(infInclThresh*totalChans)) {
              qtCurrTable$Infant_Excluded[iQT] = 1
            } else {
              qtCurrTable$Infant_Excluded[iQT] = 0
            }
            #insert SCI param
            # Get the start and end positions for the SCI parameter
            startIndex <- gregexpr("_SCI", qtSnrFiles[iFile])[[1]][1] + nchar("_SCI")
            endIndex <- gregexpr("_PSP", qtSnrFiles[iFile])[[1]][1] - 1
            # Extract the SCI parameter
            sciParam <- substr(qtSnrFiles[iFile], startIndex, endIndex)
            # Replace the character at position 2 with a decimal point
            dotPos <- 2
            part1 <- substr(sciParam, 1, dotPos - 1)
            part2 <- substr(sciParam, dotPos, nchar(sciParam))
            sciParam <- paste0(part1, ".", part2)

            #insert PSP param
            # Get the start and end positions for the PSP parameter
            startIndex <- gregexpr("PSP", qtSnrFiles[iFile])[[1]][1] + nchar("PSP")
            endIndex <- gregexpr(".mat", qtSnrFiles[iFile])[[1]] # chooses all instances (there are multiple)
            endIndex <- tail(endIndex, n = 1) - 1 #chooses the last one
            # Extract the PSP parameter
            pspParam <- substr(qtSnrFiles[iFile], startIndex, endIndex)
            # Replace the character at position 2 with a decimal point
            dotPos <- 2
            part1 <- substr(pspParam, 1, dotPos - 1)
            part2 <- substr(pspParam, dotPos, nchar(pspParam))
            pspParam <- paste0(part1, ".", part2)

            # Convert to numeric
            qtCurrTable$SCI_Param[iQT] <- sciParam
            qtCurrTable$PSP_Param[iQT] <- pspParam

            #increment counter
            iQT = iQT + 1
          }
        }
        #append current table to overall table
        qtTable <- rbind(qtTable, qtCurrTable)
      }
    }
  }
}

#### ---- Check for negative values and NaNs/NAs/Infs ####
# Check for negative values in cvTable
negative_cvTable <- cvTable %>%
  filter(ROI_SNR < 0)

# Check for negative values in sciTable
negative_sciTable <- sciTable %>%
  filter(ROI_SNR < 0)

# Check for negative values in qtTable
negative_qtTable <- qtTable %>%
  filter(ROI_SNR < 0)

# Print results
cat("Negative values in cvTable:\n")
print(negative_cvTable)

cat("Negative values in sciTable:\n")
print(negative_sciTable)

cat("Negative values in qtTable:\n")
print(negative_qtTable)

# Define a function to identify rows with NA, NaN, or Inf values
check_na_nan_inf <- function(table, table_name) {
  # # Check for NA, NaN, or Inf values
  # problematic_rows <- table %>%
  #   filter_all(any_vars(is.na(.) | is.nan(.) | is.infinite(.)))
  
  # Check for NA, NaN, or Inf values
  problematic_rows <- table %>%
    filter_all(any_vars(is.infinite(.)))
  
  # Display the problematic rows, if any
  if (nrow(problematic_rows) > 0) {
    cat("\nProblematic rows in", table_name, ":\n")
    print(problematic_rows)
  } else {
    cat("\nNo NA, NaN, or Inf values found in", table_name, "\n")
  }
}

# Check each table
check_na_nan_inf(cvTable, "cvTable")
check_na_nan_inf(sciTable, "sciTable")
check_na_nan_inf(qtTable, "qtTable")


# # Manual check of tables - change parameters where necessary
# # params
# ageCheck = "05"
# taskCheck = "hand"
# cohortCheck = "gm"
# #load data
# cvSnrFilePath <- paste(fullPath, cohortCheck, paste(ageCheck, "mo", sep = ""), taskCheck, "prune", paste(taskCheck, "FullPruneSNR_CV_20.mat", sep=""), sep = "/")
# cvSnrDataCheck <- readMat(cvSnrFilePath)
# #check data
# View(cvSnrDataCheck$snrMat)



##### CALCULATE CHANNEL EXCLUSION AND TOTAL INFANT EXCLUSION #########
# Find unique combinations of age, task, and cohort across the tables
uniqueCombinations <- cvTable %>%
  select(Age, Task, Cohort) %>%
  distinct()

# Initialize a list to store the results
results <- list()

# Loop through each unique combination
for (i in 1:nrow(uniqueCombinations)) {
  # Extract current combination of age, task, and cohort
  age <- uniqueCombinations$Age[i]
  task <- uniqueCombinations$Task[i]
  cohort <- uniqueCombinations$Cohort[i]
  
  # Filter the cvTable for the current combination
  cvSubset <- cvTable %>%
    filter(Age == age, Task == task, Cohort == cohort)
  
  # (i) Calculate mean channels and total infantExcluded in cvTable
  cvMeanChannels <- mean(cvSubset$Channels, na.rm = TRUE)
  cvTotalInfantExcluded <- sum(cvSubset$Infant_Excluded, na.rm = TRUE)
  
  # (ii) Calculate averages and totals by sciParam in sciTable
  sciSubset <- sciTable %>%
    filter(Age == age, Task == task, Cohort == cohort) %>%
    group_by(SCI_Param) %>%
    summarise(
      meanChannels = mean(Channels, na.rm = TRUE),
      totalInfantExcluded = sum(Infant_Excluded, na.rm = TRUE)
    )
  
  # (iii) Calculate averages and totals by sciParam and pspParam in qtTable
  qtSubset <- qtTable %>%
    filter(Age == age, Task == task, Cohort == cohort) %>%
    group_by(SCI_Param, PSP_Param) %>%
    summarise(
      meanChannels = mean(Channels, na.rm = TRUE),
      totalInfantExcluded = sum(Infant_Excluded, na.rm = TRUE)
    )
  
  # Store the results for this combination
  results[[i]] <- list(
    age = age,
    task = task,
    cohort = cohort,
    cvMeanChannels = cvMeanChannels,
    cvTotalInfantExcluded = cvTotalInfantExcluded,
    sciResults = sciSubset,
    qtResults = qtSubset
  )
}

# Display results
results



##### RANK PARAMS ##################
# BASED ON SIMILARITY TO CHANNEL AND INFANT EXCLUSION

# Loop through each result to calculate rankings
for (i in 1:length(results)) {
  # Extract necessary data for the current combination
  cvMeanChannels <- results[[i]]$cvMeanChannels
  cvTotalInfantExcluded <- results[[i]]$cvTotalInfantExcluded
  sciResults <- results[[i]]$sciResults
  qtResults <- results[[i]]$qtResults
  
  # Step 1: Rank SCI_Params based on similarity to cvTable metrics
  
  # Calculate differences in mean Channels and total Infant_Excluded for SCI_Params
  sciResults <- sciResults %>%
    mutate(
      diffChannels = abs(meanChannels - cvMeanChannels),
      diffInfantExcluded = abs(totalInfantExcluded - cvTotalInfantExcluded)
    ) %>%
    # Rank based on similarity (smallest difference) for both metrics
    mutate(
      rankChannels = rank(diffChannels),
      rankInfantExcluded = rank(diffInfantExcluded)
    ) %>%
    # Combine rankings to get an overall ranking
    mutate(
      combinedRank = rankChannels + rankInfantExcluded
    ) %>%
    # Arrange by combined rank, and if tied, by SCI_Param
    arrange(combinedRank, SCI_Param) %>%
    # Select the best SCI_Param based on combined rank
    slice(1)
  
  # Store the best SCI_Param in results
  results[[i]]$bestSCIParam <- sciResults$SCI_Param

}

# Loop through each result to calculate rankings for parameter combinations
for (i in 1:length(results)) {
  # Extract necessary data for the current combination
  cvMeanChannels <- results[[i]]$cvMeanChannels
  cvTotalInfantExcluded <- results[[i]]$cvTotalInfantExcluded
  qtResults <- results[[i]]$qtResults
  
  # First ranking step: Calculate initial differences and rankings
  qtResults <- qtResults %>%
    mutate(
      diffChannels = abs(meanChannels - cvMeanChannels),
      diffInfantExcluded = abs(totalInfantExcluded - cvTotalInfantExcluded)
    ) %>%
    mutate(
      rankChannels = rank(diffChannels),
      rankInfantExcluded = rank(diffInfantExcluded)
    ) %>%
    mutate(
      combinedRank = rankChannels + rankInfantExcluded
    ) %>%
    arrange(combinedRank, SCI_Param, PSP_Param)
  
  # Store all candidate pairs that achieved the minimum combined rank
  candidatePairs <- qtResults %>%
    filter(combinedRank == min(combinedRank))
  
  # Second ranking step: Rank among the candidate pairs
  bestPair <- candidatePairs %>%
    mutate(
      # Create new rankings among just the candidates
      finalRankChannels = rank(diffChannels),
      finalRankInfantExcluded = rank(diffInfantExcluded),
      finalCombinedRank = finalRankChannels + finalRankInfantExcluded
    ) %>%
    arrange(finalCombinedRank, SCI_Param, PSP_Param) %>%
    slice(1)
  
  # Store the single best parameter combination in results
  results[[i]]$bestParamCombo <- list(
    SCI_Param = bestPair$SCI_Param,
    PSP_Param = bestPair$PSP_Param,
    combinedRank = bestPair$finalCombinedRank,
    diffChannels = bestPair$diffChannels,
    diffInfantExcluded = bestPair$diffInfantExcluded
  )
}

# Loop through each result to calculate rankings for parameter combinations
for (i in 1:length(results)) {
  # Create a data frame from the arrays in bestParamCombo
  paramComboDF <- data.frame(
    SCI_Param = results[[i]]$bestParamCombo$SCI_Param,
    PSP_Param = results[[i]]$bestParamCombo$PSP_Param,
    diffChannels = results[[i]]$bestParamCombo$diffChannels,
    diffInfantExcluded = results[[i]]$bestParamCombo$diffInfantExcluded
  )
  
  # Rank the parameter combinations
  bestPair <- paramComboDF %>%
    mutate(
      rankChannels = rank(diffChannels),
      rankInfantExcluded = rank(diffInfantExcluded)
    ) %>%
    mutate(
      finalCombinedRank = rankChannels + rankInfantExcluded
    ) %>%
    arrange(finalCombinedRank, diffChannels) %>%  # If tied on rank, prefer lower diffChannels
    slice(1)
  
  # Store only the single best parameter combination
  results[[i]]$finalBestParams <- list(
    SCI_Param = as.character(bestPair$SCI_Param),  # Single value
    PSP_Param = as.character(bestPair$PSP_Param),  # Single value
    diffChannels = bestPair$diffChannels,          # Single value
    diffInfantExcluded = bestPair$diffInfantExcluded,  # Single value
    finalCombinedRank = bestPair$finalCombinedRank     # Single value
  )
}

##### FIT MODELS AND DEAL WITH NON_NORMAL RESIDUALS ########
##### ----- Fit models #####
# Initialize lists to store statistical test results and summary table data
# Initialize list to store models
mlmResults <- list()

# Loop through each unique combination of Age, Task, and Cohort
for (i in 1:nrow(uniqueCombinations)) {
  age <- uniqueCombinations$Age[i]
  task <- uniqueCombinations$Task[i]
  cohort <- uniqueCombinations$Cohort[i]
  
  # Get the best parameters for this combination
  bestSCIParam <- results[[i]]$bestSCIParam
  bestParamCombo <- results[[i]]$finalBestParams
  
  # Filter cvTable, sciTable, and qtTable for the current combination
  cvSubset <- cvTable %>%
    filter(Age == age, Task == task, Cohort == cohort)
  sciSubset <- sciTable %>%
    filter(Age == age, Task == task, Cohort == cohort, SCI_Param == bestSCIParam)
  qtSubset <- qtTable %>%
    filter(Age == age, Task == task, Cohort == cohort,
           SCI_Param == bestParamCombo$SCI_Param,
           PSP_Param == bestParamCombo$PSP_Param)
  
  # Extract and rename columns for merging
  cvData <- cvSubset[, c("Participant", "ROI_SNR")]
  sciData <- sciSubset[, c("Participant", "ROI_SNR")]
  qtData <- qtSubset[, c("Participant", "ROI_SNR")]
  
  names(cvData)[names(cvData) == "ROI_SNR"] <- "CV_SNR"
  names(sciData)[names(sciData) == "ROI_SNR"] <- "SCI_SNR"
  names(qtData)[names(qtData) == "ROI_SNR"] <- "QT_SNR"
  
  # Merge tables by Participants
  mergedWide <- merge(merge(cvData, sciData, by = "Participant", all = TRUE), qtData, by = "Participant", all = TRUE)
  
  # Convert to long format
  mergedLong <- pivot_longer(mergedWide, cols = c("CV_SNR", "SCI_SNR", "QT_SNR"),
                             names_to = "Source", values_to = "ROI_SNR")
  mergedLong$Source <- gsub("_SNR", "", mergedLong$Source)
  mergedLong$Source <- factor(mergedLong$Source)
  
  # Check if mergedLong has all levels of Source
  if (!all(c("CV", "SCI", "QT") %in% levels(mergedLong$Source))) {
    cat("Skipping Age:", age, "Task:", task, "Cohort:", cohort, "- missing Source levels\n")
    next
  }
  
  # Check if there's enough data in each Source level
  if (any(table(mergedLong$Source) < 2)) {
    cat("Insufficient data for Age:", age, "Task:", task, "Cohort:", cohort, "- not enough data in each Source level\n")
    next
  }
  
  # Fit mixed-effects model with lmerTest
  model <- tryCatch({
    lmer(ROI_SNR ~ Source + (1 | Participant), data = mergedLong)
  }, error = function(e) {
    cat("Error fitting model for Age:", age, "Task:", task, "Cohort:", cohort, "\n")
    return(NULL)
  })
  
  # Skip if model failed to fit
  if (is.null(model)) next
  
  # Store model
  mlmResults[[i]] <- list(age = age, task = task, cohort = cohort, mlm = model)
  
  # Plot residuals for examining model fit
  residuals <- resid(model)
  
  # QQ-plot for residuals
  qqnorm(residuals, main = paste("QQ-Plot for Residuals - Age:", age, "Task:", task, "Cohort:", cohort))
  qqline(residuals, col = "red", lwd = 2)
  
  # Histogram for residuals
  hist(residuals, breaks = 40, main = paste("Histogram of Residuals - Age:", age, "Task:", task, "Cohort:", cohort), 
       xlab = "Residuals", col = "lightblue")
}



##### ----- Investigate outliers, influential and high leverage points #####
for (i in 1:nrow(uniqueCombinations)) {
  age <- uniqueCombinations$Age[i]
  task <- uniqueCombinations$Task[i]
  cohort <- uniqueCombinations$Cohort[i]
  
  # Get the best parameters for this combination
  bestSCIParam <- results[[i]]$bestSCIParam
  bestParamCombo <- results[[i]]$finalBestParams
  
  # Filter cvTable, sciTable, and qtTable for the current combination
  cvSubset <- cvTable %>%
    filter(Age == age, Task == task, Cohort == cohort)
  
  sciSubset <- sciTable %>%
    filter(Age == age, Task == task, Cohort == cohort, SCI_Param == bestSCIParam)
  
  qtSubset <- qtTable %>%
    filter(Age == age, Task == task, Cohort == cohort,
           SCI_Param == bestParamCombo$SCI_Param,
           PSP_Param == bestParamCombo$PSP_Param)
  
  # Extract columns and rename for merging
  cvData <- cvSubset[, c("Participant", "ROI_SNR")]
  sciData <- sciSubset[, c("Participant", "ROI_SNR")]
  qtData <- qtSubset[, c("Participant", "ROI_SNR")]
  
  names(cvData)[names(cvData) == "ROI_SNR"] <- "CV_SNR"
  names(sciData)[names(sciData) == "ROI_SNR"] <- "SCI_SNR"
  names(qtData)[names(qtData) == "ROI_SNR"] <- "QT_SNR"
  
  # Merge tables by Participants
  mergedWide <- merge(merge(cvData, sciData, by = "Participant", all = TRUE), qtData, by = "Participant", all = TRUE)
  
  # Convert to long format
  mergedLong <- pivot_longer(mergedWide, cols = c("CV_SNR", "SCI_SNR", "QT_SNR"),
                             names_to = "Source", values_to = "ROI_SNR")
  mergedLong$Source <- gsub("_SNR", "", mergedLong$Source)
  mergedLong$Source <- factor(mergedLong$Source)
  
  # Check if mergedLong has all levels of Source
  if (!all(c("CV", "SCI", "QT") %in% levels(mergedLong$Source))) {
    cat("Skipping Age:", age, "Task:", task, "Cohort:", cohort, "- missing Source levels\n")
    next
  }
  
  # Check if there's enough data in each Source level
  if (any(table(mergedLong$Source) < 2)) {
    cat("Insufficient data for Age:", age, "Task:", task, "Cohort:", cohort, "- not enough data in each Source level\n")
    next
  }
  
  # Fit mixed-effects model with lmerTest
  model <- tryCatch({
    lmer(ROI_SNR ~ Source + (1 | Participant), data = mergedLong)
  }, error = function(e) {
    cat("Error fitting model for Age:", age, "Task:", task, "Cohort:", cohort, "\n")
    return(NULL)
  })
  
  # Skip if model failed to fit
  if (is.null(model)) next
  
  # Extract residuals for overall best fitting model(s)
  residuals <- resid(model)
  hatValues <- hatvalues(model)  # Leverage values
  
  # Cook's Distance (Influential Points)
  cooksD <- cooks.distance(model)
  influentialPoints <- which(cooksD > (4 / length(cooksD)))  # Influential points threshold
  influentialPointsPercentage <- length(influentialPoints) / length(cooksD) * 100
  
  # Identify outliers: Residuals greater than 3 standard deviations
  outliers <- which(abs(residuals) > 3 * sd(residuals))
  outliersPercentage <- length(outliers) / length(residuals) * 100
  
  # Identify high leverage points: Leverage greater than 2 * mean leverage
  highLeverage <- which(hatValues > 2 * mean(hatValues))
  highLeveragePercentage <- length(highLeverage) / length(hatValues) * 100
  
  # Output percentage of influential points, outliers, and high leverage points
  cat("Age:", age, "Task:", task, "Cohort:", cohort, "\n")
  cat("Influential Points (Cook's Distance) Percentage:", round(influentialPointsPercentage, 2), "%\n")
  cat("Outliers Percentage:", round(outliersPercentage, 2), "%\n")
  cat("High Leverage Points Percentage:", round(highLeveragePercentage, 2), "%\n\n")
  
}

# Nothing which is unduly concerning, so try transformation of outcome

##### ----- Transformation of outcome #####
# Try Yeo-Johnson transformation, as log, sqrt, BoxCox not suitable for 
# negative (CV) SNR values
# Loop through each unique combination of Age, Task, and Cohort
for (i in 1:nrow(uniqueCombinations)) {
  age <- uniqueCombinations$Age[i]
  task <- uniqueCombinations$Task[i]
  cohort <- uniqueCombinations$Cohort[i]
  
  # Get the best parameters for this combination
  bestSCIParam <- results[[i]]$bestSCIParam
  bestParamCombo <- results[[i]]$finalBestParams
  
  # Filter cvTable, sciTable, and qtTable for the current combination
  cvSubset <- cvTable %>%
    filter(Age == age, Task == task, Cohort == cohort)
  
  sciSubset <- sciTable %>%
    filter(Age == age, Task == task, Cohort == cohort, SCI_Param == bestSCIParam)
  
  qtSubset <- qtTable %>%
    filter(Age == age, Task == task, Cohort == cohort,
           SCI_Param == bestParamCombo$SCI_Param,
           PSP_Param == bestParamCombo$PSP_Param)
  
  # Extract columns and rename for merging
  cvData <- cvSubset[, c("Participant", "ROI_SNR")]
  sciData <- sciSubset[, c("Participant", "ROI_SNR")]
  qtData <- qtSubset[, c("Participant", "ROI_SNR")]
  
  names(cvData)[names(cvData) == "ROI_SNR"] <- "CV_SNR"
  names(sciData)[names(sciData) == "ROI_SNR"] <- "SCI_SNR"
  names(qtData)[names(qtData) == "ROI_SNR"] <- "QT_SNR"
  
  # Merge tables by Participants
  mergedWide <- merge(merge(cvData, sciData, by = "Participant", all = TRUE), qtData, by = "Participant", all = TRUE)
  
  # Convert to long format
  mergedLong <- pivot_longer(mergedWide, cols = c("CV_SNR", "SCI_SNR", "QT_SNR"),
                             names_to = "Source", values_to = "ROI_SNR")
  mergedLong$Source <- gsub("_SNR", "", mergedLong$Source)
  mergedLong$Source <- factor(mergedLong$Source)
  
  # Check if mergedLong has all levels of Source
  if (!all(c("CV", "SCI", "QT") %in% levels(mergedLong$Source))) {
    cat("Skipping Age:", age, "Task:", task, "Cohort:", cohort, "- missing Source levels\n")
    next
  }
  
  # Check if there's enough data in each Source level
  if (any(table(mergedLong$Source) < 2)) {
    cat("Insufficient data for Age:", age, "Task:", task, "Cohort:", cohort, "- not enough data in each Source level\n")
    next
  }
  
  # Apply the Yeo-Johnson transformation using car::powerTransform
  yj_transform <- powerTransform(mergedLong$ROI_SNR, family = "yjPower")
  transformedData <- yjPower(mergedLong$ROI_SNR, lambda = yj_transform$lambda)
  
  # Fit the mixed-effects model with the transformed data
  mergedLong$TransformedROI_SNR <- transformedData
  model <- tryCatch({
    lmer(TransformedROI_SNR ~ Source + (1 | Participant), data = mergedLong)
  }, error = function(e) {
    cat("Error fitting model for Age:", age, "Task:", task, "Cohort:", cohort, "\n")
    return(NULL)
  })
  
  # Skip if model failed to fit
  if (is.null(model)) next
  
  # Extract residuals for the model
  residuals <- resid(model)
  
  # QQ-plot for residuals to examine normality
  par(mfrow = c(1, 1))  # single plot layout
  qqnorm(residuals, main = "QQ-Plot for Residuals - Yeo-Johnson Transformation")
  qqline(residuals, col = "red", lwd = 2)
  
  # Histogram for residuals
  par(mfrow = c(1, 1))
  hist(residuals, breaks = 40, main = "Histogram of Residuals - Yeo-Johnson Transformation",
       xlab = "Residuals", col = "lightblue")
  
  # Store the results for further use
  mlmResults[[i]] <- list(age = age, task = task, cohort = cohort, mlm = model, transformation = "Yeo-Johnson")
}












######### STATISTICALLY TEST ROI_SNR FOR 3 PARAM CHOICES ########
# Set up parallel processing
numCores <- detectCores() - 1
registerDoParallel(numCores)

# Set number of bootstrap iterations
nBoot <- 10000
nBoot <- 20 #for testing code - comment out!
set.seed(42)

# Initialize lists to store bootstrapped summary results for each model
mlmBootResults <- list()

# Loop through each unique combination of Age, Task, and Cohort for bootstrapping
for (i in 1:nrow(uniqueCombinations)) {
  age <- uniqueCombinations$Age[i]
  task <- uniqueCombinations$Task[i]
  cohort <- uniqueCombinations$Cohort[i]
  
  # Get the best parameters for this combination
  bestSCIParam <- results[[i]]$bestSCIParam
  bestParamCombo <- results[[i]]$finalBestParams
  
  # Filter cvTable, sciTable, and qtTable for the current combination
  cvSubset <- cvTable %>%
    filter(Age == age, Task == task, Cohort == cohort)
  
  sciSubset <- sciTable %>%
    filter(Age == age, Task == task, Cohort == cohort, SCI_Param == bestSCIParam)
  
  qtSubset <- qtTable %>%
    filter(Age == age, Task == task, Cohort == cohort,
           SCI_Param == bestParamCombo$SCI_Param,
           PSP_Param == bestParamCombo$PSP_Param)
  
  # Extract columns and rename for merging
  cvData <- cvSubset[, c("Participant", "ROI_SNR")]
  sciData <- sciSubset[, c("Participant", "ROI_SNR")]
  qtData <- qtSubset[, c("Participant", "ROI_SNR")]
  
  names(cvData)[names(cvData) == "ROI_SNR"] <- "CV_SNR"
  names(sciData)[names(sciData) == "ROI_SNR"] <- "SCI_SNR"
  names(qtData)[names(qtData) == "ROI_SNR"] <- "QT_SNR"
  
  # Merge tables by Participants
  mergedWide <- merge(merge(cvData, sciData, by = "Participant", all = TRUE), qtData, by = "Participant", all = TRUE)
  
  # Convert to long format
  mergedLong <- pivot_longer(mergedWide, cols = c("CV_SNR", "SCI_SNR", "QT_SNR"),
                             names_to = "Source", values_to = "ROI_SNR")
  mergedLong$Source <- gsub("_SNR", "", mergedLong$Source)
  mergedLong$Source <- factor(mergedLong$Source)
  
  # Check if mergedLong has all levels of Source
  if (!all(c("CV", "SCI", "QT") %in% levels(mergedLong$Source))) {
    cat("Skipping Age:", age, "Task:", task, "Cohort:", cohort, "- missing Source levels\n")
    next
  }
  
  # Check if there's enough data in each Source level
  if (any(table(mergedLong$Source) < 2)) {
    cat("Insufficient data for Age:", age, "Task:", task, "Cohort:", cohort, "- not enough data in each Source level\n")
    next
  }
  
  # Fit the mixed-effects model
  model <- tryCatch({
    lmer(ROI_SNR ~ Source + (1 | Participant), data = mergedLong)
  }, error = function(e) {
    cat("Error fitting model for Age:", age, "Task:", task, "Cohort:", cohort, "\n")
    cat("Error message:", e$message, "\n")  # Print the detailed error message
    return(NULL)
  })
  
  # Skip if model failed to fit
  if (is.null(model)) next
  
  # Bootstrapping for this model
  bootResults <- foreach(j = 1:nBoot, .combine = rbind) %dopar% {
    # Resample data with replacement
    indices <- sample(nrow(mergedLong), replace = TRUE)
    bootData <- mergedLong[indices, ]
    
    # Fit the model to the bootstrapped data
    bootModel <- try(lmer(ROI_SNR ~ Source + (1 | Participant), data = bootData), silent = TRUE)
    
    # Extract fixed effects if model fit was successful
    if (!inherits(bootModel, "try-error")) {
      fixef(bootModel)
    } else {
      rep(NA, length(fixef(model)))  # Return NA if model fit failed
    }
  }
  
  # Calculate confidence intervals and standard errors for the fixed effects
  nCoef <- length(fixef(model))
  coefNames <- names(fixef(model))
  coefCIs <- t(sapply(1:nCoef, function(k) {
    quantile(bootResults[, k], probs = c(0.025, 0.975), na.rm = TRUE)
  }))
  rownames(coefCIs) <- coefNames
  coefSEs <- apply(bootResults[, 1:nCoef], 2, sd, na.rm = TRUE)
  
  # Store bootstrapped summary for this model in mlmBootResults
  mlmBootResults[[i]] <- data.frame(
    Estimate = fixef(model),
    CI_lower = coefCIs[, "2.5%"],
    CI_upper = coefCIs[, "97.5%"],
    SE = coefSEs,
    Age = age,
    Task = task,
    Cohort = cohort
  )
}

# Combine all bootstrapped summaries into a single data frame for easy viewing
combinedBootResults <- do.call(rbind, mlmBootResults)

# Display the combined bootstrapped results
print(combinedBootResults)

# Stop parallel processing
stopImplicitCluster()

# Define significance threshold and Bonferroni correction factor
pThreshold <- 0.05
bonferroniFactor <- 3 * nrow(uniqueCombinations)  # 3 contrasts per combination (SCI vs CV, QT vs CV, SCI vs QT)

# Initialize the summary table to store significance results
summaryTable <- data.frame(
  Age = integer(),
  Task = character(),
  Cohort = character(),
  sciVsCvPValue = numeric(),
  sciVsCvSignificant = logical(),
  qtVsCvPValue = numeric(),
  qtVsCvSignificant = logical(),
  sciVsQtPValue = numeric(),
  sciVsQtSignificant = logical(),
  stringsAsFactors = FALSE
)

# Loop through each model result in combinedBootResults to calculate contrasts and significance
for (i in seq_len(nrow(combinedBootResults))) {
  age <- combinedBootResults$Age[i]
  task <- combinedBootResults$Task[i]
  cohort <- combinedBootResults$Cohort[i]
  
  # Extract estimates for each Source (CV, SCI, QT)
  coefTable <- subset(combinedBootResults, Age == age & Task == task & Cohort == cohort)
  
  # Initialize p-values and significance flags
  sciVsCv <- qtVsCv <- sciVsQt <- 1
  sciVsCvSig <- qtVsCvSig <- sciVsQtSig <- FALSE
  
  # Calculate contrasts and apply Bonferroni correction
  if ("SourceSCI" %in% rownames(coefTable)) {
    sciVsCv <- coefTable["SourceSCI", "Pr(>|t|)"] * bonferroniFactor
    sciVsCvSig <- sciVsCv < pThreshold
  }
  if ("SourceQT" %in% rownames(coefTable)) {
    qtVsCv <- coefTable["SourceQT", "Pr(>|t|)"] * bonferroniFactor
    qtVsCvSig <- qtVsCv < pThreshold
  }
  if ("SourceSCI" %in% rownames(coefTable) && "SourceQT" %in% rownames(coefTable)) {
    sciVsQt <- abs(coefTable["SourceSCI", "Estimate"] - coefTable["SourceQT", "Estimate"]) * bonferroniFactor
    sciVsQtSig <- sciVsQt < pThreshold
  }
  
  # Ensure p-values are bounded at 1
  sciVsCv <- min(sciVsCv, 1)
  qtVsCv <- min(qtVsCv, 1)
  sciVsQt <- min(sciVsQt, 1)
  
  # Append the results for this Age, Task, Cohort combination to the summary table
  summaryTable <- rbind(summaryTable, data.frame(
    Age = age, 
    Task = task, 
    Cohort = cohort,
    sciVsCvPValue = sciVsCv,
    sciVsCvSignificant = sciVsCvSig,
    qtVsCvPValue = qtVsCv,
    qtVsCvSignificant = qtVsCvSig,
    sciVsQtPValue = sciVsQt,
    sciVsQtSignificant = sciVsQtSig
  ))
}

# Display the summary table of significant contrasts after Bonferroni correction
print(summaryTable)

############# PLOT ROI_SNR FOR SIG. RESULTS ############
# Check for any significant differences in the summary table
significantResults <- summaryTable %>%
  filter(sciVsCvSignificant | qtVsCvSignificant | sciVsQtSignificant)

# Generate violin plots for all three groups with significance bars if any significant differences are detected
if (nrow(significantResults) > 0) {
  cat("Significant differences detected. Generating violin plots...\n")
  
  # Loop through each significant result to create individual plots for each combination
  for (i in 1:nrow(significantResults)) {
    age <- significantResults$Age[i]
    task <- significantResults$Task[i]
    cohort <- significantResults$Cohort[i]
    
    # Filter the merged data for the current combination of Age, Task, and Cohort
    cvSubset <- cvTable %>%
      filter(Age == age, Task == task, Cohort == cohort) %>%
      mutate(Source = "CV", ROI_SNR = ROI_SNR)
    sciSubset <- sciTable %>%
      filter(Age == age, Task == task, Cohort == cohort) %>%
      mutate(Source = "SCI", ROI_SNR = ROI_SNR)
    qtSubset <- qtTable %>%
      filter(Age == age, Task == task, Cohort == cohort) %>%
      mutate(Source = "QT", ROI_SNR = ROI_SNR)
    
    # Combine data into one data frame for plotting
    plotData <- bind_rows(cvSubset, sciSubset, qtSubset) %>%
      mutate(Source = factor(Source, levels = c("CV", "SCI", "QT")))  # Set the order of the factors
    
    # Calculate adaptive yPositions based on maximum y-values in plotData
    yMaxValues <- plotData %>%
      group_by(Source) %>%
      summarize(y_max = max(ROI_SNR, na.rm = TRUE)) %>%
      pull(y_max)
    
    # Define significance levels and corresponding adaptive yPositions
    yPositions <- max(yMaxValues) + c(6, 9, 7)  # Adjust offset as needed
    
    # Generate violin plot
    p <- ggplot(plotData, aes(x = Source, y = ROI_SNR, fill = Source)) +
      geom_violin(trim = FALSE) +
      geom_boxplot(width = 0.1, outlier.shape = NA, 
                   position = position_dodge(0.9), color = "orange") +
      labs(title = paste("ROI SNR for Age:", age, "Task:", task, "Cohort:", cohort),
           x = "Source", y = "ROI SNR") +
      theme_minimal() +
      scale_fill_viridis_d()  # colourblind friendly
    
    # Set up comparisons and apply the appropriate y_position
    plotPvals <- list()
    yPositionValues <- c()
    if (significantResults$sciVsCvSignificant[i]) {
      p <- p + 
        geom_segment(aes(x = 1, xend = 1.99, y = yPositions[1], yend = yPositions[1]), 
                     size = 0.3,  # Thinner line between group A and B
                     color = "black") + 
        geom_text(aes(x = 1.5, y = yPositions[1]+0.5, label = assignStars(significantResults$sciVsCvPValue[i])), size = 5, color = "black")
    }
    if (significantResults$qtVsCvSignificant[i]) {
      p <- p + 
        geom_segment(aes(x = 1, xend = 3, y = yPositions[2], yend = yPositions[2]), 
                   size = 0.3,  # Thinner line between group A and B
                   color = "black") + 
        geom_text(aes(x = 2, y = yPositions[2]+0.5, label = assignStars(significantResults$qtVsCvPValue[i])), size = 5, color = "black")
    }
    if (significantResults$sciVsQtSignificant[i]) {
      p <- p + 
        geom_segment(aes(x = 2.01, xend = 3, y = yPositions[3], yend = yPositions[3]), 
                     size = 0.3,  # Thinner line between group A and B
                     color = "black") + 
        geom_text(aes(x = 2.5, y = yPositions[3]+0.5, label = assignStars(significantResults$sciVsQtPValue[i])), size = 5, color = "black")
    }
    
    # Print the plot for each combination
    print(p)
  }
} else {
  cat("No significant differences detected.\n")
}

############ SUMMARY STATISTICS FOR SIG. RESULTS ############
# After the main loop, outside of the analysis loop

# Initialize a new data frame to store summary statistics
summaryStats <- data.frame(Age = integer(), Cohort = character(), Task = character(),
                           Method = character(), SCI_Param = character(),
                           PSP_Param = character(), Mean = numeric(), IQR = numeric(), SD = numeric())

# Loop through each significant result to calculate summary statistics
for (i in 1:nrow(significantResults)) {
  age <- significantResults$Age[i]
  task <- significantResults$Task[i]
  cohort <- significantResults$Cohort[i]
  
  # Filter the merged data for the current combination of Age, Task, and Cohort
  cvSubset <- cvTable %>%
    filter(Age == age, Task == task, Cohort == cohort) %>%
    mutate(Method = "CV")
  
  sciSubset <- sciTable %>%
    filter(Age == age, Task == task, Cohort == cohort) %>%
    mutate(Method = "SCI")
  
  qtSubset <- qtTable %>%
    filter(Age == age, Task == task, Cohort == cohort) %>%
    mutate(Method = "QT")
  
  # Combine data into one data frame for calculating statistics
  combinedStatsData <- bind_rows(cvSubset, sciSubset, qtSubset)
  
  # Calculate summary statistics for each Method
  for (method in c("CV", "SCI", "QT")) {
    methodData <- combinedStatsData %>% filter(Method == method)
    
    if (nrow(methodData) > 0) {
      meanValue <- mean(methodData$ROI_SNR, na.rm = TRUE)
      iqrValue <- IQR(methodData$ROI_SNR, na.rm = TRUE)
      sdValue <- sd(methodData$ROI_SNR, na.rm = TRUE)
      
      # Determine SCI_Param and PSP_Param based on the Method
      if (method == "CV") {
        sciParamValue <- "-"
        pspParamValue <- "-"
      } else if (method == "SCI") {
        sciParamValue <- results[[i]]$bestSCIParam  # Access bestSCIParam
        pspParamValue <- "-"
      } else if (method == "QT") {
        sciParamValue <- results[[i]]$finalBestParams$SCI_Param  # Access finalBestParams
        pspParamValue <- results[[i]]$finalBestParams$PSP_Param;
      }
      
      # Add the calculated statistics to the summaryStats data frame
      summaryStats <- rbind(summaryStats, data.frame(Age = age, Cohort = cohort, Task = task,
                                                     Method = method, SCI_Param = sciParamValue, 
                                                     PSP_Param = pspParamValue, Mean = meanValue, 
                                                     IQR = iqrValue, SD = sdValue))
    }
  }
}

# Display the summary statistics table
print(summaryStats)


