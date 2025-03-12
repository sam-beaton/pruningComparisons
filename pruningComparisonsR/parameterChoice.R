  ######################
  library(dplyr)
  library(tidyr)
  
  # Choosing best parameters
  
  #load the data, just in case
  data <- read.csv("[path-to...]/pruningComparisonsMatlab/stats/[project]/overall/pruneMLMInputTable.csv")
  
  head(data)
  
  # Define threshold for Infant_Excluded
  thresholdInfant <- 0.01  # Change this to adjust the threshold (e.g., 0.05 for 5%)
  
  # Create an empty data frame to store results
  grid_results <- expand.grid(SCI = unique(data$SCI), PSP = unique(data$PSP), Age = unique(data$Age))
  
  # Function to calculate ROI_SNR and Channels_Retained for each combination
  for (i in 1:nrow(grid_results)) {
    subset_data <- data %>%
      filter(SCI == grid_results$SCI[i], PSP == grid_results$PSP[i], Age == grid_results$Age[i])
    
    grid_results$ROI_SNR[i] <- mean(subset_data$ROI_SNR, na.rm = TRUE)
    grid_results$Channels_Retained[i] <- mean(subset_data$Channels_Retained, na.rm = TRUE)
  }
  
  grid_results <- as.data.frame(grid_results)
  
  # Fit a linear model for each age group and calculate perpendicular distances
  grid_results <- grid_results %>%
    group_by(Age) %>%
    do({
      model <- lm(ROI_SNR ~ Channels_Retained, data = .)
      
      # Create a data frame to store model predictions and distances
      resultsParamSelect <- data.frame(
        Channels_Retained = .$Channels_Retained,
        ROI_SNR = .$ROI_SNR,
        predicted = predict(model, newdata = .),
        residuals = residuals(model),
        Age = .$Age,
        SCI = .$SCI,  # Include SCI and PSP for the best point later
        PSP = .$PSP
      )
      
      # Calculate perpendicular distance
      resultsParamSelect$perpendicular_distance <- abs(resultsParamSelect$residuals) / sqrt(1 + (coef(model)[2])^2)
      resultsParamSelect
    }) %>%
    ungroup()
  
  # Find the points with the greatest perpendicular distance AND higher values than predicted
  max_distance_points <- grid_results %>%
    filter(perpendicular_distance > 0) %>%
    group_by(Age) %>%
    filter(ROI_SNR > predicted & Channels_Retained > predicted) %>%  # Only points above/right of line
    slice(which.max(perpendicular_distance)) %>%  # Maximize perpendicular distance
    ungroup()
  
  # Print the results for the selected best points along with SCI and PSP
  max_distance_points <- max_distance_points %>%
    dplyr::select(Age, Channels_Retained, ROI_SNR, SCI, PSP, perpendicular_distance)
  
  print(max_distance_points)
  
  # Plotting the decision grid with fitted lines and highlighted best points
  ggplot() +
    geom_point(data = grid_results, aes(x = Channels_Retained, y = ROI_SNR, color = as.factor(Age), shape = as.factor(Age)), size = 3) +  # Data points
    geom_smooth(data = grid_results, aes(x = Channels_Retained, y = ROI_SNR, group = Age), method = "lm", se = FALSE, color = "black", size = 0.5) +  # All lines black
    geom_point(data = max_distance_points, aes(x = Channels_Retained, y = ROI_SNR), color = "red", size = 4, shape = 8) +  # Highlight best points
    labs(title = "Average SNR v Channel Retention across all ages", x = "Channels Retained", y = "ROI SNR", color = "Age", shape = "Age") +
    theme_minimal() +
    scale_shape_manual(values = c(16, 17, 18, 19, 15, 8)) +  # Different shapes for each age
    scale_color_manual(values = RColorBrewer::brewer.pal(6, "Set1"))  # Color palette for age groups
  
  # Plot each age separately
  # Unique ages
  unique_ages <- unique(grid_results$Age)
  
  # Define a custom shape scale for PSP
  custom_shapes <- c(16, 17, 18, 19, 15, 8, 13, 14, 10, 11, 7, 12, 9, 4, 3, 1, 2, 5, 6, 20)  # 18 shapes, adjust as needed
  
  # Loop through each unique age group
  for (age_group in unique_ages) {
    # Filter the data for the current age group
    age_data <- grid_results %>% filter(Age == age_group)
    best_point <- max_distance_points %>% filter(Age == age_group)
    
    ### Plot grouped by SCI ###
    
    # Create the plot, grouped by SCI using a gradient color scale
    p_sci <- ggplot() +
      # Data points grouped by SCI with a color gradient
      geom_point(data = age_data, aes(x = Channels_Retained, y = ROI_SNR, color = SCI, shape = as.factor(SCI)), size = 3) +
      
      # Fitted line for the current age group
      geom_smooth(data = age_data, aes(x = Channels_Retained, y = ROI_SNR), method = "lm", se = FALSE, color = "black", size = 0.5) +
      
      # Highlight the best point with a larger size, keeping the same shape and color as its group
      geom_point(data = best_point, aes(x = Channels_Retained, y = ROI_SNR, color = SCI, shape = as.factor(SCI)), size = 5, stroke = 1.5) +
      
      # Add a black outline around the best point to make it stand out
      geom_point(data = best_point, aes(x = Channels_Retained, y = ROI_SNR), size = 7, shape = 21, fill = NA, color = "black", stroke = 1.5) +
      
      labs(title = paste("Average SNR v Channel Retention for Age", age_group, "Grouped by SCI Threshold"), 
           x = "Channels Retained", y = "Average SNR", color = "SCI", shape = "SCI") +
      
      theme_minimal() +
      
      # Use a viridis color gradient for SCI
      scale_color_viridis_c(option = "D") +  # Change option as desired for different gradients
      scale_shape_manual(values = custom_shapes)  # Different shapes for each SCI
    
    # Display the SCI plot
    print(p_sci)
    
    ### Plot grouped by PSP ###
    
    # Create the plot, grouped by PSP using a gradient color scale
    p_psp <- ggplot() +
      # Data points grouped by PSP with a color gradient
      geom_point(data = age_data, aes(x = Channels_Retained, y = ROI_SNR, color = PSP, shape = as.factor(PSP)), size = 3) +
      
      # Fitted line for the current age group
      geom_smooth(data = age_data, aes(x = Channels_Retained, y = ROI_SNR), method = "lm", se = FALSE, color = "black", size = 0.5) +
      
      # Highlight the best point with a larger size, keeping the same shape and color as its group
      geom_point(data = best_point, aes(x = Channels_Retained, y = ROI_SNR, color = PSP, shape = as.factor(PSP)), size = 5, stroke = 1.5) +
      
      # Add a black outline around the best point to make it stand out
      geom_point(data = best_point, aes(x = Channels_Retained, y = ROI_SNR), size = 7, shape = 21, fill = NA, color = "black", stroke = 1.5) +
      
      labs(title = paste("Average SNR v Channel Retention for Age", age_group, "Grouped by PSP Threshold"), 
           x = "Channels Retained", y = "Average SNR", color = "PSP", shape = "PSP") +
      
      theme_minimal() +
      
      scale_shape_manual(values = custom_shapes) +  # Use custom shapes for PSP
      scale_color_viridis_c(option = "D")  # Use viridis gradient for PSP
    
    # Display the PSP plot
    print(p_psp)
  }
  
  ##### Heatmaps
  # Create heatmaps for each age, ranking parameter combinations while accounting 
  # for variable relationships across ages with SCI and PSP parameters
  unique_ages <- unique(grid_results$Age)
  
  for (age in unique_ages) {
    heatmap_data <- grid_results %>%
      filter(Age == age) %>%
      dplyr::select(SCI, PSP, perpendicular_distance) %>%
      mutate(rank = rank(-perpendicular_distance))  # Rank based on perpendicular_distance (larger values get lower ranks)
    
    # Ensure SCI and PSP are factors for proper heatmap layout
    heatmap_data$SCI <- factor(heatmap_data$SCI)
    heatmap_data$PSP <- factor(heatmap_data$PSP)
    
    # Create the plot with custom color gradient
    heatmap_plot <- ggplot(heatmap_data, aes(x = SCI, y = PSP)) +
      geom_tile(aes(fill = rank), color = "white") +  # Fill by rank instead of perpendicular distance
      scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                           midpoint = median(heatmap_data$rank, na.rm = TRUE), 
                           name = "Rank") +  # Custom gradient for rank
      labs(title = paste("Parameter Combination Rankings for Data Quality v Retention (", age, "months)"),
           x = "SCI",
           y = "PSP") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_text(aes(label = rank), color = "black", size = 3)  # Add ranking as text labels
    
    # Display the plot
    print(heatmap_plot)
  }
  
  
  
  
  
  
  
  ### Look at % of participants excluded by parameter
  # Flexible proportion
  proportion <- 0.6  # Change this to adjust the proportion (e.g., 0.7, 0.65)
  
  # Calculate threshold based on the total number of channels (68) and the selected proportion
  thresholdChannel <- floor(68 * proportion)
  
  ### Look at % of participants excluded by parameter
  # Function: percentage of participants with Channels_Retained <= threshold for each value of SCI & PSP
  plot_percentage_by_sci_psp <- function(df, age_group, thresholdChannel) {
    # Filter by age group
    df_age <- df %>% filter(Age == age_group)
    
    # Calculate percentage of participants with Channels_Retained <= threshold by SCI
    sci_group <- df_age %>%
      group_by(SCI) %>%
      summarise(Percentage = mean(Channels_Retained <= thresholdChannel) * 100)
    
    # Calculate percentage of participants with Channels_Retained <= threshold by PSP
    psp_group <- df_age %>%
      group_by(PSP) %>%
      summarise(Percentage = mean(Channels_Retained <= thresholdChannel) * 100)
    
    # Plot for SCI
    p1 <- ggplot(sci_group, aes(x = SCI, y = Percentage)) +
      geom_line(color = "blue") +
      geom_point(color = "blue") +
      labs(title = paste("Percentage of Infants with Channels Retained <= ", thresholdChannel, " by SCI (Age:", age_group, "months)"),
           x = "SCI", y = "Percentage of Participants Excluded") +
      theme_minimal()
    
    # Plot for PSP
    p2 <- ggplot(psp_group, aes(x = PSP, y = Percentage)) +
      geom_line(color = "green") +
      geom_point(color = "green") +
      labs(title = paste("Percentage of Infants with Channels Retained <= ", thresholdChannel, " by PSP (Age:", age_group, "months)"),
           x = "PSP", y = "Percentage of Participants Excluded") +
      theme_minimal()
    
    # Print both plots
    gridExtra::grid.arrange(p1, p2, ncol = 2)
  }
  
  # Heatmap of SCI vs PSP with percentage <= thresholdChannel
  plot_heatmap <- function(df, age_group, thresholdChannel) {
    # Check if df is a data frame and Age exists
    if(!is.data.frame(df)) stop("Input df is not a data frame")
    if(!"Age" %in% colnames(df)) stop("Column 'Age' not found in the data frame")
    
    # Filter by age group
    df_age <- df %>% filter(Age == age_group)
    
    # Create a summary table with the percentage of Channels_Retained <= thresholdChannel for each combination of SCI and PSP
    summary_table <- df_age %>%
      group_by(SCI, PSP) %>%
      summarise(Percentage = mean(Channels_Retained <= thresholdChannel) * 100, .groups = 'drop')
    
    # Pivot to long format to make it suitable for heatmap plotting
    summary_table_long <- summary_table %>%
      pivot_longer(cols = Percentage, names_to = "Metric", values_to = "Value")
    
    # Plot heatmap with PSP on the y-axis and SCI on the x-axis
    ggplot(summary_table_long, aes(x = SCI, y = PSP, fill = Value)) +
      geom_tile(color = "white") +
      scale_fill_gradient(low = "blue", high = "red", na.value = "white", name = "Percentage") +
      labs(title = paste("Heatmap of Channels Retained <= ", thresholdChannel, " by SCI and PSP (Age:", age_group, "months)"),
           x = "SCI", y = "PSP") +
      theme_minimal()
  }
  
  
  # Plot
  for (age in unique(data$Age)) {
    plot_percentage_by_sci_psp(data, age, thresholdChannel)
    print(plot_heatmap(data, age, thresholdChannel))
  }
  
  # Plot all ages
  plot_percentage_by_sci_psp_all_ages <- function(df, thresholdChannel) {
    # Calculate percentage of participants with Channels_Retained <= thresholdChannel by SCI for all ages
    sci_group <- df %>%
      group_by(Age, SCI) %>%
      summarise(Percentage = mean(Channels_Retained <= thresholdChannel) * 100)
    
    # Calculate percentage of participants with Channels_Retained <= thresholdChannel by PSP for all ages
    psp_group <- df %>%
      group_by(Age, PSP) %>%
      summarise(Percentage = mean(Channels_Retained <= thresholdChannel) * 100)
    
    # Plot for SCI across all ages with different shapes and colors
    p1 <- ggplot(sci_group, aes(x = SCI, y = Percentage, shape = factor(Age), color = factor(Age), group = Age)) +
      geom_line() +
      geom_point(size = 3) +
      scale_shape_manual(values = 1:length(unique(df$Age))) +  # Different shapes for each age
      scale_color_manual(values = rainbow(length(unique(df$Age)))) +  # Different colors for each age
      labs(title = paste("Percentage of Infants excluded by SCI parameter value (all ages)"),
           x = "SCI", y = "Percentage of Participants Excluded", shape = "Age Group", color = "Age Group") +
      theme_minimal()
    
    # Plot for PSP across all ages with different shapes and colors
    p2 <- ggplot(psp_group, aes(x = PSP, y = Percentage, shape = factor(Age), color = factor(Age), group = Age)) +
      geom_line() +
      geom_point(size = 3) +
      scale_shape_manual(values = 1:length(unique(df$Age))) +  # Different shapes for each age
      scale_color_manual(values = rainbow(length(unique(df$Age)))) +  # Different colors for each age
      labs(title = paste("Percentage of Infants excluded by PSP parameter value (all ages)"),
           x = "PSP", y = "Percentage of Participants Excluded", shape = "Age Group", color = "Age Group") +
      theme_minimal()
    
    # Print both plots
    gridExtra::grid.arrange(p1, p2, ncol = 2)
  }
  
  # Use function to plot
  plot_percentage_by_sci_psp_all_ages(data, thresholdChannel)
  
  
  
  
  
  
  
  
  
  
  
  
  
  #### Best parameter choice restricted to those combinations which do not exclude
  #### too many infants
  
  # Generate combinations of SCI and PSP
  SCI_vals <- c(0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05)
  PSP_vals <- c(0.1, 0.095, 0.09, 0.085, 0.08, 0.075, 0.07, 0.065, 0.06, 0.055, 0.05, 0.045, 0.04, 0.035, 0.03, 0.025, 0.02, 0.015, 0.01, 0.05)
  ages <- unique(data$Age)
  
  # Create an empty data frame to store results
  grid_results <- expand.grid(SCI = SCI_vals, PSP = PSP_vals, Age = ages)
  
  # Function to calculate ROI_SNR and Channels_Retained for each combination
  for (i in 1:nrow(grid_results)) {
    subset_data <- data %>%
      filter(SCI == grid_results$SCI[i], PSP == grid_results$PSP[i], Age == grid_results$Age[i])
    
    grid_results$ROI_SNR[i] <- mean(subset_data$ROI_SNR, na.rm = TRUE)
    grid_results$Channels_Retained[i] <- mean(subset_data$Channels_Retained, na.rm = TRUE)
    grid_results$Infant_Excluded_Rate[i] <- mean(subset_data$Infant_Excluded, na.rm = TRUE)  # Calculate exclusion rate
  }
  
  # Fit a linear model for each age group to calculate predictions for all combinations
  overall_results <- grid_results %>%
    group_by(Age) %>%
    do({
      model <- lm(ROI_SNR ~ Channels_Retained, data = .)
      
      # Create a data frame to store model predictions and distances
      resultsParamSelect <- data.frame(
        Channels_Retained = .$Channels_Retained,
        ROI_SNR = .$ROI_SNR,
        predicted = predict(model, newdata = .),
        residuals = residuals(model),
        Age = .$Age,
        SCI = .$SCI,  
        PSP = .$PSP,
        Infant_Excluded_Rate = .$Infant_Excluded_Rate
      )
      
      # Calculate perpendicular distance
      resultsParamSelect$perpendicular_distance <- abs(resultsParamSelect$residuals) / sqrt(1 + (coef(model)[2])^2)
      resultsParamSelect
    }) %>%
    ungroup()
  
  # Filter grid_results based on the Infant_Excluded threshold for best parameter selection
  filtered_grid_results <- grid_results %>%
    filter(Infant_Excluded_Rate <= thresholdInfant)  # Keep only those with exclusion rate below threshold
  
  # Fit a linear model for filtered grid_results to get predictions and distances
  filtered_results <- filtered_grid_results %>%
    group_by(Age) %>%
    do({
      model <- lm(ROI_SNR ~ Channels_Retained, data = .)
      
      # Create a data frame to store predictions and distances
      resultsParamSelect <- data.frame(
        Channels_Retained = .$Channels_Retained,
        ROI_SNR = .$ROI_SNR,
        predicted = predict(model, newdata = .),
        residuals = residuals(model),
        Age = .$Age,
        SCI = .$SCI,  
        PSP = .$PSP,
        Infant_Excluded_Rate = .$Infant_Excluded_Rate
      )
      
      # Calculate perpendicular distance
      resultsParamSelect$perpendicular_distance <- abs(resultsParamSelect$residuals) / sqrt(1 + (coef(model)[2])^2)
      resultsParamSelect
    }) %>%
    ungroup()
  
  # Find the points with the greatest perpendicular distance AND higher values than predicted
  max_distance_points <- filtered_results %>%
    group_by(Age) %>%
    filter(ROI_SNR > predicted & Channels_Retained > predicted) %>%  # Only points above/right of line
    slice(which.max(perpendicular_distance)) %>%  # Maximize perpendicular distance
    ungroup()
  
  # If no combinations meet the threshold for a particular age group, select the one with the smallest exclusion rate
  for (age in unique(ages)) {
    if (nrow(max_distance_points[max_distance_points$Age == age, ]) == 0) {
      best_choice <- filtered_grid_results %>%
        filter(Age == age) %>%
        arrange(Infant_Excluded_Rate) %>%
        slice(1)
      
      max_distance_points <- rbind(max_distance_points, best_choice)
    }
  }
  
  # Print the updated results for the selected best points along with SCI, PSP, and Infant_Excluded Rate
  max_distance_points <- max_distance_points %>%
    dplyr::select(Age, Channels_Retained, ROI_SNR, SCI, PSP, perpendicular_distance, Infant_Excluded_Rate) %>%
    arrange(Age)  # Sort by Age in ascending order
  
  print(max_distance_points)
  
  # # Plotting the decision grid with fitted lines and highlighted best points
  # ggplot() +
  #   geom_point(data = filtered_grid_results, aes(x = Channels_Retained, y = ROI_SNR, color = as.factor(Age), shape = as.factor(Age)), size = 3, alpha = 0.6) +  # Data points
  #   geom_smooth(data = filtered_results, aes(x = Channels_Retained, y = ROI_SNR, group = Age), method = "lm", se = FALSE, color = "black", size = 0.5) +  # All lines black
  #   geom_point(data = max_distance_points, aes(x = Channels_Retained, y = ROI_SNR), color = "red", size = 4, shape = 8) +  # Highlight best points
  #   labs(
  #     title = paste("Average SNR v Channel Retention with Exclusion Rate <", threshold * 100, "%"),
  #     x = "Channels Retained", 
  #     y = "ROI SNR", 
  #     color = "Age", 
  #     shape = "Age"
  #   )
  #   theme_minimal() +
  #   scale_shape_manual(values = c(16, 17, 18, 19, 15, 8)) +  # Different shapes for each age
  #   scale_color_manual(values = RColorBrewer::brewer.pal(6, "Set1"))  # Color palette for age groups
  
  # Plot each age separately
  # Unique ages
  unique_ages <- unique(filtered_results$Age)
  
  # Loop through each unique age group
  for (age_group in unique_ages) {
    # Filter the data for the current age group
    age_data <- filtered_results %>% filter(Age == age_group)
    best_point <- max_distance_points %>% filter(Age == age_group)
    
    ### Plot grouped by SCI ###
    
    # Create the plot, grouped by SCI using a gradient color scale
    p_sci <- ggplot() +
      geom_point(data = age_data, aes(x = Channels_Retained, y = ROI_SNR, color = SCI, shape = as.factor(SCI)), size = 3) +
      geom_smooth(data = age_data, aes(x = Channels_Retained, y = ROI_SNR), method = "lm", se = FALSE, color = "black", size = 0.5) +
      geom_point(data = best_point, aes(x = Channels_Retained, y = ROI_SNR, color = SCI, shape = as.factor(SCI)), size = 5, stroke = 1.5) +
      geom_point(data = best_point, aes(x = Channels_Retained, y = ROI_SNR), size = 7, shape = 21, fill = NA, color = "black", stroke = 1.5) +
      labs(title = paste("Average SNR v Channel Retention for Age", age_group, "Grouped by SCI"), 
           x = "Channels Retained", y = "Average SNR", color = "SCI", shape = "SCI") +
      theme_minimal() +
      scale_color_viridis_c(option = "D") +  # Change option as desired for different gradients
      scale_shape_manual(values = custom_shapes)  # Different shapes for each SCI
    
    # Display the SCI plot
    print(p_sci)
    
    ### Plot grouped by PSP ###
    
    # Create the plot, grouped by PSP using a gradient color scale
    p_psp <- ggplot() +
      geom_point(data = age_data, aes(x = Channels_Retained, y = ROI_SNR, color = PSP, shape = as.factor(PSP)), size = 3) +
      geom_smooth(data = age_data, aes(x = Channels_Retained, y = ROI_SNR), method = "lm", se = FALSE, color = "black", size = 0.5) +
      geom_point(data = best_point, aes(x = Channels_Retained, y = ROI_SNR, color = PSP, shape = as.factor(PSP)), size = 5, stroke = 1.5) +
      geom_point(data = best_point, aes(x = Channels_Retained, y = ROI_SNR), size = 7, shape = 21, fill = NA, color = "black", stroke = 1.5) +
      labs(title = paste("Average SNR v Channel Retention for Age", age_group, "Grouped by PSP"), 
           x = "Channels Retained", y = "Average SNR", color = "PSP", shape = "PSP") +
      theme_minimal() +
      scale_shape_manual(values = custom_shapes) +  # Use custom shapes for PSP
      scale_color_viridis_c(option = "D")  # Use viridis gradient for PSP
    
    # Display the PSP plot
    print(p_psp)
  }
  
  
  ### Heatmaps
  
  # Generate combinations of SCI and PSP
  SCI_vals <- c(0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05)
  PSP_vals <- c(0.1, 0.095, 0.09, 0.085, 0.08, 0.075, 0.07, 0.065, 0.06, 0.055, 0.05, 0.045, 0.04, 0.035, 0.03, 0.025, 0.02, 0.015, 0.01, 0.05)
  ages <- unique(data$Age)
  
  # Create all possible combinations
  all_combinations <- expand.grid(SCI = SCI_vals, PSP = PSP_vals, Age = ages)
  
  # Create an empty data frame to store results
  all_results <- all_combinations %>%
    rowwise() %>%
    mutate(
      ROI_SNR = mean(data$ROI_SNR[data$SCI == SCI & data$PSP == PSP & data$Age == Age], na.rm = TRUE),
      Channels_Retained = mean(data$Channels_Retained[data$SCI == SCI & data$PSP == PSP & data$Age == Age], na.rm = TRUE),
      Infant_Excluded_Rate = mean(data$Infant_Excluded[data$SCI == SCI & data$PSP == PSP & data$Age == Age], na.rm = TRUE)
    ) %>%
    ungroup()
  
  # Fit a linear model for all combinations to calculate predictions and distances
  overall_results <- all_results %>%
    group_by(Age) %>%
    do({
      model <- lm(ROI_SNR ~ Channels_Retained, data = .)
      
      resultsParamSelect <- data.frame(
        Channels_Retained = .$Channels_Retained,
        ROI_SNR = .$ROI_SNR,
        predicted = predict(model, newdata = .),
        residuals = residuals(model),
        Age = .$Age,
        SCI = .$SCI,
        PSP = .$PSP,
        Infant_Excluded_Rate = .$Infant_Excluded_Rate
      )
      
      # Calculate perpendicular distance
      resultsParamSelect$perpendicular_distance <- abs(resultsParamSelect$residuals) / sqrt(1 + (coef(model)[2])^2)
      resultsParamSelect
    }) %>%
    ungroup()
  
  # Create heatmaps for each age, ranking parameter combinations
  unique_ages <- unique(overall_results$Age)
  
  for (age in unique_ages) {
    heatmap_data <- overall_results %>%
      filter(Age == age) %>%
      dplyr::select(SCI, PSP, perpendicular_distance, Infant_Excluded_Rate) %>%
      mutate(rank = rank(-perpendicular_distance))  # Rank based on perpendicular_distance
    
    # Ensure SCI and PSP are factors for proper heatmap layout
    heatmap_data$SCI <- factor(heatmap_data$SCI)
    heatmap_data$PSP <- factor(heatmap_data$PSP)
    
    # Identify the best rank (lowest rank) among acceptable infant exclusion rates
    best_rank <- min(heatmap_data$rank[heatmap_data$Infant_Excluded_Rate <= thresholdInfant], na.rm = TRUE)
    
    # Create the plot with custom color gradient
    heatmap_plot <- ggplot(heatmap_data, aes(x = SCI, y = PSP)) +
      geom_tile(aes(fill = rank), color = "white", 
                alpha = ifelse(heatmap_data$Infant_Excluded_Rate > thresholdInfant, 0.3, 1)) +  #change opacity based on exclusion rate
      scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                           midpoint = median(heatmap_data$rank, na.rm = TRUE), 
                           name = "Rank") +  # Custom gradient for rank
      labs(title = paste("Parameter Choice for Data Quality vs Retention (", age, "months)"),
           x = "SCI",
           y = "PSP") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_text(aes(label = rank), color = "black", 
                size = ifelse(heatmap_data$Infant_Excluded_Rate > thresholdInfant, 2.5, 3)) +  # change text size based on exclusion rate
      # Highlight the cell with the best rank among acceptable rates
      geom_tile(data = heatmap_data %>% filter(rank == best_rank & Infant_Excluded_Rate <= thresholdInfant),
                color = "black", size = 1.5, fill = NA)  # Draw border around the best rank cell
    
    # Display the plot
    print(heatmap_plot)
  }
