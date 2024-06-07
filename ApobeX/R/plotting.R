# Plotting functions

####--------------------------------------------------------------------------------####
# time ordered plot of signature refitting after bootstrap correction

#' Master Function to Calculate, Sort, and Replot Results
#'
#' This function serves as a master routine to calculate, sort, and replot results based on mutation data.
#'
#' @param mut_data The mutation data used for analysis.
#' @param sample_names A character vector specifying the sample names.
#' @param prior_knowledge_subset A subset of prior knowledge to be used in the analysis.
#' @details This function performs the following steps:
#' 1. Calculates the best subset fitting based on mutation data and prior knowledge subset.
#' 2. Reorders the resulting dataframe based on the specified sample names.
#' 3. Plots the reordered contribution of each subset.
#' @examples
#' \dontrun{
#' # Example usage:
#' # Calculate, sort, and replot results based on mutation data:
#' best_subset_plot(mut_data, sample_names, prior_knowledge_subset)
#' }
#' @export
best_subset_plot = function(mut_data, sample_names, prior_knowledge_subset){
  
  result_df = best_subset_fitting(mut_data, prior_knowledge_subset)
  result_df_reordered = reorder_dataframe(result_df, sample_names)
  plot_reordered_contribution(result_df_reordered)
}


#' Best Subset Fitting
#'
#' This function fits mutation data to the best subset based on prior knowledge.
#'
#' @param mut_data The mutation data to be analyzed.
#' @param prior_knowledge_subset A subset of prior knowledge to be used in the fitting.
#' @return A tibble containing the contribution data of the best subset fitting.
#' @details This function fits mutation data to the best subset based on the provided prior knowledge subset. 
#' It uses the method of fitting to signatures strictly with a maximum delta of 0.0001 and the backward method. 
#' The function then generates a dummy plot to pull out the total contribution data, which is used to reorder bars 
#' in the plot. The contribution data is returned as a tibble.
#' @examples
#' \dontrun{
#' # Example usage:
#' # Fit mutation data to the best subset based on prior knowledge:
#' contribution_data = best_subset_fitting(mut_data, prior_knowledge_subset)
#' }
#' @export
best_subset_fitting = function(mut_data, prior_knowledge_subset){
  
  best_subset_refit = fit_to_signatures_strict(mut_data, prior_knowledge_subset, max_delta = 0.0001, method = "backwards")
  result = best_subset_refit$fit_res$contribution
  
  # Generate dummy plot
  p_tmp = plot_contribution(result,
                            coord_flip = FALSE,
                            mode = "absolute"
  )
  # Pull out total contribution data to reorder bars
  bar_data = p_tmp$data
  
  return(tibble(bar_data))
}

#' Reorder Dataframe Based on Sample Names
#'
#' This function reorders a dataframe based on specified sample names.
#'
#' @param df The dataframe to be reordered.
#' @param sample_names A character vector specifying the desired order of sample names.
#' @return A reordered dataframe based on sample names.
#' @details This function takes a dataframe and a vector of sample names as input. It converts the 'Sample' 
#' column of the dataframe to a factor with the specified order of sample names and removes duplicate rows. 
#' The dataframe is then reordered based on the 'Signature' column and the new order of 'Sample'. The reordered 
#' dataframe is returned as the result.
#' @examples
#' \dontrun{
#' # Example usage:
#' # Reorder a dataframe based on sample names:
#' reordered_df = reorder_dataframe(df, sample_names)
#' }
#' @export
reorder_dataframe = function(df, sample_names) {
  
  # Convert Sample and Signature columns to factors with specified order
  df$Sample = factor(df$Sample, levels = sample_names)
  df$Signature = factor(df$Signature, levels = unique(df$Signature))
  
  # Remove duplicate rows
  df = distinct(df)
  
  # Reorder the dataframe based on Signature and the new order of Sample
  result_df = df %>% arrange(Signature, Sample)
  
  return(result_df)
}


#' Plot Reordered Contribution
#'
#' This function plots the reordered contribution data.
#'
#' @param reordered_contribution A dataframe containing reordered contribution data.
#' @return A ggplot object representing the contribution by signatures.
#' @details This function takes a dataframe containing reordered contribution data as input. It generates 
#' a bar plot using ggplot2 library, where the x-axis represents the samples, the y-axis represents the 
#' contribution, and the bars are filled by signature types. The colors of signature types are manually specified. 
#' The function returns the ggplot object representing the contribution by signatures.
#' @examples
#' \dontrun{
#' # Example usage:
#' # Plot reordered contribution data:
#' plot_reordered_contribution(reordered_contribution)
#' }
#' @export
plot_reordered_contribution = function(reordered_contribution){
  ggplot(reordered_contribution, aes(x = Sample, y = Contribution, fill = Signature)) +
    geom_bar(stat = "identity") +
    labs(title = "Contribution by Signatures",
         x = "Sample",
         y = "Contribution") +
    scale_fill_manual(values = c("SBS1" = "#e67474", "SBS2" = "#d1af49", "SBS5" = "#9abf32",
                                  "SBS6" = "#3acf65", "SBS13" = "#3ccbb3", "SBS15" = "#44b9e6", "SBS18" = "#c970ec",
                                  "SBS40" = "#eb70df" )) +
    theme_minimal()
}

####--------------------------------------------------------------------------------####

#' Make Zoomed-in Heatmap of APOBEC-Related Signatures
#'
#' This function creates a zoomed-in heatmap of APOBEC-related signatures.
#'
#' @param refit_model The refitted model containing contribution data.
#' @param sample_selection A vector specifying the samples to be included in the heatmap.
#' @param cell_line A character string specifying the cell line for the plot title.
#' @return A ggplot object representing the zoomed-in heatmap of APOBEC-related signatures.
#' @details This function takes a refitted model containing contribution data, a vector specifying 
#' the samples to be included in the heatmap, and a character string specifying the cell line for 
#' the plot title. It extracts contribution data for APOBEC-related signatures, creates a melted 
#' dataframe, and generates a heatmap using ggplot2 library. The colors of the heatmap are based on 
#' a gradient color palette. The function returns the ggplot object representing the zoomed-in heatmap 
#' of APOBEC-related signatures.
#' @examples
#' \dontrun{
#' # Example usage:
#' # Create a zoomed-in heatmap of APOBEC-related signatures:
#' plot_apobec_heatmap(refit_model, sample_selection, cell_line)
#' }
#' @export
plot_apobec_heatmap = function(refit_model, sample_selection, cell_line){
  contribution_df = as.data.frame(refit_model$contribution)
  contribution_df_apo2 = contribution_df[c(2,6),sample_selection]
  contribution_df_apo2 = t(contribution_df_apo2)
  melted_data_apo2 = melt(contribution_df_apo2, varnames = c("Sample", "Signature"), value.name = "Contribution")

  num_colors = 100
  color_palette = colorRampPalette(brewer.pal(9, "Blues"))(num_colors)

  # Create the heatmap
  ggplot(melted_data_apo2, aes(x = Sample, y = Signature, fill = Contribution)) +
    geom_tile(color = "white")+
    scale_fill_gradientn(colors = color_palette)+
    labs(title = paste(cell_line))
}

####--------------------------------------------------------------------------------####

#' Visualize Specific Sequence Context by Genotype
#'
#' This function reshapes and plots data to visualize a specific sequence context by genotype.
#'
#' @param input_data The input data to be reshaped and plotted.
#' @param cell_line A character string specifying the cell line for the plot title.
#' @param context_name A character string specifying the name of the specific sequence context.
#' @return A ggplot object representing the visualization of the specific sequence context by genotype.
#' @details This function takes input data, melts it to a long format using the 'melt' function from the 
#' 'reshape2' package, and then creates a bar plot using ggplot2 library. Each bar represents the percentage 
#' of mutation contexts for a specific genotype. The x-axis represents the samples, the y-axis represents 
#' the percentage, and bars are filled by variable. The function returns the ggplot object representing 
#' the visualization of the specific sequence context by genotype.
#' @examples
#' \dontrun{
#' # Example usage:
#' # Visualize SBS2 specific sequence context by genotype:
#' plot_specific_context(input_data, cell_line, context_name)
#' }
#' @export
plot_specific_context = function(input_data, cell_line, context_name){
  input_data_melted = melt(input_data, id.vars = "Sample")
  input_data_melted$value = as.numeric(input_data_melted$value)

  ggplot(input_data_melted, aes(x = Sample, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(x = "Sample", y = "Percentage", title =  paste("Percentage mutation contexts for genotype", context_name, cell_line))
}