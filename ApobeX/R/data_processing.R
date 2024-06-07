#' Retrieve Sample Names Based on Metadata Conditions
#'
#' This function retrieves a vector of sample names based on specified conditions to filter the metadata.
#'
#' @param input_data A data frame containing metadata information, where each row represents a sample.
#' @param genotype A character vector specifying the genotype condition to filter.
#' @param culture_days An integer vector specifying the culture days condition to filter.
#' @param clone A character vector specifying the clone condition to filter.
#' @return A character vector containing the sample names that meet the specified conditions.
#' @details This function filters the input data frame based on the provided genotype, culture days, 
#' and clone conditions. If no conditions are specified, all sample names are returned. If multiple 
#' conditions are provided, the function applies them cumulatively (i.e., as an AND operation).
#' @examples
#' \dontrun{
#' # Example usage:
#' # Retrieve sample names for samples with genotype "A", culture days 7, and clone "X":
#' sample_names = retrieve_sample_names(metadata_df, genotype = "A", culture_days = 7, clone = "X")
#' }
#' @export
retrieve_sample_names = function(input_data, genotype = NULL, culture_days = NULL, clone = NULL) {
  if (!is.data.frame(input_data)) {
    stop("Input must be a dataframe.")
  }

  subset_condition = rep(TRUE, nrow(input_data))

  # Filter based on genotype
  if (!is.null(genotype)) {
    subset_condition = subset_condition & (input_data$E3_genotype == genotype)
  }

  # Filter based on culture days
  if (!is.null(culture_days)) {
    subset_condition = subset_condition & (input_data$days_in_culture == culture_days)
  }

  #Filter based on clone
  if (!is.null(clone)) {
      subset_condition = subset_condition & (input_data$clone == clone)
    }

  # Apply the subset condition
  subset_data = input_data[subset_condition, , drop = FALSE]

  sample_ids_vector = unlist(as.character(subset_data$sample_ID))

  return(sample_ids_vector)
}


#' Reorder Sample Names Based on Sorting Index
#'
#' This function reorders sample names based on a sorting index. The sorting index specifies the desired order 
#' of sample names, typically based on genotype, day, and clone as serialized in the sample_mapping.csv file.
#'
#' @param sorting_index A data frame containing a column 'sample_ID' specifying the desired order of sample names.
#' @param sample_names A character vector containing the original sample names to be reordered.
#' @return A character vector containing the sample names reordered according to the sorting index.
#' @details This function reorders the sample names according to the specified sorting index. It matches the 
#' sample names in the sorting index with the original sample names and reorders them accordingly. Sample names 
#' not found in the sorting index are omitted from the result.
#' @examples
#' \dontrun{
#' # Example usage:
#' # Reorder sample names based on a sorting index:
#' sorted_names = order_sample_names(sorting_index_df, sample_names_vector)
#' }
#' @export
order_sample_names = function(sorting_index, sample_names){
  
  ordered_sample_names = sample_names[match(sorting_index$sample_ID, sample_names)]
  ordered_sample_names = ordered_sample_names[!is.na(ordered_sample_names)]
  
  return(ordered_sample_names)
}


#' Reorder Matrix Columns According to Sample Name Vector
#'
#' This function reorders the columns of a matrix according to a given sample name vector.
#'
#' @param matrix_data The matrix whose columns are to be reordered.
#' @param sample_vector A character vector specifying the desired order of sample names.
#' @return A matrix with columns reordered according to the provided sample name vector.
#' @details This function takes a matrix and a sample name vector as input. It matches the sample names 
#' in the sample vector with the column names of the matrix and reorders the columns accordingly. If a 
#' sample name in the vector does not match any column name in the matrix, it is omitted from the result.
#' @examples
#' \dontrun{
#' # Example usage:
#' # Reorder matrix columns according to a sample name vector:
#' reordered_matrix = reorder_matrix_columns(matrix_data, sample_vector)
#' }
#' @export
reorder_matrix_columns = function(matrix_data, sample_vector) {
  # Get the column indices based on the order of the sample vector
  column_indices = match(sample_vector, colnames(matrix_data))
  
  # Reorder the columns of the matrix
  reordered_matrix = matrix_data[, column_indices]
  
  return(reordered_matrix)
}

#' Subset Mutational Matrix by Sample Names
#'
#' This function subsets a mutational matrix based on a vector of sample names.
#'
#' @param original_matrix The original mutational matrix to be subsetted.
#' @param sample_names A character vector specifying the sample names to be included in the subset.
#' @return A subset of the original mutational matrix containing only the specified sample names.
#' @details This function takes an original mutational matrix and a vector of sample names as input. 
#' It subsets the original matrix to include only the columns corresponding to the specified sample names. 
#' If a sample name in the vector does not match any column name in the matrix, it is omitted from the result.
#' @examples
#' \dontrun{
#' # Example usage:
#' # Subset a mutational matrix by a vector of sample names:
#' subsetted_matrix = subset_matrix_by_samples(original_matrix, sample_names)
#' }
#' @export
subset_matrix_by_samples = function(original_matrix, sample_names) {
  selected_columns = colnames(original_matrix) %in% sample_names
  return(original_matrix[, selected_columns, drop = FALSE])
}

#' Rescale Mutational Matrix to Have the Same Total Number of Mutations
#'
#' This function rescales a mutational matrix to have the same total number of mutations across all samples.
#'
#' @param mat The original mutational matrix to be rescaled.
#' @return A rescaled mutational matrix where each sample has the same total number of mutations.
#' @details This function takes an original mutational matrix as input. It calculates the total number of 
#' mutations for each sample and rescales the mutation counts such that each sample has the same total number 
#' of mutations. This is achieved by adjusting the mutation counts in each column proportionally. If a sample 
#' already has the minimum total number of mutations, it remains unchanged.
#' @examples
#' \dontrun{
#' # Example usage:
#' # Rescale a mutational matrix to have the same total number of mutations:
#' rescaled_matrix = scale_matrix(original_matrix)
#' }
#' @export
scale_matrix = function(mat) {
  col_sums = colSums(mat)
  min_col_sum = min(col_sums)
  
  scaled_matrix = mat
  
  for (col_index in seq_along(col_sums)) {
    if (col_sums[col_index] != min_col_sum) {
      scaling_factor = min_col_sum / col_sums[col_index]
      scaled_matrix[, col_index] = round(mat[, col_index] * scaling_factor)
    }
  }
  
  return(scaled_matrix)
}

#' Check if Column Contains All Zeros
#'
#' This function checks if a specified column in a matrix contains all zeros.
#'
#' @param matrix The matrix in which to check the column.
#' @param column_index The index of the column to be checked for all zeros.
#' @return TRUE if the specified column contains all zeros, FALSE otherwise.
#' @details This function takes a matrix and a column index as input. It checks if the column specified 
#' by the index contains all zeros. If all elements in the column are zero, the function returns TRUE; otherwise, 
#' it returns FALSE.
#' @examples
#' \dontrun{
#' # Example usage:
#' # Check if a column contains all zeros:
#' is_all_zeros = check_all_zeros(matrix_data, 2)
#' }
#' @export
check_all_zeros = function(matrix, column_index) {
  all(matrix[, column_index] == 0)
}

#' Eliminate All-Zero Columns from Matrix
#'
#' This function eliminates all-zero columns from a matrix.
#'
#' @param matrix The original matrix from which to eliminate all-zero columns.
#' @return A list containing:
#' \item{non_zero_matrix}{The matrix with all-zero columns removed.}
#' \item{name_of_zero_columns}{A character vector containing the names of the removed all-zero columns.}
#' @details This function takes a matrix as input and identifies columns containing all zeros. It then removes 
#' these columns from the matrix and returns the modified matrix along with the names of the removed columns. 
#' If there are no all-zero columns, the original matrix is returned along with an empty character vector.
#' @examples
#' \dontrun{
#' # Example usage:
#' # Eliminate all-zero columns from a matrix:
#' result = eliminate_all_zero_columns(matrix_data)
#' non_zero_matrix = result[[1]]
#' names_of_zero_columns = result[[2]]
#' }
#' @export
eliminate_all_zero_columns = function(matrix){

  all_zeros_columns = c()

  # Loop over all columns and check for all zeros
  for (col_index in seq_len(ncol(matrix))) {
    if (check_all_zeros(matrix, col_index)) {
      all_zeros_columns = c(all_zeros_columns, col_index)
    }
  }
  
  original_column_names = colnames(matrix)
  
  if(is.null(all_zeros_columns)){
    non_zero_matrix = matrix
    name_of_zero_columns = c()
  }
  else{
    non_zero_matrix = matrix[,-all_zeros_columns]
    name_of_zero_columns = original_column_names[all_zeros_columns]
  }
  
 return(list(non_zero_matrix, name_of_zero_columns))
}

#' Check if Sum of Column is Less Than Threshold
#'
#' This function checks if the sum of a specified column in a matrix is less than a given threshold.
#'
#' @param matrix The matrix in which to check the column sum.
#' @param column_index The index of the column to be checked for its sum.
#' @param threshold The threshold value to compare against the sum of the column.
#' @return TRUE if the sum of the specified column is less than the threshold, FALSE otherwise.
#' @details This function takes a matrix, a column index, and a threshold value as input. It calculates 
#' the sum of the column specified by the index and checks if it is less than the threshold. If the sum 
#' is less than the threshold, the function returns TRUE; otherwise, it returns FALSE.
#' @examples
#' \dontrun{
#' # Example usage:
#' # Check if the sum of a column is less than a threshold:
#' is_sum_less_than_threshold = check_sum_less_than_TH(matrix_data, 2, 1000)
#' }
#' @export
check_sum_less_than_TH = function(matrix, column_index, threshold) {
  sum(matrix[, column_index]) < threshold
}


#' Eliminate Columns with Sum Less Than Threshold from Matrix
#'
#' This function eliminates columns from a matrix where the sum of values is less than a specified threshold.
#'
#' @param matrix The original matrix from which to eliminate columns.
#' @param threshold The threshold value for the sum of column values.
#' @return A list containing:
#' \item{filtered_matrix}{The matrix with columns removed if their sum is less than the threshold.}
#' \item{filtered_columns}{A character vector containing the names of removed columns.}
#' @details This function takes a matrix and a threshold value as input. It checks the sum of values in each column 
#' of the matrix and removes columns where the sum is less than the specified threshold. The function returns the modified 
#' matrix with columns removed along with a character vector containing the names of the removed columns.
#' @examples
#' \dontrun{
#' # Example usage:
#' # Eliminate columns with sum less than a threshold from a matrix:
#' result = eliminate_low_mut_samples(matrix_data, 1000)
#' filtered_matrix = result[[1]]
#' removed_columns = result[[2]]
#' }
#' @export
eliminate_low_mut_samples = function(matrix, threshold) {
  selected_columns = c()

  # Loop over all columns and check the sum
  for (col_index in seq_len(ncol(matrix))) {
    if (check_sum_less_than_TH(matrix, col_index, threshold)) {
      selected_columns = c(selected_columns, col_index)
    }
  }

  original_column_names = colnames(matrix)

  if (length(selected_columns) == 0) {
    filtered_matrix = matrix
    filtered_columns = character(0)
  } else {
    filtered_matrix = matrix[, -selected_columns, drop = FALSE]
    filtered_columns = original_column_names[selected_columns]
  }

  return(list(filtered_matrix, filtered_columns))
}

####--------------------------------------------------------------------------------####

#' Aggregate Refitting Results (Beta)
#'
#' This function aggregates refitting results by calculating row means for each set of 3 columns.
#'
#' @param data The data to be aggregated, typically refitting results.
#' @param column_names A character vector specifying the column names for the result matrix.
#' @return A matrix containing aggregated refitting results with row means for each set of 3 columns.
#' @details This function takes input data, typically refitting results, and aggregates them by calculating 
#' row means for each set of 3 columns. It checks if the number of columns is divisible by 3, creates an empty 
#' matrix to store the aggregated results, and then loops through each set of 3 columns to calculate row means. 
#' The function returns the matrix containing aggregated refitting results with row means for each set of 3 columns.
#' @examples
#' \dontrun{
#' # Example usage:
#' # Aggregate refitting results:
#' aggregated_results = aggregate_refitting_results_beta(data, column_names)
#' }
#' @export
aggregate_refitting_results = function(data, column_names) {
  num_rows = nrow(data)
  num_columns = ncol(data)
  row_names = rownames(data)
  
  # Check if the number of columns is divisible by 3
  if (num_columns %% 3 != 0) {
    stop("Number of columns is not divisible by 3")
  }
  
  # Create an empty matrix to store the results
  result = matrix(NA, nrow = num_rows, ncol = num_columns / 3)
  
  # Loop through each set of 3 columns
  for (i in 1:(num_columns / 3)) {
    # Calculate row means for each set of 3 columns
    result[, i] = rowMeans(data[, ((i - 1) * 3 + 1):(i * 3), drop = FALSE])
  }
  
  # Reassign row names
  rownames(result) = row_names
  
  # Assign column names
  colnames(result) = column_names
  # Return the result
  return(result)
}

