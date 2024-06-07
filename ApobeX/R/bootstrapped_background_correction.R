# Functions for bootstrapping routine

#' Master Bootstrap Routine
#'
#' This function serves as a master routine for bootstrapping a given input matrix.
#'
#' @param input_matrix The input matrix to be bootstrapped.
#' @return A hash containing bootstrapped matrices obtained from the input matrix.
#' @details This function takes an input matrix and performs the following steps:
#' 1. Extracts the input matrix into a hash.
#' 2. Prepares bootstrap iterations in the hash.
#' 3. Performs bootstrapping using the prepared iterations.
#' The function then returns a hash containing bootstrapped matrices derived from the input matrix.
#' @examples
#' \dontrun{
#' # Example usage:
#' # Perform bootstrapping routine on an input matrix:
#' bootstrapped_results = bootstrap_routine(input_matrix)
#' }
#' @export
bootstrap_routine = function(input_matrix){

    # Extract input matrix into hash
    input_matrix_hash = matrix_to_hash(input_matrix)
    # Prepare bootstrap iterations in hash
    input_matrix_hash_10k = gather_replicates_into_hash(input_matrix_hash, 10000)
    # Perform bootstrapping
    input_matrix_bootstrapped = bootstrap(input_matrix_hash_10k, seed=1234)

    return(input_matrix_bootstrapped)
}

###############################################################################################
############################### BOOSTRAPPING FUNCTIONS ########################################
###############################################################################################

#' Convert Matrix Columns to Hash
#'
#' This function takes a matrix as input and stores each column of the matrix 
#' in a hash, using the column names as keys.
#'
#' @param input_matrix A matrix with named columns.
#' @return A hash object where each key is a column name and each value is a vector 
#' corresponding to the matrix column.
#' @details This function is useful for separating out all the clones in a matrix 
#' and storing each with its own key.
#' @examples
#' \dontrun{
#' library(hash)
#' mat = matrix(1:9, nrow = 3, dimnames = list(NULL, c("A", "B", "C")))
#' result_hash = matrix_to_hash(mat)
#' print(result_hash)
#' }
#' 
matrix_to_hash = function(input_matrix)
{
  output_hash = hash()
  hash_keys = colnames(input_matrix)
  for (i in seq(1:length(hash_keys)))
  {
    output_hash[[hash_keys[i]]] = input_matrix[,i]
  }
  
  return(output_hash)
}



#' Replicate Columns in a Matrix
#'
#' This helper function replicates the columns of an input matrix a specified number of times.
#' It is useful for preparing matrices for bootstrapping.
#'
#' @param input_matrix A matrix whose columns will be replicated.
#' @param num_replica An integer specifying the number of times to replicate the columns of the input matrix.
#' @return A matrix with replicated columns.
#' @details The function takes an input matrix and replicates its columns, creating a larger matrix
#' with repeated columns. This is useful for certain statistical techniques such as bootstrapping.
#' @examples
#' \dontrun{
#' mat = matrix(1:6, nrow = 2, dimnames = list(NULL, c("A", "B", "C")))
#' boot_matrix = prepare_bootstrap_iterations(mat, 3)
#' print(boot_matrix)
#' }
#' 
prepare_bootstrap_iterations = function(input_matrix, num_replica)
{
  boot_matrix = input_matrix
  for (i in seq(2:num_replica))
  {
    boot_matrix = cbind(boot_matrix, input_matrix)
  }
  
  return(boot_matrix)
}



#' Gather Replicate Matrices into a Hash for Bootstrapping
#'
#' This function creates a hash of replicate matrices, ready for bootstrapping. 
#' The number of bootstrapped samples is equal to the specified number of replications.
#'
#' @param mutation_hash A hash object where each key contains a matrix to be replicated.
#' @param num_replica An integer specifying the number of replications for each matrix in the hash.
#' @return A hash object where each key contains a replicated matrix.
#' @details This function iterates over each matrix in the input hash, replicates its columns 
#' the specified number of times, and stores the resulting matrix in a new hash. This is useful 
#' for preparing data for bootstrapping.
#' @examples
#' \dontrun{
#' library(hash)
#' mat1 = matrix(1:6, nrow = 2, dimnames = list(NULL, c("A", "B", "C")))
#' mat2 = matrix(7:12, nrow = 2, dimnames = list(NULL, c("D", "E", "F")))
#' mutation_hash = hash()
#' mutation_hash[["clone1"]] = mat1
#' mutation_hash[["clone2"]] = mat2
#' replicated_hash = gather_replicates_into_hash(mutation_hash, 3)
#' print(replicated_hash)
#' }
#' 
gather_replicates_into_hash = function(mutation_hash, num_replica)
{
  output_hash = hash()
  for(i in ls(mutation_hash))
  {
    output_hash[[i]] = prepare_bootstrap_iterations(mutation_hash[[i]], num_replica)
  }
  
  return(output_hash)
}


#' Bootstrap Matrices from Hash
#'
#' This function pulls out individual sample matrices from a hash, performs a single bootstrap 
#' iteration, and stores the result as a list of bootstrapped matrices in a new hash. Each matrix
#' in the hash is one iteration of the bootstrapping. The seed 
#' can be set to obtain repeatable results.
#'
#' @param input_hash_of_matrices A hash object where each key contains a matrix to be bootstrapped.
#' @param seed An integer value to set the seed for repeatable results.
#' @return A hash object where each key contains a list of bootstrapped matrices.
#' @details This function iterates over each matrix in the input hash, performs a single bootstrap 
#' iteration using the `mutSignatures:::bootstrapCancerGenomes` function, and stores the resulting 
#' bootstrapped matrices in a new hash. Each bootstrapped matrix is one bootstrapped resampling of the
#' input data.
#' @examples
#' \dontrun{
#' library(hash)
#' library(mutSignatures)
#' mat1 = matrix(1:6, nrow = 2, dimnames = list(NULL, c("A", "B", "C")))
#' mat2 = matrix(7:12, nrow = 2, dimnames = list(NULL, c("D", "E", "F")))
#' input_hash = hash()
#' input_hash[["sample1"]] = mat1
#' input_hash[["sample2"]] = mat2
#' seed = 123
#' bootstrapped_hash = bootstrap(input_hash, seed)
#' print(bootstrapped_hash)
#' }
#' 
bootstrap = function(input_hash_of_matrices, seed)
{
  bootstrap_results = hash()
  for (i in ls(input_hash_of_matrices))
  {
    bootstrap_results[[i]] = mutSignatures:::bootstrapCancerGenomes(input_hash_of_matrices[[i]], seed)
  }
  
  return(bootstrap_results)
}


###############################################################################################
############################### AGGREGATION FUNCTIONS #########################################
###############################################################################################

#' Compute Median of Bootstrapped Samples
#'
#' This function pulls out the median (a measure of central tendency) of bootstrapped results 
#' stored in a hash.
#'
#' @param hash_of_bootstrap_results A hash object where each key contains a list of bootstrapped matrices.
#' @return A hash object where each key contains the median values of the bootstrapped matrices.
#' @details This function iterates over each list of bootstrapped matrices in the input hash, 
#' computes the median for each row across the matrices, and stores the resulting medians in a new hash.
#' @examples
#' \dontrun{
#' library(hash)
#' bootstrapped_results = hash()
#' bootstrapped_results[["sample1"]] = matrix(runif(20), nrow=5)
#' bootstrapped_results[["sample2"]] = matrix(runif(20), nrow=5)
#' median_results = median_bootstrapped_samples(bootstrapped_results)
#' print(median_results)
#' }
#' @export
median_bootstrapped_samples = function(hash_of_bootstrap_results) {
  median_clones_hash = hash()
  
  for (i in ls(hash_of_bootstrap_results)) {
    median_clones_hash[[i]] = apply(hash_of_bootstrap_results[[i]], 1, median)
  }
  
  return(median_clones_hash)
}



#' Compute Average of Bootstrapped Samples
#'
#' This function computes the row-wise average of bootstrapped matrices stored in a hash 
#' and stores the results in a new hash.
#'
#' @param hash_of_bootstrap_results A hash object where each key contains a list of bootstrapped matrices.
#' @return A hash object where each key contains the average values of the bootstrapped matrices.
#' @details This function iterates over each list of bootstrapped matrices in the input hash, 
#' computes the row-wise average for each matrix, and stores the resulting averages in a new hash.
#' @examples
#' \dontrun{
#' library(hash)
#' bootstrapped_results = hash()
#' bootstrapped_results[["sample1"]] = matrix(runif(20), nrow = 5)
#' bootstrapped_results[["sample2"]] = matrix(runif(20), nrow = 5)
#' average_results = average_bootstrapped_samples(bootstrapped_results)
#' print(average_results)
#' }
#' @export
average_bootstrapped_samples = function(hash_of_bootstrap_results)
{
  averaged_samples_hash = hash()
  
  for (i in ls(hash_of_bootstrap_results))
  {
    averaged_samples_hash[[i]] = rowMeans(hash_of_bootstrap_results[[i]])
  }
  
  return(averaged_samples_hash)
}


#' Convert Hash of Averaged Results to Matrix
#'
#' This function converts a hash of averaged results into a matrix.
#'
#' @param averaged_samples_hash A hash object where each key contains a single column of averaged results (row-wise average).
#' @return A matrix with averaged results, where each column corresponds to a key in the input hash.
#' @details This function iterates over each item in the input hash, constructs a data frame from the averaged results, 
#' and then converts this data frame into a matrix. The row names of the resulting matrix are set to match 
#' the row names of the original mutation matrix (`mut_mat`).
#' @examples
#' \dontrun{
#' library(hash)
#' library(S4Vectors)
#' averaged_results = hash()
#' averaged_results[["sample1"]] = runif(5)
#' averaged_results[["sample2"]] = runif(5)
#' mut_mat = matrix(1:10, nrow = 5)  # Example mutation matrix to provide row names
#' average_matrix = hash_to_average_matrix(averaged_results)
#' print(average_matrix)
#' }
#' @export
hash_to_average_matrix = function(averaged_samples_hash)
{
  tmp_df = DataFrame()
  for(i in ls(averaged_samples_hash))
  {
    tmp_df[i] = averaged_samples_hash[[i]] 
  }
  
  averaged_boot_matrix = as.matrix(tmp_df)
  row.names(averaged_boot_matrix) = row.names(mut_mat)
  
  return(averaged_boot_matrix)
}



#' Convert Hash of Bootstrapped Results to Matrix
#'
#' This function converts a hash of bootstrapped results into a matrix for error bars.
#'
#' @param input_hash A hash object where each key contains a list of bootstrapped matrices.
#' @return A matrix containing all bootstrapped samples concatenated column-wise.
#' @details This function iterates over each list of bootstrapped matrices in the input hash, 
#' combines them into a single data frame, and converts this data frame into a matrix.
#' @examples
#' \dontrun{
#' library(hash)
#' bootstrapped_results = hash()
#' bootstrapped_results[["sample1"]] = matrix(runif(20), nrow = 5)
#' bootstrapped_results[["sample2"]] = matrix(runif(20), nrow = 5)
#' bootstrap_matrix = hash_to_matrix(bootstrapped_results)
#' print(bootstrap_matrix)
#' }
#' @export
hash_to_matrix = function(input_hash)
{
  tmp_df = DataFrame()
  for(i in ls(input_hash))
  {
    for (j in ls(input_hash[[i]])){
      tmp_df = cbind(tmp_df, input_hash[[,j]])
    }
  }
  
  bootstrap_matrix = as.matrix(tmp_df)
  
  return(bootstrap_matrix)
}