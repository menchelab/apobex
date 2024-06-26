% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrapped_background_correction.R
\name{bootstrap}
\alias{bootstrap}
\title{Bootstrap Matrices from Hash}
\usage{
bootstrap(input_hash_of_matrices, seed)
}
\arguments{
\item{input_hash_of_matrices}{A hash object where each key contains a matrix to be bootstrapped.}

\item{seed}{An integer value to set the seed for repeatable results.}
}
\value{
A hash object where each key contains a list of bootstrapped matrices.
}
\description{
This function pulls out individual sample matrices from a hash, performs a single bootstrap
iteration, and stores the result as a list of bootstrapped matrices in a new hash. Each matrix
in the hash is one iteration of the bootstrapping. The seed
can be set to obtain repeatable results.
}
\details{
This function iterates over each matrix in the input hash, performs a single bootstrap
iteration using the \code{mutSignatures:::bootstrapCancerGenomes} function, and stores the resulting
bootstrapped matrices in a new hash. Each bootstrapped matrix is one bootstrapped resampling of the
input data.
}
\examples{
\dontrun{
library(hash)
library(mutSignatures)
mat1 = matrix(1:6, nrow = 2, dimnames = list(NULL, c("A", "B", "C")))
mat2 = matrix(7:12, nrow = 2, dimnames = list(NULL, c("D", "E", "F")))
input_hash = hash()
input_hash[["sample1"]] = mat1
input_hash[["sample2"]] = mat2
seed = 123
bootstrapped_hash = bootstrap(input_hash, seed)
print(bootstrapped_hash)
}

}
