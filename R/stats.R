#' Compute Multiple Knockoff Statistics (Kappa and Tau)
#'
#' Calculates the feature importance statistics for multiple knockoffs.
#' Following the methodology in Gimenez and Zou, this function computes \eqn{\kappa}
#' (the index of the feature with the largest importance score) and \eqn{\tau}
#' (the difference between the largest score and the median of the remaining scores).
#'
#' @param T_0 A numeric vector or a one-column matrix representing the importance scores 
#'   of the **original** features (variables).
#' @param T_k A numeric matrix representing the importance scores of the **knockoff** features.
#'   It must have the same number of rows as \code{T_0} and \code{M} columns 
#'   (where \code{M} is the number of knockoff replicates).
#' @param method A character string specifying the aggregation method for \eqn{\tau}. 
#'   Currently supports 'median'. Default is 'median'.
#'
#' @return A matrix with two columns:
#' \describe{
#'   \item{kappa}{The index of the original (denoted as 0) or knockoff feature (1 to M) 
#'     that has the largest importance score. Formally: \eqn{\kappa = \text{argmax}_{0 \le m \le M} T^{(m)}}.}
#'   \item{tau}{The difference between the largest importance score and the median of the 
#'     remaining importance scores.}
#' }
#' 
#' @references 
#' Gimenez, J. R., & Zou, J. (2019). Improving the stability of the knockoff procedure: 
#' Multiple knockoffs and pantheon.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate dummy importance scores for 100 variables
#' # Original scores (T_0)
#' T0 <- abs(rnorm(100)) 
#' # Knockoff scores (T_k) with M=5 knockoffs
#' Tk <- matrix(abs(rnorm(100 * 5)), ncol = 5)
#' 
#' # Calculate statistics
#' stats <- MK.statistic(T0, Tk)
#' head(stats)
#' }
MK.statistic <- function(T_0, T_k, method = 'median') {
  T_0 <- as.matrix(T_0)
  T_k <- as.matrix(T_k)
  
  T.temp <- cbind(T_0, T_k)
  T.temp[is.na(T.temp)] <- 0
  
  kappa <- apply(T.temp, 1, which.max) - 1
  
  if (method == 'median') {
    Get.OtherMedian <- function(x) {
      rem <- x[-which.max(x)]
      if(length(rem) == 0) return(0) # 防止报错
      median(rem)
    }
    tau <- apply(T.temp, 1, max) - apply(T.temp, 1, Get.OtherMedian)
  }
  
  return(cbind(kappa, tau))
}
