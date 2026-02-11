#' Example Summary Statistics and LD Matrices for KAMA
#'
#' Summary statistics and LD matrices for 233 variants across four independent populations
#' (AFR, EAS, EUR, SAS).
#'
#' @description
#' This dataset contains summary statistics (Z-scores, Betas, SEs) and reference LD correlation
#' matrices. It serves as input for generating Knockoff statistics (e.g., using GhostKnockoff)
#' and performing subsequent meta-analysis.
#'
#' @format A named list of 4 elements (AFR, EAS, EUR, SAS). Each element contains:
#' \describe{
#'   \item{sum_stats}{A data frame with 233 rows and 9 columns:
#'     \describe{
#'       \item{chr}{Integer. The chromosome number.}
#'       \item{SNP}{Character. The variant identifier (e.g., rsID).}
#'       \item{dis}{Integer. Genetic distance (cM).}
#'       \item{pos}{Integer. Physical base pair position.}
#'       \item{A1, A2}{Character. The alternative and reference alleles.}
#'       \item{Beta}{Numeric. The effect size estimate.}
#'       \item{SE}{Numeric. The standard error of the effect size.}
#'       \item{Zscore}{Numeric. The Z-statistic (Beta / SE).}
#'     }
#'   }
#'   \item{cor_mat}{A 233 x 233 symmetric numeric matrix representing the LD correlation structure.}
#' }
#'
#' @usage data(pops.sum_res)
#'
#' @examples
#' \dontrun{
#' # 1. Load data
#' data(pops.sum_res)
#'
#' # Note: This example requires the 'GhostKnockoff' package to generate statistics.
#' library(GhostKnockoff)
#'
#' # Parameters
#' M <- 5 # number of knockoffs
#' n.study <- 2000 # sample size per population
#' kappa_list <- list()
#' tau_list <- list()
#'
#' # 2. Process each population
#' for (pop in names(pops.sum_res)) {
#'   sum_stats <- pops.sum_res[[pop]]$sum_stats
#'   cor_mat   <- pops.sum_res[[pop]]$cor_mat
#'
#'   # Generate Knockoffs (Set seed for reproducibility)
#'   set.seed(99)
#'   fit.prelim <- GhostKnockoff.prelim(cor_mat, M = M, method = 'asdp', max.size = 500)
#'
#'   # Fit Knockoff statistics
#'   GK.stat <- GhostKnockoff.fit(
#'     Zscore_0 = as.matrix(sum_stats$Zscore),
#'     n.study = n.study,
#'     fit.prelim = fit.prelim,
#'     gamma = 1
#'   )
#'
#'   # Calculate KAMA statistics (Kappa and Tau)
#'   stats <- MK.statistic(GK.stat$T_0, GK.stat$T_k)
#'   kappa_list[[pop]] <- stats[, "kappa"]
#'   tau_list[[pop]]   <- stats[, "tau"]
#' }
#'
#' # 3. Prepare input matrices
#' kappa_set <- do.call(cbind, kappa_list)
#' tau_set   <- do.call(cbind, tau_list)
#'
#' # Add prefixes to distinguish columns
#' colnames(kappa_set) <- paste0("kappa_", names(pops.sum_res))
#' colnames(tau_set)   <- paste0("tau_", names(pops.sum_res))
#'
#' # Extract genetic information from the first study
#' genetic_info <- pops.sum_res[[1]]$sum_stats[, c("chr", "SNP", "pos", "A1", "A2")]
#'
#' # 4. Run Meta-analysis and Localization
#' meta_res <- KAMA_meta_analysis(kappa_set, tau_set, M = M, genetic_info = genetic_info)
#'
#' loc_res  <- KAMA_localization(kappa_set, tau_set, M = M, q = 0.1,
#'                               genetic_info = genetic_info,
#'                               study_label = names(pops.sum_res))
#'
#' print(loc_res)
#' }
"pops.sum_res"

#' Example Association Statistics for KAMA
#'
#' Association statistics (P-values) from four independent studies, including results
#' for both original variables and Knockoff replicates.
#'
#' @description
#' This dataset contains P-values regarding genomic regions across four studies.
#' It serves as input for calculating Knockoff statistics (Kappa and Tau) and performing
#' subsequent meta-analysis and localization.
#'
#' @format A named list of 4 data frames. Each data frame corresponds to a study and contains:
#' \describe{
#'   \item{chr}{Integer. The chromosome number.}
#'   \item{start, end}{Integer. The start and end position of the genomic region.}
#'   \item{actual_start, actual_end}{Integer. The position of the first and last variant in the region.}
#'   \item{p_value}{Numeric. The P-value from the original association test.}
#'   \item{p_value_kc1 ... p_value_kc5}{Numeric. P-values from 5 Knockoff replicates.}
#' }
#'
#' @usage data(study.assoc_stats)
#'
#' @examples
#' # 1. Load data
#' data(study.assoc_stats)
#'
#' # Parameters
#' M <- 5 # number of knockoffs
#' kappa_list <- list()
#' tau_list <- list()
#'
#' # 2. Process each study
#' for (study in names(study.assoc_stats)) {
#'   df <- study.assoc_stats[[study]]
#'
#'   # Extract P-values
#'   P_0 <- df$p_value
#'   P_k <- df[, paste0("p_value_kc", 1:M)]
#'
#'   # Transform P-values to importance scores (larger is better)
#'   # Using -log10(P) transformation
#'   T_0 <- -log10(P_0)
#'   T_k <- -log10(as.matrix(P_k))
#'
#'   # Calculate KAMA statistics
#'   stats <- MK.statistic(T_0, T_k)
#'   kappa_list[[study]] <- stats[, "kappa"]
#'   tau_list[[study]]   <- stats[, "tau"]
#' }
#'
#' # 3. Prepare input matrices
#' kappa_set <- do.call(cbind, kappa_list)
#' tau_set   <- do.call(cbind, tau_list)
#'
#' # Add prefixes to distinguish columns
#' colnames(kappa_set) <- paste0("kappa_", names(study.assoc_stats))
#' colnames(tau_set)   <- paste0("tau_", names(study.assoc_stats))
#'
#' # Extract genetic coordinates from the first study
#' genetic_info <- study.assoc_stats[[1]][, c("chr", "start", "end", "actual_start", "actual_end")]
#'
#' # 4. Run Meta-analysis and Localization
#' meta_res <- KAMA_meta_analysis(kappa_set, tau_set, M = M, genetic_info = genetic_info)
#'
#' loc_res  <- KAMA_localization(kappa_set, tau_set, M = M, q = 0.1,
#'                               genetic_info = genetic_info,
#'                               study_label = names(study.assoc_stats))
#'
#' print(loc_res)
"study.assoc_stats"
