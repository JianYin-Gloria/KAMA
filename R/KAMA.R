#' KAMA Meta-Analysis
#'
#' @param kappa_set A numeric matrix of kappa statistics. Rows correspond to variables,
#'   and columns correspond to studies.
#'   **Note:** The order of variables (rows) and studies (columns) must strictly match that of \code{tau_set}.
#' @param tau_set A numeric matrix of tau statistics. Rows correspond to variables,
#'   and columns correspond to studies.
#'   **Note:** The order of variables (rows) and studies (columns) must strictly match that of \code{kappa_set}.
#' @param M Integer. The number of knockoffs used. Must be >= max(kappa_set).
#' @param genetic_info Optional data frame or matrix containing additional information (e.g., SNP IDs, positions) to be appended to the results.
#' @param seed Integer. Random seed for reproducibility (default 111).
#'
#' @return A data frame containing:
#' \itemize{
#'   \item \code{genetic_info} (if provided)
#'   \item \code{kappa_set} columns
#'   \item \code{tau_set} columns
#'   \item \code{qvalue}: KAMA's q-value for meta-analysis FDR control.
#' }
#' @export
KAMA_meta_analysis <- function(kappa_set, tau_set, M, genetic_info = NULL, seed = 111){

  # Input validation: check if tau_set and kappa_set are provided
  if (missing(tau_set) || missing(kappa_set)) {
    stop("tau_set and kappa_set are required arguments.")
  }

  # Convert to matrix if necessary
  tau_set <- as.matrix(tau_set)
  kappa_set <- as.matrix(kappa_set)

  # Input validation: check dimensions of tau_set and kappa_set
  if (!all(dim(tau_set) == dim(kappa_set))) {
    stop("tau_set and kappa_set must have the same dimensions. ",
         "tau_set: ", paste(dim(tau_set), collapse = "x"),
         ", kappa_set: ", paste(dim(kappa_set), collapse = "x"))
  }

  # Input validation: check if there are at least 1 study
  if (ncol(tau_set) < 1 || nrow(tau_set) < 1) {
    stop("tau_set and kappa_set must have at least 1 row and 1 column. ",
         "Current dimensions: ", paste(dim(tau_set), collapse = "x"))
  }

  num_study <- ncol(kappa_set)
  if (is.null(colnames(kappa_set))) {
    colnames(kappa_set) <- paste0("kappa", 1:num_study)
  }
  if (is.null(colnames(tau_set))) {
    colnames(tau_set) <- paste0("tau", 1:num_study)
  }

  # Input validation: check for NA, NaN, or Inf values
  if (anyNA(tau_set) || any(!is.finite(tau_set))) {
    stop("tau_set contains NA, NaN, or Inf values.")
  }
  if (anyNA(kappa_set) || any(!is.finite(kappa_set))) {
    stop("kappa_set contains NA, NaN, or Inf values.")
  }

  # Input validation: check M
  if (missing(M)) {
    stop("M is a required argument.")
  }
  if (!is.numeric(M) || length(M) != 1 || M <= 0 || M != as.integer(M)) {
    stop("M must be a single positive integer. ",
         "Current value: ", M)
  }

  # Input validation: check M is >= max(kappa_set)
  max_kappa <- max(kappa_set, na.rm = TRUE)
  if (M < max_kappa) {
    stop("M must be greater than or equal to the maximum value in kappa_set. ",
         "M: ", M, ", max(kappa_set): ", max_kappa)
  }

  # Input validation: check seed
  if (!is.numeric(seed) || length(seed) != 1) {
    stop("seed must be a single numeric value. ",
         "Current value: ", seed)
  }

  # Input validation: check genetic_info dimensions if provided
  if (!is.null(genetic_info)) {
    if (!is.data.frame(genetic_info) && !is.matrix(genetic_info)) {
      stop("genetic_info must be a data frame or matrix.")
    }
    if (nrow(genetic_info) != nrow(tau_set)) {
      stop("genetic_info must have the same number of rows as tau_set. ",
           "genetic_info rows: ", nrow(genetic_info),
           ", tau_set rows: ", nrow(tau_set))
    }
  }

  n <- nrow(tau_set)

  sign <- as.matrix(kappa_set)
  W_mag <- apply(tau_set, 1, combine_mag, r = 1)

  set.seed(seed + 1)
  pvals <- apply(
    sign,
    MARGIN = 1,
    FUN = calc_p_value,
    r = 1,
    M = M
  )
  ord <- order(-W_mag, pvals)
  pvals_sorted <- pvals[ord]

  ### ForwardStop
  hfun_FS <- create_ForwardStop_function()
  FDPest_order <- AccumulationTest(pvals_sorted, hfun_FS, output_type = 'FDPest')

  FDPest <- numeric(length(pvals))
  FDPest[ord] <- FDPest_order
  FDPest[FDPest >= 1] <- 1

  all_kappa_positive <- apply(kappa_set, 1, function(x) all(x > 0))
  FDPest[all_kappa_positive] <- 1

  # Build result dataframe
  res <- data.frame(kappa_set, tau_set, qvalue = FDPest)

  # Append genetic_info if provided
  if (!is.null(genetic_info)) {
    res <- cbind(genetic_info, res)
  }

  return(res)
}

#' KAMA Meta-Analysis Based Variable Selection
#'
#' Performs KAMA meta-analysis and selects variables based on FDR control and stability thresholds.
#'
#' @param kappa_set A numeric matrix of kappa statistics. Rows correspond to variables,
#'   and columns correspond to studies.
#'   **Note:** The order of variables (rows) and studies (columns) must strictly match that of \code{tau_set}.
#' @param tau_set A numeric matrix of tau statistics. Rows correspond to variables,
#'   and columns correspond to studies.
#'   **Note:** The order of variables (rows) and studies (columns) must strictly match that of \code{kappa_set}.
#' @param M Integer. The number of knockoffs used. Must be >= max(kappa_set).
#' @param q Numeric. The FDR control level (default 0.1).
#' @param genetic_info Optional data frame or matrix containing additional information (e.g., SNP IDs, positions) to be appended to the results.
#' @param rep Integer. Number of repetitions for stability selection (default 1).
#' @param Uj_thre Numeric. Stability threshold for selection (default 0.51).
#' @param seed Integer. Random seed for reproducibility (default 111).
#'
#' @return A data frame of the **selected** variables containing:
#' \itemize{
#'   \item \code{genetic_info}: (If provided) The appended metadata.
#'   \item \code{kappa_set}: The kappa statistics for the selected variables.
#'   \item \code{tau_set}: The tau statistics for the selected variables.
#'   \item \code{qvalue}: KAMA's q-value for meta-analysis FDR control.
#' }
#' Returns an empty data frame if no variables are selected.
#' @export
KAMA_meta_analysis_select <- function(kappa_set, tau_set, M, q = 0.1, genetic_info = NULL, rep = 1, Uj_thre=0.51, seed = 111){

  # Convert to matrix if necessary
  tau_set <- as.matrix(tau_set)
  kappa_set <- as.matrix(kappa_set)

  # Input validation: check dimensions of tau_set and kappa_set
  if (!all(dim(tau_set) == dim(kappa_set))) {
    stop("tau_set and kappa_set must have the same dimensions. ",
         "tau_set: ", paste(dim(tau_set), collapse = "x"),
         ", kappa_set: ", paste(dim(kappa_set), collapse = "x"))
  }

  # Input validation: check if there are at least 1 study
  if (ncol(tau_set) < 1 || nrow(tau_set) < 1) {
    stop("tau_set and kappa_set must have at least 1 row and 1 column. ",
         "Current dimensions: ", paste(dim(tau_set), collapse = "x"))
  }

  num_study <- ncol(kappa_set)
  if (is.null(colnames(kappa_set))) {
    colnames(kappa_set) <- paste0("kappa", 1:num_study)
  }
  if (is.null(colnames(tau_set))) {
    colnames(tau_set) <- paste0("tau", 1:num_study)
  }

  # Input validation: check for NA, NaN, or Inf values
  if (anyNA(tau_set) || any(!is.finite(tau_set))) {
    stop("tau_set contains NA, NaN, or Inf values.")
  }
  if (anyNA(kappa_set) || any(!is.finite(kappa_set))) {
    stop("kappa_set contains NA, NaN, or Inf values.")
  }

  # Input validation: check M
  if (missing(M)) {
    stop("M is a required argument.")
  }
  if (!is.numeric(M) || length(M) != 1 || M <= 0 || M != as.integer(M)) {
    stop("M must be a single positive integer. ",
         "Current value: ", M)
  }

  # Input validation: check q is in valid range
  if (!is.numeric(q) || length(q) != 1 || q < 0 || q > 1) {
    stop("q must be a single numeric value between 0 and 1. ",
         "Current value: ", q)
  }

  # Input validation: check M is >= max(kappa_set)
  max_kappa <- max(kappa_set, na.rm = TRUE)
  if (M < max_kappa) {
    stop("M must be greater than or equal to the maximum value in kappa_set. ",
         "M: ", M, ", max(kappa_set): ", max_kappa)
  }

  # Input validation: check seed
  if (!is.numeric(seed) || length(seed) != 1) {
    stop("seed must be a single numeric value. ",
         "Current value: ", seed)
  }

  # Input validation: check genetic_info dimensions if provided
  if (!is.null(genetic_info)) {
    if (!is.data.frame(genetic_info) && !is.matrix(genetic_info)) {
      stop("genetic_info must be a data frame or matrix.")
    }
    if (nrow(genetic_info) != nrow(tau_set)) {
      stop("genetic_info must have the same number of rows as tau_set. ",
           "genetic_info rows: ", nrow(genetic_info),
           ", tau_set rows: ", nrow(tau_set))
    }
  }

  n <- nrow(tau_set)

  sign <- as.matrix(kappa_set)

  select_all <- c()
  W_mag = apply(tau_set, 1, combine_mag, r = 1)

  for (sd in 1:rep){

    set.seed(seed+sd)
    pvals = apply(sign, MARGIN = 1, FUN = calc_p_value, r = 1, M = M)

    ### improve the p-value
    pvals_sorted = pvals[order(-W_mag)]

    ### ForwardStop
    hfun_FS <- create_ForwardStop_function()
    select_num = AccumulationTest(pvals_sorted, hfun_FS, alpha=q)
    if (select_num==0){
      select = NULL
    }else{
      select = 1:select_num
    }

    select_all <- c(select_all,method_selection(select,W_mag))
  }

  ### select variables
  select_idx <- sort(final_select(select_all,rep,Uj_thre))

  if (length(select_idx) > 0) {
    selected_kappa_rows <- kappa_set[select_idx, , drop = FALSE]
    has_zeros <- rowSums(selected_kappa_rows == 0) > 0
    select_idx <- select_idx[has_zeros]
  }

  if (!is.null(genetic_info)){
    all_info <- data.frame(genetic_info, kappa_set, tau_set)
    res <- all_info[select_idx,]
  }else{
    all_info <- data.frame(kappa_set, tau_set)
    res <- all_info[select_idx,]
  }

  return(res)
}


#' KAMA Localization Analysis
#'
#' The main pipeline for KAMA to identify study-specific associations.
#' It performs meta-analysis and localizes the signal to specific studies/cohorts
#'
#' @param kappa_set A numeric matrix of kappa statistics. Rows correspond to variables,
#'   and columns correspond to studies.
#'   **Note:** The order of variables (rows) and studies (columns) must strictly match that of \code{tau_set}.
#' @param tau_set A numeric matrix of tau statistics. Rows correspond to variables,
#'   and columns correspond to studies.
#'   **Note:** The order of variables (rows) and studies (columns) must strictly match that of \code{kappa_set}.
#' @param M Integer. The number of knockoffs used.
#' @param q Numeric. The FDR control level (default 0.1).
#' @param genetic_info Optional data frame containing variable metadata.
#' @param study_label Character vector. Custom labels for the studies (columns). Defaults to "study1", "study2", etc.
#' @param seed Integer. Random seed.
#'
#' @return A data frame of the significant variables containing:
#' \itemize{
#'   \item \code{genetic_info}: (If provided) The appended metadata.
#'   \item \code{kappa_set}: The kappa statistics.
#'   \item \code{tau_set}: The tau statistics.
#'   \item \code{qvalue}: KAMA's q-value for meta-analysis FDR control.
#'   \item \code{localization}: A character column specifying the studies/cohorts where the variable is significant.
#' }
#' @export
KAMA_localization <- function(kappa_set, tau_set, M, q = 0.1, genetic_info = NULL, study_label = NULL, seed = 111){

  # Input validation: check if tau_set and kappa_set are provided
  if (missing(tau_set) || missing(kappa_set)) {
    stop("tau_set and kappa_set are required arguments.")
  }

  # 转换为矩阵（data.frame会被转换为matrix）
  tau_set <- as.matrix(tau_set)
  kappa_set <- as.matrix(kappa_set)

  # Input validation: check dimensions of tau_set and kappa_set
  if (!all(dim(tau_set) == dim(kappa_set))) {
    stop("tau_set and kappa_set must have the same dimensions. ",
         "tau_set: ", paste(dim(tau_set), collapse = "x"),
         ", kappa_set: ", paste(dim(kappa_set), collapse = "x"))
  }

  # Input validation: check if there are at least 1 study
  if (ncol(tau_set) < 1 || nrow(tau_set) < 1) {
    stop("tau_set and kappa_set must have at least 1 row and 1 column. ",
         "Current dimensions: ", paste(dim(tau_set), collapse = "x"))
  }

  num_study <- ncol(kappa_set)
  if (is.null(colnames(kappa_set))) {
    colnames(kappa_set) <- paste0("kappa", 1:num_study)
  }
  if (is.null(colnames(tau_set))) {
    colnames(tau_set) <- paste0("tau", 1:num_study)
  }

  # Input validation: check for NA, NaN, or Inf values
  if (anyNA(tau_set) || any(!is.finite(tau_set))) {
    stop("tau_set contains NA, NaN, or Inf values.")
  }
  if (anyNA(kappa_set) || any(!is.finite(kappa_set))) {
    stop("kappa_set contains NA, NaN, or Inf values.")
  }

  # Input validation: check M
  if (missing(M)) {
    stop("M is a required argument.")
  }
  if (!is.numeric(M) || length(M) != 1 || M <= 0 || M != as.integer(M)) {
    stop("M must be a single positive integer. ",
         "Current value: ", M)
  }

  # Input validation: check q is in valid range
  if (!is.numeric(q) || length(q) != 1 || q < 0 || q > 1) {
    stop("q must be a single numeric value between 0 and 1. ",
         "Current value: ", q)
  }

  # Input validation: check M is >= max(kappa_set)
  max_kappa <- max(kappa_set, na.rm = TRUE)
  if (M < max_kappa) {
    stop("M must be greater than or equal to the maximum value in kappa_set. ",
         "M: ", M, ", max(kappa_set): ", max_kappa)
  }

  # Input validation: check seed
  if (!is.numeric(seed) || length(seed) != 1) {
    stop("seed must be a single numeric value. ",
         "Current value: ", seed)
  }

  # Input validation: check genetic_info dimensions if provided
  if (!is.null(genetic_info)) {
    if (!is.data.frame(genetic_info) && !is.matrix(genetic_info)) {
      stop("genetic_info must be a data frame or matrix.")
    }
    if (nrow(genetic_info) != nrow(tau_set)) {
      stop("genetic_info must have the same number of rows as tau_set. ",
           "genetic_info rows: ", nrow(genetic_info),
           ", tau_set rows: ", nrow(tau_set))
    }
  }

  # Input validation and default generation for study_label
  if (is.null(study_label)) {
    study_label <- paste0("study", 1:num_study)
  } else {
    if (length(study_label) != num_study) {
      stop("study_label must have the same length as the number of studies (columns in kappa_set/tau_set). ",
           "study_label length: ", length(study_label),
           ", number of studies: ", num_study)
    }
    # Check for duplicate study labels
    if (any(duplicated(study_label))) {
      stop("study_label contains duplicate values. Each study must have a unique label.")
    }
  }

  meta_res <- KAMA_meta_analysis(kappa_set, tau_set, M, genetic_info)

  selected_sig <- lapply(1:num_study, function(i)
    KAMA_select(kappa_set, tau_set, M, r = i, q = q)
  )

  KAMA_env_diff_snps <- env_setdiff_sig(selected_sig)

  KAMA_env_diff_detailed <- lapply(KAMA_env_diff_snps, function(rows) {
    if (length(rows) == 0) {
      meta_res[0, , drop = FALSE]
    } else {
      meta_res[rows, , drop = FALSE]
    }
  })

  KAMA_env_spec_res <- lapply(1:num_study, function(i) {
    multi_env_spec_func(kappa_set, tau_set, KAMA_env_diff_snps[[i]], study_label, i)
  })

  KAMA_env_spec_res_total <- do.call(rbind, KAMA_env_spec_res)

  ind_match <- as.numeric(KAMA_env_spec_res_total[, "row_idx"])
  info <- meta_res[ind_match, ,drop = FALSE]
  loc_vec <- KAMA_env_spec_res_total[, "localization"]

  KAMA_env_spec_res_pos <- data.frame(
    info,
    row_idx = ind_match,
    localization = loc_vec,
    stringsAsFactors = FALSE
  )

  KAMA_env_spec_res_pos <- KAMA_env_spec_res_pos[order(KAMA_env_spec_res_pos[, "row_idx"]), ]
  KAMA_local_res <- KAMA_env_spec_res_pos[KAMA_env_spec_res_pos[, "qvalue"] <= q, ]

  rownames(KAMA_local_res) <- NULL
  return(KAMA_local_res)
}
