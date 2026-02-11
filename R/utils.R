##--------- MEKF -----------------

#' @noRd
create_ForwardStop_function = function(){
  function(x){log(1/(1-x))}
}

#' @noRd
ForwardStop = function(pvals,alpha=0.2,output_type='khat'){
  AccumulationTest(pvals,create_ForwardStop_function(),alpha=alpha,output_type=output_type,check_integrate_to_one=FALSE)
}

#' @noRd
AccumulationTest = function(pvals,hfun,alpha=0.2,numerator_plus=0,denominator_plus=0,output_type='khat',check_integrate_to_one=TRUE){

  # check for valid arguments
  check_inputs = CheckInputs_AccumulationTest(pvals,hfun,alpha,numerator_plus,denominator_plus,output_type,check_integrate_to_one)
  if(length(check_inputs)>0){
    stop(check_inputs)
  }

  # perform the test
  n=length(pvals)
  FDPest=(numerator_plus+cumsum(unlist(lapply(pvals,hfun))))/(denominator_plus+1:n)
  FDPest_vs_alpha=(FDPest%*%t(rep(1,length(alpha)))<=rep(1,n)%*%t(alpha))
  findlast=function(x){max(c(0,which(x)))}
  khat=apply(FDPest_vs_alpha,2,findlast)
  if(output_type=='khat'){
    return(khat)
  }else{if(output_type=='FDPest'){
    return(FDPest)
  }else{
    output=list()
    output$FDPest=FDPest
    output$khat=khat
    return(output)
  }}
}

#' @noRd
CheckInputs_AccumulationTest = function(pvals,hfun,alpha=0.2,numerator_plus=0,denominator_plus=0,output_type='khat', check_integrate_to_one){
  # check_integrate_to_one should be logical
  if(!is.logical(check_integrate_to_one)){
    return('check_integrate_to_one must be logical')
  }
  # check that pvals and alpha are each sequences with values in [0,1]
  if(!is.numeric(pvals) || !is.vector(pvals) || min(pvals)<0 || max(pvals)>1){
    return('pvals must be a number or numeric vector with values in [0,1]')
  }
  n=length(pvals)

  if(!is.numeric(alpha) || !is.vector(alpha) || min(alpha)<0 || max(alpha)>1){
    return('alpha must be a number or numeric vector with values in [0,1]')
  }

  # check that hfun is a function that gives a nonnegative value for each pvalue
  if(!is.function(hfun)){
    return('hfun must be a function')
  }
  if(!is.numeric(try(hfun(pvals),silent=TRUE)) || any(is.na(try(hfun(pvals),silent=TRUE)))){
    return('The function hfun must take as input a vector of p-values in [0,1], and return a vector of nonnegative numbers of the same length')
  }
  if(length(hfun(pvals))!=length(pvals)){
    return('The function hfun must take as input a vector of p-values in [0,1], and return a vector of nonnegative numbers of the same length')
  }
  if(any(hfun(pvals)<0)){
    return('The function hfun must take as input a vector of p-values in [0,1], and return a vector of nonnegative numbers of the same length')
  }
  if(check_integrate_to_one){
    if(abs(integrate(hfun,0,1)$value-1)>1e-2){
      return('The function hfun must have expectation 1 when applied to a uniform variable p~Uniform[0,1] (set check_integrate_to_one=FALSE to override)')
    }
  }

  # check that numerator_plus and denominator_plus are numbers
  if(!is.numeric(numerator_plus) || !is.numeric(denominator_plus) || length(numerator_plus)!=1 || length(denominator_plus)!=1){
    return('numerator_plus and denominator_plus must each be a scalar')
  }

  # check that output_type is in {'khat', 'FDPest', 'both'}
  if(!is.element(output_type,c('khat','FDPest','both'))){
    return('Invalid output type: choose from "khat", "FDPest", or "both"')
  }

  return(NULL)
}

#' @noRd
calc_p_value = function(vec, r, M){

  # l = length(vec)
  l = length(vec)
  kappa_null <- sum(vec > 0)

  a = l - r + 1

  if (kappa_null > a){
    kappa_null <- a
  }

  if (a > 0){
    discrete_p = pbinom(kappa_null,a,M/(M+1))
    discrete_p_l = pbinom(kappa_null - 1,a,M/(M+1))
  } else {
    discrete_p = 1
    discrete_p_l = 0
  }
  p_cont = runif(1,discrete_p_l,discrete_p)

  return(p_cont)
}

#' @noRd
combine_mag = function(vec, r){
  vec_abs = abs(vec)
  vec_sort = -sort(-vec_abs, partial = r)[1:r]
  return(prod(vec_sort))
}

#' @noRd
final_select <- function(select_all, rep, Uj_thre){
  table_tmp <- table(select_all)/rep
  selecteda <- as.numeric(names(table_tmp[table_tmp>=Uj_thre]))

  return(selecteda)
}

#' @noRd
method_selection <- function(select,W_mag){
  select_ind <- NULL
  if (length(select) > 0){
    select_ind = order(-W_mag)[select]
  }
  return(select_ind)
}

##--------- localization --------------
#' @noRd
env_setdiff_sig <- function(big_list) {
  n <- length(big_list)
  new_list <- vector("list", n)

  new_list[[n]] <- big_list[[n]]

  for (i in (n - 1):1) {

    to_subtract <- unique(unlist(big_list[(i + 1):n]))

    new_list[[i]] <- sort(setdiff(big_list[[i]], to_subtract))
  }

  return(new_list)
}

#' @noRd
multi_env_spec_func <- function(kappa_set,tau_set,diff_ele,pops_label,num_pop){

  res <- data.frame(
    row_idx = numeric(), localization = character()
  )

  if (length(diff_ele)){

    kappa_set_tmp = kappa_set[diff_ele, , drop = FALSE]
    tau_set_tmp = tau_set[diff_ele, , drop = FALSE]
    num_sig = length(diff_ele)

    pop = rep(0,num_sig)
    for (i in 1:num_sig){
      # i = 2
      ind_tmp <- order(as.numeric(tau_set_tmp[i, , drop = FALSE]), decreasing = TRUE)[1:num_pop]
      sig_pop <- which(kappa_set_tmp[i,ind_tmp]==0)
      sig_ind <- ind_tmp[sig_pop]

      if (length(sig_pop)==0){
        next
      }else{
        sig_pop_info <- data.frame(
          row_idx = diff_ele[i],
          localization = paste0(pops_label[sort(sig_ind)], collapse = " "))

        res <- rbind(res, sig_pop_info)
      }
    }
  }

  return(res)
}

#' @noRd
env_spec_reorg_format <- function(df, pops_label) {

  # Step 1: split pred.label
  df_long <- df %>%
    separate_rows(localization, sep = " ") %>%
    mutate(value = TRUE)

  # Step 2: extend to wide，fill FALSE
  effect_binary_df <- df_long %>%
    pivot_wider(
      id_cols = row_idx,
      names_from = localization,
      values_from = value,
      values_fill = FALSE
    )

  # Step 3: add missing pop
  missing_cols <- setdiff(pops_label, colnames(effect_binary_df))
  if (length(missing_cols) > 0) {
    effect_binary_df[missing_cols] <- FALSE
  }

  # Step 4: sort column
  effect_binary_df <- effect_binary_df %>%
    dplyr::select(row_idx, all_of(pops_label))

  # Step 5:
  effect_binary <- effect_binary_df %>%
    arrange(row_idx) %>%
    column_to_rownames("row_idx") %>%
    as.matrix()

  return(effect_binary)
}

#' @noRd
KAMA_select <- function(kappa_set, tau_set, r, M, q = 0.1, rep = 1, Uj_thre=0.51, seed = 111){

  n <- nrow(tau_set)

  sign <- as.matrix(kappa_set)

  select_all <- c()
  W_mag = apply(tau_set, 1, combine_mag, r = r)

  for (sd in 1:rep){

    set.seed(seed+sd)
    pvals = apply(sign, MARGIN = 1, FUN = calc_p_value, r = r, M = M)

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
  res <- sort(final_select(select_all,rep,Uj_thre))

  return(res)
}
