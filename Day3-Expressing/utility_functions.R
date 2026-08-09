
#Functions for printing out equations from contrast matrices.

neat_term <- function(coef_i, name_i, i){

  use_coef_i <- abs(coef_i)
  if(use_coef_i==1) {
    use_coef_i <- ""
  }

  if(i>1){
    if(coef_i>0) {
      tmp <- paste0("+ ", paste0(use_coef_i, name_i))
    } else {
      tmp <- paste0("- ", paste0(use_coef_i, name_i))
    }
  } else if (i==1){

    if(coef_i>0) {
      tmp <- paste0(use_coef_i, name_i)
    } else {
      tmp <- paste0("-", use_coef_i, name_i)
    }

  }

  tmp
}

read_rows <- function(contrast_mat, rowname){
  selected_row <- contrast_mat[rowname,]
  nonzero_entries <- selected_row[selected_row!=0]
  form <- paste0(rowname, " =")
  if(length(nonzero_entries)>0){
    for(i in 1:length(nonzero_entries)){
      coef_i <- nonzero_entries[i]
      name_i <- names(nonzero_entries)[i]
      term = neat_term(coef_i, name_i, i)
      form <- paste(form, term)
    }
  }
  form
}

#Functions for summarizing posterior predictives
calc_stat_postfit <- function(dat, pp, facs, stat_fn, stat_cols) {

  model_stat <- pp %>% group_by(across(all_of(c("postn", facs))))%>%
    summarise(
      poststat= stat_fn(!!!syms(stat_cols)), .groups = "drop_last") %>% group_by(across(all_of(facs))) %>%
    summarise(postmn_stat=mean(poststat), lower=quantile(poststat, prob=0.025),
              upper=quantile(poststat, prob=0.975), .groups = "drop_last")

  data_stat <- dat %>% group_by(across(all_of(facs))) %>%
    summarise(
      stat= stat_fn(!!!syms(stat_cols)), .groups = "drop_last")

  full_join(model_stat, data_stat, by=facs)
}

stack_stat_postfits <- function(dat, pp, facs, stat_fn_list, stat_col_list, stat_fn_names){

  for(i in 1:length(stat_fn_list)){
    stat <- calc_stat_postfit(dat, pp, facs=facs,
                              stat_fn= stat_fn_list[[i]],
                              stat_cols=stat_col_list[[i]]
    )
    stat$nam <- stat_fn_names[i]
    if(i==1) stats <- stat else stats <- rbind(stats, stat)
  }
  stats
}

