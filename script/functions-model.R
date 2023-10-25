#################################
# Helper functions for modeling
#
# Some of these may eventually
# make it into bbr or bbr.bayes
#################################

#' Add loo_compare() columns to a run log
#' @inheritParams bbr::run_log
#' @param ll_regex String, a pattern to match for variable(s) to pass to .fit$loo(variables)
#' @param name_col String, column name in .log_df to use as names for loo_compare output
loo_compare.bbi_log_df <- function(.log_df, ..., ll_regex = "^log_lik", name_col = "run") {
  require(loo)
  checkmate::assert_string(ll_regex)
  checkmate::assert_string(name_col)
  if (!(name_col %in% names(.log_df))) stop("name_col must be one of ", paste(names(.log_df), collapse = ", "))
  
  purrr::map(.log_df$absolute_model_path, ~{
    .fit <- try(bbr.bayes::read_fit_model(.x), silent = TRUE)
    if (inherits(.fit, "try-error")) return(NULL)
    ll_var <- grep(ll_regex, .fit$metadata()$variables, value = TRUE) %>% gsub(pattern = "\\[.+", x = ., "") %>% unique()
    if (length(ll_var)) .fit$loo(variables = ll_var)
  }) %>% 
    purrr::set_names(.log_df[[name_col]]) %>% 
    purrr::compact() %>%
    loo::loo_compare(...)
}



#' Add loo_compare() columns to a run log
#' @inheritParams bbr::run_log
#' @param ll_regex String, a pattern to match for variable(s) to pass to .fit$loo(variables)
add_loo_compare <- function(.log_df, ...) {
  require(loo)
  checkmate::assert_class(.log_df, "bbi_log_df")
  
  loo_df <- loo_compare(.log_df, ..., name_col = "absolute_model_path") %>% 
    as.data.frame() %>%
    tibble::rownames_to_column(var = "absolute_model_path")
  
  dplyr::left_join(.log_df, loo_df, by = "absolute_model_path")
}


#######################
# render diagnostic template
# from model object

model_diagnostics <- function(
    .mod, 
    .p = list(),
    template = here::here("script", "diagnostic-templates", "diagnostics-basic.Rmd"),
    ...
) {
  UseMethod("model_diagnostics")
}

#' Renders diagnostic template from bbr Stan model
model_diagnostics.bbi_stan_model <- function(
    .mod, 
    .p = list(),
    template = here::here("script", "diagnostic-templates", "stan-diagnostics-basic.Rmd")
) {
  checkmate::assert_list(.p, names = "named")
  checkmate::assert_string(template)
  
  .p$mod <- .mod
  
  ##### experimental idea:
  ##### pull any param that doesn't have `[\\d]` in the name for plotting
  .r <- read_fit_model(.mod)
  model_params <- .r$metadata()$variables
  if (is.null(.p$pars)) {
    .p$pars <- model_params %>%
      stringr::str_subset("(lp__)|(\\[\\d+\\]$)", negate = TRUE)
  } else if (checkmate::test_string(.p$pars, pattern = "^all$")) {
    .p$pars <- model_params
  } else {
    wrong_params <- !(.p$pars %in% model_params)
    if (any(wrong_params)) {
      stop(paste(
        "Requested the following param(s) not present in the model:",
        paste(.p$pars[wrong_params], collapse = ", ")
      ))
    }
  }
  #####
  #TODO: add pars_regex, to be parsed here too ^
  
  out_dir <- get_output_dir(.mod)
  out_file <- paste(get_model_id(.mod), gsub('\\.[Rr][Mm][Dd]$', '.html', basename(template)), sep = "-")
  rmarkdown::render(
    template,
    params = .p,
    output_dir = out_dir,
    output_file = out_file
  )
  
  out_path <- file.path(out_dir, out_file)
  message(paste("Diagnostic HTML saved to", out_path))
  
  return(invisible(out_path)) # return path so user can pipe to browseURL()
}


#' Alias for model_diagnostics.bbi_stan_model that calls MCMC template
mcmc_diagnostics <- function(
    .mod, 
    .p = list(),
    template = here::here("script", "diagnostic-templates", "mcmc-diagnostics.Rmd")
) {
  checkmate::assert_class(.mod, "bbi_stan_model")
  model_diagnostics(.mod, .p = .p, template = template)
}