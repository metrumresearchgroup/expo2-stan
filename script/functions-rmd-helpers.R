#these are helper functions for stepping through the demo rmd files 
#these are not functions you would use in a real project
library(fs)

# message to run models at top of each script. this_model should be full path to model you want to check
check_model <- function(this_model){
  tryCatch({
    read_fit_model(this_model)
  },
  error = function(.e) {
    cli::cli({
      cli::cli_h1("Warning") 
      cli::cli_text(glue::glue("{basename(this_model)} has not been run. 
                             Please run the following: \n
                             {basename(this_model)} <- read_model({glue::double_quote(this_model)}) \n
                             submit_model({basename(this_model)})"))
    })
  })
}

#this code reverts model files in expo folder back to what they were before `copy_model` or `new_model` was run. This does not copy back output folder
revert_model <- function(.mod, .expo_model_dir, .demo_model_dir){
                               
  all_files <- dir_ls(.demo_model_dir, recurse = 0, type = "file")
  stan_files <- all_files[grepl(paste0(.mod, ".stan"), basename(all_files))]
  files_we_want <- all_files[grepl(paste0(.mod, "-"), basename(all_files))]
  files_we_want <- files_we_want[!(grepl(".csv", files_we_want))]
  files_we_want <- all_files[all_files %in% c(files_we_want, stan_files)]
  purrr::walk(files_we_want, ~ file_copy(.x, here(.expo_model_dir, .mod), overwrite = TRUE))
  
}

