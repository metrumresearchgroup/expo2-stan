# Create Stan data
#
# This function must return the list that will be passed to `data` argument
#   of `cmdstanr::sample()`
#
# The `.dir` argument represents the absolute path to the directory containing
#   this file. This is useful for building file paths to the input files you will
#   load. Note: you _don't_ need to pass anything to this argument, you only use
#   it within the function. `bbr` will pass in the correct path when it calls
#   `make_standata()` under the hood.
make_standata <- function(.dir) {
  # read in any input data
  in_data <- readr::read_csv(
    file.path(.dir, "..", "..", "..", "..", "data/derived/exdata3.csv"),
    show_col_types = FALSE
  )
  
  in_data <- in_data %>% mutate(row=1:n())
  
  list(
    N = nrow(in_data),
    n_id = length(unique(in_data$ID)),
    K = 1, 
    DV = in_data$fxa.inh,
    X_e0 = in_data %>% distinct(ID,sex) %>% mutate(sex=sex-1) %>% pull(sex) %>% as.matrix(),
    X_ec50 = in_data %>% distinct(ID,sex) %>% mutate(sex=sex-1) %>% pull(sex) %>% as.matrix(),
    Conc = in_data$cobs,
    ID = in_data$ID,
    id_start = in_data %>% group_by(ID) %>% summarise(min=min(row)) %>% pull(min),
    id_end = in_data %>% group_by(ID) %>% summarise(max=max(row)) %>% pull(max),
    prior_only=0,
    n_sim = 1000
  )
}
