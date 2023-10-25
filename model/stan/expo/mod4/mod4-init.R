# Create Stan initial values
#
# This function must return something that can be passed to the `init` argument
#   of `cmdstanr::sample()`. There are several options; see `?cmdstanr::sample`
#   for details.
#
# `.data` represents the list returned from `make_standata()` for this model.
#   This is provided in case any of your initial values are dependent on some
#   aspect of the data (e.g. the number of rows).
#
# `.args` represents the list of attached arguments that will be passed through to
#   cmdstanr::sample(). This is provided in case any of your initial values are
#   dependent on any of these arguments (e.g. the number of chains).
#
# Note: you _don't_ need to pass anything to either of these arguments, you only
#   use it within the function. `bbr` will pass in the correct objects when it calls
#   `make_init()` under the hood.
#
make_init <- function(.data, .args) {
  # returning NULL causes cmdstanr::sample() to use the default initial values
  return(
    function() {
      list(tv_e0 = rnorm(1),
           emax = runif(1,90,100),
           beta_e0 = as.array(rnorm(.data$K,0,1)),
           beta_log_ec50 = as.array(rnorm(.data$K,0,1)),
           tv_log_ec50 = rnorm(1,log(100),1),
           log_gamma = rnorm(1,0,.1),
           etas = matrix(rnorm(2*.data$n_id, 0, 1),nrow=2),
           omega_sds = abs(rnorm(2)),
           L_Omega_corr = diag(2),
           sigma = abs(rnorm(1,10,2)))
    }
  )
}
