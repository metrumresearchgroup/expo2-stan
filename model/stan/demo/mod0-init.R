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
           tv_log_ec50 = rnorm(1,log(100),1),
           log_gamma = rnorm(1,0,.1),
           eta_e0 = rnorm(.data$n_id),
           eta_log_ec50 = rnorm(.data$n_id),
           omega_e0 = abs(rnorm(1)),
           omega_log_ec50 = abs(rnorm(1)),
           sigma = abs(rnorm(1,10,2)))
    }
  )
}
