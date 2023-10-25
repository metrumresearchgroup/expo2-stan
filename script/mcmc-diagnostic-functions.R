###
### Bill Gillespie's MCMC diagnostic functions
### as of 26 January 2023
###

mcmc_history <- function(posterior, pars = dimnames(posterior)[[3]], 
                         nParPerPage = 5, myTheme = NULL, np = NULL){
  ## Create MCMC history plots
  ## posterior = 3-D array of MCMC samples. Dims = {iterations, chains, parameters}. 
  ## pars = names of parameters to plot
  ## nParPerPage = maximum number of parameters to plot per page
  ## myTheme = ggplot2 theme
  
  require(bayesplot)
  allPars <- dimnames(posterior)[[3]]
  pars <- unname(unlist(sapply(pars,
                               function(x){
                                 allPars[startsWith(allPars, x)]
                               })))
  posterior <- posterior[, , pars]
  ##  posterior <- posterior[, , unlist(sapply(c(paste("^", pars, "$",sep=""),
  ##                                             paste("^", pars, "\\[",sep="")),
  ##                                           grep, 
  ##                                           x = dimnames(posterior)[[3]]))]
  ##  pars <- dimnames(posterior)[[3]]
  nPars <- length(pars)
  nPost <- dim(posterior)[1]
  nPages <- ceiling(nPars / nParPerPage)
  parameters <- data.frame(parameter = pars,
                           page = sort(rep(1:nPages, length = nPars)),
                           stringsAsFactors = FALSE)
  
  lapply(1:nPages,
         function(i){
           postPage <- posterior[, , with(parameters, pars[page == i])]
           res <- mcmc_trace(postPage,
                             facet_args = list(ncol = 1, strip.position = "left"),
                             np = np) +
             myTheme +
             scale_x_continuous(breaks = seq(0, nPost, len = 5))
           ## Code to remedy an apparent bug in Bayesplot
           res$data$Value <- res$data$value
           res
         }
  )
}

mcmc_density <- function(posterior, 
                         pars = dimnames(posterior)[[length(dim(posterior))]],
                         byChain = FALSE, nParPerPage = 16, myTheme = NULL, prior = NULL){
  ## Create density plots for marginal distributions of MCMC samples
  ## posterior = 3-D array of MCMC samples. Dims = {iterations, chains, parameters}
  ## pars = names of parameters to plot
  ## byChain = logical indicating whether to plot sensities by chain
  ## nParPerPage = maximum number of parameters to plot per page
  ## myTheme = ggplot2 theme
  ## prior = data.frame with columns value, density and parameter
  require(bayesplot)
  ndim <- length(dim(posterior))
  if(ndim == 2){
    allPars <- dimnames(posterior)[[2]]
    pars <- unname(unlist(sapply(pars,
                                 function(x){
                                   allPars[startsWith(allPars, x)]
                                 })))
    posterior <- posterior[, pars]
    ##    posterior <- posterior[, unlist(sapply(c(paste("^", pars, "$",sep=""),
    ##                                               paste("^", pars, "\\[",sep="")),
    ##                                             grep, 
    ##                                             x = dimnames(posterior)[[2]]))]
  }
  if(ndim == 3){
    allPars <- dimnames(posterior)[[3]]
    pars <- unname(unlist(sapply(pars,
                                 function(x){
                                   allPars[startsWith(allPars, x)]
                                 })))
    posterior <- posterior[, , pars]
    ##    posterior <- posterior[, , unlist(sapply(c(paste("^", pars, "$",sep=""),
    ##                                               paste("^", pars, "\\[",sep="")),
    ##                                             grep,
    ##                                             x = dimnames(posterior)[[3]]))]
  }
  ##  pars <- dimnames(posterior)[[ndim]]
  nPars <- length(pars)
  nPages <- ceiling(nPars / nParPerPage)
  parameters <- data.frame(parameter = pars,
                           page = sort(rep(1:nPages, length = nPars)),
                           stringsAsFactors = FALSE)
  
  if(!is.null(prior)) prior <- prior %>% left_join(parameters)
  
  lapply(1:nPages,
         function(i){
           if(ndim ==2)
             postPage <- posterior[, with(parameters, pars[page == i])]
           else
             postPage <- posterior[, , with(parameters, pars[page == i])]
           
           if(byChain){
             p1 <- mcmc_dens_overlay(postPage)
           }else{
             p1 <- mcmc_dens(postPage)
           }
           if(!is.null(prior))
             p1 <- p1 + geom_line(data = subset(prior, page == i), 
                                  aes(x = value, y = density),
                                  color = "red")
           p1 <- p1 + myTheme
           ## Code to remedy an apparent bug in Bayesplot
           #p1$data$Value <- p1$data$value
           p1
         })
}

mcmc_ptable <- function(posterior, pars = dimnames(posterior)[[3]]){
  ## Make table summarizing parameter estimates and sampling performance
  ## posterior = 3-D array of MCMC samples. Dims = {iterations, chains, parameters}
  ## pars = names of parameters to plot
  
  ndim <- length(dim(posterior))
  if(ndim == 2){
    allPars <- dimnames(posterior)[[2]]
    pars <- unname(unlist(sapply(pars,
                                 function(x){
                                   allPars[startsWith(allPars, x)]
                                 })))
    posterior <- posterior[, pars]
    ##    posterior <- posterior[, unlist(sapply(c(paste("^", pars, "$",sep=""),
    ##                                               paste("^", pars, "\\[",sep="")),
    ##                                             grep, 
    ##                                             x = dimnames(posterior)[[2]]))]
  }
  if(ndim == 3){
    allPars <- dimnames(posterior)[[3]]
    pars <- unname(unlist(sapply(pars,
                                 function(x){
                                   allPars[startsWith(allPars, x)]
                                 })))
    posterior <- posterior[, , pars]
    ##    posterior <- posterior[, , unlist(sapply(c(paste("^", pars, "$",sep=""),
    ##                                               paste("^", pars, "\\[",sep="")),
    ##                                             grep,
    ##                                             x = dimnames(posterior)[[3]]))]
  }
  
  ptable <- summarize_draws(posterior) %>%
    mutate_if(is.numeric, ~formatC(., 3)) %>% 
    rename(parameter = variable) %>%
    mutate("90% CI" = paste("(", q5, ", ", q95, ")", 
                            sep = "")) %>%
    select(parameter, mean, median, sd, mad, "90% CI", ess_bulk, ess_tail, rhat)
  ptable
}

mcmc_plots <- function(posterior, pars = dimnames(posterior)[[3]], 
                       np = NULL, myTheme = NULL){
  
  allPars <- dimnames(posterior)[[3]]
  pars <- unname(unlist(sapply(pars,
                               function(x){
                                 allPars[startsWith(allPars, x)]
                               })))
  posterior <- posterior[, , pars]
  ##  posterior <- posterior[, , unlist(sapply(c(paste("^", pars, "$",sep=""),
  ##                                             paste("^", pars, "\\[",sep="")),
  ##                                           grep,
  ##                                           x = dimnames(posterior)[[3]]))]
  nSamp <- dim(posterior)[1]
  nChains <- dim(posterior)[2]
  
  diags <- summarize_draws(posterior, "rhat", "ess_bulk", "ess_tail")
  
  rhats <- diags$rhat
  names(rhats) <- diags$variable
  rhat_plot <- mcmc_rhat(rhats) + yaxis_text() + myTheme
  
  ess_bulk <- diags$ess_bulk
  names(ess_bulk) <- diags$variable
  ess_bulk_ratios <- ess_bulk / (nSamp * nChains)
  ess_bulk_plot <- mcmc_neff(ess_bulk_ratios) + yaxis_text() + myTheme +
    labs(title = "Bulk ESS ratios")
  
  ess_tail <- diags$ess_tail
  names(ess_tail) <- diags$variable
  ess_bulk_ratios <- ess_bulk / (nSamp * nChains)
  ess_tail_plot <- mcmc_neff(ess_bulk_ratios) + yaxis_text() + myTheme +
    labs(title = "Bulk ESS ratios")
  
  mcmc_history_plots <- mcmc_history(posterior, pars = parametersToPlot, 
                                     nParPerPage = 4, myTheme = myTheme, np = np)
  mcmc_density_chain_plots <- mcmc_density(posterior, pars = parametersToPlot, 
                                           nParPerPage = 16, byChain = TRUE, 
                                           myTheme = theme(text = element_text(size = 12),
                                                           axis.text = element_text(size = 10)))
  mcmc_density_plots <- mcmc_density(posterior, pars = parametersToPlot, 
                                     nParPerPage = 16, 
                                     myTheme = theme(text = element_text(size = 12),
                                                     axis.text = element_text(size = 10)))
  
  list(rhat_plot = rhat_plot,
       ess_bulk_plot = ess_bulk_plot,
       ess_tail_plot = ess_tail_plot,
       mcmc_history_plots = mcmc_history_plots,
       mcmc_density_chain_plots = mcmc_density_chain_plots,
       mcmc_density_plots = mcmc_density_plots)
}




mcmc_plots_v2 <- function(.mod, 
                          .pars = NULL, #dimnames(posterior)[[3]], 
                          .n_per_page = 4,
                          np = NULL, 
                          myTheme = NULL){
  
  ###--->> Add pairs plots?
  ###--->> Add regx pairs like bayesplot allows
  
  # Check that either .mod or .posterior supplied
  # If both, use .posterior
  # If neither, error.
  
   .posterior <- as_draws(.mod0)
   
  # if ("bbi_stan_model" %in% class(.mod)) {
  #   .posterior <- read_fit_model(.mod) %>% as_draws()
  # } else {
  #   if ("CmdStanMCMC" %in% class(.mod)) {
  #     # A fit object was used
  #     .posterior <- as_draws(.mod)
  #   }
  #   else {
  #     stop("Must submit either a bbr.tan model or fit object.")
  #   }
  # }

  allPars <- dimnames(.posterior)[[3]]
  
  if (!is.null(.pars)) {
  .pars <- unname(unlist(sapply(.pars,
                               function(x){
                                 allPars[startsWith(allPars, x)]
                               })))
  } else {
    .pars = allPars
  }
  
  .posterior <- .posterior[, , .pars]
  ##  posterior <- posterior[, , unlist(sapply(c(paste("^", pars, "$",sep=""),
  ##                                             paste("^", pars, "\\[",sep="")),
  ##                                           grep,
  ##                                           x = dimnames(posterior)[[3]]))]
  nSamp <- dim(.posterior)[1]
  nChains <- dim(.posterior)[2]
  
  diags <- summarize_draws(.posterior, "rhat", "ess_bulk", "ess_tail")
  
  rhats <- diags$rhat
  names(rhats) <- diags$variable
  rhat_plot <- mcmc_rhat(rhats) + yaxis_text() + myTheme
  
  ess_bulk <- diags$ess_bulk
  names(ess_bulk) <- diags$variable
  ess_bulk_ratios <- ess_bulk / (nSamp * nChains)
  ess_bulk_plot <- mcmc_neff(ess_bulk_ratios) + yaxis_text() + myTheme +
    labs(title = "Bulk ESS ratios")
  
  ess_tail <- diags$ess_tail
  names(ess_tail) <- diags$variable
  ess_tail_ratios <- ess_tail / (nSamp * nChains)
  ess_tail_plot <- mcmc_neff(ess_tail_ratios) + yaxis_text() + myTheme +
    labs(title = "Tail ESS ratios")
  
  mcmc_history_plots <- mcmc_history(.posterior, pars = .pars, 
                                     nParPerPage = .n_per_page, myTheme = myTheme, np = np)
  mcmc_density_chain_plots <- mcmc_density(.posterior, pars = .pars, 
                                           nParPerPage = .n_per_page, byChain = TRUE, 
                                           myTheme = theme(text = element_text(size = 12),
                                                           axis.text = element_text(size = 10)))
  mcmc_density_plots <- mcmc_density(.posterior, pars = .pars, 
                                     nParPerPage = .n_per_page, 
                                     myTheme = theme(text = element_text(size = 12),
                                                     axis.text = element_text(size = 10)))
  
  list(rhat_plot = rhat_plot,
       ess_bulk_plot = ess_bulk_plot,
       ess_tail_plot = ess_tail_plot,
       mcmc_history_plots = mcmc_history_plots,
       mcmc_density_chain_plots = mcmc_density_chain_plots,
       mcmc_density_plots = mcmc_density_plots)
}

