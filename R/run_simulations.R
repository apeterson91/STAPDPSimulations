#' Runs Distance Distribution Simulations using data in Data/
#'
#' @param n_sim number of simulations to run
#' @param fix_alpha whether to fix_alpha in simulations
#' @param iter_max maximum number of iterations to run Gibbs Sampler
#' @param burn_in number of iterations to burn_in
#' @param K Number of components to truncate DP to
#' @param chains number of MCMC chains to run per simulated dataset
#' @export
run_dd_simulations <- function(n_sim = 25,
                            fix_alpha = FALSE,
                            iter_max = 4E3,
                            burn_in = 2E3,
                            K = 20L,
                            chains = 1){


  dta_names <- stringr::str_replace(list.files("Data/"),".rda","")
  dta_names <- dta_names[stringr::str_detect(dta_names,"ef",negate = TRUE)]
  mdta <- stringr::str_split_fixed(dta_names,"_", n = 3)
  get_lbls <- function(x){
    out <- dplyr::case_when(x == "cskew" ~ "CA Skew",
                      (x == "lskew" | x== "skew") ~ "Light Skew",
                      x == "u" ~ "Uniform")
    return(out)
  }

  ddst_one <- get_lbls(mdta[,1])
  ddst_two <- get_lbls(mdta[,2])
  num <- dplyr::case_when(mdta[,3] == "ff" ~ 15,
                   mdta[,3] == "tt" ~ 30,
                   mdta[,3] == "ffff" ~ 45)

  dataset_ix <- 1:length(dta_names)
  sim_ix <- 1:n_sim
  gd <- expand.grid(dataset_ix=dataset_ix,
              sim_ix = sim_ix)

  simdf <- purrr::map2_dfr(gd$dataset_ix,gd$sim_ix,function(data_ix,sim_ix){
    dta <- get(dta_names[data_ix])
    set.seed(sim_ix + data_ix + 34134L)
    dta$benvo$subject_data$y <- dta$benvo$subject_data$y + rnorm(length(dta$benvo$subject_data$y),sd=1)
    cat("Distance Dist: ",ddst_one[data_ix],", ", ddst_two[data_ix], "# BEF: ",num[data_ix])
    cat("; Simulation #",sim_ix)
    cat("\n")

    sink("tempfile.txt")
    fit <- rstapDP::fdp_staplm(y ~ foo_1 + sap(BEF), benvo = dta$benvo,
                               fix_alpha = fix_alpha,
                               iter_max = iter_max,
                               burn_in = burn_in,
                               K = K,
                               tau_a = 1,tau_b = 1,
                               chains = chains,
                               seed = 34134L + data_ix + sim_ix)
    sink()
    file.remove("tempfile.txt")
    nums <- apply(fit$cmat,1,function(x) length(unique(x)))
    lnum <- quantile(nums,0.025)
    mnum <- median(nums)
    unum <- quantile(nums,0.975)
    top_five <- sort(apply(fit$probs,2,median),decreasing=TRUE)[1:5]
    names(top_five) <- paste0("pi_",1:5)
    md <- rstapDP::green_loss(fit,truth=dta$pmat)
    loss <- md$loss
    md <- md$mode
    xt <- xtabs(~dta$kcat + md)
    pcc <- apply(xt,1,max)/rowSums(xt)
    out <- dplyr::tibble(sim_id = sim_ix,
                         Distance_Distribution_One = ddst_one[data_ix],
                         Distance_Distribution_Two = ddst_two[data_ix],
                         Num_BEF = num[data_ix],
                         mn_loss = mean(loss),
                         best_loss = min(loss),
                         lnum = lnum,
                         mnum = mnum,
                         unum  = unum,
                         pi_one = top_five[1],
                         pi_two = top_five[2],
                         pi_three = top_five[3],
                         pi_four = top_five[4],
                         pi_five = top_five[5],
                         pcc_one = pcc[1],
                         pcc_two = pcc[2],
                         rhat_sigma = summary(fit)["sigma","Rhat"])
    return(out)


  })

  return(simdf)


}

#' Runs Effect Size Simulations using data in Data/
#'
#' @param n_sim number of simulations to run
#' @param fix_alpha whether to fix_alpha in simulations
#' @param iter_max maximum number of iterations to run Gibbs Sampler
#' @param burn_in number of iterations to burn_in
#' @param K Number of components to truncate DP to
#' @param chains number of MCMC chains to run per simulated dataset
#' @export
run_es_simulations <- function(n_sim=25,
                               fix_alpha = FALSE,
                               iter_max = 4E3,
                               burn_in = 2E3,
                               K = 20L,
                               chains = 1){


  dta_names <- stringr::str_replace(list.files("Data/"),".rda","")
  dta_names <- dta_names[stringr::str_detect(dta_names,"ef")]

  effect_size <- as.integer(stringr::str_extract(dta_names,"[0-9][0-9]?"))
  dataset_ix <- 1:length(dta_names)
  sim_ix <- 1:n_sim
  gd <- expand.grid(dataset_ix = dataset_ix,
                    sim_ix = sim_ix)

  simdf <- purrr::map2_dfr(gd$dataset_ix,gd$sim_ix,function(data_ix,sim_ix){
    dta <- get(dta_names[data_ix])
    set.seed(sim_ix + data_ix + 34134L)
    dta$benvo$subject_data$y <- dta$benvo$subject_data$y + rnorm(length(dta$benvo$subject_data$y),sd=1)
    cat("Effect Size: ",effect_size[data_ix])
    cat("; Simulation #",sim_ix)
    cat("\n")

    sink("tempfile.txt")
    fit <- rstapDP::fdp_staplm(y ~ foo_1 + sap(BEF), benvo = dta$benvo,
                               fix_alpha = fix_alpha,
                               iter_max = iter_max,
                               burn_in = burn_in,
                               K = K,
                               tau_a = 1,tau_b = 1,
                               chains = chains,
                               seed = 34134L + data_ix + sim_ix)
    sink()
    file.remove("tempfile.txt")
    nums <- apply(fit$cmat,1,function(x) length(unique(x)))
    lnum <- quantile(nums,0.025)
    mnum <- median(nums)
    unum <- quantile(nums,0.975)
    top_five <- sort(apply(fit$probs,2,median),decreasing=TRUE)[1:5]
    names(top_five) <- paste0("pi_",1:5)
    md <- rstapDP::green_loss(fit,truth=dta$pmat)
    loss <- md$loss
    md <- md$mode
    xt <- xtabs(~dta$kcat + md)
    out <- dplyr::tibble(sim_id = sim_ix,
                         EffectSize = effect_size[data_ix],
                         mn_loss = mean(loss),
                         best_loss = min(loss),
                         lnum = lnum,
                         mnum = mnum,
                         unum  = unum,
                         pi_one = top_five[1],
                         pi_two = top_five[2],
                         pi_three = top_five[3],
                         pi_four = top_five[4],
                         pi_five = top_five[5],
                         rhat_sigma = summary(fit)["sigma","Rhat"])
    return(out)


  })

  return(simdf)

}
