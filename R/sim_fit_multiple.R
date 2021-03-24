
#' Simulation Study
#'
#' @param num_sim number of simulations to run
#' @param .... other arguments for \code{\link{generate_dataset}}
#'
#' @export
sim_fit_multiple <- function(num_sim = 10L,
                             K = 50,
                             fix_alpha  = FALSE,
                             iter_max = 5E3,
                             burn_in = 3E3,
                             chains = 4,
                             dataset_function){

  sim_and_loss <- function(sim_id,dataset_function){
    dta <- dataset_function(3434314 + sim_id)
    cat("Simulation #",sim_id)
    cat("\n")
    sink("tempfile.txt")
    fit <- rstapDP::fdp_staplm(y ~ foo_1 + sap(BEF), benvo = dta$benvo,
                               fix_alpha = fix_alpha,
                               iter_max = iter_max,
                               burn_in = burn_in,
                               K = K,
                               chains = chains,
                               seed = sim_id)
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
    irmse <- calculate_irmse(fit,dta$functions)
    out <- dplyr::tibble(sim_id = sim_id,
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
                  irmse_one = irmse$iRMSE[1],
                  irmse_two = irmse$iRMSE[2],
                  pcc_one = pcc[1],
                  pcc_two = pcc[2],
                  rhat = summary(fit)["sigma","Rhat"])
    return(out)
  }
  out <- purrr::map_dfr(1L:num_sim,function(x) sim_and_loss(x,dataset_function))

}


calculate_irmse <- function(fit,fs){

  ddf <- data.frame(Distance=seq(from = 0, to = 1, by = 0.01))
  sobj <- fit$spec$smooth_objs[[1]]
  dgrid <- mgcv::Predict.matrix(sobj,ddf)
  beta_one <- fit$beta[,,1]
  beta_two <- fit$beta[,,2]
  true_fs <- purrr::map2_dfr(names(fs),fs,function(x,y) dplyr::tibble(Distance = ddf$Distance,
                                                                     f_x = y(ddf$Distance),
                                                                     f_name =x))

  f_hats <- dplyr::tibble(Distance = ddf$Distance,
                          f_name = names(fs)[1],
                          fhat_x = apply(tcrossprod(dgrid,beta_one),1,median)
                          ) %>% rbind(.,
    dplyr::tibble(Distance = ddf$Distance,
                  f_name = names(fs[2]),
                  fhat_x = apply(tcrossprod(dgrid,beta_two),1,median)))

  out <- dplyr::inner_join(true_fs,f_hats,by=c("Distance","f_name")) %>%
    dplyr::group_by(f_name,Distance) %>% dplyr::mutate(RMSE = (f_x - fhat_x)^2) %>%
    dplyr::group_by(f_name) %>% dplyr::summarise(iRMSE = sum(RMSE))

  return(out)
}
