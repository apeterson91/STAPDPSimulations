#' Generate Weak Dataset
#'
#' Generates clustered nonlinear spatial exposure effects in a linear model where distance distribution and count is correlated with the cluster function
#'
#' @param n number of observations
#' @param prop_cluster vector of cluster proportions
#' @param cluster_functions list of length \code{num_clusters}
#' @param mean_num_distances mean number of simulated distances to include
#' @param prop_exposed proportion of observations with exposure to spatial features
#' @param Z matrix of standard regression coefficients - assumes first column is an intercept
#' @param delta vector of regression coefficients
#' @param seed random number generator seed
#'
#' @export
generate_weak_dataset <- function(n = 100,
                             prop_cluster = c(0.5,.5),
                             cluster_functions = list(fun_one =function(x) pweibull(x,shape = 5,scale = .5,lower.tail = FALSE),
                                                      fun_two = function(x) 0),
                             num_distances = list(Rural=15,Urban=15),
                             rdist_generator = list(Rural = function(x) rbeta(x,1,1),
                                                    Urban = function(x) rbeta(x,1,1)),
                             prop_exposed = .95,
                             Z = cbind(1,rbinom(n = n, 1,.5)),
                             delta = c(26,-.5),
                             seed = NULL
){

  stopifnot(sum(prop_cluster)==1)
  stopifnot(prop_exposed<1 && prop_exposed >0)
  if(is.null(seed))
    seed <- 3431343L
  set.seed(seed)

  K <- length(prop_cluster)
  cmat <- rmultinom(n,size=1,prob = prop_cluster)
  which_cluster <- apply(cmat,2,function(x) which(x==1))
  pmat <- create_pmat(t(cmat))

  has_exp <- rbinom(n = n,size=1,prob=prop_exposed)
  cnt <- purrr::map_dbl(which_cluster,function(x) {
    num <- num_distances[[x]]
    cntr <- sample((num-5):(num+5),size=1,replace=F)
    return(cntr)
  })*has_exp

  distances <- purrr::map2(cnt,which_cluster,function(x,y) rdist_generator[[y]](x))

  exposure <- purrr::map2_dbl(which_cluster,distances,function(x,y) sum(cluster_functions[[x]](y)) )

  y <- Z %*% delta + exposure

  sdf <- dplyr::tibble(id=1:n,
                       y = y)

  colnames(Z) <- paste0("foo_",0:(ncol(Z)-1))
  Z_ <- Z[,2:ncol(Z), drop = F]
  sdf <- cbind(sdf,Z_)

  ddf <- purrr::map2_dfr(1:length(distances),distances,function(x,y) dplyr::tibble(id=x,Distance=y))

  bdf <- rbenvo::benvo(subject_data = sdf,
                       sub_bef_data = list(BEF=ddf),
                       by='id')

  out <- list(benvo=bdf,
              kcat = which_cluster,
              functions = cluster_functions,
              pmat = pmat)

}
