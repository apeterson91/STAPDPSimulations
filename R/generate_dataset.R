#' Generate Dataset
#'
#' Generates clustered nonlinear spatial exposure effects in a linear model
#'
#' @param n number of observations
#' @param prop_cluster vector of cluster proportions
#' @param cluster_functions list of length \code{num_clusters}
#' @param mean_num_distances mean number of simulated distances to include
#' @param prop_exposed proportion of observations with exposure to spatial features
#' @param Z matrix of standard regression coefficients - assumes first column is an intercept
#' @param delta vector of regression coefficients
#' @param sigma residual variability standard deviation
#' @param seed random number generator seed
#'
#' @export
generate_dataset <- function(n = 200,
                             prop_cluster = c(0.5,.5),
                             cluster_functions = list(fun_one =function(x) pweibull(x,shape = 5,scale = .4,lower.tail = FALSE),
                                                      fun_two = function(x) pweibull(x,shape = 5,scale=.4,lower.tail = FALSE)),
                             mean_num_distances = 50,
                             rdist_generator = runif,
                             prop_exposed = .95,
                             Z = cbind(1,rbinom(n = n,1,.5)),
                             delta = c(26,-.5),
                             sigma = 1,
                             seed = NULL
                             ){

  stopifnot(sum(prop_cluster)==1)
  stopifnot(prop_exposed<1 && prop_exposed >0)
  if(is.null(seed))
    seed <- 3431343L
  set.seed(seed)

  has_exp <- rbinom(n,size=1,p=prop_exposed)
  cnt <- rpois(n, mean_num_distances)*has_exp
  distances <- lapply(cnt,function(x) rdist_generator(x))

  K <- length(prop_cluster)
  cmat <- rmultinom(n,size=1,prob = prop_cluster)
  which_cluster <- apply(cmat,2,function(x) which(x==1))
  pmat <- create_pmat(t(cmat))

  exposure <- purrr::map2_dbl(which_cluster,distances,function(x,y) sum(cluster_functions[[x]](y)) )

  y <- Z %*% delta + exposure + rnorm(n,sd=sigma)

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

create_pmat <- function(cmat){

  K <- ncol(cmat)
  mat <- matrix(0,nrow=nrow(cmat),ncol=nrow(cmat))

  for(k in 1:K){
    ics <- which(cmat[,k]==1)
    mat[ics,ics] <- 1
  }
  return(mat)
}
