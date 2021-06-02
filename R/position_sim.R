


position_sim <- function(n = 200,
                         num_sims = 15,
                         effect_functions = list(lowrisk = function(x) .25*pweibull(x,shape = 5, scale = .5,lower.tail = FALSE),
                                          highrisk = function(x) pweibull(x,shape = 5, scale = .5,lower.tail = FALSE)),
                         num_BEFs = c(15,30,45)){

  subj <- cbind(x=runif(200,0,2),
                y=runif(200,0,2))

  bef <- cbind(x=runif(30,0,1),
               y=runif(30,0,1))

  kcat <- sample(c(0,1),size=200,replace=T)

  sdf <- dplyr::tibble(subj_ix = 1:200,
                       Cluster = kcat,
                       foo = sample(c(0,1),size=200,replace=T))

  dmat <- fields::rdist(subj,bef)
  colnames(dmat) <- stringr::str_c("BEF_",1:ncol(dmat))
  ddf <- dplyr::as_tibble(dmat) %>%
    dplyr::mutate(subj_ix = 1:n()) %>%
    tidyr::gather(everything(),-subj_ix,key = "BEF",value="Distance") %>%
    dplyr::filter(Distance<=1)


  mdf <- ddf %>%
    dplyr::left_join(sdf %>% dplyr::select(-foo)) %>%
    split(.$Cluster) %>%
    purrr::map2_dfr(.,effect_functions,function(df,f){
      df %>%
        dplyr::group_by(subj_ix) %>%
        dplyr::summarise(Exposure = sum(f(Distance)))
    }) %>% dplyr::arrange(subj_ix) %>%
    dplyr::right_join(sdf) %>%
    dplyr::mutate(y = 26 + -2.2*foo + replace_na(Exposure,0) + rnorm(200))


  bdf <- rbenvo::benvo(subject_data = mdf,sub_bef_data = list(FFR=ddf %>% filter(!is.na(Distance))))

  fit <- rstapDP::fdp_staplm(y ~ foo + sap(FFR),
                             benvo = bdf)

  fit2 <- rsstap::sstap_glm(y ~ foo + sap(FFR),benvo=bdf)

}
