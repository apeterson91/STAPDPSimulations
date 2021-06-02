

CA_sims_effectsize <- function(n = 200,
                               num_sim = 10,
                               functions = list(lowrisk = function(x) .25*pweibull(x,shape = 2.3, scale = 2.3,lower.tail = FALSE),
                                                highrisk = function(x) pweibull(x,shape = 2.3, scale = 2.3,lower.tail = FALSE)),
                               sampled_BEFs = c(15,30,45)
){

  f <- function(x,pos) x %>% dplyr::distinct(cdscode,Urbanicity) %>% dplyr::filter(Urbanicity %in% c("Urban","Sub-Urban"))
  sdf <- readr::read_csv_chunked("~/../../Volumes/Biostat/Brisa/B4SI/atpvyc/CA_datasets_creation/filtered_grouped_child_data/student_df.csv",
                          readr::DataFrameCallback$new(f),chunk_size = 1000)

  smp_n <- floor(n/2)
  sdf <- sdf %>% dplyr::group_by(Urbanicity) %>% dplyr::sample_n(size = smp_n,replace=F) %>%
    dplyr::ungroup()

  cds <- sdf$cdscode
  f <- function(x,pos) {
    x %>% dplyr::filter(calendar_year == 2001,
                        cdscode %in% cds) %>%
    dplyr::select(cdscode,Distance)
  }
  ddf <- readr::read_csv_chunked("~/../../Volumes/Biostat/Brisa/B4SI/atpvyc/CA_datasets_creation/grouped_child_CA/FFR_CA.csv",
                          readr::DataFrameCallback$new(f),chunk_size = 1E3)


  pseudo_coef <- sample(c(0,1), size = n,replace = T)

  ddf_ <- ddf %>%
    dplyr::filter(Distance<=5) %>%
    dplyr::group_by(cdscode) %>%
    dplyr::sample_n(15)

  X <- ddf_ %>%
    right_join(sdf) %>%
    split(.$Urbanicity) %>% map2_dfr(.,functions,function(x,f){
      x %>%
        dplyr::group_by(cdscode) %>%
        dplyr::summarise(Exposure = sum(f(Distance))) %>%
        dplyr::mutate_at(vars(Exposure),function(x) replace_na(x,0))
    })

  sdf_ <- dplyr::tibble(y = 26 + pseudo_coef*-2.2 + X$Exposure + rnorm(n,sd=.5),
                        foo = pseudo_coef,
                        cdscode = ddf_ %>% dplyr::right_join(sdf) %>%
                          dplyr::distinct(cdscode) %>% dplyr::pull(cdscode))

  bdf <- rbenvo::benvo(subject_data = sdf_ %>%
                         dplyr::mutate(foo=pseudo_coef,y=y),
                       sub_bef_data = list(FFR = ddf_))

  fit <- rstapDP::fdp_staplm(y ~ foo + sap(FFR),
                             benvo = bdf,
                             K = 20,
                             iter_max = 2E3,
                             burn_in = 1E3)

}
