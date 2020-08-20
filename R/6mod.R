####### MODEL 0

model0 <- lapply(c("LSO", "MOZ", "NAM", "UGA", "ZMB", "ETH", "TZA", "MWI", "ZWE"), function(iso3_current) {
  
  iso3_current <- "ZWE"
  list2env(make_areas_population(iso3_current, naomi_data_path, full = FALSE), globalenv())
  exclude_districts= ""
  
  asfr <- get_asfr_pred_df(iso3_current, area_level = "naomi", areas_long, project = FALSE)
  
  if(filter(mics_key, iso3 == iso3_current)$mics) {
    mics_asfr <- readRDS(here(grep("mics", list.files(paste0("countries/", iso3_current, "/data"), full.names = TRUE), value=TRUE)))
  } else {
    mics_asfr <- NULL
  }
  
  mf <- make_model_frames(iso3_current, population, asfr, mics_asfr, exclude_districts, project=FALSE)
  
  Z <- list()
  Z$Z_spatial <- sparse.model.matrix(~0 + area_id, mf$mf_model)
  Z$Z_age <- sparse.model.matrix(~0 + age_group, mf$mf_model)
  Z$Z_period <- sparse.model.matrix(~0 + period, mf$mf_model)
  
  M_obs <- sparse.model.matrix(~0 + idx, mf$dist$obs) 
  Z$Z_tips <- sparse.model.matrix(~0 + tips_f, mf$dist$obs)
  Z$Z_tips_dhs <- sparse.model.matrix(~0 + tips_f, mf$dist$obs %>% filter(ais_dummy ==0))
  Z$Z_tips_ais <- sparse.model.matrix(~0 + tips_f, mf$dist$obs %>% filter(ais_dummy ==1))
  X_tips_dummy <- model.matrix(~0 + tips_dummy, mf$dist$obs)
  
  ais_join <- mf$dist$obs %>% 
    mutate(col_idx = row_number()) %>%
    select(col_idx, ais_dummy) %>%
    filter(ais_dummy == 1) %>%
    mutate(row_idx = row_number(),
           x=1)
  
  X_extract_ais <- spMatrix(nrow(ais_join), nrow(mf$dist$obs), i=ais_join$row_idx, j=ais_join$col_idx, x=ais_join$x)
  
  dhs_join <- mf$dist$obs %>% 
    mutate(col_idx = row_number()) %>%
    select(col_idx, ais_dummy) %>%
    filter(ais_dummy == 0) %>%
    mutate(row_idx = row_number(),
           x=1)
  
  X_extract_dhs <- spMatrix(nrow(dhs_join), nrow(mf$dist$obs), i=dhs_join$row_idx, j=dhs_join$col_idx, x=dhs_join$x)
  
  if(mf$mics_toggle) {
    M_obs_mics <- sparse.model.matrix(~0 + idx, mf$mics$obs) 
    Z$Z_tips_mics <- sparse.model.matrix(~0 + tips_f, mf$mics$obs)
    X_tips_dummy_mics <- model.matrix(~0 + tips_dummy, mf$mics$obs)
  }
  
  R <- list()
  R$R_spatial <- make_adjacency_matrix(iso3_current, areas_long, boundaries, exclude_districts, level= unique(asfr %>% left_join(areas_long) %>% .$area_level))
  R$R_tips <- make_rw_structure_matrix(ncol(Z$Z_tips), 1, adjust_diagonal = TRUE)
  R$R_age <- make_rw_structure_matrix(ncol(Z$Z_age), 1, adjust_diagonal = TRUE)
  R$R_period <- make_rw_structure_matrix(ncol(Z$Z_period), 2, adjust_diagonal = TRUE)
  
  # dyn.unload(dynlib(here("tmb/fertility_tmb_dev")))
  compile(here("tmb/model_0.cpp"))               # Compile the C++ file
  dyn.load(dynlib(here("tmb/model_0")))
  
  tmb_int <- list()
  
  tmb_int$data <- list(M_obs = M_obs,
                       X_tips_dummy = X_tips_dummy,
                       X_extract_dhs = X_extract_dhs,
                       X_extract_ais = X_extract_ais,
                       Z_tips = Z$Z_tips,
                       Z_tips_dhs = Z$Z_tips_dhs,
                       Z_tips_ais = Z$Z_tips_ais,
                       Z_age = Z$Z_age,
                       Z_period = Z$Z_period,
                       Z_spatial = Z$Z_spatial,
                       Z_interaction1 = sparse.model.matrix(~0 + id.interaction1, mf$mf_model),
                       Z_interaction2 = sparse.model.matrix(~0 + id.interaction2, mf$mf_model),
                       Z_interaction3 = sparse.model.matrix(~0 + id.interaction3, mf$mf_model),
                       R_tips = R$R_tips,
                       R_age = R$R_age,
                       R_period = R$R_period,
                       R_spatial = R$R_spatial,
                       rankdef_R_spatial = 1,
                       log_offset = log(mf$dist$obs$pys),
                       births_obs = mf$dist$obs$births,
                       log_offset_dhs = log(filter(mf$dist$obs, ais_dummy ==0)$pys),
                       births_obs_dhs = filter(mf$dist$obs, ais_dummy ==0)$births,
                       log_offset_ais = log(filter(mf$dist$obs, ais_dummy ==1)$pys),
                       births_obs_ais = filter(mf$dist$obs, ais_dummy ==1)$births,
                       pop = mf$mf_model$population,
                       A_out = mf$out$A_out,
                       mics_toggle = mf$mics_toggle,
                       out_toggle = mf$out_toggle
  )
  
  tmb_int$par <- list(
    beta_0 = 0,
    
    beta_tips_dummy = rep(0, ncol(X_tips_dummy)),
    u_tips = rep(0, ncol(Z$Z_tips)),
    log_prec_rw_tips = 4,
    
    u_age = rep(0, ncol(Z$Z_age)),
    log_prec_rw_age = 4,
    
    u_period = rep(0, ncol(Z$Z_period)),
    log_prec_rw_period = 4,
    
    u_spatial_str = rep(0, ncol(Z$Z_spatial)),
    log_prec_spatial = 0,
    
    eta1 = array(0, c(ncol(Z$Z_period), ncol(Z$Z_age))),
    log_prec_eta1 = 4,
    lag_logit_eta1_phi_age = 0,
    lag_logit_eta1_phi_period = 0,
    
    eta2 = array(0, c(ncol(Z$Z_spatial), ncol(Z$Z_period))),
    log_prec_eta2 = 4,
    lag_logit_eta2_phi_period = 0,
    #
    eta3 = array(0, c(ncol(Z$Z_spatial), ncol(Z$Z_age))),
    log_prec_eta3 = 4,
    lag_logit_eta3_phi_age = 0
  )
  
  tmb_int$random <- c("beta_0", "u_spatial_str", "u_age", "u_period", "beta_tips_dummy", "u_tips", "eta1", "eta2", "eta3")
  
  if(mf$mics_toggle) {
    tmb_int$data <- c(tmb_int$data, "M_obs_mics" = M_obs_mics,
                      "X_tips_dummy_mics" = list(X_tips_dummy_mics),
                      "Z_tips_mics" = Z$Z_tips_mics,
                      "births_obs_mics" = list(mf$mics$obs$births),
                      "log_offset_mics" = list(log(mf$mics$obs$pys)),
                      "A_mics" = mf$mics$A_mics)
    # tmb_int$par <- c(tmb_int$par, 
                     # "beta_tips_dummy_mics" = rep(0, ncol(X_tips_dummy_mics)),
                     # "u_tips_mics" = list(rep(0, ncol(Z$Z_tips_mics))),
                     # "log_prec_rw_tips_mics" = 0)
    # tmb_int$random <- c(tmb_int$random, "beta_tips_dummy_mics")
  }
  
  f <- mcparallel({TMB::MakeADFun(data = tmb_int$data,
                                  parameters = tmb_int$par,
                                  DLL = "model_0",
                                  silent=0,
                                  checkParameterOrder=FALSE)
  })
  
  mccollect(f)
  
  obj <-  MakeADFun(data = tmb_int$data,
                    parameters = tmb_int$par,
                    DLL = "model_0",
                    random = tmb_int$random,
                    hessian = FALSE)
  
  f <- nlminb(obj$par, obj$fn, obj$gr)
  f$par.fixed <- f$par
  f$par.full <- obj$env$last.par
  
  fit <- c(f, obj = list(obj))
  fit$sdreport <- sdreport(fit$obj, fit$par)
  
  class(fit) <- "naomi_fit"  # this is hacky...
  fit <- sample_tmb(fit, random_only=FALSE)
  
  saveRDS(fit, paste0("~/Downloads/model0/", iso3_current, ".rds"))
  
  return(fit)
  
})

####### MODEL 1 #####

model1 <- lapply(c("LSO", "MOZ", "NAM", "UGA", "ZMB", "ETH", "TZA", "MWI", "ZWE"), function(iso3_current) {
  
  list2env(make_areas_population(iso3_current, naomi_data_path, full = FALSE), globalenv())
  exclude_districts= ""
  
  asfr <- get_asfr_pred_df(iso3_current, area_level = "naomi", areas_long, project = FALSE)
  
  if(filter(mics_key, iso3 == iso3_current)$mics) {
    mics_asfr <- readRDS(here(grep("mics", list.files(paste0("countries/", iso3_current, "/data"), full.names = TRUE), value=TRUE)))
  } else {
    mics_asfr <- NULL
  }
  
  mf <- make_model_frames(iso3_current, population, asfr, mics_asfr, exclude_districts, project=FALSE)
  
  Z <- list()
  Z$Z_spatial <- sparse.model.matrix(~0 + area_id, mf$mf_model)
  Z$Z_age <- sparse.model.matrix(~0 + age_group, mf$mf_model)
  Z$Z_period <- sparse.model.matrix(~0 + period, mf$mf_model)
  
  M_obs <- sparse.model.matrix(~0 + idx, mf$dist$obs) 
  Z$Z_tips <- sparse.model.matrix(~0 + tips_f, mf$dist$obs)
  Z$Z_tips_dhs <- sparse.model.matrix(~0 + tips_f, mf$dist$obs %>% filter(ais_dummy ==0))
  Z$Z_tips_ais <- sparse.model.matrix(~0 + tips_f, mf$dist$obs %>% filter(ais_dummy ==1))
  X_tips_dummy <- model.matrix(~0 + tips_dummy, mf$dist$obs)
  
  ais_join <- mf$dist$obs %>% 
    mutate(col_idx = row_number()) %>%
    select(col_idx, ais_dummy) %>%
    filter(ais_dummy == 1) %>%
    mutate(row_idx = row_number(),
           x=1)
  
  X_extract_ais <- spMatrix(nrow(ais_join), nrow(mf$dist$obs), i=ais_join$row_idx, j=ais_join$col_idx, x=ais_join$x)
  
  dhs_join <- mf$dist$obs %>% 
    mutate(col_idx = row_number()) %>%
    select(col_idx, ais_dummy) %>%
    filter(ais_dummy == 0) %>%
    mutate(row_idx = row_number(),
           x=1)
  
  X_extract_dhs <- spMatrix(nrow(dhs_join), nrow(mf$dist$obs), i=dhs_join$row_idx, j=dhs_join$col_idx, x=dhs_join$x)
  
  if(mf$mics_toggle) {
    M_obs_mics <- sparse.model.matrix(~0 + idx, mf$mics$obs) 
    Z$Z_tips_mics <- sparse.model.matrix(~0 + tips_f, mf$mics$obs)
    X_tips_dummy_mics <- model.matrix(~0 + tips_dummy, mf$mics$obs)
  }
  
  R <- list()
  R$R_spatial <- make_adjacency_matrix(iso3_current, areas_long, boundaries, exclude_districts, level= unique(asfr %>% left_join(areas_long) %>% .$area_level))
  R$R_tips <- make_rw_structure_matrix(ncol(Z$Z_tips), 1, adjust_diagonal = TRUE)
  R$R_age <- make_rw_structure_matrix(ncol(Z$Z_age), 1, adjust_diagonal = TRUE)
  R$R_period <- make_rw_structure_matrix(ncol(Z$Z_period), 2, adjust_diagonal = TRUE)
  
  # dyn.unload(dynlib(here("tmb/fertility_tmb_dev")))
  compile(here("tmb/model_1.cpp"))               # Compile the C++ file
  dyn.load(dynlib(here("tmb/model_1")))
  
  tmb_int <- list()
  
  tmb_int$data <- list(M_obs = M_obs,
                       X_tips_dummy = X_tips_dummy,
                       X_extract_dhs = X_extract_dhs,
                       X_extract_ais = X_extract_ais,
                       Z_tips = Z$Z_tips,
                       Z_tips_dhs = Z$Z_tips_dhs,
                       Z_tips_ais = Z$Z_tips_ais,
                       Z_age = Z$Z_age,
                       Z_period = Z$Z_period,
                       Z_spatial = Z$Z_spatial,
                       Z_interaction1 = sparse.model.matrix(~0 + id.interaction1, mf$mf_model),
                       Z_interaction2 = sparse.model.matrix(~0 + id.interaction2, mf$mf_model),
                       Z_interaction3 = sparse.model.matrix(~0 + id.interaction3, mf$mf_model),
                       R_tips = R$R_tips,
                       R_age = R$R_age,
                       R_period = R$R_period,
                       R_spatial = R$R_spatial,
                       rankdef_R_spatial = 1,
                       log_offset = log(mf$dist$obs$pys),
                       births_obs = mf$dist$obs$births,
                       log_offset_dhs = log(filter(mf$dist$obs, ais_dummy ==0)$pys),
                       births_obs_dhs = filter(mf$dist$obs, ais_dummy ==0)$births,
                       log_offset_ais = log(filter(mf$dist$obs, ais_dummy ==1)$pys),
                       births_obs_ais = filter(mf$dist$obs, ais_dummy ==1)$births,
                       pop = mf$mf_model$population,
                       A_out = mf$out$A_out,
                       mics_toggle = mf$mics_toggle,
                       out_toggle = mf$out_toggle
  )
  
  tmb_int$par <- list(
    beta_0 = 0,
    
    beta_tips_dummy = rep(0, ncol(X_tips_dummy)),
    u_tips = rep(0, ncol(Z$Z_tips)),
    log_prec_rw_tips = 4,
    
    u_age = rep(0, ncol(Z$Z_age)),
    log_prec_rw_age = 4,
    
    u_period = rep(0, ncol(Z$Z_period)),
    log_prec_rw_period = 4,
    
    u_spatial_str = rep(0, ncol(Z$Z_spatial)),
    log_prec_spatial = 0,
    
    eta1 = array(0, c(ncol(Z$Z_period), ncol(Z$Z_age))),
    log_prec_eta1 = 4,
    lag_logit_eta1_phi_age = 0,
    lag_logit_eta1_phi_period = 0,
    
    eta2 = array(0, c(ncol(Z$Z_spatial), ncol(Z$Z_period))),
    log_prec_eta2 = 4,
    lag_logit_eta2_phi_period = 0,
    #
    eta3 = array(0, c(ncol(Z$Z_spatial), ncol(Z$Z_age))),
    log_prec_eta3 = 4,
    lag_logit_eta3_phi_age = 0
  )
  
  tmb_int$random <- c("beta_0", "u_spatial_str", "u_age", "u_period", "beta_tips_dummy", "u_tips", "eta1", "eta2", "eta3")
  
  if(mf$mics_toggle) {
    tmb_int$data <- c(tmb_int$data, "M_obs_mics" = M_obs_mics,
                      "X_tips_dummy_mics" = list(X_tips_dummy_mics),
                      "Z_tips_mics" = Z$Z_tips_mics,
                      "births_obs_mics" = list(mf$mics$obs$births),
                      "log_offset_mics" = list(log(mf$mics$obs$pys)),
                      "A_mics" = mf$mics$A_mics)
    # tmb_int$par <- c(tmb_int$par, 
    # "beta_tips_dummy_mics" = rep(0, ncol(X_tips_dummy_mics)),
    # "u_tips_mics" = list(rep(0, ncol(Z$Z_tips_mics))),
    # "log_prec_rw_tips_mics" = 0)
    # tmb_int$random <- c(tmb_int$random, "beta_tips_dummy_mics")
  }
  
  obj <-  MakeADFun(data = tmb_int$data,
                    parameters = tmb_int$par,
                    DLL = "model_1",
                    random = tmb_int$random,
                    hessian = FALSE)
  
  f <- nlminb(obj$par, obj$fn, obj$gr)
  f$par.fixed <- f$par
  f$par.full <- obj$env$last.par
  
  fit <- c(f, obj = list(obj))
  fit$sdreport <- sdreport(fit$obj, fit$par)
  
  class(fit) <- "naomi_fit"  # this is hacky...
  fit <- sample_tmb(fit, random_only=FALSE)
  
  saveRDS(fit, paste0("~/Downloads/model1/", iso3_current, ".rds"))
  
  return(fit)
  
})

###### MODEL 2 ######

model2 <- lapply(c("LSO", "MOZ", "NAM", "UGA", "ZMB", "ETH", "TZA", "MWI", "ZWE"), function(iso3_current) {
  
  list2env(make_areas_population(iso3_current, naomi_data_path, full = FALSE), globalenv())
  exclude_districts= ""
  
  asfr <- get_asfr_pred_df(iso3_current, area_level = "naomi", areas_long, project = FALSE)
  
  if(filter(mics_key, iso3 == iso3_current)$mics) {
    mics_asfr <- readRDS(here(grep("mics", list.files(paste0("countries/", iso3_current, "/data"), full.names = TRUE), value=TRUE)))
  } else {
    mics_asfr <- NULL
  }
  
  mf <- make_model_frames(iso3_current, population, asfr, mics_asfr, exclude_districts, project=FALSE)
  
  Z <- list()
  Z$Z_spatial <- sparse.model.matrix(~0 + area_id, mf$mf_model)
  Z$Z_age <- sparse.model.matrix(~0 + age_group, mf$mf_model)
  Z$Z_period <- sparse.model.matrix(~0 + period, mf$mf_model)
  
  M_obs <- sparse.model.matrix(~0 + idx, mf$dist$obs) 
  Z$Z_tips <- sparse.model.matrix(~0 + tips_f, mf$dist$obs)
  Z$Z_tips_dhs <- sparse.model.matrix(~0 + tips_f, mf$dist$obs %>% filter(ais_dummy ==0))
  Z$Z_tips_ais <- sparse.model.matrix(~0 + tips_f, mf$dist$obs %>% filter(ais_dummy ==1))
  X_tips_dummy <- model.matrix(~0 + tips_dummy, mf$dist$obs %>% filter(ais_dummy == 0))
  
  ais_join <- mf$dist$obs %>% 
    mutate(col_idx = row_number()) %>%
    select(col_idx, ais_dummy) %>%
    filter(ais_dummy == 1) %>%
    mutate(row_idx = row_number(),
           x=1)
  
  X_extract_ais <- spMatrix(nrow(ais_join), nrow(mf$dist$obs), i=ais_join$row_idx, j=ais_join$col_idx, x=ais_join$x)
  
  dhs_join <- mf$dist$obs %>% 
    mutate(col_idx = row_number()) %>%
    select(col_idx, ais_dummy) %>%
    filter(ais_dummy == 0) %>%
    mutate(row_idx = row_number(),
           x=1)
  
  X_extract_dhs <- spMatrix(nrow(dhs_join), nrow(mf$dist$obs), i=dhs_join$row_idx, j=dhs_join$col_idx, x=dhs_join$x)
  
  if(mf$mics_toggle) {
    M_obs_mics <- sparse.model.matrix(~0 + idx, mf$mics$obs) 
    Z$Z_tips_mics <- sparse.model.matrix(~0 + tips_f, mf$mics$obs)
    X_tips_dummy_mics <- model.matrix(~0 + tips_dummy, mf$mics$obs)
  }
  
  R <- list()
  R$R_spatial <- make_adjacency_matrix(iso3_current, areas_long, boundaries, exclude_districts, level= unique(asfr %>% left_join(areas_long) %>% .$area_level))
  R$R_tips <- make_rw_structure_matrix(ncol(Z$Z_tips), 1, adjust_diagonal = TRUE)
  R$R_age <- make_rw_structure_matrix(ncol(Z$Z_age), 1, adjust_diagonal = TRUE)
  R$R_period <- make_rw_structure_matrix(ncol(Z$Z_period), 2, adjust_diagonal = TRUE)
  
  # dyn.unload(dynlib(here("tmb/fertility_tmb_dev")))
  compile(here("tmb/model_2.cpp"))               # Compile the C++ file
  dyn.load(dynlib(here("tmb/model_2")))
  
  tmb_int <- list()
  
  tmb_int$data <- list(M_obs = M_obs,
                       X_tips_dummy = X_tips_dummy,
                       X_extract_dhs = X_extract_dhs,
                       X_extract_ais = X_extract_ais,
                       Z_tips = Z$Z_tips,
                       Z_tips_dhs = Z$Z_tips_dhs,
                       Z_tips_ais = Z$Z_tips_ais,
                       Z_age = Z$Z_age,
                       Z_period = Z$Z_period,
                       Z_spatial = Z$Z_spatial,
                       Z_interaction1 = sparse.model.matrix(~0 + id.interaction1, mf$mf_model),
                       Z_interaction2 = sparse.model.matrix(~0 + id.interaction2, mf$mf_model),
                       Z_interaction3 = sparse.model.matrix(~0 + id.interaction3, mf$mf_model),
                       R_tips = R$R_tips,
                       R_age = R$R_age,
                       R_period = R$R_period,
                       R_spatial = R$R_spatial,
                       rankdef_R_spatial = 1,
                       log_offset = log(mf$dist$obs$pys),
                       births_obs = mf$dist$obs$births,
                       log_offset_dhs = log(filter(mf$dist$obs, ais_dummy ==0)$pys),
                       births_obs_dhs = filter(mf$dist$obs, ais_dummy ==0)$births,
                       log_offset_ais = log(filter(mf$dist$obs, ais_dummy ==1)$pys),
                       births_obs_ais = filter(mf$dist$obs, ais_dummy ==1)$births,
                       pop = mf$mf_model$population,
                       A_out = mf$out$A_out,
                       mics_toggle = mf$mics_toggle,
                       out_toggle = mf$out_toggle
  )
  
  tmb_int$par <- list(
    beta_0 = 0,
    
    beta_tips_dummy = rep(0, ncol(X_tips_dummy)),
    u_tips = rep(0, ncol(Z$Z_tips)),
    log_prec_rw_tips = 4,
    
    u_age = rep(0, ncol(Z$Z_age)),
    log_prec_rw_age = 4,
    
    u_period = rep(0, ncol(Z$Z_period)),
    log_prec_rw_period = 4,
    
    u_spatial_str = rep(0, ncol(Z$Z_spatial)),
    log_prec_spatial = 0,
    
    eta1 = array(0, c(ncol(Z$Z_period), ncol(Z$Z_age))),
    log_prec_eta1 = 4,
    lag_logit_eta1_phi_age = 0,
    lag_logit_eta1_phi_period = 0,
    
    eta2 = array(0, c(ncol(Z$Z_spatial), ncol(Z$Z_period))),
    log_prec_eta2 = 4,
    lag_logit_eta2_phi_period = 0,
    #
    eta3 = array(0, c(ncol(Z$Z_spatial), ncol(Z$Z_age))),
    log_prec_eta3 = 4,
    lag_logit_eta3_phi_age = 0
  )
  
  tmb_int$random <- c("beta_0", "u_spatial_str", "u_age", "u_period", "beta_tips_dummy", "u_tips", "eta1", "eta2", "eta3")
  
  if(mf$mics_toggle) {
    tmb_int$data <- c(tmb_int$data, "M_obs_mics" = M_obs_mics,
                      "X_tips_dummy_mics" = list(X_tips_dummy_mics),
                      "Z_tips_mics" = Z$Z_tips_mics,
                      "births_obs_mics" = list(mf$mics$obs$births),
                      "log_offset_mics" = list(log(mf$mics$obs$pys)),
                      "A_mics" = mf$mics$A_mics)
    # tmb_int$par <- c(tmb_int$par, 
    # "beta_tips_dummy_mics" = rep(0, ncol(X_tips_dummy_mics)),
    # "u_tips_mics" = list(rep(0, ncol(Z$Z_tips_mics))),
    # "log_prec_rw_tips_mics" = 0)
    # tmb_int$random <- c(tmb_int$random, "beta_tips_dummy_mics")
  }
  
  obj <-  MakeADFun(data = tmb_int$data,
                    parameters = tmb_int$par,
                    DLL = "model_2",
                    random = tmb_int$random,
                    hessian = FALSE)
  
  f <- nlminb(obj$par, obj$fn, obj$gr)
  f$par.fixed <- f$par
  f$par.full <- obj$env$last.par
  
  fit <- c(f, obj = list(obj))
  fit$sdreport <- sdreport(fit$obj, fit$par)
  
  class(fit) <- "naomi_fit"  # this is hacky...
  fit <- sample_tmb(fit, random_only=FALSE)
  
  saveRDS(fit, paste0("~/Downloads/model2/", iso3_current, ".rds"))
  
  return(fit)
  
})

#### MODEL 3 ####

model3 <- lapply(c("LSO", "MOZ", "MWI", "ZWE"), function(iso3_current) {
  
  list2env(make_areas_population(iso3_current, naomi_data_path, full = FALSE), globalenv())
  exclude_districts= ""
  
  asfr <- get_asfr_pred_df(iso3_current, area_level = "naomi", areas_long, project = FALSE)
  
  if(filter(mics_key, iso3 == iso3_current)$mics) {
    mics_asfr <- readRDS(here(grep("mics", list.files(paste0("countries/", iso3_current, "/data"), full.names = TRUE), value=TRUE)))
  } else {
    mics_asfr <- NULL
  }
  
  mf <- make_model_frames(iso3_current, population, asfr, mics_asfr, exclude_districts, project=FALSE)
  
  Z <- list()
  Z$Z_spatial <- sparse.model.matrix(~0 + area_id, mf$mf_model)
  Z$Z_age <- sparse.model.matrix(~0 + age_group, mf$mf_model)
  Z$Z_period <- sparse.model.matrix(~0 + period, mf$mf_model)
  
  M_obs <- sparse.model.matrix(~0 + idx, mf$dist$obs) 
  Z$Z_tips <- sparse.model.matrix(~0 + tips_f, mf$dist$obs)
  Z$Z_tips_dhs <- sparse.model.matrix(~0 + tips_f, mf$dist$obs %>% filter(ais_dummy ==0))
  Z$Z_tips_ais <- sparse.model.matrix(~0 + tips_f, mf$dist$obs %>% filter(ais_dummy ==1))
  X_tips_dummy <- model.matrix(~0 + tips_dummy, mf$dist$obs %>% filter(ais_dummy == 0))
  
  ais_join <- mf$dist$obs %>% 
    mutate(col_idx = row_number()) %>%
    select(col_idx, ais_dummy) %>%
    filter(ais_dummy == 1) %>%
    mutate(row_idx = row_number(),
           x=1)
  
  X_extract_ais <- spMatrix(nrow(ais_join), nrow(mf$dist$obs), i=ais_join$row_idx, j=ais_join$col_idx, x=ais_join$x)
  
  dhs_join <- mf$dist$obs %>% 
    mutate(col_idx = row_number()) %>%
    select(col_idx, ais_dummy) %>%
    filter(ais_dummy == 0) %>%
    mutate(row_idx = row_number(),
           x=1)
  
  X_extract_dhs <- spMatrix(nrow(dhs_join), nrow(mf$dist$obs), i=dhs_join$row_idx, j=dhs_join$col_idx, x=dhs_join$x)
  
  if(mf$mics_toggle) {
    M_obs_mics <- sparse.model.matrix(~0 + idx, mf$mics$obs) 
    Z$Z_tips_mics <- sparse.model.matrix(~0 + tips_f, mf$mics$obs)
    X_tips_dummy_mics <- model.matrix(~0 + tips_dummy, mf$mics$obs)
  }
  
  R <- list()
  R$R_spatial <- make_adjacency_matrix(iso3_current, areas_long, boundaries, exclude_districts, level= unique(asfr %>% left_join(areas_long) %>% .$area_level))
  R$R_tips <- make_rw_structure_matrix(ncol(Z$Z_tips), 1, adjust_diagonal = TRUE)
  R$R_age <- make_rw_structure_matrix(ncol(Z$Z_age), 1, adjust_diagonal = TRUE)
  R$R_period <- make_rw_structure_matrix(ncol(Z$Z_period), 2, adjust_diagonal = TRUE)
  
  # dyn.unload(dynlib(here("tmb/fertility_tmb_dev")))
  compile(here("tmb/model_3.cpp"))               # Compile the C++ file
  dyn.load(dynlib(here("tmb/model_3")))
  
  tmb_int <- list()
  
  tmb_int$data <- list(M_obs = M_obs,
                       X_tips_dummy = X_tips_dummy,
                       X_extract_dhs = X_extract_dhs,
                       X_extract_ais = X_extract_ais,
                       Z_tips = Z$Z_tips,
                       Z_tips_dhs = Z$Z_tips_dhs,
                       Z_tips_ais = Z$Z_tips_ais,
                       Z_age = Z$Z_age,
                       Z_period = Z$Z_period,
                       Z_spatial = Z$Z_spatial,
                       Z_interaction1 = sparse.model.matrix(~0 + id.interaction1, mf$mf_model),
                       Z_interaction2 = sparse.model.matrix(~0 + id.interaction2, mf$mf_model),
                       Z_interaction3 = sparse.model.matrix(~0 + id.interaction3, mf$mf_model),
                       R_tips = R$R_tips,
                       R_age = R$R_age,
                       R_period = R$R_period,
                       R_spatial = R$R_spatial,
                       rankdef_R_spatial = 1,
                       log_offset = log(mf$dist$obs$pys),
                       births_obs = mf$dist$obs$births,
                       log_offset_dhs = log(filter(mf$dist$obs, ais_dummy ==0)$pys),
                       births_obs_dhs = filter(mf$dist$obs, ais_dummy ==0)$births,
                       log_offset_ais = log(filter(mf$dist$obs, ais_dummy ==1)$pys),
                       births_obs_ais = filter(mf$dist$obs, ais_dummy ==1)$births,
                       pop = mf$mf_model$population,
                       A_out = mf$out$A_out,
                       mics_toggle = mf$mics_toggle,
                       out_toggle = mf$out_toggle
  )
  
  tmb_int$par <- list(
    beta_0 = 0,
    
    beta_tips_dummy = rep(0, ncol(X_tips_dummy)),
    u_tips = rep(0, ncol(Z$Z_tips)),
    log_prec_rw_tips = 4,
    
    u_age = rep(0, ncol(Z$Z_age)),
    log_prec_rw_age = 4,
    
    u_period = rep(0, ncol(Z$Z_period)),
    log_prec_rw_period = 4,
    
    u_spatial_str = rep(0, ncol(Z$Z_spatial)),
    log_prec_spatial = 0,
    
    eta1 = array(0, c(ncol(Z$Z_period), ncol(Z$Z_age))),
    log_prec_eta1 = 4,
    lag_logit_eta1_phi_age = 0,
    lag_logit_eta1_phi_period = 0,
    
    eta2 = array(0, c(ncol(Z$Z_spatial), ncol(Z$Z_period))),
    log_prec_eta2 = 4,
    lag_logit_eta2_phi_period = 0,
    #
    eta3 = array(0, c(ncol(Z$Z_spatial), ncol(Z$Z_age))),
    log_prec_eta3 = 4,
    lag_logit_eta3_phi_age = 0
  )
  
  tmb_int$random <- c("beta_0", "u_spatial_str", "u_age", "u_period", "beta_tips_dummy", "u_tips", "eta1", "eta2", "eta3")
  
  if(mf$mics_toggle) {
    tmb_int$data <- c(tmb_int$data, "M_obs_mics" = M_obs_mics,
                      "X_tips_dummy_mics" = list(X_tips_dummy_mics),
                      "Z_tips_mics" = Z$Z_tips_mics,
                      "births_obs_mics" = list(mf$mics$obs$births),
                      "log_offset_mics" = list(log(mf$mics$obs$pys)),
                      "A_mics" = mf$mics$A_mics)
    tmb_int$par <- c(tmb_int$par,
    # "beta_tips_dummy_mics" = rep(0, ncol(X_tips_dummy_mics)),
    "u_tips_mics" = list(rep(0, ncol(Z$Z_tips_mics))))
    # "log_prec_rw_tips_mics" = 0)
    # tmb_int$random <- c(tmb_int$random, "beta_tips_dummy_mics")
  }
  
  obj <-  MakeADFun(data = tmb_int$data,
                    parameters = tmb_int$par,
                    DLL = "model_3",
                    random = tmb_int$random,
                    hessian = FALSE)
  
  f <- nlminb(obj$par, obj$fn, obj$gr)
  f$par.fixed <- f$par
  f$par.full <- obj$env$last.par
  
  fit <- c(f, obj = list(obj))
  fit$sdreport <- sdreport(fit$obj, fit$par)
  
  class(fit) <- "naomi_fit"  # this is hacky...
  fit <- sample_tmb(fit, random_only=FALSE)
  
  saveRDS(fit, paste0("~/Downloads/model3/", iso3_current, ".rds"))
  
  return(fit)
  
})

#### MODEL 4 ####

model4 <- lapply(c("LSO", "MOZ", "MWI", "ZWE"), function(iso3_current) {
  
  list2env(make_areas_population(iso3_current, naomi_data_path, full = FALSE), globalenv())
  exclude_districts= ""
  
  asfr <- get_asfr_pred_df(iso3_current, area_level = "naomi", areas_long, project = FALSE)
  
  if(filter(mics_key, iso3 == iso3_current)$mics) {
    mics_asfr <- readRDS(here(grep("mics", list.files(paste0("countries/", iso3_current, "/data"), full.names = TRUE), value=TRUE)))
  } else {
    mics_asfr <- NULL
  }
  
  mf <- make_model_frames(iso3_current, population, asfr, mics_asfr, exclude_districts, project=FALSE)
  
  Z <- list()
  Z$Z_spatial <- sparse.model.matrix(~0 + area_id, mf$mf_model)
  Z$Z_age <- sparse.model.matrix(~0 + age_group, mf$mf_model)
  Z$Z_period <- sparse.model.matrix(~0 + period, mf$mf_model)
  
  M_obs <- sparse.model.matrix(~0 + idx, mf$dist$obs) 
  Z$Z_tips <- sparse.model.matrix(~0 + tips_f, mf$dist$obs)
  Z$Z_tips_dhs <- sparse.model.matrix(~0 + tips_f, mf$dist$obs %>% filter(ais_dummy ==0))
  Z$Z_tips_ais <- sparse.model.matrix(~0 + tips_f, mf$dist$obs %>% filter(ais_dummy ==1))
  X_tips_dummy <- model.matrix(~0 + tips_dummy, mf$dist$obs %>% filter(ais_dummy == 0))
  
  ais_join <- mf$dist$obs %>% 
    mutate(col_idx = row_number()) %>%
    select(col_idx, ais_dummy) %>%
    filter(ais_dummy == 1) %>%
    mutate(row_idx = row_number(),
           x=1)
  
  X_extract_ais <- spMatrix(nrow(ais_join), nrow(mf$dist$obs), i=ais_join$row_idx, j=ais_join$col_idx, x=ais_join$x)
  
  dhs_join <- mf$dist$obs %>% 
    mutate(col_idx = row_number()) %>%
    select(col_idx, ais_dummy) %>%
    filter(ais_dummy == 0) %>%
    mutate(row_idx = row_number(),
           x=1)
  
  X_extract_dhs <- spMatrix(nrow(dhs_join), nrow(mf$dist$obs), i=dhs_join$row_idx, j=dhs_join$col_idx, x=dhs_join$x)
  
  if(mf$mics_toggle) {
    M_obs_mics <- sparse.model.matrix(~0 + idx, mf$mics$obs) 
    Z$Z_tips_mics <- sparse.model.matrix(~0 + tips_f, mf$mics$obs)
    X_tips_dummy_mics <- model.matrix(~0 + tips_dummy, mf$mics$obs)
  }
  
  R <- list()
  R$R_spatial <- make_adjacency_matrix(iso3_current, areas_long, boundaries, exclude_districts, level= unique(asfr %>% left_join(areas_long) %>% .$area_level))
  R$R_tips <- make_rw_structure_matrix(ncol(Z$Z_tips), 1, adjust_diagonal = TRUE)
  R$R_age <- make_rw_structure_matrix(ncol(Z$Z_age), 1, adjust_diagonal = TRUE)
  R$R_period <- make_rw_structure_matrix(ncol(Z$Z_period), 2, adjust_diagonal = TRUE)
  
  # dyn.unload(dynlib(here("tmb/fertility_tmb_dev")))
  compile(here("tmb/model_4.cpp"))               # Compile the C++ file
  dyn.load(dynlib(here("tmb/model_4")))
  
  tmb_int <- list()
  
  tmb_int$data <- list(M_obs = M_obs,
                       X_tips_dummy = X_tips_dummy,
                       X_extract_dhs = X_extract_dhs,
                       X_extract_ais = X_extract_ais,
                       Z_tips = Z$Z_tips,
                       Z_tips_dhs = Z$Z_tips_dhs,
                       Z_tips_ais = Z$Z_tips_ais,
                       Z_age = Z$Z_age,
                       Z_period = Z$Z_period,
                       Z_spatial = Z$Z_spatial,
                       Z_interaction1 = sparse.model.matrix(~0 + id.interaction1, mf$mf_model),
                       Z_interaction2 = sparse.model.matrix(~0 + id.interaction2, mf$mf_model),
                       Z_interaction3 = sparse.model.matrix(~0 + id.interaction3, mf$mf_model),
                       R_tips = R$R_tips,
                       R_age = R$R_age,
                       R_period = R$R_period,
                       R_spatial = R$R_spatial,
                       rankdef_R_spatial = 1,
                       log_offset = log(mf$dist$obs$pys),
                       births_obs = mf$dist$obs$births,
                       log_offset_dhs = log(filter(mf$dist$obs, ais_dummy ==0)$pys),
                       births_obs_dhs = filter(mf$dist$obs, ais_dummy ==0)$births,
                       log_offset_ais = log(filter(mf$dist$obs, ais_dummy ==1)$pys),
                       births_obs_ais = filter(mf$dist$obs, ais_dummy ==1)$births,
                       pop = mf$mf_model$population,
                       A_out = mf$out$A_out,
                       mics_toggle = mf$mics_toggle,
                       out_toggle = mf$out_toggle
  )
  
  tmb_int$par <- list(
    beta_0 = 0,
    
    beta_tips_dummy = rep(0, ncol(X_tips_dummy)),
    u_tips = rep(0, ncol(Z$Z_tips)),
    log_prec_rw_tips = 4,
    
    u_age = rep(0, ncol(Z$Z_age)),
    log_prec_rw_age = 4,
    
    u_period = rep(0, ncol(Z$Z_period)),
    log_prec_rw_period = 4,
    
    u_spatial_str = rep(0, ncol(Z$Z_spatial)),
    log_prec_spatial = 0,
    
    eta1 = array(0, c(ncol(Z$Z_period), ncol(Z$Z_age))),
    log_prec_eta1 = 4,
    lag_logit_eta1_phi_age = 0,
    lag_logit_eta1_phi_period = 0,
    
    eta2 = array(0, c(ncol(Z$Z_spatial), ncol(Z$Z_period))),
    log_prec_eta2 = 4,
    lag_logit_eta2_phi_period = 0,
    #
    eta3 = array(0, c(ncol(Z$Z_spatial), ncol(Z$Z_age))),
    log_prec_eta3 = 4,
    lag_logit_eta3_phi_age = 0
  )
  
  tmb_int$random <- c("beta_0", "u_spatial_str", "u_age", "u_period", "beta_tips_dummy", "u_tips", "eta1", "eta2", "eta3")
  
  if(mf$mics_toggle) {
    tmb_int$data <- c(tmb_int$data, "M_obs_mics" = M_obs_mics,
                      "X_tips_dummy_mics" = list(X_tips_dummy_mics),
                      "Z_tips_mics" = Z$Z_tips_mics,
                      "births_obs_mics" = list(mf$mics$obs$births),
                      "log_offset_mics" = list(log(mf$mics$obs$pys)),
                      "A_mics" = mf$mics$A_mics)
    tmb_int$par <- c(tmb_int$par,
                     # "beta_tips_dummy_mics" = rep(0, ncol(X_tips_dummy_mics)),
                     "u_tips_mics" = list(rep(0, ncol(Z$Z_tips_mics))))
    # "log_prec_rw_tips_mics" = 0)
    # tmb_int$random <- c(tmb_int$random, "beta_tips_dummy_mics")
  }
  
  obj <-  MakeADFun(data = tmb_int$data,
                    parameters = tmb_int$par,
                    DLL = "model_4",
                    random = tmb_int$random,
                    hessian = FALSE)
  
  f <- nlminb(obj$par, obj$fn, obj$gr)
  f$par.fixed <- f$par
  f$par.full <- obj$env$last.par
  
  fit <- c(f, obj = list(obj))
  fit$sdreport <- sdreport(fit$obj, fit$par)
  
  class(fit) <- "naomi_fit"  # this is hacky...
  fit <- sample_tmb(fit, random_only=FALSE)
  
  saveRDS(fit, paste0("~/Downloads/model4/", iso3_current, ".rds"))
  
  return(fit)
  
})

#### MODEL 5 ####

model5 <- lapply(c("NAM", "UGA", "ZMB", "ETH", "MWI", "ZWE"), function(iso3_current) {
  
  list2env(make_areas_population(iso3_current, naomi_data_path, full = FALSE), globalenv())
  exclude_districts= ""
  
  asfr <- get_asfr_pred_df(iso3_current, area_level = "naomi", areas_long, project = FALSE)
  
  if(filter(mics_key, iso3 == iso3_current)$mics) {
    mics_asfr <- readRDS(here(grep("mics", list.files(paste0("countries/", iso3_current, "/data"), full.names = TRUE), value=TRUE)))
  } else {
    mics_asfr <- NULL
  }
  
  mf <- make_model_frames(iso3_current, population, asfr, mics_asfr, exclude_districts, project=FALSE)
  
  Z <- list()
  Z$Z_spatial <- sparse.model.matrix(~0 + area_id, mf$mf_model)
  Z$Z_age <- sparse.model.matrix(~0 + age_group, mf$mf_model)
  Z$Z_period <- sparse.model.matrix(~0 + period, mf$mf_model)
  
  M_obs <- sparse.model.matrix(~0 + idx, mf$dist$obs) 
  Z$Z_tips <- sparse.model.matrix(~0 + tips_f, mf$dist$obs)
  Z$Z_tips_dhs <- sparse.model.matrix(~0 + tips_f, mf$dist$obs %>% filter(ais_dummy ==0))
  Z$Z_tips_ais <- sparse.model.matrix(~0 + tips_f, mf$dist$obs %>% filter(ais_dummy ==1))
  X_tips_dummy <- model.matrix(~0 + tips_dummy, mf$dist$obs %>% filter(ais_dummy == 0))
  
  ais_join <- mf$dist$obs %>% 
    mutate(col_idx = row_number()) %>%
    select(col_idx, ais_dummy) %>%
    filter(ais_dummy == 1) %>%
    mutate(row_idx = row_number(),
           x=1)
  
  X_extract_ais <- spMatrix(nrow(ais_join), nrow(mf$dist$obs), i=ais_join$row_idx, j=ais_join$col_idx, x=ais_join$x)
  
  dhs_join <- mf$dist$obs %>% 
    mutate(col_idx = row_number()) %>%
    select(col_idx, ais_dummy) %>%
    filter(ais_dummy == 0) %>%
    mutate(row_idx = row_number(),
           x=1)
  
  X_extract_dhs <- spMatrix(nrow(dhs_join), nrow(mf$dist$obs), i=dhs_join$row_idx, j=dhs_join$col_idx, x=dhs_join$x)
  
  if(mf$mics_toggle) {
    M_obs_mics <- sparse.model.matrix(~0 + idx, mf$mics$obs) 
    Z$Z_tips_mics <- sparse.model.matrix(~0 + tips_f, mf$mics$obs)
    X_tips_dummy_mics <- model.matrix(~0 + tips_dummy, mf$mics$obs)
  }
  
  R <- list()
  R$R_spatial <- make_adjacency_matrix(iso3_current, areas_long, boundaries, exclude_districts, level= unique(asfr %>% left_join(areas_long) %>% .$area_level))
  R$R_tips <- make_rw_structure_matrix(ncol(Z$Z_tips), 1, adjust_diagonal = TRUE)
  R$R_age <- make_rw_structure_matrix(ncol(Z$Z_age), 1, adjust_diagonal = TRUE)
  R$R_period <- make_rw_structure_matrix(ncol(Z$Z_period), 2, adjust_diagonal = TRUE)
  
  # dyn.unload(dynlib(here("tmb/fertility_tmb_dev")))
  compile(here("tmb/model_5.cpp"))               # Compile the C++ file
  dyn.load(dynlib(here("tmb/model_5")))
  
  tmb_int <- list()
  
  tmb_int$data <- list(M_obs = M_obs,
                       X_tips_dummy = X_tips_dummy,
                       X_extract_dhs = X_extract_dhs,
                       X_extract_ais = X_extract_ais,
                       Z_tips = Z$Z_tips,
                       Z_tips_dhs = Z$Z_tips_dhs,
                       Z_tips_ais = Z$Z_tips_ais,
                       Z_age = Z$Z_age,
                       Z_period = Z$Z_period,
                       Z_spatial = Z$Z_spatial,
                       Z_interaction1 = sparse.model.matrix(~0 + id.interaction1, mf$mf_model),
                       Z_interaction2 = sparse.model.matrix(~0 + id.interaction2, mf$mf_model),
                       Z_interaction3 = sparse.model.matrix(~0 + id.interaction3, mf$mf_model),
                       R_tips = R$R_tips,
                       R_age = R$R_age,
                       R_period = R$R_period,
                       R_spatial = R$R_spatial,
                       rankdef_R_spatial = 1,
                       log_offset = log(mf$dist$obs$pys),
                       births_obs = mf$dist$obs$births,
                       log_offset_dhs = log(filter(mf$dist$obs, ais_dummy ==0)$pys),
                       births_obs_dhs = filter(mf$dist$obs, ais_dummy ==0)$births,
                       log_offset_ais = log(filter(mf$dist$obs, ais_dummy ==1)$pys),
                       births_obs_ais = filter(mf$dist$obs, ais_dummy ==1)$births,
                       pop = mf$mf_model$population,
                       A_out = mf$out$A_out,
                       mics_toggle = mf$mics_toggle,
                       out_toggle = mf$out_toggle
  )
  
  tmb_int$par <- list(
    beta_0 = 0,
    
    beta_tips_dummy = rep(0, ncol(X_tips_dummy)),
    u_tips = rep(0, ncol(Z$Z_tips)),
    u_tips_ais = rep(0, ncol(Z$Z_tips_ais)),
    log_prec_rw_tips = 4,
    
    u_age = rep(0, ncol(Z$Z_age)),
    log_prec_rw_age = 4,
    
    u_period = rep(0, ncol(Z$Z_period)),
    log_prec_rw_period = 4,
    
    u_spatial_str = rep(0, ncol(Z$Z_spatial)),
    log_prec_spatial = 0,
    
    eta1 = array(0, c(ncol(Z$Z_period), ncol(Z$Z_age))),
    log_prec_eta1 = 4,
    lag_logit_eta1_phi_age = 0,
    lag_logit_eta1_phi_period = 0,
    
    eta2 = array(0, c(ncol(Z$Z_spatial), ncol(Z$Z_period))),
    log_prec_eta2 = 4,
    lag_logit_eta2_phi_period = 0,
    #
    eta3 = array(0, c(ncol(Z$Z_spatial), ncol(Z$Z_age))),
    log_prec_eta3 = 4,
    lag_logit_eta3_phi_age = 0
  )
  
  tmb_int$random <- c("beta_0", "u_spatial_str", "u_age", "u_period", "beta_tips_dummy", "u_tips", "eta1", "eta2", "eta3")
  
  if(mf$mics_toggle) {
    tmb_int$data <- c(tmb_int$data, "M_obs_mics" = M_obs_mics,
                      "X_tips_dummy_mics" = list(X_tips_dummy_mics),
                      "Z_tips_mics" = Z$Z_tips_mics,
                      "births_obs_mics" = list(mf$mics$obs$births),
                      "log_offset_mics" = list(log(mf$mics$obs$pys)),
                      "A_mics" = mf$mics$A_mics)
    tmb_int$par <- c(tmb_int$par,
                     # "beta_tips_dummy_mics" = rep(0, ncol(X_tips_dummy_mics)),
                     "u_tips_mics" = list(rep(0, ncol(Z$Z_tips_mics))))
    # "log_prec_rw_tips_mics" = 0)
    # tmb_int$random <- c(tmb_int$random, "beta_tips_dummy_mics")
  }
  
  obj <-  MakeADFun(data = tmb_int$data,
                    parameters = tmb_int$par,
                    DLL = "model_5",
                    random = tmb_int$random,
                    hessian = FALSE)
  
  f <- nlminb(obj$par, obj$fn, obj$gr)
  f$par.fixed <- f$par
  f$par.full <- obj$env$last.par
  
  fit <- c(f, obj = list(obj))
  fit$sdreport <- sdreport(fit$obj, fit$par)
  
  class(fit) <- "naomi_fit"  # this is hacky...
  fit <- sample_tmb(fit, random_only=FALSE)
  
  saveRDS(fit, paste0("~/Downloads/model5/", iso3_current, ".rds"))
  
  return(fit)
  
})

#### MODEL 6 ####

model6 <- lapply(c("LSO", "MOZ", "MWI", "ZWE"), function(iso3_current) {
  
  list2env(make_areas_population(iso3_current, naomi_data_path, full = FALSE), globalenv())
  exclude_districts= ""
  
  asfr <- get_asfr_pred_df(iso3_current, area_level = "naomi", areas_long, project = FALSE)
  
  if(filter(mics_key, iso3 == iso3_current)$mics) {
    mics_asfr <- readRDS(here(grep("mics", list.files(paste0("countries/", iso3_current, "/data"), full.names = TRUE), value=TRUE)))
  } else {
    mics_asfr <- NULL
  }
  
  mf <- make_model_frames(iso3_current, population, asfr, mics_asfr, exclude_districts, project=FALSE)
  
  Z <- list()
  Z$Z_spatial <- sparse.model.matrix(~0 + area_id, mf$mf_model)
  Z$Z_age <- sparse.model.matrix(~0 + age_group, mf$mf_model)
  Z$Z_period <- sparse.model.matrix(~0 + period, mf$mf_model)
  
  M_obs <- sparse.model.matrix(~0 + idx, mf$dist$obs) 
  Z$Z_tips <- sparse.model.matrix(~0 + tips_f, mf$dist$obs)
  Z$Z_tips_dhs <- sparse.model.matrix(~0 + tips_f, mf$dist$obs %>% filter(ais_dummy ==0))
  Z$Z_tips_ais <- sparse.model.matrix(~0 + tips_f, mf$dist$obs %>% filter(ais_dummy ==1))
  X_tips_dummy <- model.matrix(~0 + tips_dummy, mf$dist$obs %>% filter(ais_dummy == 0))
  
  ais_join <- mf$dist$obs %>% 
    mutate(col_idx = row_number()) %>%
    select(col_idx, ais_dummy) %>%
    filter(ais_dummy == 1) %>%
    mutate(row_idx = row_number(),
           x=1)
  
  X_extract_ais <- spMatrix(nrow(ais_join), nrow(mf$dist$obs), i=ais_join$row_idx, j=ais_join$col_idx, x=ais_join$x)
  
  dhs_join <- mf$dist$obs %>% 
    mutate(col_idx = row_number()) %>%
    select(col_idx, ais_dummy) %>%
    filter(ais_dummy == 0) %>%
    mutate(row_idx = row_number(),
           x=1)
  
  X_extract_dhs <- spMatrix(nrow(dhs_join), nrow(mf$dist$obs), i=dhs_join$row_idx, j=dhs_join$col_idx, x=dhs_join$x)
  
  if(mf$mics_toggle) {
    M_obs_mics <- sparse.model.matrix(~0 + idx, mf$mics$obs) 
    Z$Z_tips_mics <- sparse.model.matrix(~0 + tips_f, mf$mics$obs)
    X_tips_dummy_mics <- model.matrix(~0 + tips_dummy, mf$mics$obs)
  }
  
  R <- list()
  R$R_spatial <- make_adjacency_matrix(iso3_current, areas_long, boundaries, exclude_districts, level= unique(asfr %>% left_join(areas_long) %>% .$area_level))
  R$R_tips <- make_rw_structure_matrix(ncol(Z$Z_tips), 1, adjust_diagonal = TRUE)
  R$R_age <- make_rw_structure_matrix(ncol(Z$Z_age), 1, adjust_diagonal = TRUE)
  R$R_period <- make_rw_structure_matrix(ncol(Z$Z_period), 2, adjust_diagonal = TRUE)
  
  # dyn.unload(dynlib(here("tmb/fertility_tmb_dev")))
  compile(here("tmb/model_6.cpp"))               # Compile the C++ file
  dyn.load(dynlib(here("tmb/model_6")))
  
  tmb_int <- list()
  
  tmb_int$data <- list(M_obs = M_obs,
                       X_tips_dummy = X_tips_dummy,
                       X_extract_dhs = X_extract_dhs,
                       X_extract_ais = X_extract_ais,
                       Z_tips = Z$Z_tips,
                       Z_tips_dhs = Z$Z_tips_dhs,
                       Z_tips_ais = Z$Z_tips_ais,
                       Z_age = Z$Z_age,
                       Z_period = Z$Z_period,
                       Z_spatial = Z$Z_spatial,
                       Z_interaction1 = sparse.model.matrix(~0 + id.interaction1, mf$mf_model),
                       Z_interaction2 = sparse.model.matrix(~0 + id.interaction2, mf$mf_model),
                       Z_interaction3 = sparse.model.matrix(~0 + id.interaction3, mf$mf_model),
                       R_tips = R$R_tips,
                       R_age = R$R_age,
                       R_period = R$R_period,
                       R_spatial = R$R_spatial,
                       rankdef_R_spatial = 1,
                       log_offset = log(mf$dist$obs$pys),
                       births_obs = mf$dist$obs$births,
                       log_offset_dhs = log(filter(mf$dist$obs, ais_dummy ==0)$pys),
                       births_obs_dhs = filter(mf$dist$obs, ais_dummy ==0)$births,
                       log_offset_ais = log(filter(mf$dist$obs, ais_dummy ==1)$pys),
                       births_obs_ais = filter(mf$dist$obs, ais_dummy ==1)$births,
                       pop = mf$mf_model$population,
                       A_out = mf$out$A_out,
                       mics_toggle = mf$mics_toggle,
                       out_toggle = mf$out_toggle
  )
  
  tmb_int$par <- list(
    beta_0 = 0,
    
    beta_tips_dummy = rep(0, ncol(X_tips_dummy)),
    u_tips = rep(0, ncol(Z$Z_tips)),
    log_prec_rw_tips = 4,
    
    u_age = rep(0, ncol(Z$Z_age)),
    log_prec_rw_age = 4,
    
    u_period = rep(0, ncol(Z$Z_period)),
    log_prec_rw_period = 4,
    
    u_spatial_str = rep(0, ncol(Z$Z_spatial)),
    log_prec_spatial = 0,
    
    eta1 = array(0, c(ncol(Z$Z_period), ncol(Z$Z_age))),
    log_prec_eta1 = 4,
    lag_logit_eta1_phi_age = 0,
    lag_logit_eta1_phi_period = 0,
    
    eta2 = array(0, c(ncol(Z$Z_spatial), ncol(Z$Z_period))),
    log_prec_eta2 = 4,
    lag_logit_eta2_phi_period = 0,
    #
    eta3 = array(0, c(ncol(Z$Z_spatial), ncol(Z$Z_age))),
    log_prec_eta3 = 4,
    lag_logit_eta3_phi_age = 0
  )
  
  tmb_int$random <- c("beta_0", "u_spatial_str", "u_age", "u_period", "beta_tips_dummy", "u_tips", "eta1", "eta2", "eta3")
  
  if(mf$mics_toggle) {
    tmb_int$data <- c(tmb_int$data, "M_obs_mics" = M_obs_mics,
                      "X_tips_dummy_mics" = list(X_tips_dummy_mics),
                      "Z_tips_mics" = Z$Z_tips_mics,
                      "births_obs_mics" = list(mf$mics$obs$births),
                      "log_offset_mics" = list(log(mf$mics$obs$pys)),
                      "A_mics" = mf$mics$A_mics)
    tmb_int$par <- c(tmb_int$par,
                     "beta_tips_dummy_mics" = rep(0, ncol(X_tips_dummy_mics)),
                     "u_tips_mics" = list(rep(0, ncol(Z$Z_tips_mics))))
    # "log_prec_rw_tips_mics" = 0)
    tmb_int$random <- c(tmb_int$random, "beta_tips_dummy_mics")
  }
  
  obj <-  MakeADFun(data = tmb_int$data,
                    parameters = tmb_int$par,
                    DLL = "model_6",
                    random = tmb_int$random,
                    hessian = FALSE)
  
  f <- nlminb(obj$par, obj$fn, obj$gr)
  f$par.fixed <- f$par
  f$par.full <- obj$env$last.par
  
  fit <- c(f, obj = list(obj))
  fit$sdreport <- sdreport(fit$obj, fit$par)
  
  class(fit) <- "naomi_fit"  # this is hacky...
  fit <- sample_tmb(fit, random_only=FALSE)
  
  saveRDS(fit, paste0("~/Downloads/model6/", iso3_current, ".rds"))
  
  return(fit)
  
})

saveRDS(model6, "~/Downloads/model6.rds")
