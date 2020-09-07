library(TMB)
library(INLA)
library(tidyverse)
library(abind)
library(sf)
library(spdep)
library(Matrix)
library(countrycode)
library(haven)
library(survival)
library(parallel)
library(demogsurv)
devtools::load_all("~/Documents/GitHub/naomi")
# library(naomi)
library(here)

naomi_data_path <- "~/Imperial College London/HIV Inference Group - WP - Documents/Analytical datasets/naomi-data"
mics_key <- read.csv(here("countries/mics_data_key.csv"))
# naomi_data_path <- "~/naomi-data"

source(here("R/inputs.R"))
source(here("R/fertility_funs.R"))

# iso3 <- c("LSO", "MOZ", "NAM", "UGA", "ZMB", "ETH", "TZA", "MWI")
# iso3 <- c("NAM", "UGA", "ZMB", "TZA")

iso3_current <-  c("ZMB", "ZWE", "MOZ", "MWI", "SWZ", "TZA")

lvl_df <- data.frame("iso3" = rep(iso3_current, times=2), "area_level_name" = c(rep("province", times=6), rep("district", times=6)), "area_level_id" = c(1,1,1,1,1,2,2,2,2,5,2,3))

list2env(make_areas_population(iso3_current, naomi_data_path, full = FALSE, return_list = FALSE), globalenv())

# exclude_districts <- areas_wide$area_id[areas_wide$area_id1 == "TZA_1_2"]
# exclude_districts <- c("UGA_3_029", "UGA_3_046")
# exclude_districts <- c("MOZ_2_0107", "MOZ_2_1009")
# exclude_districts <- "MWI_5_07"
exclude_districts= ""

asfr <- Map(function(iso3_current, level) {
  get_asfr_pred_df(iso3_current, area_level = level, areas_long, project = FALSE)
}, iso3_current = filter(lvl_df, area_level_name == "district")$iso3, level = filter(lvl_df, area_level_name == "district")$area_level_id)

mics_asfr <- lapply(iso3_current, function(iso3_current) {
  if(filter(mics_key, iso3 == iso3_current)$mics) {
    readRDS(here(grep("mics", list.files(paste0("countries/", iso3_current, "/data"), full.names = TRUE), value=TRUE)))
  } else {
    mics_asfr <- NULL
  }
})

# names(mics_asfr) <- iso3_current
# 
# asfr[["ZWE"]] <- asfr[["ZWE"]] %>%
#   bind_rows(mics_asfr[["ZWE"]])
# 
# asfr[["SWZ"]] <- asfr[["SWZ"]] %>%
#   bind_rows(mics_asfr[["SWZ"]])
# 
# asfr[["MOZ"]] <- asfr[["MOZ"]] %>%
#   bind_rows(mics_asfr[["MOZ"]])

# mics_dat <- readRDS("input_data/mics_extract.rds")
# # #
# mics_asfr <- Map(calc_asfr_mics, mics_dat$wm[c(4,5)], y=list(1),
#                  by = list(~area_id + survey_id),
#                  tips = list(c(0,15)),
#                  agegr= list(3:10*5),
#                  period = list(1995:2019),
#                  counts = TRUE,
#                  bhdata = mics_dat$bh_df[c(4,5)]) %>%
#   bind_rows %>%
#   type.convert() %>%
#   separate(col=survey_id, into=c(NA, "survyear", NA), sep=c(3,7), remove = FALSE, convert = TRUE) %>%
#   filter(period <= survyear) %>%
#   rename(age_group = agegr)

mf <- make_model_frames(iso3_current, population, asfr, mics_asfr,
                        exclude_districts,
                        project=FALSE,
                        mics_flag = TRUE,
                        level = "naomi")

# saveRDS(mf, here(paste0("countries/", iso3_current, "/mods/", iso3_current, "_mf.rds")))

Z <- list()
Z$Z_spatial <- sparse.model.matrix(~0 + area_id, mf$mf_model)
Z$Z_age <- sparse.model.matrix(~0 + age_group, mf$mf_model)
Z$Z_period <- sparse.model.matrix(~0 + period, mf$mf_model)
Z$Z_country <- sparse.model.matrix(~0 + iso3, mf$mf_model)

M_obs <- sparse.model.matrix(~0 + idx, mf$dist$obs) 
Z$Z_tips <- sparse.model.matrix(~0 + tips_f, mf$dist$obs)
Z$Z_tips_dhs <- sparse.model.matrix(~0 + tips_f, mf$dist$obs %>% filter(ais_dummy ==0))
# Z$Z_tips_ais <- sparse.model.matrix(~0 + tips_f, mf$dist$obs %>% filter(ais_dummy ==1))
# Z_tips[which(mf$dist$obs$survtype != "DHS"), ] <- 0
X_tips_dummy <- model.matrix(~0 + tips_dummy, mf$dist$obs %>% filter(ais_dummy == 0))
# X_urban_dummy <- model.matrix(~0 + urban_dummy, mf$dist$obs)

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
  # X_tips_dummy_mics <- model.matrix(~0 + tips_dummy, mf$mics$obs)
}

R <- list()
R$R_spatial <- make_adjacency_matrix(iso3_current, areas_long, boundaries, exclude_districts, level= "naomi")
# unique(asfr %>% left_join(areas_long) %>% .$area_level)
R$R_tips <- make_rw_structure_matrix(ncol(Z$Z_tips), 1, adjust_diagonal = TRUE)
R$R_tips_mics <- make_rw_structure_matrix(ncol(Z$Z_tips_mics), 1, adjust_diagonal = TRUE)
R$R_age <- make_rw_structure_matrix(ncol(Z$Z_age), 1, adjust_diagonal = TRUE)
R$R_period <- make_rw_structure_matrix(ncol(Z$Z_period), 2, adjust_diagonal = TRUE)
R$R_country <- as(diag(1, length(iso3_current)), "dgTMatrix")

# dyn.unload(dynlib(here("tmb/fertility_tmb_dev")))
compile(here("tmb/multi.cpp"))               # Compile the C++ file
dyn.load(dynlib(here("tmb/multi")))

tmb_int <- list()

tmb_int$data <- list(M_obs = M_obs,
             X_tips_dummy = X_tips_dummy,
             # X_urban_dummy = X_urban_dummy,
             X_extract_dhs = X_extract_dhs,
             X_extract_ais = X_extract_ais,
             Z_tips = Z$Z_tips,
             Z_tips_dhs = Z$Z_tips_dhs,
             # Z_tips_ais = Z$Z_tips_ais,
             Z_age = Z$Z_age,
             Z_period = Z$Z_period,
             Z_spatial = Z$Z_spatial,
             Z_interaction1 = sparse.model.matrix(~0 + id.interaction1, mf$mf_model),
             Z_interaction2 = sparse.model.matrix(~0 + id.interaction2, mf$mf_model),
             Z_interaction3 = sparse.model.matrix(~0 + id.interaction3, mf$mf_model),
             Z_omega1 = sparse.model.matrix(~0 + id.omega1, mf$mf_model),
             Z_omega2 = sparse.model.matrix(~0 + id.omega2, mf$mf_model),
             R_tips = R$R_tips,
             R_age = R$R_age,
             R_period = R$R_period,
             R_spatial = R$R_spatial,
             R_country = R$R_country,
             rankdef_R_spatial = 1,
             # log_offset = log(mf$dist$obs$pys),
             # births_obs = mf$dist$obs$births,
             log_offset_dhs = log(filter(mf$dist$obs, ais_dummy ==0)$pys),
             births_obs_dhs = filter(mf$dist$obs, ais_dummy ==0)$births,
             log_offset_ais = log(filter(mf$dist$obs, ais_dummy ==1)$pys),
             births_obs_ais = filter(mf$dist$obs, ais_dummy ==1)$births,
             pop = mf$mf_model$population,
             A_out = mf$out$A_out,
             mics_toggle = mf$mics_toggle,
             out_toggle = mf$out_toggle
             
             # beta_tips_dummy = 0.05
             # log_prec_rw_tips = 7
             
             # log_prec_eta1 = 0.7717237,
             # lag_logit_eta1_phi_age = 1.8265262,
             # lag_logit_eta1_phi_period  = 4.9531259,
             # log_prec_eta2 = 3.0783250,
             # lag_logit_eta2_phi_period = 1.9670780,
             # log_prec_eta3 = 3.4021674,
             # lag_logit_eta3_phi_age = 2.3139407

)

tmb_int$par <- list(
  beta_0 = 0,

  beta_tips_dummy = rep(0, ncol(X_tips_dummy)),
  # beta_urban_dummy = rep(0, ncol(X_urban_dummy)),
  u_tips = rep(0, ncol(Z$Z_tips)),
  log_prec_rw_tips = 4,

  u_age = rep(0, ncol(Z$Z_age)),
  log_prec_rw_age = 4,
  
  omega1 = array(0, c(ncol(Z$Z_country), ncol(Z$Z_age))),
  log_prec_omega1 = 4,
  lag_logit_omega1_phi_age = 0,
  
  omega2 = array(0, c(ncol(Z$Z_country), ncol(Z$Z_period))),
  log_prec_omega2 = 4,
  lag_logit_omega2_phi_period = 0,

  u_period = rep(0, ncol(Z$Z_period)),
  log_prec_rw_period = 4,

  u_spatial_str = rep(0, ncol(Z$Z_spatial)),
  log_prec_spatial = 0,

  # u_spatial_iid = rep(0, ncol(Z_spatial)),
  # logit_spatial_rho = 0,

  eta1 = array(0, c(ncol(Z$Z_country), ncol(Z$Z_period), ncol(Z$Z_age))),
  log_prec_eta1 = 4,
  lag_logit_eta1_phi_age = 0,
  lag_logit_eta1_phi_period = 0,
  #
  eta2 = array(0, c(ncol(Z$Z_spatial), ncol(Z$Z_period))),
  log_prec_eta2 = 4,
  lag_logit_eta2_phi_period = 0,
  # #
  eta3 = array(0, c(ncol(Z$Z_spatial), ncol(Z$Z_age))),
  log_prec_eta3 = 4,
  lag_logit_eta3_phi_age = 0
)


# "u_spatial_str", "u_spatial_iid", "eta1" , "eta1" "beta_tips_dummy",,"beta_tips_dummy",  "u_tips" "eta1""eta1", "u_tips",  "eta1", "eta2", "eta3""beta_tips_dummy", , "u_spatial_iid", "eta3" , 
tmb_int$random <- c("beta_0", "u_spatial_str", "u_age", "u_period", "beta_tips_dummy", "u_tips", "omega1", "omega2", "eta1", "eta2", "eta3")

if(mf$mics_toggle) {
  tmb_int$data <- c(tmb_int$data, 
            "M_obs_mics" = M_obs_mics,
            # "X_tips_dummy_mics" = list(X_tips_dummy_mics),
            "Z_tips_mics" = Z$Z_tips_mics,
            "R_tips_mics" = R$R_tips_mics,
            "births_obs_mics" = list(mf$mics$obs$births),
            "log_offset_mics" = list(log(mf$mics$obs$pys)),
            "A_mics" = mf$mics$A_mics)
  tmb_int$par <- c(tmb_int$par,
                   "u_tips_mics" = list(rep(0, ncol(Z$Z_tips_mics)))
                   )
}


f <- mcparallel({TMB::MakeADFun(data = tmb_int$data,
                                parameters = tmb_int$par,
                                DLL = "multi",
                                silent=0,
                                checkParameterOrder=FALSE)
})

mccollect(f)

obj <-  MakeADFun(data = tmb_int$data,
                  parameters = tmb_int$par,
                  DLL = "multi",
                  random = tmb_int$random,
                  hessian = FALSE)

f <- nlminb(obj$par, obj$fn, obj$gr)
f$par.fixed <- f$par
f$par.full <- obj$env$last.par

fit <- c(f, obj = list(obj))
fit$sdreport <- sdreport(fit$obj, fit$par)

class(fit) <- "naomi_fit"  # this is hacky...
fit <- sample_tmb(fit, random_only=FALSE)

# saveRDS(fit, paste0("countries/", iso3_current, "/mods/", iso3_current, "_latest_tmb_mod.rds"))



qtls1 <- apply(fit$sample$lambda, 1, quantile, c(0.025, 0.5, 0.975))

asfr_plot <- readRDS(here("countries/NAM/data/NAM_asfr_plot.rds"))
tfr_plot <- readRDS(here("countries/NAM/data/NAM_tfr_plot.rds"))

mf$out$mf_out %>%
# mf$mf_model %>%
  mutate(lower = qtls1[1,],
         median = qtls1[2,],
         upper = qtls1[3,],
         source = "tmb") %>%
  type.convert() %>%
  left_join(areas_long) %>%
  filter(area_level == 1) %>%
  ggplot(aes(x=period, y=median)) +
  geom_line() +
  # geom_point(data=asfr_plot %>% left_join(areas_long) %>% filter(area_level == 1, asfr<2), aes(y=asfr, group=survtype, color=survtype)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=age_group), alpha = 0.3) +
  facet_grid(age_group~area_name) +
  ylim(0,0.5)

# mf$out$mf_out %>%
mf$mf_model %>%
  mutate(lower = qtls1[1,],
         median = qtls1[2,],
         upper = qtls1[3,],
         source = "tmb") %>%
  type.convert() %>%
  left_join(areas_long) %>%
  filter(area_level == 1) %>%
  ggplot(aes(x=period, y=median, group=age_group)) +
  geom_line(aes(color=age_group)) +
  # geom_point(data=asfr_plot %>% left_join(areas_long) %>% filter(survtype == "DHS"), aes(y=asfr, group=survtype, color=survtype)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=age_group), alpha = 0.3) +
  facet_wrap(~area_id)

  
p <- mf$mf_model %>%
    select(period, age_group, area_id) %>%
  cbind(fit$sample$lambda) %>%
  type.convert() %>%
  group_by(area_id, period) %>%
  # summarise_at(.vars = vars(-c(age_group, out_idx)), sum) %>%
  summarise_at(.vars = vars(-c(age_group)), sum) %>%
  ungroup %>%
  mutate_at(.vars = vars(-c(area_id, period)), function(x) 5*x) %>%
  mutate(
    lower = select(., -c(area_id, period)) %>% apply(1, quantile, 0.025),
    median = select(., -c(area_id, period)) %>% apply(1, quantile, 0.5),
    upper = select(., -c(area_id, period)) %>% apply(1, quantile, 0.975)
  )

p %>%
  select(area_id, period, lower, median, upper) %>%
  left_join(areas_long) %>%
  type.convert() %>%
  # filter(area_level == 4) %>%
  ggplot(aes(x=period, y=median)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  geom_line() +
  # geom_point(data=tfr_plot %>% left_join(areas_long) %>% filter(area_level == 1), aes(y=tfr, group=survtype, color=survtype)) +
  facet_wrap(~area_id)

# saveRDS(mod[[2]], "countries/MOZ/mods/MOZ_no18MIS.rds")

tips <- lapply(mods, function(mods) {
  
  int <- data.frame(summary(mods$sdreport))
  
  int %>% 
    mutate(id = rownames(.),
           idx = row_number()) %>%
    filter(idx %in% c(grep("u_tips", .data$id), 
                      grep("sigma_rw_tips", .data$id),
                      grep("beta_tips", .data$id)
    )
    ) %>%
    mutate(log_sigma_rw_tips = Estimate[id=="log_sigma_rw_tips"],
           beta_tips_dummy = Estimate[id == "beta_tips_dummy"]) %>%
    filter(idx > 200) %>%
    mutate(idx = row_number() - 1,
           dummy = ifelse(idx<6, 0, 1),
           beta_tips_dummy = beta_tips_dummy * dummy,
           est = (Estimate*exp(log_sigma_rw_tips)) + beta_tips_dummy)
})

names(tips) <- iso3

tips <- tips %>%
  bind_rows(.id="iso3")


hyper %>%
  pivot_longer(-c(hyper, mod, mics_data)) %>%
  ggplot(aes(y=hyper, x=value)) +
    geom_point(aes(color=mics_data)) +
    facet_wrap(~name)

tfr_plot <- readRDS(here("countries/LSO/data/LSO_tfr_plot.rds"))
  
  
mf$out$mf_out %>%
  cbind(fit$sample$lambda_out) %>%
  # select(-c(population, id.interaction1, id.interaction2, id.interaction3)) %>%
  type.convert() %>%
  group_by(area_id, period) %>%
  summarise_at(.vars = vars(-c(age_group, out_idx)), sum) %>%
  ungroup %>%
  mutate_at(.vars = vars(-c(area_id, period)), function(x) 5*x) %>%
  mutate(
    lower = select(., -c(area_id, period)) %>% apply(1, quantile, 0.025),
    median = select(., -c(area_id, period)) %>% apply(1, quantile, 0.5),
    upper = select(., -c(area_id, period)) %>% apply(1, quantile, 0.975)
         ) %>%
  select(area_id, period, lower, median, upper) %>%
  left_join(areas_long) %>%
  type.convert() %>%
  filter(area_level == 1) %>%
  ggplot(aes(x=period, y=median)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
    geom_line() +
    geom_point(data=tfr_plot %>% left_join(areas_long) %>% filter(area_level == 1), aes(y=tfr, group=survtype, color=survtype)) +
    facet_wrap(~area_name)



hyper <- hyper %>%
  bind_rows(data.frame("x" = mics_leave_out[[1]]$sample$log_prec_rw_age, "source" = "1995-1999", "hyper" = "log_prec_rw_age")) %>%
  bind_rows(data.frame("x" = mics_leave_out[[1]]$sample$log_prec_rw_period, "source" = "1995-1999", "hyper" = "log_prec_rw_period")) %>%
  bind_rows(data.frame("x" = mics_leave_out[[1]]$sample$log_prec_rw_tips, "source" = "1995-1999", "hyper" = "log_prec_rw_tips")) %>%
  bind_rows(data.frame("x" = mics_leave_out[[1]]$sample$beta_tips_dummy, "source" = "1995-1999", "hyper" = "beta_tips")) %>%
  bind_rows(data.frame("x" = mics_leave_out[[1]]$sample$log_prec_eta1, "source" = "1995-1999", "hyper" = "log_prec_eta1")) %>%
  bind_rows(data.frame("x" = mics_leave_out[[1]]$sample$eta1_phi_age, "source" = "1995-1999", "hyper" = "eta1_phi_age")) %>%
  bind_rows(data.frame("x" = mics_leave_out[[1]]$sample$eta1_phi_period, "source" = "1995-1999", "hyper" = "eta1_phi_period")) %>%
  bind_rows(data.frame("x" = mics_leave_out[[1]]$sample$log_prec_spatial, "source" = "1995-1999", "hyper" = "log_prec_spatial")) %>%
  bind_rows(data.frame("x" = mics_leave_out[[1]]$sample$log_prec_eta2, "source" = "1995-1999", "hyper" = "log_prec_eta2")) %>%
  bind_rows(data.frame("x" = mics_leave_out[[1]]$sample$eta2_phi_period, "source" = "1995-1999", "hyper" = "eta2_phi_period")) %>%
  bind_rows(data.frame("x" = mics_leave_out[[1]]$sample$log_prec_eta3, "source" = "1995-1999", "hyper" = "log_prec_eta3")) %>%
  bind_rows(data.frame("x" = mics_leave_out[[1]]$sample$eta3_phi_age, "source" = "1995-1999", "hyper" = "eta3_phi_age"))

# p1, p2, p3, p3a, p4, p5, p6, p7, p8, 
gridExtra::grid.arrange(p1, p2, p3, p3a, p4, p5, p6, p7, p9, p10, p11, p12,  ncol=3)
gridExtra::grid.arrange(p1, p2, p7, p9, p10, p11, p12, ncol=3)
gridExtra::grid.arrange(p1, p7, p8, p11, p12, ncol=3)


