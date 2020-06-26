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

naomi_data_path <- "~/Imperial College London/HIV Inference Group - Documents/Analytical datasets/naomi-data"
mics_key <- read.csv(here("countries/mics_data_key.csv"))
# naomi_data_path <- "~/naomi-data"

source(here("R/inputs.R"))
source(here("R/fertility_funs.R"))

iso3_current <-  "ETH"
list2env(make_areas_population(iso3_current, naomi_data_path, full = FALSE), globalenv())
exclude_districts=""

asfr <- get_asfr_pred_df(iso3_current, area_level = 2, areas_long, project = FALSE)

mf <- make_model_frames(iso3_current, population, asfr, mics_asfr = NULL, exclude_districts, project=FALSE)

Z_spatial <- sparse.model.matrix(~0 + area_id, mf$mf_model)
Z_age <- sparse.model.matrix(~0 + age_group, mf$mf_model)
Z_period <- sparse.model.matrix(~0 + period, mf$mf_model)

M_obs <- sparse.model.matrix(~0 + idx, mf$dist$obs) 
Z_tips <- sparse.model.matrix(~0 + tips_f, mf$dist$obs)
X_tips_dummy <- model.matrix(~0 + tips_dummy, mf$dist$obs)

R_spatial <- make_adjacency_matrix(iso3_current, areas_long, boundaries, exclude_districts, level= 2)
R_tips <- make_rw_structure_matrix(ncol(Z_tips), 1, adjust_diagonal = TRUE)
R_age <- make_rw_structure_matrix(ncol(Z_age), 1, adjust_diagonal = TRUE)
R_period <- make_rw_structure_matrix(ncol(Z_period), 2, adjust_diagonal = TRUE)


# dyn.unload(dynlib(here("tmb/fertility_tmb_dev")))
compile(here("tmb/besag.cpp"))               # Compile the C++ file
dyn.load(dynlib(here("tmb/besag")))

tmb_int <- list()

tmb_int$data <- list(M_obs = M_obs,
             X_tips_dummy = X_tips_dummy,
             # X_urban_dummy = X_urban_dummy,
             Z_tips = Z_tips,
             Z_age = Z_age,
             Z_period = Z_period,
             Z_spatial = Z_spatial,
             Z_interaction1 = sparse.model.matrix(~0 + id.interaction1, mf$mf_model),
             Z_interaction2 = sparse.model.matrix(~0 + id.interaction2, mf$mf_model),
             Z_interaction3 = sparse.model.matrix(~0 + id.interaction3, mf$mf_model),
             R_tips = R_tips,
             R_age = R_age,
             R_period = R_period,
             R_spatial = R_spatial,
             log_offset = log(mf$dist$obs$pys),
             births_obs = mf$dist$obs$births,
             # pop = mf$mf_model$population,
             # A_out = mf$out$A_out,
             # A_out_restype = mf$out$A_out_restype,
             mics_toggle = mf$mics_toggle,
             out_toggle = mf$out_toggle
)

tmb_int$par <- list(
  beta_0 = 0,
  
  # beta_tips_dummy = rep(0, ncol(X_tips_dummy)),
  # u_tips = rep(0, ncol(Z_tips)),
  # log_prec_rw_tips = 4,
  
  
  # u_age = rep(0, ncol(Z_age)),
  # log_prec_rw_age = 4,

  u_period = rep(0, ncol(Z_period)),
  log_prec_rw_period = 4,

  
  u_spatial_str = rep(0, ncol(Z_spatial)),
  log_prec_spatial = 0,
  
  # u_spatial_iid = rep(0, ncol(Z_spatial)),
  # logit_spatial_rho = 0,
  
  # eta1 = array(0, c(ncol(Z_period), ncol(Z_age))),
  # log_prec_eta1 = 4,
  # lag_logit_eta1_phi_age = 0,
  # lag_logit_eta1_phi_period = 0,
  
  eta2 = array(0, c(ncol(Z_spatial), ncol(Z_period))),
  log_prec_eta2 = 4,
  lag_logit_eta2_phi_period = 0
  
  # eta3 = array(0, c(ncol(Z_spatial), ncol(Z_age))),
  # log_prec_eta3 = 4,
  # lag_logit_eta3_phi_age = 0
)

# "u_spatial_str", "u_spatial_iid", "eta1" , "eta1" "beta_tips_dummy",, "eta1""eta1", "u_tips", "beta_tips_dummy", , "u_spatial_iid", "eta3""u_age", "u_period",
tmb_int$random <- c("beta_0", "u_spatial_str", "u_period", "eta2")

if(mf$mics_toggle) {
  tmb_int$data <- c(tmb_int$data, "M_obs_mics" = M_obs_mics,
            "X_tips_dummy_mics" = list(X_tips_dummy_mics),
            "Z_tips_mics" = Z_tips_mics,
            "births_obs_mics" = list(mf$mics$obs$births),
            "log_offset_mics" = list(log(mf$mics$obs$pys)),
            "A_mics" = mf$mics$A_mics)
  tmb_int$par <- c(tmb_int$par, "beta_tips_dummy_mics" = rep(0, ncol(X_tips_dummy_mics)))
  tmb_int$random <- c(tmb_int$random, "beta_tips_dummy_mics")
}


f <- mcparallel({TMB::MakeADFun(data = tmb_int$data,
                                parameters = tmb_int$par,
                                DLL = "besag",
                                silent=0,
                                checkParameterOrder=FALSE)
})

mccollect(f)

obj <-  MakeADFun(data = tmb_int$data,
                  parameters = tmb_int$par,
                  DLL = "besag",
                  random = tmb_int$random,
                  hessian = FALSE)

f <- nlminb(obj$par, obj$fn, obj$gr)
f$par.fixed <- f$par
f$par.full <- obj$env$last.par

fit <- c(f, obj = list(obj))
fit$sdreport <- sdreport(fit$obj, fit$par)

class(fit) <- "naomi_fit"  # this is hacky...
fit <- sample_tmb(fit, random_only=FALSE)

qtls1 <- apply(fit$sample$lambda, 1, quantile, c(0.025, 0.5, 0.975))

mf$mf_model %>%
  mutate(lower = qtls1[1,],
         median = qtls1[2,],
         upper = qtls1[3,],
         source = "tmb") %>%
  type.convert() %>%
  left_join(areas_long) %>%
  ggplot(aes(x=period, y=median, group=age_group, color=age_group)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.3) +
  facet_wrap(~area_name)