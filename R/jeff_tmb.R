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
## devtools::load_all("~/Documents/GitHub/naomi")
library(naomi)
library(here)

naomi_data_path <- "~/Documents/GitHub/naomi-data"
# naomi_data_path <- "~/naomi-data"

source(here("R/inputs.R"))
source(here("R/fertility_funs.R"))

iso3_current <- "ZWE"
exc  <- ""

##sorry..
list2env(make_areas_population(iso3_current, naomi_data_path, full = FALSE), globalenv())

asfr <- get_asfr_pred_df(iso3_current, 2, project = FALSE)

mics_data <- read_mics(iso3_current)
mics_asfr <- Map(calc_asfr_mics, mics_data$wm, y=list(1),
                 by = list(~area_id + survyear + surveyid + survtype),
                 tips = list(c(0:15)),
                 agegr= list(3:10*5),
                 period = list(1995:2019),
                 counts = TRUE,
                 bhdata = mics_data$bh_df) %>%
  bind_rows %>%
  type.convert() %>%
  filter(period <= survyear) %>%
  rename(age_group = agegr)


mf <- make_model_frames(iso3_current, population, asfr, mics_asfr = mics_asfr, exclude_districts = "", project=FALSE)

X_mf <- model.matrix(~1, mf$mf_model)

Z_spatial <- sparse.model.matrix(~0 + area_id, mf$mf_model)
Z_age <- sparse.model.matrix(~0 + age_group, mf$mf_model)
Z_period <- sparse.model.matrix(~0 + period, mf$mf_model)

M_obs <- sparse.model.matrix(~0 + idx, mf$dist$obs) 
Z_tips <- sparse.model.matrix(~0 + tips_f, mf$dist$obs)
X_tips_dummy <- model.matrix(~0 + tips_dummy, mf$dist$obs)

M_obs_mics <- sparse.model.matrix(~0 + idx, mf$mics$obs) 
Z_tips_mics <- sparse.model.matrix(~0 + tips_f, mf$mics$obs)
X_tips_dummy_mics <- model.matrix(~0 + tips_dummy, mf$mics$obs)

R_spatial <- make_adjacency_matrix(iso3_current, areas_long, boundaries, exclude_districts = "", level = 2)
R_tips <- make_rw_structure_matrix(ncol(Z_tips), 1, TRUE)
R_age <- make_rw_structure_matrix(ncol(Z_age), 1, TRUE)
R_period <- make_rw_structure_matrix(ncol(Z_period), 2, TRUE)


dyn.unload(dynlib(here("tmb/fertility_tmb_dev")))
compile(here("tmb/fertility_tmb_dev.cpp"))               # Compile the C++ file
dyn.load(dynlib(here("tmb/fertility_tmb_dev")))

ar1_phi_age <- 0.99
ar1_phi_period <- 0.99

data <- list(X_mf = X_mf,
             M_obs = M_obs,
             M_obs_mics = M_obs_mics,
             X_tips_dummy_mics = X_tips_dummy_mics,
             Z_tips_mics = Z_tips_mics,
             births_obs_mics = mf$mics$obs$births,
             log_offset_mics = log(mf$mics$obs$pys),
             A_mics = mf$mics$A_mics,
             X_tips_dummy = X_tips_dummy,
             Z_tips = Z_tips,
             Z_age = Z_age,
             Z_period = Z_period,
             Z_spatial = Z_spatial,
             # Z_interaction = sparse.model.matrix(~0 + id.interaction, mf$mf_model),
             Z_interaction1 = sparse.model.matrix(~0 + id.interaction1, mf$mf_model),
             # Z_interaction2 = sparse.model.matrix(~0 + id.interaction2, mf$mf_model),
             # Z_interaction3 = sparse.model.matrix(~0 + id.interaction3, mf$mf_model),
             R_tips = R_tips,
             R_age = R_age,
             R_period = R_period,
             R_spatial = R_spatial,
             ar1_phi_age = ar1_phi_age,
             ar1_phi_period = ar1_phi_period,
             log_offset = log(mf$dist$obs$pys),
             births_obs = mf$dist$obs$births,
             pop = mf$mf_model$population,
             A_out = mf$out$A_out
)

par <- list(
  # beta_mf = rep(0, ncol(X_mf)),
  beta_0 = 0,
  beta_tips_dummy = rep(0, ncol(X_tips_dummy)),
  beta_tips_dummy_mics = rep(0, ncol(X_tips_dummy_mics)),
  u_tips = rep(0, ncol(Z_tips)),
  u_age = rep(0, ncol(Z_age)),
  u_period = rep(0, ncol(Z_period)),
  u_spatial_str = rep(0, ncol(Z_spatial)),
  u_spatial_iid = rep(0, ncol(Z_spatial)),
  # eta = array(0, c(ncol(Z_spatial), ncol(Z_age), ncol(Z_period))),
  eta1 = array(0, c(ncol(Z_period), ncol(Z_age))),
  log_sigma_eta1 = log(2.5),
  # eta2 = array(0, c(ncol(Z_spatial), ncol(Z_period))),
  # log_sigma_eta2 = log(2.5),
  # eta3 = array(0, c(ncol(Z_spatial), ncol(Z_age))),
  # log_sigma_eta3 = log(2.5),
  log_sigma_rw_period = log(2.5),
  log_sigma_rw_age = log(2.5),
  log_sigma_rw_tips = log(2.5),
  log_sigma_spatial = log(2.5),
  logit_spatial_rho = 0
)


f <- mcparallel({TMB::MakeADFun(data = data,
                                parameters = par,
                                DLL = "fertility_tmb_dev",
                                silent=0,
                                checkParameterOrder=FALSE)
})

mccollect(f)

obj <-  MakeADFun(data = data,
                  parameters = par,
                  DLL = "fertility_tmb_dev",
                  #  random = c("beta_mf", "beta_tips_dummy", "u_tips", "u_age", "u_period", "u_spatial_str", "u_spatial_iid", "eta1", "eta2", "eta3"),
                  random = c("beta_0", "beta_tips_dummy", "beta_tips_dummy_mics", "u_tips", "u_age", "u_period", "u_spatial_str", "u_spatial_iid", "eta1"),
                  hessian = FALSE)

f <- nlminb(obj$par, obj$fn, obj$gr)
f$par.fixed <- f$par
f$par.full <- obj$env$last.par

fit <- c(f, obj = list(obj))
fit$sdreport <- sdreport(fit$obj, fit$par)

class(fit) <- "naomi_fit"  # this is hacky...
fit <- sample_tmb(fit)

qtls <- apply(fit$sample$lambda, 1, quantile, c(0.025, 0.5, 0.975))

mf$out$mf_out %>%
  mutate(lower = qtls[1,],
         median = qtls[2,],
         upper = qtls[3,],
         source = "tmb") %>%
  type.convert() %>%
  left_join(areas_long) %>%
  filter(area_level == 0) %>%
  # bind_rows(inla_res %>% mutate(source = "inla")) %>%
  ggplot(aes(x=period, y=median, group=age_group)) +
    geom_line(aes(color=age_group)) +
    # geom_point(data=mics_plot %>% bind_rows(asfr_admin1_plot) %>% filter(age_group == "20-24"), aes(y=asfr, group=survtype, color=survtype)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill=age_group), alpha=0.3) +
    facet_wrap(~area_id)

mics_plot <- Map(calc_asfr_mics, mics_data$wm, y=list(1),
                 by = list(~area_id + survyear + surveyid + survtype),
                 tips = list(c(0,15)),
                 agegr= list(3:10*5),
                 period = list(1995:2019),
                 counts = TRUE,
                 bhdata = mics_data$bh_df) %>%
  bind_rows %>%
  type.convert() %>%
  filter(period <= survyear) %>%
  rename(age_group = agegr)
