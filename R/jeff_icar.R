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
# devtools::load_all("~/Documents/GitHub/naomi")
library(naomi)
library(here)

naomi_data_path <- "~/Imperial College London/HIV Inference Group - Documents/Analytical datasets/naomi-data"

source(here("R/inputs.R"))
source(here("R/fertility_funs.R"))

iso3_current <-  "ZWE"
list2env(make_areas_population(iso3_current, naomi_data_path, full = FALSE), globalenv())
exclude_districts=""

asfr <- get_asfr_pred_df(iso3_current, area_level = 1, areas_long, project = FALSE)
mics_asfr <- NULL

mf <- make_model_frames(iso3_current, population, asfr, mics_asfr, exclude_districts, project=FALSE)

Z_spatial <- sparse.model.matrix(~0 + area_id, mf$mf_model)
Z_age <- sparse.model.matrix(~0 + age_group, mf$mf_model)
Z_period <- sparse.model.matrix(~0 + period, mf$mf_model)

M_obs <- sparse.model.matrix(~0 + idx, mf$dist$obs) 
Z_tips <- sparse.model.matrix(~0 + tips_f, mf$dist$obs)
X_tips_dummy <- model.matrix(~0 + tips_dummy, mf$dist$obs)

R_spatial <- make_adjacency_matrix(iso3_current, areas_long, boundaries, exclude_districts, level= 1)
R_tips <- make_rw_structure_matrix(ncol(Z_tips), 1, TRUE)
R_age <- make_rw_structure_matrix(ncol(Z_age), 1, TRUE)
R_period <- make_rw_structure_matrix(ncol(Z_period), 2, TRUE)

compile(here("tmb/jeff_icar.cpp"))               # Compile the C++ file
dyn.load(dynlib(here("tmb/jeff_icar")))

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
  
  u_age = rep(0, ncol(Z_age)),
  log_prec_rw_age = 4,

  u_spatial_str = rep(0, ncol(Z_spatial)),
  log_prec_spatial = 0
)

tmb_int$random <- c("beta_0", "u_spatial_str", "u_age")

f <- mcparallel({TMB::MakeADFun(data = tmb_int$data,
                                parameters = tmb_int$par,
                                DLL = "jeff_icar",
                                silent=0,
                                checkParameterOrder=FALSE)
})

mccollect(f)

obj <-  MakeADFun(data = tmb_int$data,
                  parameters = tmb_int$par,
                  DLL = "jeff_icar",
                  random = tmb_int$random,
                  hessian = FALSE)

f <- nlminb(obj$par, obj$fn, obj$gr)
f$par.fixed <- f$par
f$par.full <- obj$env$last.par

fit <- c(f, obj = list(obj))
fit$sdreport <- sdreport(fit$obj, fit$par)

class(fit) <- "naomi_fit"  # this is hacky...
fit <- sample_tmb(fit, random_only=FALSE)

zwe.m <- inla(
  births ~ f(id.district, model = "besag", graph = here("countries/ZWE/adj/ZWE_admin1.adj"), scale.model = TRUE) +
    f(id.age_group, model="rw1"),
  family="xpoisson", data=asfr, E=pys, 
  control.family=list(link='log'),
  control.predictor=list(compute=TRUE, link=1),
  control.inla = list(strategy = "gaussian", int.strategy = "eb"),
  control.compute=list(config = TRUE, dic= FALSE, cpo=FALSE),
  verbose=TRUE)

inla_log_prec_rw_age <- data.frame(zwe.m$internal.marginals.hyperpar$`Log precision for id.age_group`)
inla_log_prec_spatial <- data.frame(zwe.m$internal.marginals.hyperpar$`Log precision for id.district`)

p1 <- 
  bind_rows(data.frame("x" = fit$sample$log_prec_rw_age, "source" = "tmb")) %>%
  ggplot(aes(x=x)) +
  geom_density(aes(group=source, fill=source), alpha=0.7) +
  geom_point(data=inla_log_prec_rw_age, aes(y=y))+
  labs(title = "log_prec_rw_age")

p7 <- 
  bind_rows(data.frame("x" = fit$sample$log_prec_spatial, "source" = "tmb")) %>%
  ggplot(aes(x=x)) +
  geom_density(aes(group=source, fill=source), alpha=0.7) +
  geom_point(data=inla_log_prec_spatial, aes(y=y))+
  labs(title = "log_prec_spatial")

gridExtra::grid.arrange(p1, p7, ncol=2)