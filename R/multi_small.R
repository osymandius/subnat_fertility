library(TMB)
library(INLA)
library(tidyverse)
library(sf)
library(spdep)
library(Matrix)
library(countrycode)
# devtools::load_all("~/Documents/GitHub/naomi")
library(naomi)
library(here)
library(MASS)

naomi_data_path <- "~/Imperial College London/HIV Inference Group - WP - Documents/Analytical datasets/naomi-data"
mics_key <- read.csv(here("countries/mics_data_key.csv"))
# naomi_data_path <- "~/naomi-data"

source(here("R/inputs.R"))
source(here("R/fertility_funs.R"))

iso3_current <-  sort(c("ZMB", "ZWE", "MOZ", "MWI", "SWZ", "TZA"))

lvl_df <- read.csv(here("input_data/lvl_df.csv")) %>%
  filter(iso3 %in% iso3_current)

list2env(make_areas_population(iso3_current, naomi_data_path, full = FALSE, return_list = FALSE), globalenv())

exclude_districts= ""

asfr <- Map(function(iso3_current, level) {
  get_asfr_pred_df(iso3_current, area_level = level, areas_long, project = FALSE)
}, iso3_current = filter(lvl_df, area_level_name == "province")$iso3, level = 0)
# filter(lvl_df, area_level_name == "province")$area_level_id)



mf <- make_model_frames(iso3_current, population, asfr, mics_asfr,
                        exclude_districts,
                        project=FALSE,
                        mics_flag = FALSE,
                        level = "national")



Z <- list()
Z$Z_country <- sparse.model.matrix(~0 + iso3, mf$mf_model)

M_obs <- sparse.model.matrix(~0 + idx, mf$dist$obs) 

R <- list()
R$R_country <- as(diag(1, length(iso3_current)), "dgTMatrix")

# dyn.unload(dynlib(here("tmb/fertility_tmb_dev")))
compile(here("tmb/multi_small_example.cpp"))               # Compile the C++ file
dyn.load(dynlib(here("tmb/multi_small_example")))

tmb_int <- list()

tmb_int$data <- list(M_obs = M_obs,
             
             Z_country = Z$Z_country,
             R_country = R$R_country,
             log_offset = log(mf$dist$obs$pys),
             births_obs = mf$dist$obs$births,
             pop = mf$mf_model$population
)

tmb_int$par <- list(
  beta_0 = 0,
  u_country = rep(0, ncol(Z$Z_country)),
  log_prec_country = 0
)

tmb_int$random <- c("beta_0", "u_country")



f <- mcparallel({TMB::MakeADFun(data = tmb_int$data,
                                parameters = tmb_int$par,
                                DLL = "multi_small_example",
                                silent=0,
                                checkParameterOrder=FALSE)
})

mccollect(f)

obj <-  MakeADFun(data = tmb_int$data,
                  parameters = tmb_int$par,
                  DLL = "multi_small_example",
                  random = tmb_int$random,
                  hessian = FALSE)

f <- nlminb(obj$par, obj$fn, obj$gr)
f$par.fixed <- f$par
f$par.full <- obj$env$last.par

fit <- c(f, obj = list(obj))
fit$sdreport <- sdreport(fit$obj, fit$par)

class(fit) <- "naomi_fit"  # this is hacky...
fit <- sample_tmb(fit, random_only=FALSE)

fit_prior <- fitdistr(fit$sample$log_prec_country, "normal")

data.frame(x = as.numeric(fit$sample$log_prec_country)) %>%
  ggplot(aes(x=x)) +
    geom_density(aes(color="red"), show.legend=FALSE) +
    geom_density(data = data.frame(x = rnorm(1000, fit_prior$estimate[[1]], fit_prior$estimate[[2]]))) +
  labs(title = "Red = sampled precision. Black = fit to sample")

## Mean = 3.91
# SD = 0.554

### Go down to a single country, using informed prior:

iso3_current <- "ZWE"

list2env(make_areas_population(iso3_current, naomi_data_path, full = FALSE, return_list = FALSE), globalenv())

exclude_districts= ""

asfr <- get_asfr_pred_df(iso3_current, area_level = 0, areas_long, project = FALSE)

mf <- make_model_frames(iso3_current, population, asfr, mics_asfr=NULL,
                        exclude_districts,
                        project=FALSE,
                        mics_flag = FALSE,
                        level = "national")

Z <- list()
Z$Z_country <- sparse.model.matrix(~0 + iso3, mf$mf_model)

# ?????

# M_obs <- sparse.model.matrix(~0 + idx, mf$dist$obs) 
# 
# R <- list()
# R$R_country <- as(diag(1, length(iso3_current)), "dgTMatrix")
# 
# # dyn.unload(dynlib(here("tmb/fertility_tmb_dev")))
# compile(here("tmb/multi_small_example.cpp"))               # Compile the C++ file
# dyn.load(dynlib(here("tmb/multi_small_example")))
# 
# tmb_int <- list()
# 
# tmb_int$data <- list(M_obs = M_obs,
#                      
#                      Z_country = Z$Z_country,
#                      R_country = R$R_country,
#                      log_offset = log(mf$dist$obs$pys),
#                      births_obs = mf$dist$obs$births,
#                      pop = mf$mf_model$population
# )
# 
# tmb_int$par <- list(
#   beta_0 = 0,
#   u_country = rep(0, ncol(Z$Z_country)),
#   log_prec_country = 0
# )
# 
# tmb_int$random <- c("beta_0", "u_country")
# 
# 
# 
# f <- mcparallel({TMB::MakeADFun(data = tmb_int$data,
#                                 parameters = tmb_int$par,
#                                 DLL = "multi_small_example",
#                                 silent=0,
#                                 checkParameterOrder=FALSE)
# })
# 
# mccollect(f)
# 
# obj <-  MakeADFun(data = tmb_int$data,
#                   parameters = tmb_int$par,
#                   DLL = "multi_small_example",
#                   random = tmb_int$random,
#                   hessian = FALSE)
# 
# f <- nlminb(obj$par, obj$fn, obj$gr)
# f$par.fixed <- f$par
# f$par.full <- obj$env$last.par
# 
# fit <- c(f, obj = list(obj))
# fit$sdreport <- sdreport(fit$obj, fit$par)
# 
# class(fit) <- "naomi_fit"  # this is hacky...
# fit <- sample_tmb(fit, random_only=FALSE)