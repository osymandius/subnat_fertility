library(TMB)
library(tidyverse)
library(abind)
library(sf)
library(spdep)
library(naomi)
library(Matrix)

setwd("~/Documents/GitHub/subnat_fertility/")

asfr <- readRDS("asfr_pred_subnat_15_2020_NEW.rds")[["MWI"]]
boundaries <- readRDS("input_data/area_boundaries.rds")

iso3_codes <- c("LSO", "MOZ", "MWI", "NAM", "SWZ", "TZA", "UGA", "ZMB", "ZWE")

paths <- paste0("~/Documents/GitHub/naomi-data/", iso3_codes, "/data/", tolower(iso3_codes), "_areas.geojson")
area_cols <- c("iso3", "area_id", "area_name", "area_level", "parent_area_id", "naomi_level")

areas_long <- lapply(paths, read_sf) %>% 
  lapply(function(x) {
  iso3_code <- x %>%
    filter(area_level == 0) %>%
    select(area_id) %>%
    unique %>%
    .$area_id
  
  x <- x %>%
    mutate(iso3 = iso3_code) %>%
    st_drop_geometry() %>%
    select(area_cols)
  
  return(x)
}) %>% 
  bind_rows

mf <- crossing(period = factor(1995:2016),
               agegr = factor(c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49")),
               area_id = filter(areas_long, iso3 == "MWI", area_level == 5)$area_id
               ) %>%
  mutate(row_index = factor(row_number()))


obs <- asfr %>%
  filter(!is.na(surveyid)) %>%
  select(area_id, period, agegr, tips, births, pys) %>%
  mutate(agegr = factor(agegr, levels(mf$agegr)),
         period = factor(period))

obs <- obs %>%
  left_join(mf, by = c("area_id", "period", "agegr"))

#' TIPS binary for +/- 5 years
obs <- obs %>%
  mutate(tips_dummy = as.integer(tips > 5),
         tips_f = factor(tips))

X_mf <- model.matrix(~1 + agegr, mf)

#' This has dimensions (number of observations) x (number of rows in model frame (i.e crossing of age x time x space))
#' Many more rows than mf because observations has things we need to adjust for bias (e.g. tips), but not required in model frame. 
#' Column index (idx) has been joined onto obs dataframe from mf. So now each observation is labelled with the appropriate index in the mf
M_mf_obs <- Matrix::sparse.model.matrix(~0 + row_index, obs) 
 
### TIPS RANDOM WALK

Z_tips <- Matrix::sparse.model.matrix(~0 + tips_f, obs)
#' Create precision matrix for RW1
D_tips <- diff(diag(ncol(Z_tips)), differences = 1)
Q_tips <- as(t(D_tips) %*% D_tips, "dgCMatrix")

### AGE RANDOM WALK

Z_age <- sparse.model.matrix(~0 + agegr, obs)
Q_age <- as(INLA:::inla.rw(ncol(Z_age), 1), "dgCMatrix")

### TIME RANDOM WALK

Z_period <- sparse.model.matrix(~0 + period, obs)
Q_period <- as(INLA:::inla.rw(ncol(Z_period), 2), "dgCMatrix")

### TIPS FIXED EFFECT

X_tips_dummy <- model.matrix(~0 + tips_dummy, obs)

### ICAR

Z_spatial <- sparse.model.matrix(~0 + area_id, obs)

sh <- areas_long %>%
  filter(iso3 == "MWI", naomi_level) %>%
  mutate(area_idx = row_number())

#' Neighbor list
nb <- sh %>%
  left_join(boundaries) %>%
  st_as_sf %>%
  as("Spatial") %>%
  spdep::poly2nb() %>%
  `names<-`(sh$area_idx)

adj <- nb2mat(nb, zero.policy=TRUE, style="B")
Q_spatial <- INLA::inla.scale.model(diag(rowSums(adj)) - adj,
                            constr = list(A = matrix(1, 1, nrow(adj)), e = 0))


compile("tmb/fertility_tmb_dev.cpp")               # Compile the C++ file
dyn.load(dynlib("tmb/fertility_tmb_dev"))

data <- list(X_mf = X_mf,
             M_all_observations = M_mf_obs,
             X_tips_dummy = X_tips_dummy,
             Z_tips = Z_tips,
             Z_age = Z_age,
             Z_period = Z_period,
             Z_spatial = Z_spatial,
             Q_tips = Q_tips,
             Q_age = Q_age,
             Q_period = Q_period,
             Q_spatial = Q_spatial,
             log_offset = log(obs$pys),
             births_obs = obs$births)

par <- list(beta_mf = rep(0, ncol(X_mf)),
            beta_tips_dummy = rep(0, ncol(X_tips_dummy)),
            u_tips = rep(0, ncol(Z_tips)),
            u_age = rep(0, ncol(Z_age)),
            u_period = rep(0, ncol(Z_period)),
            u_spatial_str = rep(0, ncol(Z_spatial)),
            u_spatial_iid = rep(0, ncol(Z_spatial)),
            log_sigma_rw_tips = log(2.5),
            log_sigma_rw_age = log(2.5),
            log_sigma_rw_period = log(2.5),
            log_sigma_spatial = log(2.5),
            logit_spatial_rho = 0
            )
             
f <-  MakeADFun(data = data,
                parameters = par,
                DLL = "fertility_tmb_dev",
                random = c("beta_mf", "beta_tips_dummy", "u_tips", "u_age", "u_period", "u_spatial_str", "u_spatial_iid"),
                hessian = TRUE,
                checkParameterOrder=FALSE)

# f$env$tracepar <- TRUE
# f$report()

fit <- nlminb(f$par, f$fn, f$gr)

exp(f$env$last.par) %>%
  split(., names(.))
