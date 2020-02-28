library(TMB)
library(tidyverse)
library(abind)
library(sf)
library(spdep)
library(Matrix)
devtools::load_all("~/Documents/GitHub/naomi")

setwd("~/Documents/GitHub/subnat_fertility/")

asfr <- readRDS("asfr_pred_subnat_15_2020_NEW.rds")[["ZWE"]]
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

mf <- crossing(period = factor(1995:2015),
               agegr = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49"),
               area_id = filter(areas_long, iso3 == "ZWE", naomi_level)$area_id) %>%
  left_join(read_csv("~/Documents/GitHub/naomi-data/ZWE/data/zwe_population_nso.csv") %>%
              mutate(year = year_labels(calendar_quarter_to_quarter_id(calendar_quarter))) %>%
              filter(sex == "female", year == 2019) %>%
              select(area_id, age_group, population) %>%
              rename(agegr = age_group)) %>%
  mutate(area_id = factor(area_id, levels = filter(areas_long, iso3 == "ZWE", naomi_level)$area_id),
         agegr = factor(agegr, levels = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49"))
         ) %>%
  arrange(period, area_id, agegr) %>%
  mutate(idx = factor(row_number()),
         # id.interaction = group_indices(., agegr, period, area_id)
         id.interaction = factor(group_indices(., agegr, period, area_id)))

str(mf)


obs <- asfr %>%
  filter(!is.na(surveyid)) %>%
  select(area_id, period, agegr, tips, births, pys) %>%
  mutate(agegr = factor(agegr, levels(mf$agegr)),
         period = factor(period)) %>%
  left_join(mf, by = c("area_id", "period", "agegr")) %>%
  mutate(tips_dummy = as.integer(tips > 5),
         tips_f = factor(tips),
         area_id = factor(area_id, levels(mf$area_id))) 

X_mf <- model.matrix(~1 + agegr + period + area_id, mf)

#' This has dimensions (number of observations) x (number of rows in model frame (i.e crossing of age x time x space))
#' Many more rows than mf because observations has things we need to adjust for bias (e.g. tips), but not required in model frame. 
#' Column index (idx) has been joined onto obs dataframe from mf. So now each observation is labelled with the appropriate index in the mf
M_mf_obs <- Matrix::sparse.model.matrix(~0 + idx, obs) 
 
### TIPS RANDOM WALK

Z_tips <- Matrix::sparse.model.matrix(~0 + tips_f, obs)
#' Create precision matrix for RW1
D_tips <- diff(diag(ncol(Z_tips)), differences = 1)
Q_tips <- as(t(D_tips) %*% D_tips, "dgCMatrix")

### AGE RANDOM WALK

Z_age <- sparse.model.matrix(~0 + agegr, obs)
D_age <- diff(diag(ncol(Z_age)), differences = 1)
Q_age <- t(D_age) %*% D_age
diag(Q_age) <- diag(Q_age) + 1E-6
Q_age <- as(Q_age, "dgCMatrix")

### TIME RANDOM WALK

Z_period <- sparse.model.matrix(~0 + period, obs)
D_period <- diff(diag(ncol(Z_period)), differences = 2)
Q_period <- t(D_period) %*% D_period
diag(Q_period) <- diag(Q_period) + 1E-6               
Q_period <- as(Q_period, "dgCMatrix")

interaction_idx <- as.numeric(sort(unique(obs$id.interaction)))

### TIPS FIXED EFFECT

X_tips_dummy <- model.matrix(~0 + tips_dummy, obs)

### ICAR

Z_spatial <- sparse.model.matrix(~0 + area_id, obs)

sh <- areas_long %>%
  filter(iso3 == "ZWE", naomi_level) %>%
  mutate(area_idx = row_number())

#' Neighbor list
nb <- sh %>%
  left_join(boundaries) %>%
  st_as_sf %>%
  as("Spatial") %>%
  spdep::poly2nb() %>%
  `names<-`(sh$area_idx)

adj <- nb2mat(nb, zero.policy=TRUE, style="B")
Q_spatial <- INLA::inla.scale.model(diag(rowSums(adj)) - 0.99*adj,
                            constr = list(A = matrix(1, 1, nrow(adj)), e = 0))

## Outputs
area_merged <-  st_read("~/Documents/GitHub/naomi-data/ZWE/data/zwe_areas.geojson")
areas <- create_areas(area_merged = area_merged)
area_aggregation <- create_area_aggregation(area_merged$area_id[area_merged$naomi_level], areas)

mf_out <- crossing(
  area_id = area_aggregation$area_id,
  agegr = unique(mf$agegr),
  period = unique(mf$period)
) %>%
  mutate(out_idx = row_number())

join_out <- crossing(area_aggregation, 
         agegr = unique(mf$agegr),
         period = unique(mf$period)) %>%
  full_join(mf %>%
              select(area_id, agegr, period, idx), by = c("model_area_id" = "area_id", "agegr", "period")) %>%
  full_join(mf_out) %>%
  mutate(x=1)

A_out <- spMatrix(nrow(mf_out), nrow(mf), join_out$out_idx, as.integer(join_out$idx), join_out$x)

compile("tmb/fertility_tmb_dev.cpp")               # Compile the C++ file
dyn.load(dynlib("tmb/fertility_tmb_dev"))

data <- list(X_mf = X_mf,
             M_all_observations = M_mf_obs,
             X_tips_dummy = X_tips_dummy,
             Z_tips = Z_tips,
             Z_age = Z_age,
             Z_period = Z_period,
             Z_spatial = Z_spatial,
             # Z_interaction = sparse.model.matrix(~0 + id.interaction, obs),
             # interaction_idx = interaction_idx,
             Q_tips = Q_tips,
             Q_age = Q_age,
             Q_period = Q_period,
             Q_spatial = Q_spatial,
             log_offset = log(obs$pys),
             births_obs = obs$births,
             pop = mf$population,
             A_out = A_out
             )

par <- list(beta_mf = rep(0, ncol(X_mf)),
            beta_tips_dummy = rep(0, ncol(X_tips_dummy)),
            u_tips = rep(0, ncol(Z_tips)),
            u_age = rep(0, ncol(Z_age)),
            u_period = rep(0, ncol(Z_period)),
            u_spatial_str = rep(0, ncol(Z_spatial)),
            u_spatial_iid = rep(0, ncol(Z_spatial)),
            # eta = array(0, c(ncol(Z_spatial), ncol(Z_age), ncol(Z_period))),
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
                hessian = FALSE,
                checkParameterOrder=FALSE)

# f$env$tracepar <- TRUE
# f$report()

fit <- nlminb(f$par, f$fn, f$gr)
rep <- sdreport(f)
tmb_res <- summary(rep)

exp(f$env$last.par) %>%
  split(., names(.))

mf %>% 
  cbind(data.frame(val = tmb_res[,1][rownames(tmb_res) == "omega"])) %>%
  ggplot(aes(x=period, y=val, group=agegr, color=agegr)) +
  geom_line() +
  facet_wrap(~area_id)

obs %>%
  type.convert %>%
  left_join(
    data.frame(id.interaction = 1:147, val = tmb_res[,1][rownames(tmb_res) == "eta"])
  ) %>%
  ggplot(aes(x=period, y=val, color=agegr, group=agegr)) +
  geom_line()

