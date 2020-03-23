library(TMB)
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
library(naomi)

library(here)

naomi_data_path <- "~/naomi-data"
## naomi_data_path <- "~/Documents/GitHub/naomi-data/"

source(here("R/inputs.R"))
source(here("R/fertility_funs.R"))

# lapply(list.files("~/Documents/GitHub/naomi/R", full.names = TRUE), source)

iso3_current <- "ZWE"

##sorry..
list2env(make_areas_population("ZWE", naomi_data_path, full = FALSE), globalenv())

asfr <- get_asfr_pred_df("ZWE", 2, project = FALSE)

mics_data <- read_mics(iso3_current)
mics_asfr <- Map(calc_asfr_mics, mics_data$wm, y=list(1),
                 by = list(~area_id + survyear + surveyid + survtype),
                 tips = list(c(0:5)),
                 agegr= list(3:10*5),
                 period = list(1995:2019),
                 counts = TRUE,
                 bhdata = mics_data$bh_df) %>%
  bind_rows %>%
  type.convert() %>%
  filter(period <= survyear) %>%
  rename(age_group = agegr)

population <- population %>%
  filter(period == min(period))

area_merged <- st_read(file.path(naomi_data_path, "ZWE/data/zwe_areas.geojson"))
areas <- create_areas(area_merged = area_merged)
area_aggregation <- create_area_aggregation(area_merged$area_id[area_merged$naomi_level], areas)


## Make model frame.
mf <- crossing(period = factor(1995:2018),
               age_group = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49"),
               area_id = filter(areas_long, iso3 == iso3_current, area_level == 2)$area_id) %>%
  left_join(population %>%
              filter(sex == "female") %>%
              select(area_id, age_group, population)
  ) %>%
  mutate(area_id = factor(area_id, levels = filter(areas_long, iso3 == iso3_current,  area_level == 2)$area_id),
         age_group = factor(age_group, levels = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49"))
         ) %>%
  arrange(period, area_id, age_group) %>%
  mutate(idx = factor(row_number()),
         id.interaction_3d = factor(group_indices(., age_group, period, area_id)),
         # id.interaction_age_time = factor(group_indices(., age_group, period)),
         id.interaction1 = factor(group_indices(., age_group, period)),
         id.interaction2 = factor(group_indices(., period, area_id)),
         id.interaction3 = factor(group_indices(., age_group, area_id)),
  )
# 
# ## Existing incorrect method that produces 63 replicates in the obs dataframe for each national row of data.
# obs <- asfr %>%
#   bind_rows(mics_asfr) %>%
#   type.convert() %>%
#   filter(!is.na(surveyid)) %>%
#   select(area_id, period, age_group, tips, births, pys) %>%
#   left_join(area_aggregation) %>%
#   mutate(idx_obs = row_number()) %>%
#   left_join(mf %>% type.convert, by = c("model_area_id" = "area_id", "period", "age_group")) %>%
#   # arrange(period, area_id, age_group) %>%
#   mutate(tips_dummy = as.integer(tips > 5),
#          tips_f = factor(tips),
#          model_area_id = factor(model_area_id, levels(mf$area_id)),
#          age_group = factor(age_group, levels(mf$age_group)),
#          period = factor(period),
#          id.interaction = factor(id.interaction, levels(mf$id.interaction)),
#          # id.interaction1 = factor(group_indices(., age_group, period)),
#          # id.interaction2 = factor(group_indices(., period, model_area_id)),
#          # id.interaction3 = factor(group_indices(., age_group, model_area_id)),
#          x=1
#          )
# 


obs <- asfr %>%
  # bind_rows(mics_asfr) %>%
  mutate(period = factor(period, levels(mf$period))) %>%
  filter(!is.na(surveyid)) %>%
  select(area_id, period, age_group, tips, births, pys) %>%
  left_join(mf) %>%
  mutate(tips_dummy = as.integer(tips > 5),
                  tips_f = factor(tips),
                  age_group = factor(age_group, levels(mf$age_group)),
                  area_id = factor(area_id, levels(mf$area_id)),
                  period = factor(period, levels(mf$period)),
                  )

# join_aggr_births <- obs %>%
#   mutate(idx_obs = row_number(),
#          x=1) %>%
#   left_join(area_aggregation) %>%
#   left_join(mf, by=c("period", "age_group", "model_area_id" = "area_id")) %>%
#   type.convert()
# 
# A_aggr_births <- sparseMatrix(i = join_aggr_births$idx_obs, j=join_aggr_births$idx, x=join_aggr_births$x, use.last.ij = TRUE)

# join_aggr_coefs <- obs %>%
#   mutate(idx_obs = row_number()) %>%
#   left_join(area_aggregation) %>%
#   mutate(
#          id.age = group_indices(., age_group),
#          id.period = max(id.age) + group_indices(., period),
#          id.period = ifelse(id.period == min(id.period), 1, id.period-1),
#          id.area_id = max(id.period) + group_indices(., model_area_id),
#          id.area_id = ifelse(id.area_id == min(id.area_id), 1, id.area_id-1),
#          x=1) %>%
#   select(idx_obs, id.age, id.period, id.area_id, x) %>%
#   pivot_longer(-c(idx_obs, x))
# 
# A_aggr_coefs <- sparseMatrix(i = join_aggr_coefs$idx_obs, j=join_aggr_coefs$value, x=join_aggr_coefs$x, use.last.ij = TRUE)

join_nat <- crossing(area_id = "ZWE",
                     period = unique(mf$period),
                     age_group = unique(mf$age_group)
  ) %>%
  mutate(idx_out = row_number()) %>%
  left_join(area_aggregation %>% left_join(areas_long) %>% filter(area_level !=1) %>% select(area_id, model_area_id)) %>%
  left_join(obs) %>%
  mutate(idx = row_number(),
         x=1)

A_nat <- sparseMatrix(i = join_nat$idx_out, j=join_nat$idx, x=join_nat$x, use.last.ij = TRUE)
dim(A_nat)

mf_nat <- crossing(area_id = "ZWE",
                    period = unique(mf$period),
                    age_group = unique(mf$age_group)
 ) %>%
  mutate(idx = factor(row_number()))

obs_nat <- mics_asfr %>%
  mutate(period = factor(period, levels(mf$period))) %>%
  left_join(mf_nat) %>%
  select(area_id, period, age_group, tips, births, pys, idx) %>%
  mutate(tips_dummy = as.integer(tips > 5),
         tips_f = factor(tips, levels(obs$tips_f)),
         age_group = factor(age_group, levels(mf$age_group)),
         idx =factor(idx, levels(mf_nat$idx))
  )

X_mf <- model.matrix(~1 + age_group + period + area_id, mf)
# X_mf <- model.matrix(~1 + age_group + period, mf)

## Model frame design matrices

Z_spatial_mf <- sparse.model.matrix(~0 + area_id, mf)
Z_age_mf <- sparse.model.matrix(~0 + age_group, mf)
Z_period_mf <- sparse.model.matrix(~0 + period, mf)

#' This has dimensions (number of observations) x (number of rows in model frame (i.e crossing of age x time x space))
#' Many more rows than mf because observations has things we need to adjust for bias (e.g. tips), but not required in model frame. 
#' Column index (idx) has been joined onto obs dataframe from mf. So now each observation is labelled with the appropriate index in the mf
M_obs <- sparse.model.matrix(~0 + idx, obs) 
M_obs_nat <- sparse.model.matrix(~0 + idx, obs_nat)

### SPATIAL MODEL

# Z_spatial <- sparse.model.matrix(~0 + area_id, obs)

sh <- areas_long %>%
  filter(iso3 == iso3_current, naomi_level) %>%
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

#####
 
### TIPS RANDOM WALK

Z_tips <- sparse.model.matrix(~0 + tips_f, obs)
Z_tips_nat <- sparse.model.matrix(~0 + tips_f, obs_nat)

#' Create precision matrix for RW1
D_tips <- diff(diag(ncol(Z_tips)), differences = 1)
Q_tips <- as(t(D_tips) %*% D_tips, "dgCMatrix")

### AGE RANDOM WALK

Z_age <- sparse.model.matrix(~0 + age_group, obs)
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
# Q_period <- INLA:::inla.rw(ncol(Z_period), 2, scale.model = F)

# interaction_idx <- as.numeric(sort(unique(obs$id.interaction)))

### TIPS FIXED EFFECT

X_tips_dummy <- model.matrix(~0 + tips_dummy, obs)
X_tips_dummy_nat <- model.matrix(~0 + tips_dummy, obs_nat)

## Outputs

mf_out <- crossing(
    area_id = area_aggregation$area_id,
    age_group = unique(mf$age_group),
    period = unique(mf$period)
  ) %>%
  arrange(area_id, age_group, period) %>%
  mutate(out_idx = row_number())

join_out <- crossing(area_aggregation, 
         age_group = unique(mf$age_group),
         period = unique(mf$period)) %>%
  full_join(mf %>%
              select(area_id, age_group, period, idx), by = c("model_area_id" = "area_id", 
                                                          "age_group", 
                                                          "period")
            ) %>%
  full_join(mf_out) %>%
  # full_join(mf_out, by=c("area_id" = "area_id",
  #                        "period",
  #                        "age_group_out" = "age_group")
  #           ) %>%
  mutate(x=1)

A_out <- spMatrix(nrow(mf_out), nrow(mf), join_out$out_idx, as.integer(join_out$idx), join_out$x)

compile(here("tmb/fertility_tmb_dev.cpp"))               # Compile the C++ file
dyn.load(dynlib(here("tmb/fertility_tmb_dev")))

data <- list(X_mf = X_mf,
             M_obs = M_obs,
             M_obs_nat = M_obs_nat,
             X_tips_dummy = X_tips_dummy,
             X_tips_dummy_nat = X_tips_dummy_nat,
             Z_tips = Z_tips,
             Z_tips_nat = Z_tips_nat,
             # Z_age = Z_age,
             # Z_period = Z_period,
             # Z_spatial = Z_spatial,
             Z_spatial_mf = Z_spatial_mf,
             Z_age_mf = Z_age_mf,
             Z_period_mf = Z_period_mf,
             Z_interaction = sparse.model.matrix(~0 + id.interaction_3d, mf),
             # Z_interaction1 = sparse.model.matrix(~0 + id.interaction1, obs),
             # Z_interaction2 = sparse.model.matrix(~0 + id.interaction2, obs),
             # Z_interaction3 = sparse.model.matrix(~0 + id.interaction3, obs),
             # interaction_idx = interaction_idx,
             Q_tips = Q_tips,
             Q_age = Q_age,
             Q_period = Q_period,
             Q_spatial = Q_spatial,
             # A_national = A_national,
             log_offset = log(obs$pys),
             log_offset_nat = log(obs_nat$pys),
             births_obs = obs$births,
             births_obs_nat = obs_nat$births,
             pop = mf$population,
             A_out = A_out,
             A_nat = A_nat
             )

par <- list(beta_mf = rep(0, ncol(X_mf)),
            beta_tips_dummy = rep(0, ncol(X_tips_dummy)),
            u_tips = rep(0, ncol(Z_tips)),
            u_age = rep(0, ncol(Z_age_mf)),
            u_period = rep(0, ncol(Z_period_mf)),
            u_spatial_str = rep(0, ncol(Z_spatial_mf)),
            u_spatial_iid = rep(0, ncol(Z_spatial_mf)),
            eta = array(0, c(ncol(Z_spatial_mf), ncol(Z_age_mf), ncol(Z_period_mf))),
            # eta1 = array(0, c(ncol(Z_age), ncol(Z_period))),
            # eta2 = array(0, c(ncol(Z_spatial), ncol(Z_period))),
            # eta3 = array(0, c(ncol(Z_spatial), ncol(Z_age))),
            log_sigma_rw_tips = log(2.5),
            log_sigma_rw_age = log(2.5),
            log_sigma_rw_period = log(2.5),
            log_sigma_spatial = log(2.5),
            logit_spatial_rho = 0
            )


f <- mcparallel({TMB::MakeADFun(data = data,
                                      parameters = par,
                                      DLL = "fertility_tmb_dev",
                                      random = c("beta_mf",  "u_age", "beta_tips_dummy", "u_period", "u_spatial_str", "u_spatial_iid"),
                                      silent=0,
                                      checkParameterOrder=FALSE)
  })

mccollect(f)
             
f <-  MakeADFun(data = data,
                parameters = par,
                DLL = "fertility_tmb_dev",
                #random = c("beta_mf", "beta_tips_dummy", "u_tips", "u_age", "u_period", "u_spatial_str", "u_spatial_iid", "eta1", "eta2", "eta3"),
                random = c("beta_mf", "beta_tips_dummy", "u_tips", "u_age", "u_period", "u_spatial_str", "u_spatial_iid", "eta"),
                hessian = FALSE,
                checkParameterOrder=FALSE)

# f$env$tracepar <- TRUE
# f$report()

f$report()

fit <- nlminb(f$par, f$fn, f$gr)
rep <- sdreport(f)
tmb_res <- summary(rep)

mf %>% 
  # left_join(areas_long) %>%
  cbind(data.frame(val = f$report()$omega)) %>%
  type.convert() %>%
  mutate(source = "tmb") %>%
  bind_rows(inla_res) %>%
  filter(age_group == "20-24", area_id %in% mf$area_id) %>%
  ggplot(aes(x=period, y=val)) +
    geom_line(aes(group=source, color=source)) +
    # geom_point(data = filter(readRDS("countries/ZWE/data/ZWE_asfr_admin0.rds"), age_group == "20-24", !is.na(surveyid)) %>% rename(val = asfr), aes(color=surveyid), alpha=0.6) +
    facet_wrap(~area_id)

mf %>% 
  # left_join(areas_long) %>%
  cbind(data.frame(val = f$report()$omega)) %>%
  type.convert() %>%
  ggplot(aes(x=period, y=val, group=age_group, color=age_group)) +
  geom_line() +
  facet_wrap(~area_id)

mics_plot <- Map(calc_asfr_mics, mics_data$wm, y=list(1),
                 by = list(~area_id + survyear + surveyid + survtype),
                 tips = list(c(0,5)),
                 agegr= list(3:10*5),
                 period = list(1995:2019),
                 counts = TRUE,
                 bhdata = mics_data$bh_df) %>%
  bind_rows %>%
  type.convert() %>%
  filter(period <= survyear) %>%
  rename(age_group = agegr)

mf_out %>% 
  # left_join(areas_long) %>%
  cbind(data.frame(val = int$report()$lambda_out)) %>%
  type.convert() %>%
  left_join(areas_long) %>%
  filter(area_level ==0) %>%
  ggplot(aes(x=period, y=val, group=age_group, color=age_group)) +
  geom_line() +
  # geom_point(data=mics_plot %>% bind_rows(asfr_nat), aes(y=asfr, color=survtype)) +
  facet_wrap(~area_id)
