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
devtools::load_all("~/Documents/GitHub/naomi")

library(here)

naomi_data_path <- "~/Documents/GitHub/naomi-data"

source(here("R/inputs.R"))
source(here("R/fertility_funs.R"))

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
  filter(period == min(period), sex=="female")

# population <- population %>%
#   left_join(areas_long) %>%
#   filter(period == min(period), sex == "female") %>%
#   group_by(parent_area_id, age_group) %>%
#   summarise(population = sum(population)) %>%
#   ungroup %>%
#   rename(area_id = parent_area_id)

area_merged <- st_read(file.path(naomi_data_path, "ZWE/data/zwe_areas.geojson"))
areas <- create_areas(area_merged = area_merged)
area_aggregation <- create_area_aggregation(area_merged$area_id[area_merged$naomi_level], areas)


## Make model frame.
mf <- crossing(period = factor(1995:2015),
               age_group = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49"),
               area_id = filter(areas_long, iso3 == iso3_current, area_level == 2)$area_id) %>%
  left_join(population %>%
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
         id.interaction3 = factor(group_indices(., age_group, area_id))
  )

obs <- asfr %>%
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

mf_nat <- crossing(area_id = unique(mics_asfr$area_id),
                   period = unique(mf$period),
                   age_group = unique(mf$age_group)
) %>%
  mutate(idx = factor(row_number()))

join_nat <- mf_nat %>%
  rename(idx_row = idx) %>%
  left_join(area_aggregation) %>%
  left_join(mf, by=c("age_group", "period", "model_area_id" = "area_id")) %>%
  mutate(idx_col = row_number(),
         x=1) %>%
  type.convert()

A_nat <- sparseMatrix(i = join_nat$idx_row, j=join_nat$idx_col, x=join_nat$x, use.last.ij = TRUE)

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

Z_spatial <- sparse.model.matrix(~0 + area_id, mf)
Z_age <- sparse.model.matrix(~0 + age_group, mf)
Z_period <- sparse.model.matrix(~0 + period, mf)

#' This has dimensions (number of observations) x (number of rows in model frame (i.e crossing of age x time x space))
#' Many more rows than mf because observations has things we need to adjust for bias (e.g. tips), but not required in model frame. 
#' Column index (idx) has been joined onto obs dataframe from mf. So now each observation is labelled with the appropriate index in the mf
M_obs <- sparse.model.matrix(~0 + idx, obs) 
M_obs_nat <- sparse.model.matrix(~0 + idx, obs_nat)

### SPATIAL MODEL

sh <- areas_long %>%
  filter(iso3 == iso3_current, area_level ==2) %>%
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

D_age <- diff(diag(ncol(Z_age)), differences = 1)
Q_age <- t(D_age) %*% D_age
diag(Q_age) <- diag(Q_age) + 1E-6
Q_age <- as(Q_age, "dgCMatrix")

### TIME RANDOM WALK

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

age_group_out <- c(as.character(unique(mf$age_group)), "15-49")

age_group_join <- get_age_groups() %>%
  filter(age_group %in% age_group_out) %>%
  setNames(paste0(names(.), "_out")) %>%
  crossing(get_age_groups() %>%
                    filter(age_group %in% unique(mf$age_group))) %>%
  filter(age_group_start_out <= age_group_start,
                age_group_span_out == Inf |
                  (age_group_start + age_group_span) <=
                  (age_group_start_out + age_group_span_out)) %>%
  select(age_group_out, age_group)

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
  # full_join(mf_out, by=c("area_id",
  #                        "period",
  #                        "age_group_out" = "age_group")
  #           ) %>%
  mutate(x=1)

A_out <- spMatrix(nrow(mf_out), nrow(mf), join_out$out_idx, as.integer(join_out$idx), join_out$x)

join_rate <- crossing(area_id = unique(join_out$area_id),
                      age_group_join,
                      period = unique(mf$period)) %>%
  mutate(idx = group_indices(., area_id, age_group_out, period)) %>%
  left_join(join_out %>% select(area_id, age_group, period, out_idx) %>% distinct %>% rename(idx_col = out_idx)) %>%
  mutate(x=ifelse(age_group_out == "15-49", 5, 1))
  
A_rate <- spMatrix(nrow(join_rate), nrow(A_out), join_rate$idx, as.integer(join_rate$idx_col), join_rate$x)

dyn.unload(dynlib(here("tmb/fertility_tmb_dev")))
compile(here("tmb/fertility_tmb_dev.cpp"))               # Compile the C++ file
dyn.load(dynlib(here("tmb/fertility_tmb_dev")))

data <- list(X_mf = X_mf,
             M_obs = M_obs,
             # M_obs_nat = M_obs_nat,
             X_tips_dummy = X_tips_dummy,
             # X_tips_dummy_nat = X_tips_dummy_nat,
             Z_tips = Z_tips,
             # Z_tips_nat = Z_tips_nat,
             Z_age = Z_age,
             Z_period = Z_period,
             Z_spatial = Z_spatial,
             # Z_interaction = sparse.model.matrix(~0 + id.interaction, mf),
             # Z_interaction1 = sparse.model.matrix(~0 + id.interaction1, mf),
             # Z_interaction2 = sparse.model.matrix(~0 + id.interaction2, mf),
             # Z_interaction3 = sparse.model.matrix(~0 + id.interaction3, mf),
             Q_tips = Q_tips,
             Q_age = Q_age,
             Q_period = Q_period,
             Q_spatial = Q_spatial,
             log_offset = log(obs$pys),
             # log_offset_nat = log(obs_nat$pys),
             births_obs = obs$births,
             # births_obs_nat = obs_nat$births,
             pop = mf$population,
             A_out = A_out
             # A_nat = A_nat
             # A_rate = A_rate
             )

par <- list(
            beta_mf = rep(0, ncol(X_mf)),
            # beta_0 = 0,
            beta_tips_dummy = rep(0, ncol(X_tips_dummy)),
            u_tips = rep(0, ncol(Z_tips)),
            u_age = rep(0, ncol(Z_age)),
            u_period = rep(0, ncol(Z_period)),
            u_spatial_str = rep(0, ncol(Z_spatial)),
            u_spatial_iid = rep(0, ncol(Z_spatial)),
            # eta = array(0, c(ncol(Z_spatial), ncol(Z_age), ncol(Z_period))),
            # eta1 = array(0, c(ncol(Z_period), ncol(Z_age))),
            # eta2 = array(0, c(ncol(Z_spatial), ncol(Z_period))),
            # eta3 = array(0, c(ncol(Z_spatial), ncol(Z_age))),
            # log_sigma_rw_period = log(2.5),
            # log_sigma_rw_tips = log(2.5),
            # log_sigma_rw_age = log(2.5),
            # log_sigma_eta1 = log(2.5),
            log_prec_rw_period = 0,
            log_prec_rw_tips = 0,
            log_prec_rw_age = 0,
            # log_prec_eta1 = log(2.5),
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
                random = c("beta_mf", "beta_tips_dummy", "u_tips", "u_age", "u_period", "u_spatial_str", "u_spatial_iid"),
                hessian = FALSE,
                checkParameterOrder=FALSE)

f <- nlminb(obj$par, obj$fn, obj$gr)
f$par.fixed <- f$par
f$par.full <- obj$env$last.par

fit <- c(f, obj = list(obj))
fit$sdreport <- sdreport(fit$obj, fit$par)
fit <- sample_tmb_test(fit)



summary(fit$sdreport)

qtls <- apply(fit$sample$lambda_out, 1, quantile, c(0.025, 0.5, 0.975))

mf_out %>%
  mutate(lower = qtls[1,],
         median = qtls[2,],
         upper = qtls[3,]) %>%
  left_join(areas_long) %>%
  filter(area_level ==0) %>%
  ggplot(aes(x=period, y=median, group=age_group)) +
  geom_line(aes(color=age_group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=age_group), alpha=0.3) +
  facet_wrap(~area_id) +
  ylim(0,0.5)
 
# tfr <- mf_out %>%
#   cbind(fit$sample$lambda_out) %>%
#   select(-c(age_group, out_idx)) %>%
#   group_by(area_id, period) %>%
#   summarise_at(vars(-group_cols()), function(x) 5*sum(x)) %>%
#   ungroup
# 
# tfr <- tfr %>%
#   select(area_id, period) %>%
#   cbind(tfr %>%
#           select(-c(area_id, period)) %>%
#           apply(., 1, quantile, c(0.025, 0.5, 0.975)) %>%
#           data.frame() %>%
#           mutate(var = c("lower", "median", "upper")) %>%
#           pivot_longer(-var) %>%
#           pivot_wider(names_from=var) %>%
#           select(-name)
#   )

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
  cbind(data.frame(val = f$report()$lambda_out)) %>%
  type.convert() %>%
  left_join(areas_long) %>%
  filter(area_level ==0) %>%
  ggplot(aes(x=period, y=val, group=source, color=source)) +
  geom_line() +
  # geom_point(data=mics_plot %>% bind_rows(asfr_admin1), aes(y=asfr, color=survtype)) +
  facet_wrap(~area_id)

tfr %>%
  type.convert() %>%
  mutate(source = "model") %>%
  left_join(areas_long) %>%
  filter(area_level == 1) %>%
  rename(val = median) %>%
  bind_rows(zim_province_tfr) %>%
  filter(period %in% 2000:2018) %>%
  ggplot(aes(x=period, y=val)) +
    geom_line(aes(group=source, color=source)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill=source), alpha=0.3) +
    # geom_point(data=mics_plot %>% bind_rows(asfr_admin1), aes(y=asfr, color=survtype)) +
    facet_wrap(~area_id) +
    ylim(0, 20)



df <- summary(fit$sdreport)

mf %>%
 cbind(val = fit$obj$report()$lambda) %>%
  ggplot(aes(x=period, y=val, group=age_group))
