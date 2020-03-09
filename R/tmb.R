library(TMB)
library(tidyverse)
library(abind)
library(sf)
library(spdep)
library(Matrix)

devtools::load_all("~/GitHub/naomi")

setwd("~/GitHub/subnat_fertility/")
source("R/inputs.R")

lapply(list.files("~/GitHub/naomi/R", full.names = TRUE), source)

iso3_current <- "ZWE"

list2env(make_areas_population("ZWE", "~/GitHub/naomi-data/"), globalenv())

asfr <- readRDS(paste0("countries/", iso3_current, "/data/", iso3_current, "_asfr_admin", 2, ".rds"))

population <- population %>%
  filter(period == min(period))


population <- population %>%
  group_by(period, sex, age_group) %>%
  summarise(population = sum(population)) %>%
  mutate(area_id = "ZWE") %>%
  ungroup

mf <- crossing(period = factor(1995:2015),
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
         # id.interaction = factor(group_indices(., age_group, period, area_id))
         id.interaction1 = factor(group_indices(., age_group, period)),
         id.interaction2 = factor(group_indices(., period, area_id)),
         id.interaction3 = factor(group_indices(., age_group, area_id))
         )


obs <- asfr %>%
  filter(!is.na(surveyid)) %>%
  select(area_id, period, age_group, tips, births, pys) %>%
  mutate(age_group = factor(age_group, levels(mf$age_group)),
         period = factor(period)) %>%
  left_join(mf, by = c("area_id", "period", "age_group")) %>%
  mutate(tips_dummy = as.integer(tips > 5),
         tips_f = factor(tips),
         area_id = factor(area_id, levels(mf$area_id))) 

X_mf <- model.matrix(~1 + age_group + period + area_id, mf)
#X_mf <- model.matrix(~1 + age_group + period, mf)

#' This has dimensions (number of observations) x (number of rows in model frame (i.e crossing of age x time x space))
#' Many more rows than mf because observations has things we need to adjust for bias (e.g. tips), but not required in model frame. 
#' Column index (idx) has been joined onto obs dataframe from mf. So now each observation is labelled with the appropriate index in the mf
M_mf_obs <- Matrix::sparse.model.matrix(~0 + idx, obs) 
 
### TIPS RANDOM WALK

Z_tips <- sparse.model.matrix(~0 + tips_f, obs)
#' Create precision matrix for RW1
D_tips <- diff(diag(ncol(Z_tips)), differences = 1)
Q_tips <- as(t(D_tips) %*% D_tips, "dgCMatrix")
# Q_tips <- INLA:::inla.rw(ncol(Z_tips), 1, scale.model = F)

### AGE RANDOM WALK

Z_age <- sparse.model.matrix(~0 + age_group, obs)
D_age <- diff(diag(ncol(Z_age)), differences = 1)
Q_age <- t(D_age) %*% D_age
diag(Q_age) <- diag(Q_age) + 1E-6
Q_age <- as(Q_age, "dgCMatrix")
# Q_age <- INLA:::inla.rw(ncol(Z_age), 1, scale.model = F)

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

### ICAR

Z_spatial <- sparse.model.matrix(~0 + area_id, obs)

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

## Outputs
area_merged <-  st_read(files[[iso3_current]][["areas"]])
areas <- create_areas(area_merged = area_merged)
area_aggregation <- create_area_aggregation(area_merged$area_id[area_merged$naomi_level], areas)

age_groups <- unique(mf$age_group)
# age_group_out <- c(unique(as.character(mf$age_group)), "15-49")
# 
# age_group_join <- get_age_groups() %>%
#   dplyr::filter(age_group %in% age_group_out) %>%
#   stats::setNames(paste0(names(.), "_out")) %>%
#   tidyr::crossing(get_age_groups() %>%
#                     dplyr::filter(age_group %in% age_groups)) %>%
#   dplyr::filter(age_group_start_out <= age_group_start,
#                 age_group_span_out == Inf |
#                   (age_group_start + age_group_span) <=
#                   (age_group_start_out + age_group_span_out)) %>%
#   dplyr::select(age_group_out, age_group)

mf_out <- crossing(
  area_id = area_aggregation$area_id,
  # age_group = age_group_join$age_group_out,
  age_group = age_groups,
  period = unique(mf$period)
) %>%
  # mutate(indicator = ifelse(age_group == "15-49", "tfr", "asfr")) %>%
  arrange(area_id, age_group, period) %>%
  mutate(out_idx = row_number())

join_out <- crossing(area_aggregation, 
         age_group = age_groups,
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

dyn.unload(dynlib("tmb/fertility_tmb_dev"))
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
             # Z_interaction1 = sparse.model.matrix(~0 + id.interaction1, obs),
             # Z_interaction2 = sparse.model.matrix(~0 + id.interaction2, obs),
             # Z_interaction3 = sparse.model.matrix(~0 + id.interaction3, obs),
             # interaction_idx = interaction_idx,
             Q_tips = Q_tips,
             Q_age = Q_age,
             Q_period = Q_period,
             Q_spatial = Q_spatial,
             log_offset = log(obs$pys),
             births_obs = obs$births
             # pop = mf$population,
             # A_out = A_out
             )

par <- list(beta_mf = rep(0, ncol(X_mf)),
            beta_tips_dummy = rep(0, ncol(X_tips_dummy)),
            u_tips = rep(0, ncol(Z_tips)),
            u_age = rep(0, ncol(Z_age)),
            u_period = rep(0, ncol(Z_period)),
            u_spatial_str = rep(0, ncol(Z_spatial)),
            u_spatial_iid = rep(0, ncol(Z_spatial)),
            # eta = array(0, c(ncol(Z_spatial), ncol(Z_age), ncol(Z_period))),
            # eta1 = array(0, c(ncol(Z_age), ncol(Z_period))),
            # eta2 = array(0, c(ncol(Z_spatial), ncol(Z_period))),
            # eta3 = array(0, c(ncol(Z_spatial), ncol(Z_age))),
            log_sigma_rw_tips = log(2.5),
            log_sigma_rw_age = log(2.5),
            log_sigma_rw_period = log(2.5),
            log_sigma_spatial = log(2.5),
            logit_spatial_rho = 0
            )
             
f <-  MakeADFun(data = data,
                parameters = par,
                DLL = "fertility_tmb_dev",
                # random = c("beta_mf", "beta_tips_dummy", "u_tips", "u_age", "u_period", "u_spatial_str", "u_spatial_iid", "eta1", "eta2", "eta3"),
                random = c("beta_mf", "beta_tips_dummy", "u_tips", "u_age", "u_period", "u_spatial_str", "u_spatial_iid"),
                #random = c("beta_mf", "beta_tips_dummy", "eta1"),
                hessian = FALSE,
                checkParameterOrder=FALSE)

# f$env$tracepar <- TRUE
# f$report()

fit <- nlminb(f$par, f$fn, f$gr)
rep <- sdreport(f)
tmb_res <- summary(rep)

split(exp(tmb_res[,1]), rownames(tmb_res))

inla_r <- r_list[["ZWE"]]

inla_r$id.age_group

inla_res <- get_mod_results_test(zwe.m, dat)%>%
  mutate(source = "inla") %>%
  rename(val = median)

inla_res <- res_list[["ZWE"]] %>%
  mutate(source = "inla") %>%
  rename(val = median)

mf %>% 
  # left_join(areas_long) %>%
  cbind(data.frame(val = f$report()$omega)) %>%
  type.convert() %>%
  mutate(source = "tmb") %>%
  bind_rows(inla_res) %>%
  filter(age_group == "20-24") %>%
  ggplot(aes(x=period, y=val, group=source, color=source)) +
    geom_line() +
    facet_wrap(~area_id)

data.frame(val = f$report()$u_period)


f$report()

asfr <- readRDS("~/Downloads/zwe_nat_asfr.rds")

asfr <- crossing(
  agegr = unique(asfr$agegr),
  period = 1995:2015,
  area_id = "ZWE",
  country = "Zimbabwe"
) %>%
  bind_rows(asfr) %>%
  mutate(area_id = "ZWE") %>%
  mutate(id.agegr = group_indices(., agegr), 
         id.period = group_indices(., period), 
         id.tips = ifelse(is.na(tips), NA, group_indices(., tips)), 
         tips_dummy = ifelse(tips<5, 0, 1),
         id = row_number())

national_mod <- run_mod(asfr)
debugonce(get_mod_results)
zwe_res_nat <- get_mod_results_test(national_mod, asfr)
get_mod_results(zwe_run, asfr)

int <- mf_out %>% 
  cbind(data.frame(nat = f$report()$omega)) %>%
  left_join(areas_long) %>%
  select(area_id,area_level,period,age_group, nat) %>%
  mutate(period = as.numeric(as.character(period)))%>%
  left_join(zwe_res_nat %>%
              rename(inla = median,age_group = agegr) %>%
              select(area_id,period,age_group,inla) %>%
              left_join(areas_long)
            ) %>%
  select(area_id, area_level, period, age_group, nat, inla) %>%
  filter(area_level == 0) %>%
  pivot_longer(-c(area_id,period,area_level,age_group))

int %>%
  filter(age_group == "20-24") %>%
  ggplot(aes(x=period, y=value)) +
    geom_line(aes(color=name, group=name)) +
    # geom_point(data= asfr_admin1 %>% filter(agegr == "20-24"), aes(x=period, y=asfr)) +
    facet_wrap(~area_id)

asfr_admin1 <- readRDS("~/Documents/GitHub/subnat_fertility/asfr_admin1_newh.rds")[["ZWE"]]

obs %>%
  type.convert %>%
  left_join(
    data.frame(id.interaction = 1:147, val = tmb_res[,1][rownames(tmb_res) == "eta"])
  ) %>%
  ggplot(aes(x=period, y=val, color=age_group, group=age_group)) +
  geom_line()

inla_mods <- readRDS("~/Documents/GitHub/subnat_fertility/2019_12_4_mods.rds")

zwe_inla <- national_mod

zwe_random <- lapply(zwe_inla$summary.random, function(x) select(x, ID, mean))

zwe_random$id.agegr <- zwe_random$id.agegr %>%
  rename(id.agegr = ID) %>%
  rename(age_rw1 = mean)
  
zwe_random$id.period <- zwe_random$id.period %>%
  rename(id.period = ID)%>%
  rename(period_rw2 = mean)

zwe_random$id.agegr2 <- zwe_random$id.agegr2 %>%
  mutate(id.period = rep(1:26, each=7)) %>%
  rename(id.agegr = ID, age_period = mean)

zwe_random$id.tips <- zwe_random$id.tips %>%
  rename(id.tips = ID, tips_rw = mean)

zwe_random$id.district <- zwe_random$id.district %>%
  mutate(thing = c(rep("top", times=63), rep("bottom", times="63"))) %>%
  rename(id.district = ID, spatial = mean) %>%
  filter(thing == "top")

zwe_random$id.district2 <- zwe_random$id.district2 %>% 
  mutate(thing = rep(c(rep("top", times=63), rep("bottom", times="63")), times=26),
        id.period = rep(1:26, each=126)
  ) %>%
  rename(id.district = ID, space_time = mean) %>%
  filter(thing == "top")

zwe_random$id.district3 <- zwe_random$id.district3 %>%
  mutate(thing = rep(c(rep("top", times=63), rep("bottom", times="63")), times=7),
         id.agegr = rep(1:7, each=126)
  ) %>%
  rename(id.district = ID, space_age = mean) %>%
  filter(thing == "top")

int <- f$report()

test <- mf %>%
  # mutate(id.interaction1 = as.integer(id.interaction1),
  #        id.interaction2 = as.integer(id.interaction2),
  #        id.interaction3 = as.integer(id.interaction3)) %>%
  # left_join(
  #   data.frame(age_period = int$eta1_v) %>%
  #     mutate(id.interaction1 = row_number())
  # ) %>%
  # left_join(
  #   data.frame(space_time = int$eta2_v) %>%
  #     mutate(id.interaction2 = row_number())
  # ) %>%
  # left_join(
  #   data.frame(space_age = int$eta3_v) %>%
  #     mutate(id.interaction3 = row_number())
  # ) %>%
  left_join(data.frame(age_rw1 = int$u_age) %>%
              mutate(age_group = unique(mf$age_group))
            ) %>%
  left_join(data.frame(spatial = int$spatial) %>%
              mutate(area_id = unique(mf$area_id))) %>%
  left_join(data.frame(period_rw2 = int$u_period) %>%
              mutate(period = unique(mf$period))) %>%
  select(-c(population, idx, id.interaction1, id.interaction2, id.interaction3)) %>%
  mutate(source = "TMB",
         period = as.integer(as.character(period))) %>%
  pivot_longer(-c(area_id, period, age_group, source)) %>%
  bind_rows(
    asfr %>%
      select(area_id, period, age_group, id.period, id.agegr, id.district) %>%
      distinct %>%
      left_join(zwe_random[[1]]) %>%
      left_join(zwe_random[[2]]) %>%
      left_join(zwe_random[[3]]) %>%
      select(-c(id.period, id.agegr, id.district, thing)) %>%
      mutate(source = "INLA") %>%
      filter(period < 2016) %>%
      pivot_longer(-c(area_id, period, age_group, source))
    
  )

test %>%
  filter(name == "period_rw2") %>%
  ggplot(aes(x=period, y=value, group=source, color=source)) +
    geom_line()

test %>%
  filter(name == "age_rw1") %>%
  ggplot(aes(x=age_group, y=value, group=source, color=source)) +
  geom_line()

test %>%
  filter(name == "spatial") %>%
  ggplot(aes(x=area_id, y=value, group=source, color=source)) +
  geom_line()

INLA:::inla.rw(5,2,TRUE)
