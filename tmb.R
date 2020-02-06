library(TMB)
library(tidyverse)
library(abind)
library(sf)
library(spdep)
setwd("~/Documents/GitHub/subnat_fertility/")

asfr <- readRDS("~/Documents/GitHub/subnat_fertility/asfr_pred_subnat_15_2020_NEW.rds")[["MWI"]]
boundaries <- readRDS("~/Documents/GitHub/subnat_fertility/area_boundaries.rds")

iso3_codes <- c("LSO", "MOZ", "MWI", "NAM", "SWZ", "TZA", "UGA", "ZMB", "ZWE")

paths <- paste0("~/Documents/GitHub/naomi-data/", iso3_codes, "/data/", tolower(iso3_codes), "_areas.geojson")
area_cols <- c("iso3", "area_id", "area_name", "area_level", "parent_area_id", "naomi_level")

areas_long <- lapply(paths, read_sf) %>% lapply(function(x) {
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

dat_m <- asfr %>%
  filter(!is.na(surveyid)) %>%
  select(area_id, period, agegr, births, pys) %>%
  group_by(period, agegr, area_id) %>%
  summarise_all(sum) %>%
  ungroup %>%
  select(-pys) %>%
  mutate(agegr = as.numeric(factor(agegr))) %>%
  pivot_wider(names_from = agegr, values_from = births) %>%
  complete(period, area_id, fill=map(asfr[2:9], ~NA)) %>%
  select(-period) %>%
  nest(-area_id) %>%
  mutate(data = map(data, ~as.matrix(.x))) %>%
  pull(data) %>%
  invoke(abind, ., along=3)

log_offset <- asfr %>%
  filter(!is.na(surveyid)) %>%
  select(area_id, period, agegr, births, pys) %>%
  group_by(period, agegr, area_id) %>%
  summarise_all(sum) %>%
  ungroup %>%
  select(-births) %>%
  mutate(agegr = as.numeric(factor(agegr)), pys = log(pys)) %>%
  pivot_wider(names_from = agegr, values_from = pys) %>%
  complete(period, area_id, fill=map(asfr[2:9], ~NA)) %>%
  select(-period) %>%
  nest(-area_id) %>%
  mutate(data = map(data, ~as.matrix(.x))) %>%
  pull(data) %>%
  invoke(abind, ., along=3)

# dat_m <- asfr %>% 
#   filter(!is.na(surveyid)) %>%
#   select(period, births, pys) %>%
#   group_by(period) %>%
#   summarise_all(sum) %>%
#   ungroup %>%
#   select(-pys, -period) %>%
#   as.matrix()
# 
# log_offset <- asfr %>% 
#   filter(!is.na(surveyid)) %>%
#   select(period,  births, pys) %>%
#   group_by(period) %>%
#   summarise_all(sum) %>%
#   ungroup %>%
#   mutate(pys = log(pys)) %>%
#   select(-births) %>%
#   select(-period) %>%
#   as.matrix()


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
diag <- diag(rowSums(adj))
icar_adj <- diag-adj


compile("fertility_tmb.cpp")               # Compile the C++ file
dyn.load(dynlib("fertility_tmb"))

f <-  MakeADFun(data=list(dat_m = dat_m, log_offset = log_offset, P = as(icar_adj, "dgTMatrix")), 
                parameters=list(alpha_a = rep(0,7), 
                                beta_a = rep(0, 7), 
                                rw_time = rep(0,22), 
                                rw_age = rep(0,7), 
                                # rw_interaction = matrix(0, 22, 7),
                                log_prec_epsilon =1,
                                log_prec_rw_time = 1,
                                log_prec_rw_age = 1,
                                V = rep(0, 33), 
                                W = rep(0, 33), 
                                log_sigma2_V = 0, 
                                log_sigma2_U = 0
                                # log_prec_rw_interaction = 1
                                ), 
                # random = c("V", "W"),
                # map = list(), 
                DLL = "fertility_tmb")

f$env$tracepar <- TRUE
# f$report()

fit = nlminb(f$par,f$fn,f$gr, control = list(iter.max = 100000, eval.max = 100000))
print(fit)

int <- data.frame(agegr = 0:6, intercept = fit$par[1:7], slope = fit$par[8:14])

crossing(period = 0:21, agegr = 0:6) %>%
  left_join(int) %>%
  left_join(data.frame(period = 0:21, rw_time = fit$par[names(fit$par) == "rw_time"])) %>%
  left_join(data.frame(agegr= 0:6, rw_age = fit$par[names(fit$par) == "rw_age"])) %>%
  # left_join(data.frame(period = rep(0:21, each=7), agegr = rep(0:6, times=22), interaction = fit$par[names(fit$par) == "rw_interaction"])) %>%
  mutate(point = exp(intercept + period*slope + rw_time + rw_age)) %>%
  ggplot(aes(x=period, y=point, color=agegr, group=agegr)) +
    geom_line()

data.frame(period = rep(0:21, each=7), agegr = rep(0:6, times=22), nu = fit$par[names(fit$par) == "nu"]) %>%
  left_join(int) %>%
  mutate(point = exp(nu + intercept + period*slope)) %>%
  ggplot(aes(x=period, y=point, color=agegr, group=agegr)) +
    geom_line()
  

plot(exp(fit$par[["alpha"]] + 0:21*fit$par[["beta"]]))

rp = sdreport(??, fit$par)
rd = summary(rp)

outputPlot <- matrix(exp(fit$par[1:154]), nrow=22, ncol=7) %>%
  data.frame %>%
  `colnames<-`(., unique(asfr$agegr)) %>%
  mutate(period = 1995:2016) %>%
  pivot_longer(cols = contains("-"), names_to = "agegr") %>%
  ggplot(aes(x=period, y=value, color=agegr, group=agegr)) +
    geom_line()

input_rates <- dat_m/exp(log_offset) %>%
  data.frame

inputPlot <- input_rates %>%
  `colnames<-`(., unique(asfr$agegr)) %>%
  mutate(period = 1995:2016)%>%
  pivot_longer(cols = contains("-"), names_to = "agegr") %>%
  ggplot(aes(x=period, y=value, color=agegr, group=agegr)) +
    geom_line()

gridExtra::grid.arrange(inputPlot, outputPlot)
