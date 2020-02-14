library(TMB)
library(tidyverse)
library(abind)
library(sf)
library(spdep)
library(INLA)
devtools::load_all("~/Documents/GitHub/naomi")
setwd("~/Documents/GitHub/subnat_fertility/")

asfr <- readRDS("~/Documents/GitHub/subnat_fertility/asfr_pred_subnat_15_2020_NEW.rds")[["MWI"]]
boundaries <- readRDS("~/Documents/GitHub/subnat_fertility/area_boundaries.rds")

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
  complete(period, area_id, fill=map(asfr[3:9], ~NA)) %>%
  mutate_if(is.numeric , replace_na, replace = -999) %>%
  select(-period) %>%
  nest(-area_id) %>%
  mutate(data = map(data, ~as.matrix(.x))) %>%
  pull(data) %>%
  invoke(abind, ., along=3)

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
Q <- INLA::inla.scale.model(diag(rowSums(adj)) - adj,
                                constr = list(A = matrix(1, 1, nrow(adj)), e = 0))


compile("fertility_tmb.cpp")               # Compile the C++ file
dyn.load(dynlib("fertility_tmb"))

f <-  MakeADFun(data=list(dat_m = dat_m, 
                          log_offset = log_offset, 
                          Q = Q,
                          Z_i = naomi::sparse_model_matrix(~0 + unique(area_id), asfr)), 
                parameters=list(alpha = 0,
                                rw_time = rep(0,22), 
                                rw_age = rep(0,7), 
                                log_prec_epsilon = 0,
                                log_prec_rw_time = 0,
                                log_prec_rw_age = 0,
                                logit_U_str = rep(0,33),
                                U_iid = rep(0,33), 
                                log_U_sigma = 0, # Variance
                                logit_U_rho = 0 # Share of variance between structured and unstructured components
                                ), 
                random = c("log_U_sigma", "logit_U_rho"),
                # map = list(), 
                DLL = "fertility_tmb",
                hessian = FALSE,
                checkParameterOrder=FALSE)

# f$env$tracepar <- FALSE
f$report()

fit = nlminb(f$par,f$fn,f$gr, control = list(iter.max = 100000, eval.max = 100000))

# opt <- do.call("optim", f)
rep <- sdreport(f)
tmb.res <- summary(rep)

cbind(crossing(district = 1:33, agegr = 1:7, period = 1:22), data.frame(point=exp(tmb.res[,1][row.names(tmb.res) == "log_mu"]))) %>%
  ggplot(aes(x=period, y=point, color=agegr, group=agegr)) +
  geom_line() +
  facet_wrap(~district)
