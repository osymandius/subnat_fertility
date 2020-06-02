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
devtools::load_all("~/Documents/GitHub/naomi")
# library(naomi)
library(here)

naomi_data_path <- "~/Imperial College London/HIV Inference Group - Documents/Analytical datasets/naomi-data"
mics_key <- read.csv(here("countries/mics_data_key.csv"))
# naomi_data_path <- "~/naomi-data"

source(here("R/inputs.R"))
source(here("R/fertility_funs.R"))

iso3 <- c("LSO", "MOZ", "NAM", "UGA", "ZMB", "ETH", "TZA")
iso3_current <-  "ZWE"
list2env(make_areas_population(iso3_current, naomi_data_path, full = FALSE), globalenv())

# exclude_districts <- areas_wide$area_id[areas_wide$area_id1 == "TZA_1_2"]
# exclude_districts <- c("UGA_3_029", "UGA_3_046")
# exclude_districts <- c("MOZ_2_0107", "MOZ_2_1009")
# exclude_districts <- "MWI_5_07"
exclude_districts=""

asfr <- get_asfr_pred_df(iso3_current, area_level = 1, areas_long, project = FALSE)
  # filter(survtype != "AIS")
  # filter(survtype == "DHS")

if(filter(mics_key, iso3 == iso3_current)$mics) {
  mics_asfr <- readRDS(here(grep("mics", list.files(paste0("countries/", iso3_current, "/data"), full.names = TRUE), value=TRUE)))
} else {
  mics_asfr <- NULL
}
# 
# mics_dat <- readRDS("input_data/mics_extract.rds")
# mics_asfr <- Map(calc_asfr_mics, mics_dat$wm[10], y=list(1),
#                  by = list(~area_id + survey_id),
#                  tips = list(c(0,15)),
#                  agegr= list(3:10*5),
#                  period = list(1995:2019),
#                  counts = TRUE,
#                  bhdata = mics_dat$bh_df[10]) %>%
#   bind_rows %>%
#   type.convert() %>%
#   separate(col=survey_id, into=c(NA, "survyear", NA), sep=c(3,7), remove = FALSE, convert = TRUE) %>%
#   filter(period <= survyear) %>%
#   rename(age_group = agegr)

mf <- make_model_frames(iso3_current, population, asfr, mics_asfr, exclude_districts, project=FALSE)

# saveRDS(mf, here(paste0("countries/", iso3_current, "/mods/", iso3_current, "_mf.rds")))

Z_spatial <- sparse.model.matrix(~0 + area_id, mf$mf_model)
Z_age <- sparse.model.matrix(~0 + age_group, mf$mf_model)
Z_period <- sparse.model.matrix(~0 + period, mf$mf_model)

M_obs <- sparse.model.matrix(~0 + idx, mf$dist$obs) 
Z_tips <- sparse.model.matrix(~0 + tips_f, mf$dist$obs)
X_tips_dummy <- model.matrix(~0 + tips_dummy, mf$dist$obs)
# X_urban_dummy <- model.matrix(~0 + urban_dummy, mf$dist$obs)

if(mf$mics_toggle) {
  M_obs_mics <- sparse.model.matrix(~0 + idx, mf$mics$obs) 
  Z_tips_mics <- sparse.model.matrix(~0 + tips_f, mf$mics$obs)
  X_tips_dummy_mics <- model.matrix(~0 + tips_dummy, mf$mics$obs)
}

R_spatial <- make_adjacency_matrix(iso3_current, areas_long, boundaries, exclude_districts, level= 1)
R_tips <- make_rw_structure_matrix(ncol(Z_tips), 1, TRUE)
R_age <- make_rw_structure_matrix(ncol(Z_age), 1, TRUE)
R_period <- make_rw_structure_matrix(ncol(Z_period), 2, TRUE)


# dyn.unload(dynlib(here("tmb/fertility_tmb_dev")))
compile(here("tmb/log_prec_tmb.cpp"))               # Compile the C++ file
dyn.load(dynlib(here("tmb/log_prec_tmb")))

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
  
  beta_tips_dummy = rep(0, ncol(X_tips_dummy)),
  u_tips = rep(0, ncol(Z_tips)),
  log_prec_rw_tips = 4,
  
  # beta_urban_dummy = rep(0, ncol(X_urban_dummy)),
  
  u_age = rep(0, ncol(Z_age)),
  log_prec_rw_age = 4,

  u_period = rep(0, ncol(Z_period)),
  log_prec_rw_period = 4,
  
  u_spatial_str = rep(0, ncol(Z_spatial)),
  # log_prec_spatial = 0
  u_spatial_iid = rep(0, ncol(Z_spatial)),
  log_sigma_spatial = 4,
  logit_spatial_rho = 0,
  
  eta1 = array(0, c(ncol(Z_period), ncol(Z_age))),
  log_prec_eta1 = 4,
  lag_logit_eta1_phi_age = 0,
  lag_logit_eta1_phi_period = 0,
  
  # eta2 = array(0, c(ncol(Z_spatial), ncol(Z_period))),
  # log_prec_eta2 = 4,
  # lag_logit_eta2_phi_period = 0
  
  eta3 = array(0, c(ncol(Z_spatial), ncol(Z_age))),
  log_prec_eta3 = 4,
  lag_logit_eta3_phi_age = 0
)

# "u_spatial_str", "u_spatial_iid", "eta1" , "eta1" "beta_tips_dummy",, "eta1""eta1", "u_tips", "beta_tips_dummy", , "u_spatial_iid", "eta3""u_age", "u_period",
tmb_int$random <- c("beta_0",  "u_age", "u_period", "u_spatial_str", "u_spatial_iid", "eta1", "u_tips", "beta_tips_dummy", "eta3")

if(mf$mics_toggle) {
  tmb_int$data <- c(tmb_int$data, "M_obs_mics" = M_obs_mics,
            "X_tips_dummy_mics" = list(X_tips_dummy_mics),
            "Z_tips_mics" = Z_tips_mics,
            "births_obs_mics" = list(mf$mics$obs$births),
            "log_offset_mics" = list(log(mf$mics$obs$pys)),
            "A_mics" = mf$mics$A_mics)
  tmb_int$par <- c(tmb_int$par, "beta_tips_dummy_mics" = rep(0, ncol(X_tips_dummy_mics)))
  tmb_int$random <- c(tmb_int$random, "beta_tips_dummy_mics")
}


f <- mcparallel({TMB::MakeADFun(data = tmb_int$data,
                                parameters = tmb_int$par,
                                DLL = "log_prec_tmb",
                                silent=0,
                                checkParameterOrder=FALSE)
})

mccollect(f)

obj <-  MakeADFun(data = tmb_int$data,
                  parameters = tmb_int$par,
                  DLL = "log_prec_tmb",
                  random = tmb_int$random,
                  hessian = FALSE)

f <- nlminb(obj$par, obj$fn, obj$gr)
f$par.fixed <- f$par
f$par.full <- obj$env$last.par

fit <- c(f, obj = list(obj))
fit$sdreport <- sdreport(fit$obj, fit$par)

class(fit) <- "naomi_fit"  # this is hacky...
fit <- sample_tmb(fit, random_only=FALSE)

qtls1 <- apply(fit$sample$lambda_out, 1, quantile, c(0.025, 0.5, 0.975))

mf$out$mf_out %>%
  mutate(lower = qtls1[1,],
         median = qtls1[2,],
         upper = qtls1[3,],
         source = "tmb") %>%
  type.convert() %>%
  # group_by(period, area_id) %>%
  # summarise(median = sum(median)) %>%
  left_join(areas_long) %>%
  filter(area_level == 1) %>%
  ggplot(aes(x=period, y=median)) +
  geom_line() +
  # geom_point(data=asfr_plot %>% left_join(areas_long) %>% filter(area_level == 1, !(area_id == admin1 & survtype == "MICS")), aes(y=asfr, group=survtype, color=survtype)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.3) +
  # geom_point(data = moz_anc, aes(y=anc_clients, x=year))+
  facet_grid(age_group~area_name) +
  ylim(0,0.4)

# saveRDS(fit, paste0("countries/", iso3_current, "/mods/", iso3_current, "_tmb_2020_05_11.rds"))
  


tips <- lapply(mods, function(mods) {
  
  int <- data.frame(summary(mods$sdreport))
  
  int %>% 
    mutate(id = rownames(.),
           idx = row_number()) %>%
    filter(idx %in% c(grep("u_tips", .data$id), 
                      grep("sigma_rw_tips", .data$id),
                      grep("beta_tips", .data$id)
    )
    ) %>%
    mutate(log_sigma_rw_tips = Estimate[id=="log_sigma_rw_tips"],
           beta_tips_dummy = Estimate[id == "beta_tips_dummy"]) %>%
    filter(idx > 200) %>%
    mutate(idx = row_number() - 1,
           dummy = ifelse(idx<6, 0, 1),
           beta_tips_dummy = beta_tips_dummy * dummy,
           est = (Estimate*exp(log_sigma_rw_tips)) + beta_tips_dummy)
})

names(tips) <- iso3

tips <- tips %>%
  bind_rows(.id="iso3")


hyper %>%
  pivot_longer(-c(hyper, mod, mics_data)) %>%
  ggplot(aes(y=hyper, x=value)) +
    geom_point(aes(color=mics_data)) +
    facet_wrap(~name)

tfr_plot <- readRDS(here("countries/ETH/data/ETH_tfr_plot.rds"))
  
mf$out$mf_out %>%
    mutate(lower = qtls1[1,],
           median = qtls1[2,],
           upper = qtls1[3,],
           source = "tmb") %>%
    type.convert() %>%
    # group_by(period, area_id) %>%
    # summarise(median = sum(median)) %>%
    left_join(areas_long) %>%
    filter(area_level == 3) %>%
    ggplot(aes(x=period, y=median, group=age_group)) +
    geom_line(aes(color=age_group)) +
    # geom_point(data=asfr_plot %>% filter(area_level == 1, !(area_id == "ZWE_1_19" & period > 2015 & survtype == "MICS")), aes(y=asfr, group=survtype, color=survtype)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill=age_group), alpha=0.3) +
    # geom_point(data = moz_anc, aes(y=anc_clients, x=year))+
    facet_wrap(~area_name) +
    ylim(0,0.4)
  
mf$out$mf_out %>%
  cbind(fit$sample$lambda_out) %>%
  group_by(area_id, period) %>%
  summarise_at(.vars = vars(-c(age_group, out_idx)), sum) %>%
  ungroup %>%
  mutate_at(.vars = vars(-c(area_id, period)), function(x) 5*x) %>%
  mutate(
    lower = select(., -c(area_id, period)) %>% apply(1, quantile, 0.025),
    median = select(., -c(area_id, period)) %>% apply(1, quantile, 0.5),
    upper = select(., -c(area_id, period)) %>% apply(1, quantile, 0.975)
         ) %>%
  select(area_id, period, lower, median, upper) %>%
  left_join(areas_long) %>%
  type.convert() %>%
  filter(area_level == 1) %>%
  ggplot(aes(x=period, y=median)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
    geom_line() +
    geom_point(data=tfr_plot %>% left_join(areas_long) %>% filter(area_level == 1), aes(y=tfr, group=survtype, color=survtype)) +
    facet_wrap(~area_name)

moz_anc <- read_csv("~/Downloads/Moz work/11 files/ancnaomi111219clean_update_20.01.20.csv") %>%
  left_join(areas_long) %>%
  group_by(parent_area_id, year) %>%
  summarise(anc_clients = sum(anc_clients)) %>%
  rename(area_id = parent_area_id)



int_df <- mf$mf_model %>%
  select(period, age_group, area_id) %>%
  left_join(mf$mf_model %>%
              select(period, age_group, area_id, id.interaction1) %>%
              distinct(age_group, period, id.interaction1) %>%
              arrange(id.interaction1) %>%
              mutate(eta1 = fit$obj$report()$eta1_v) %>%
              select(-id.interaction1)
  ) %>%
  left_join(mf$mf_model %>%
              select(period, age_group, area_id, id.interaction2) %>%
              distinct(area_id, period, id.interaction2) %>%
              arrange(id.interaction2) %>%
              mutate(eta2 = fit$obj$report()$eta2_v) %>%
              select(-id.interaction2)
            ) %>%
  left_join(mf$mf_model %>%
              select(period, age_group, area_id, id.interaction3) %>%
              distinct(area_id, age_group, id.interaction3) %>%
              arrange(id.interaction3) %>%
              mutate(eta3 = fit$obj$report()$eta3_v) %>%
              select(-id.interaction3)
  ) %>%
  mutate(eta_sum = eta1 + eta2 + eta3)

int_df %>% 
  filter(area_id == "ZWE_2_2") %>%
  ggplot(aes(x=period, y=eta2)) +
    geom_point()






mf$out$mf_out %>%
  mutate(lower = qtls[1,],
         median = qtls[2,],
         upper = qtls[3,],
         source = "tmb") %>%
  type.convert() %>%
  filter(period == 2017) %>%
  group_by(area_id) %>%
  summarise(tfr = 5*sum(median)) %>%
  left_join(areas_long) %>%
  filter(area_level == 1) %>%
  left_join(boundaries) %>%
  st_as_sf  %>%
  ggplot() +
    geom_sf(aes(geometry = geometry, fill=tfr)) +
    scale_fill_viridis(guide = guide_colorbar(frame.colour = "black", frame.linewidth = 1, ticks.colour = "black"))+
    coord_sf(datum=NA)+
    labs(title="", fill="TFR")

mf$out$mf_out %>%
  mutate(lower = qtls[1,],
         median = qtls[2,],
         upper = qtls[3,],
         source = "tmb") %>%
  type.convert() %>%
  group_by(area_id, period) %>%
  summarise(tfr = 5*sum(median)) %>%
  left_join(areas_long) %>%
  filter(area_level == 1) %>%
  ggplot(aes(x=period, y=tfr)) +
    geom_line() +
    geom_point(data= admin1_tfr %>% left_join(areas_long %>% select(area_id, area_name)) %>% filter(area_id != "ZWE_1_10"), aes(color=survtype)) +
    facet_wrap(~area_name)

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

admin1_tfr <- Map(calc_tfr, dat$ir,
            by = list(~country + surveyid + survtype + survyear + area_id),
            tips = dat$tips_surv,
            agegr= list(3:10*5),
            period = list(1995:2017)) %>%
  bind_rows %>%
  type.convert %>%
  filter(period<=survyear) %>%
  mutate(iso3 = countrycode(country, "country.name", "iso3c"),
         iso3 = ifelse(country == "Eswatini", "SWZ", iso3)) %>%
  select(-country)

foo <- admin1_tfr %>% bind_rows(mics_plot %>%
  group_by(area_id, period, surveyid, survtype) %>%
  summarise(tfr = 5*sum(asfr))
)

debugonce(calc_tfr)
calc_tfr(mics_data$wm[[1]], by = ~surveyid + survtype + survyear + area_id,
         tips = c(0,15),
         agegr= c(3:10*5),
         period = c(1995:2017),
         clusters=~cluster,
         id="unique_id",
         dob="wdob",
         strata=NULL,
         intv = "doi",
         weight= "wmweight",
         bvars = "cdob",
         bhdata = mics_data$bh_df[[1]]
)
