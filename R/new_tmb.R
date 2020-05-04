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
# naomi_data_path <- "~/naomi-data"

source(here("R/inputs.R"))
source(here("R/fertility_funs.R"))

iso3_current <- "ZWE"
# exc <- areas_wide$area_id[areas_wide$area_id1 == "TZA_1_2"]
# exc <- c("UGA_3_029", "UGA_3_046")
# exc <- c("MOZ_2_0107", "MOZ_2_1009")
# exc <- "MWI_5_07"
exc=""

list2env(make_areas_population(iso3_current, naomi_data_path, full = FALSE), globalenv())

asfr <- get_asfr_pred_df(iso3_current, 2, project = FALSE)
mics_asfr <- readRDS(here("countries", paste0(iso3_current, "/data/", iso3_current, "_mics_admin", 1, ".rds")))

# mics_dat <- readRDS("input_data/mics_extract.rds")
# mics_asfr <- Map(calc_asfr_mics, mics_dat$wm[10], y=list(1),
#                  by = list(~area_id + survey_id),
#                  tips = list(c(0:15)),
#                  agegr= list(3:10*5),
#                  period = list(1995:2019),
#                  counts = TRUE,
#                  bhdata = mics_dat$bh_df[10]) %>%
#   bind_rows %>%
#   type.convert() %>%
#   separate(col=survey_id, into=c(NA, "survyear", NA), sep=c(3,7), remove = FALSE, convert = TRUE) %>%
#   filter(period <= survyear) %>%
#   rename(age_group = agegr)

mf <- make_model_frames(iso3_current, population, asfr, mics_asfr, exclude_districts = exc, project=2020)

Z_spatial <- sparse.model.matrix(~0 + area_id, mf$mf_model)
Z_age <- sparse.model.matrix(~0 + age_group, mf$mf_model)
Z_period <- sparse.model.matrix(~0 + period, mf$mf_model)

M_obs <- sparse.model.matrix(~0 + idx, mf$dist$obs) 
Z_tips <- sparse.model.matrix(~0 + tips_f, mf$dist$obs)
X_tips_dummy <- model.matrix(~0 + tips_dummy, mf$dist$obs)

M_obs_mics <- sparse.model.matrix(~0 + idx, mf$mics$obs) 
Z_tips_mics <- sparse.model.matrix(~0 + tips_f, mf$mics$obs)
X_tips_dummy_mics <- model.matrix(~0 + tips_dummy, mf$mics$obs)

R_spatial <- make_adjacency_matrix(iso3_current, areas_long, boundaries, exclude_districts = exc, level=2)
R_tips <- make_rw_structure_matrix(ncol(Z_tips), 1, TRUE)
R_age <- make_rw_structure_matrix(ncol(Z_age), 1, TRUE)
R_period <- make_rw_structure_matrix(ncol(Z_period), 2, TRUE)


dyn.unload(dynlib(here("tmb/fertility_tmb_dev")))
compile(here("tmb/fertility_tmb_dev.cpp"))               # Compile the C++ file
dyn.load(dynlib(here("tmb/fertility_tmb_dev")))

ar1_phi_age <- 0.99
ar1_phi_period <- 0.99

data <- list(M_obs = M_obs,
             X_tips_dummy = X_tips_dummy,
             Z_tips = Z_tips,
             Z_age = Z_age,
             Z_period = Z_period,
             Z_spatial = Z_spatial,
             # Z_interaction = sparse.model.matrix(~0 + id.interaction, mf$mf_model),
             Z_interaction1 = sparse.model.matrix(~0 + id.interaction1, mf$mf_model),
             Z_interaction2 = sparse.model.matrix(~0 + id.interaction2, mf$mf_model),
             Z_interaction3 = sparse.model.matrix(~0 + id.interaction3, mf$mf_model),
             R_tips = R_tips,
             R_age = R_age,
             R_period = R_period,
             R_spatial = R_spatial,
             ar1_phi_age = ar1_phi_age,
             ar1_phi_period = ar1_phi_period,
             log_offset = log(mf$dist$obs$pys),
             births_obs = mf$dist$obs$births,
             pop = mf$mf_model$population,
             A_out = mf$out$A_out,
             mics_toggle = mf$mics_toggle
)

par <- list(
  # beta_mf = rep(0, ncol(X_mf)),
  beta_0 = 0,
  beta_tips_dummy = rep(0, ncol(X_tips_dummy)),
  u_tips = rep(0, ncol(Z_tips)),
  u_age = rep(0, ncol(Z_age)),
  u_period = rep(0, ncol(Z_period)),
  u_spatial_str = rep(0, ncol(Z_spatial)),
  u_spatial_iid = rep(0, ncol(Z_spatial)),
  # eta = array(0, c(ncol(Z_spatial), ncol(Z_age), ncol(Z_period))),
  eta1 = array(0, c(ncol(Z_period), ncol(Z_age))),
  log_sigma_eta1 = log(2.5),
  eta2 = array(0, c(ncol(Z_spatial), ncol(Z_period))),
  log_sigma_eta2 = log(2.5),
  eta3 = array(0, c(ncol(Z_spatial), ncol(Z_age))),
  log_sigma_eta3 = log(2.5),
  log_sigma_rw_period = log(2.5),
  log_sigma_rw_age = log(2.5),
  log_sigma_rw_tips = log(2.5),
  log_sigma_spatial = log(2.5),
  logit_spatial_rho = 0
)

random <- c("beta_0", "beta_tips_dummy", "u_tips", "u_age", "u_period", "u_spatial_str", "u_spatial_iid", "eta1", "eta2", "eta3")

if(mf$mics_toggle) {
  data <- c(data, "M_obs_mics" = M_obs_mics,
            "X_tips_dummy_mics" = list(X_tips_dummy_mics),
            "Z_tips_mics" = Z_tips_mics,
            "births_obs_mics" = list(mf$mics$obs$births),
            "log_offset_mics" = list(log(mf$mics$obs$pys)),
            "A_mics" = mf$mics$A_mics)
  par <- c(par, "beta_tips_dummy_mics" = rep(0, ncol(X_tips_dummy_mics)))
  random = c(random, "beta_tips_dummy_mics")
}


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
                  random = random,
                  hessian = FALSE)

f <- nlminb(obj$par, obj$fn, obj$gr)
f$par.fixed <- f$par
f$par.full <- obj$env$last.par

fit <- c(f, obj = list(obj))
fit$sdreport <- sdreport(fit$obj, fit$par)

class(fit) <- "naomi_fit"  # this is hacky...
fit <- sample_tmb(fit)

qtls <- apply(fit$sample$lambda_out, 1, quantile, c(0.025, 0.5, 0.975))

asfr_plot <- readRDS("countries/ZWE/data/ZWE_asfr_plot.rds")
tfr_plot <- readRDS("countries/ZWE/data/ZWE_tfr_plot.rds")

moz_anc <- read_csv("~/Downloads/Moz work/11 files/ancnaomi111219clean_update_20.01.20.csv") %>%
  left_join(areas_long) %>%
  group_by(parent_area_id, year) %>%
  summarise(anc_clients = sum(anc_clients)) %>%
  rename(area_id = parent_area_id)

mf$out$mf_out %>%
  mutate(lower = qtls[1,],
         median = qtls[2,],
         upper = qtls[3,],
         source = "tmb") %>%
  type.convert() %>%
  # group_by(period, area_id) %>%
  # summarise(median = sum(median)) %>%
  left_join(areas_long) %>%
  filter(area_level == 1, age_group == "20-24") %>%
  ggplot(aes(x=period, y=median)) +
    geom_line() +
    geom_point(data=asfr_plot %>% filter(area_level == 1, age_group == "20-24"), aes(y=asfr, group=survtype, color=survtype)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.3) +
    # geom_point(data = moz_anc, aes(y=anc_clients, x=year))+
    facet_wrap(~area_name)

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
  group_by(area_id, period) %>%
  summarise(median = 5*sum(median)) %>%
  left_join(areas_long) %>%
  filter(area_level == 1) %>%
  # bind_rows(inla_res %>% mutate(source = "inla")) %>%
  ggplot(aes(x=period, y=median)) +
  geom_line() +
  geom_point(data=tfr_plot %>% filter(area_level == 1), aes(y=tfr, group=survtype, color=survtype)) +
  facet_wrap(~area_name)

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
