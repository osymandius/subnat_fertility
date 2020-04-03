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

naomi_data_path <- "~/Documents/GitHub/naomi-data"
# naomi_data_path <- "~/GitHub/naomi-data"

source(here("R/inputs.R"))
source(here("R/fertility_funs.R"))

iso3_current <- "ZWE"

##sorry..
list2env(make_areas_population("ZWE", naomi_data_path, full = FALSE), globalenv())

asfr <- get_asfr_pred_df("ZWE", 2, project = FALSE)

# mics_data <- read_mics(iso3_current)
# mics_asfr <- Map(calc_asfr_mics, mics_data$wm, y=list(1),
#                  by = list(~area_id + survyear + surveyid + survtype),
#                  tips = list(c(0:5)),
#                  agegr= list(3:10*5),
#                  period = list(1995:2019),
#                  counts = TRUE,
#                  bhdata = mics_data$bh_df) %>%
#   bind_rows %>%
#   type.convert() %>%
#   filter(period <= survyear) %>%
#   rename(age_group = agegr)

mf <- make_model_frames(iso3_current, population, asfr, mics_asfr = NULL)

X_mf <- model.matrix(~1, mf$mf_model)

Z_spatial <- sparse.model.matrix(~0 + area_id, mf$mf_model)
Z_age <- sparse.model.matrix(~0 + age_group, mf$mf_model)
Z_period <- sparse.model.matrix(~0 + period, mf$mf_model)

M_obs <- sparse.model.matrix(~0 + idx, mf$dist$obs) 
Z_tips <- sparse.model.matrix(~0 + tips_f, mf$dist$obs)
X_tips_dummy <- model.matrix(~0 + tips_dummy, mf$dist$obs)

M_obs_mics <- sparse.model.matrix(~0 + idx, mf$mics$obs)
Z_tips_mics <- sparse.model.matrix(~0 + tips_f, mf$mics$obs)
X_tips_dummy_mics <- model.matrix(~0 + tips_dummy, mf$mics$obs)

R_spatial <- make_adjacency_matrix("ZWE", areas_long, boundaries, 2)
R_tips <- make_rw_structure_matrix(ncol(Z_tips), 1, TRUE)
R_age <- make_rw_structure_matrix(ncol(Z_age), 1, TRUE)
R_period <- make_rw_structure_matrix(ncol(Z_period), 2, TRUE)

dyn.unload(dynlib(here("tmb/fertility_tmb_dev")))
compile(here("tmb/fertility_tmb_dev.cpp"))               # Compile the C++ file
dyn.load(dynlib(here("tmb/fertility_tmb_dev")))

data <- list(X_mf = X_mf,
             M_obs = M_obs,
             
             # M_obs_mics = M_obs_mics,
             # X_tips_dummy_mics = X_tips_dummy_mics,
             # Z_tips_mics = Z_tips_mics,
             # births_obs_mics = mf$mics$obs$births,
             # log_offset_mics = log(mf$mics$obs$pys),
             # A_mics = mf$mics$A_mics,
             
             X_tips_dummy = X_tips_dummy,
             Z_tips = Z_tips,
             Z_age = Z_age,
             Z_period = Z_period,
             Z_spatial = Z_spatial,
             # Z_interaction = sparse.model.matrix(~0 + id.interaction, mf$mf_model),
             # Z_interaction1 = sparse.model.matrix(~0 + id.interaction1, mf$mf_model),
             # Z_interaction2 = sparse.model.matrix(~0 + id.interaction2, mf$mf_model),
             # Z_interaction3 = sparse.model.matrix(~0 + id.interaction3, mf$mf_model),
             R_tips = R_tips,
             R_age = R_age,
             R_period = R_period,
             R_spatial = R_spatial,
             log_offset = log(mf$dist$obs$pys),
             births_obs = mf$dist$obs$births,
             pop = mf$mf_model$population,
             A_out = mf$out$A_out
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
            # eta1 = array(0, c(ncol(Z_period), ncol(Z_age))),
            # eta2 = array(0, c(ncol(Z_spatial), ncol(Z_period))),
            # eta3 = array(0, c(ncol(Z_spatial), ncol(Z_age))),
            log_sigma_rw_period = log(2.5),
            log_sigma_rw_age = log(2.5),
            log_sigma_rw_tips = log(2.5),
            # log_sigma_eta1 = log(2.5),
            # log_prec_rw_period = 4,
            # log_prec_rw_age = 4,
            # log_prec_rw_tips = 4,
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
                random = c("beta_0", "beta_tips_dummy", "u_tips", "u_age", "u_period", "u_spatial_str", "u_spatial_iid"),
                hessian = FALSE)

f <- nlminb(obj$par, obj$fn, obj$gr)
f$par.fixed <- f$par
f$par.full <- obj$env$last.par

fit <- c(f, obj = list(obj))
fit$sdreport <- sdreport(fit$obj, fit$par)
fit <- sample_tmb_test(fit)

qtls <- apply(fit$sample$lambda_out, 1, quantile, c(0.025, 0.5, 0.975))

formula <- births ~
  f(id.age_group, model="rw1") +
  f(id.period, model="rw2") +
  f(id.district, model="bym2", graph=here("countries/ZWE/adj/ZWE_admin2.adj")) +
  f(id.tips, model="rw1") +
  tips_dummy

inla_mod2 <- inla(formula, family="poisson", data=asfr, E=pys,
                              control.family=list(link='log'),
                              control.predictor=list(compute=TRUE, link=1),
                              control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                              ## control.compute=list(config = TRUE, dic= TRUE, cpo=TRUE),
                              control.compute=list(config = TRUE),
                              verbose=TRUE)

inla_res <- get_mod_results_test(inla_mod, asfr)

mf$out$mf_out %>%
  mutate(lower = qtls[1,],
         median = qtls[2,],
         upper = qtls[3,],
         source = "tmb") %>%
  type.convert() %>%
  bind_rows(inla_res %>% mutate(source = "inla")) %>%
  left_join(areas_long) %>%
  filter(area_level ==2, age_group == "20-24") %>%
  ggplot(aes(x=period, y=median, group=source)) +
    geom_line(aes(color = source)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill=source), alpha=0.3) +
    facet_wrap(~area_id) +
    ylim(0,0.5)
 
# mics_plot <- Map(calc_asfr_mics, mics_data$wm, y=list(1),
#                  by = list(~area_id + survyear + surveyid + survtype),
#                  tips = list(c(0,5)),
#                  agegr= list(3:10*5),
#                  period = list(1995:2019),
#                  counts = TRUE,
#                  bhdata = mics_data$bh_df) %>%
#   bind_rows %>%
#   type.convert() %>%
#   filter(period <= survyear) %>%
#   rename(age_group = agegr)
# 
# mf_out %>% 
#   # left_join(areas_long) %>%
#   cbind(data.frame(val = f$report()$lambda_out)) %>%
#   type.convert() %>%
#   left_join(areas_long) %>%
#   filter(area_level ==0) %>%
#   ggplot(aes(x=period, y=val, group=source, color=source)) +
#   geom_line() +
#   # geom_point(data=mics_plot %>% bind_rows(asfr_admin1), aes(y=asfr, color=survtype)) +
#   facet_wrap(~area_id)
# 
# tfr %>%
#   type.convert() %>%
#   mutate(source = "model") %>%
#   left_join(areas_long) %>%
#   filter(area_level == 1) %>%
#   rename(val = median) %>%
#   bind_rows(zim_province_tfr) %>%
#   filter(period %in% 2000:2018) %>%
#   ggplot(aes(x=period, y=val)) +
#     geom_line(aes(group=source, color=source)) +
#     geom_ribbon(aes(ymin = lower, ymax = upper, fill=source), alpha=0.3) +
#     # geom_point(data=mics_plot %>% bind_rows(asfr_admin1), aes(y=asfr, color=survtype)) +
#     facet_wrap(~area_id) +
#     ylim(0, 20)


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


# join_rate <- crossing(area_id = unique(join_out$area_id),
#                       age_group_join,
#                       period = unique(mf$mf_model$period)) %>%
#   mutate(idx = group_indices(., area_id, age_group_out, period)) %>%
#   left_join(join_out %>% select(area_id, age_group, period, out_idx) %>% distinct %>% rename(idx_col = out_idx)) %>%
#   mutate(x=ifelse(age_group_out == "15-49", 5, 1))
# 
# A_rate <- spMatrix(nrow(join_rate), nrow(A_out), join_rate$idx, as.integer(join_rate$idx_col), join_rate$x)

hyper <- inla_mod$internal.marginals.hyperpar
hyper_default <- readRDS("~/Downloads/inla_full_internal.marginals.hyperpar.rds")

inla_transformed_hyper <- list(
  age = data.frame(x=exp(hyper[[1]][,1]), y=hyper[[1]][,2]),
  period = data.frame(x=exp(hyper[[2]][,1]), y=hyper[[2]][,2]),
  district = data.frame(x=exp(hyper[[3]][,1]), y=hyper[[3]][,2]),
  district_phi = data.frame(x=exp(hyper[[4]][,1])/(1+exp(hyper[[4]][,1])), y=hyper[[4]][,2]),
  tips = data.frame(x=exp(hyper[[5]][,1]), y=hyper[[5]][,2])
)

inla_untransformed_hyper <- list(
  age = data.frame(x=hyper[[1]][,1], y=hyper[[1]][,2]),
  period = data.frame(x=hyper[[2]][,1], y=hyper[[2]][,2]),
  district = data.frame(x=hyper[[3]][,1], y=hyper[[3]][,2]),
  district_phi = data.frame(x=hyper[[4]][,1], y=hyper[[4]][,2]),
  tips = data.frame(x=hyper[[5]][,1], y=hyper[[5]][,2])
)

inla_untransformed_hyper_default <- list(
  age = data.frame(x=hyper_default[[1]][,1], y=hyper_default[[1]][,2]),
  period = data.frame(x=hyper_default[[2]][,1], y=hyper_default[[2]][,2]),
  district = data.frame(x=hyper_default[[3]][,1], y=hyper_default[[3]][,2]),
  district_phi = data.frame(x=hyper_default[[4]][,1], y=hyper_default[[4]][,2])
)

hist(fit$sample$log_prec_rw_age, breaks = 100)
hist(fit$sample$log_prec_rw_period, breaks = 100)

hist(fit$sample$log_tau_rw_age, breaks = 100)
hist(fit$sample$log_tau_rw_period, breaks = 100)

gridExtra::grid.arrange(
data.frame(val=fit$sample$log_tau2_rw_age, source = "tmb log tau") %>%
  ggplot(aes(color=source)) +
    geom_density(aes(x=val)) +
    geom_point(data=inla_untransformed_hyper$age %>% mutate(source = "INLA log precision"), aes(x=x, y=y)) +
    labs(title="Age"),

data.frame(val=fit$sample$log_tau2_rw_period, source = "tmb log tau") %>%
  ggplot(aes(color=source)) +
  geom_density(aes(x=val)) +
  geom_point(data=inla_untransformed_hyper$period %>% mutate(source = "INLA log precision"), aes(x=x, y=y)) +
  labs(title="Period"),

data.frame(val=fit$sample$log_tau2_spatial, source = "tmb log sigma") %>%
  ggplot(aes(color=source)) +
  geom_density(aes(x=val)) +
  geom_point(data=inla_untransformed_hyper$district %>% mutate(source = "INLA log precision"), aes(x=x, y=y)) +
  labs(title="Spatial"),

data.frame(val=fit$sample$logit_spatial_rho, source = "tmb logit rho") %>%
  ggplot(aes(color=source)) +
  geom_density(aes(x=val)) +
  geom_point(data=inla_untransformed_hyper$district_phi %>% mutate(source = "INLA logit phi"), aes(x=x, y=y)) +
  labs(title="Rho (tmb), phi (INLA"),

data.frame(val=fit$sample$log_tau2_rw_tips, source = "tmb log tau") %>%
  ggplot(aes(color=source)) +
  geom_density(aes(x=val)) +
  geom_point(data=inla_untransformed_hyper$tips %>% mutate(source = "INLA log precison"), aes(x=x, y=y)) +
  labs(title="TIPS")
)

fit$sample$log_tau_rw_period

quantile(fit$sample$log_sigma_rw_age, c(0.5))
quantile(fit$sample$log_sigma_rw_period, c(0.5))

quantile(fit$sample$log_tau_rw_age, c(0.5))
quantile(fit$sample$log_tau_rw_period, c(0.5))

plot(inla_transformed_hyper$age, xlim=c(0,10), main="prec age")
plot(inla_transformed_hyper$period, xlim=c(0,1000), main="prec time")
plot(inla_transformed_hyper$district, xlim=c(0,100), main="prec space")
plot(inla_transformed_hyper$district_phi, xlim=c(0.8,1.2), main="phi")

gridExtra::grid.arrange(plot(inla_untransformed_hyper$age, xlim=c(-1,3), main="Log prec age"),
plot(inla_untransformed_hyper$period, xlim=c(0,10), main="Log prec time"),
plot(inla_untransformed_hyper$district, xlim=c(2,6), main="Log prec space"),
plot(inla_untransformed_hyper$district_phi, xlim=c(2,6), main="Logit phi"),
)

