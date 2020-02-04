library(TMB)
library(tidyverse)
setwd("~/Documents/GitHub/subnat_fertility/")

asfr <- readRDS("~/Documents/GitHub/subnat_fertility/asfr_pred_subnat_15_2020_NEW.rds")[["MWI"]]

dat_m <- asfr %>%
  filter(!is.na(surveyid)) %>%
  select(period, agegr, births, pys) %>%
  group_by(period, agegr) %>%
  summarise_all(sum) %>%
  ungroup %>%
  select(-pys) %>%
  mutate(agegr = as.numeric(factor(agegr))) %>%
  pivot_wider(names_from = agegr, values_from = births) %>%
  select(-period) %>%
  as.matrix()

log_offset <- asfr %>%
  filter(!is.na(surveyid)) %>%
  select(period, agegr, births, pys) %>%
  group_by(period, agegr) %>%
  summarise_all(sum) %>%
  ungroup %>%
  select(-births) %>%
  mutate(agegr = as.numeric(factor(agegr)), pys = log(pys)) %>%
  pivot_wider(names_from = agegr, values_from = pys) %>%
  select(-period) %>%
  as.matrix()

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


compile("fertility_tmb.cpp")               # Compile the C++ file
dyn.load(dynlib("fertility_tmb"))  

#alpha_a = rep(1,7), beta_a = rep(1, 7)
f <-  MakeADFun(data=list(dat_m = dat_m, log_offset = log_offset), 
                parameters=list(alpha_a = rep(0,7), 
                                beta_a = rep(0, 7), 
                                rw_time = rep(0,22), 
                                rw_age = rep(0,7), 
                                rw_interaction = matrix(0, 22, 7),
                                log_prec_epsilon =1,
                                log_prec_rw_time = 1,
                                log_prec_rw_age = 1,
                                log_prec_rw_interaction = 1
                                ), 
                # map = list(), 
                DLL = "fertility_tmb")

f$fn(f$par)
f$report()

fit = nlminb(f$par,f$fn,f$gr, control = list(iter.max = 10000, eval.max = 10000))
print(fit)

int <- data.frame(agegr = 0:6, intercept = fit$par[1:7], slope = fit$par[8:14])

int %>%
  bind_rows %>%
  mutate(period = 1995:2016) %>%
  melt(id="period", variable.name="agegr") %>%
  ggplot(aes(x=period, y=value, color=agegr, group=agegr)) +
    geom_line()

crossing(period = 0:21, agegr = 0:6) %>%
  left_join(int) %>%
  left_join(data.frame(period = 0:21, rw_time = fit$par[names(fit$par) == "rw_time"])) %>%
  left_join(data.frame(agegr= 0:6, rw_age = fit$par[names(fit$par) == "rw_age"])) %>%
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
