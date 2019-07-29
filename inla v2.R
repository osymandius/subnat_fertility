library(tidyverse)
library(rdhs)
library(demogsurv)
library(INLA)
library(reshape2)
library(survival)

setwd("~/Documents/GitHub/subnat_fertility")
load("~/Documents/GitHub/subnat_fertility/asfr.Rda")

set_rdhs_config(email="o.stevens@imperial.ac.uk", project="Subnational fertility", config_path = "~/.rdhs.json")

##+ datasets
surveys <- dhs_surveys(countryIds = "MW", surveyYearStart=1995)
ird <- dhs_datasets(fileType = "IR", fileFormat = "flat", surveyIds = surveys$SurveyId)
ird$path <- unlist(get_datasets(ird))


ir <- lapply(ird$path, readRDS) %>%
  Map(data.frame,
      surveyid = surveys$SurveyId,
      country = surveys$CountryName,
      survyear = surveys$SurveyYear,
      survtype = surveys$SurveyType,
      .,
      stringsAsFactors = FALSE)

tips_surv <- list("DHS" = c(0, 7), "MIS" = c(0, 5))[surveys$SurveyType]
# 
# tfr <- Map(calc_tfr, ir,
#            by = list(~surveyid + country + survyear + v025),
#            tips = tips_surv,
#            period = list(1995:2017))
# 
# tfr <- tfr %>%
#   bind_rows %>%
#   type.convert %>% 
#   mutate(v025 = factor(v025, levels= c(1,2), labels = c("Urban", "Rural"))) %>%
#   filter(period <= survyear) %>%
#   mutate(lower = tfr - qnorm(0.975) * se_tfr,
#          upper = tfr + qnorm(0.975) * se_tfr)
# 
# asfr <- Map(calc_asfr, ir,
#               by = list(~surveyid + country + survyear+ v025),
#               tips = tips_surv,
#               agegr= 15:50,
#               period = list(1995:2017),
#               counts = TRUE)
debugonce(demog_pyears)
asfr <- Map(calc_asfr, ir,
            by = list(~surveyid + country + survyear+ v025),
            tips = tips_surv,
            agegr= 15:50,
            period = list(2000:2002),
            counts = TRUE)

debugonce(calc_asfr)

asfr1 <- asfr %>%
  bind_rows %>%
  type.convert %>%
  mutate(v025 = factor(v025, levels= c(1,2), labels = c("Urban", "Rural"))) %>%
  filter(period <= survyear) %>%
  mutate(id.period = group_indices(., period),
         id.period2 = id.period,
         id.agegr = group_indices(., agegr),
         id.agegr2 = id.agegr,
         id.agegr.period = group_indices(., period, agegr),
         id = 1:nrow(.))

asfr_pred <- asfr %>%
  bind_rows %>%
  type.convert %>%
  mutate(v025 = factor(v025, levels= c(1,2), labels = c("Urban", "Rural"))) %>%
  filter(period <= survyear) %>%
  mutate(births = NA_integer_,
         pys = 1,
         asfr = NA_integer_,
         se_asfr = NA_integer_)%>%
  select("v025", "agegr", "period", "pys") %>%
  group_by_all() %>%
  summarise() %>%
  ungroup() %>%
  bind_rows(asfr %>%
              bind_rows %>%
              type.convert %>%
              mutate(v025 = factor(v025, levels= c(1,2), labels = c("Urban", "Rural"))) %>%
              filter(period <= survyear)) %>%
  mutate(id.period = group_indices(., period),
         id.period2 = id.period,
         id.agegr = group_indices(., agegr),
         id.agegr2 = id.agegr,
         id.agegr.period = group_indices(., period, agegr),
         id = 1:nrow(.))

# next things to do:
#   - Make the `agegr x period` term a RW2 on period (instead of iid). I think the way to do this with `f(id.agegr.period, model = "rw2", group = agegr)` or something like that.
# - Stratify the data by single-year of age and make a RW2 over age
# - Fit an interaction of RW2 on age and RW2 on period
# - Get the district-level data sorted out and start doing space

formula.1 <- births ~ v025 + f(id.agegr.period, model="iid") + f(id.agegr, model="rw1", group=id.period, control.group=list(model="rw2"))

mod <- inla(formula.1, family="poisson", data=asfr1, E=pys,
            control.family=list(link='log'),
            control.predictor=list(compute=TRUE, link=1),
            control.compute=list(config = TRUE, dic= TRUE, cpo=TRUE))
            # control.fixed=list(mean=0, prec=0.00001,
                              # mean.intercept=0, prec.intercept=0.00001))


mod$marginals.random$id.agegr %>%
  lapply(t_emarg)

asfr_pred %>%
  filter(id<309) %>%
  select("v025", "agegr", "period","pys", "id") %>%
  left_join(mod$summary.fitted.values[1:308, ] %>%
              mutate(id = 1:308), by="id") %>%
  arrange(period, v025, agegr) %>%
  ggplot(aes(x=period, y=`0.5quant`, ymin=`0.025quant`, ymax=`0.975quant`, group=agegr))+
  geom_line(aes(color=agegr)) +
  facet_wrap(~v025)

t_emarg <- function(x){
  inla.emarginal(exp, x)
}

t_tmarg <- function(x) {
  inla.tmarginal(exp, x)
}

t_qmarg <- function(x){
  inla.qmarginal(c(0.025, 0.975), x)
}

### Age group iid random effect

mod$marginals.random$id.agegr %>%
  lapply(t_tmarg) %>%
  melt() %>%
  mutate(Var2 = as.character(Var2)) %>%
  pivot_wider(names_from=Var2) %>%
  mutate(agegr = factor(as.numeric(factor(L1)), labels=unique(asfr1$agegr))) %>%
  ggplot(aes(x=x, y=y, group=agegr)) +
    geom_line(aes(color=agegr)) +
    xlab("ASFR") +
    ylab("")+
    xlim(0,3)

mod$marginals.random$id.agegr %>%
  lapply(t_emarg) %>%
  melt() %>%
  mutate(agegr = factor(as.numeric(factor(L1)), labels=unique(asfr1$agegr))) %>%
  left_join(mod$marginals.random$id.agegr %>%
              lapply(t_qmarg) %>%
              lapply(exp) %>%
              bind_rows() %>%
              t() %>%
              `colnames<-`(value=c("lower", "upper")) %>%
              data.frame() %>%
              mutate(agegr = unique(asfr1$agegr)), by = "agegr") %>%
  ggplot(aes(x=agegr, y=value, ymin=lower, ymax=upper)) +
  geom_point() +
  geom_errorbar(width=0.3) +
  labs(title="Age group iid random effect", y="ASFR relative to 15-19", x="Age group")

### Attempt at agegr.period interaction - what is this doing? No idea.

mod$marginals.random$id.agegr.period %>%
  lapply(t_emarg) %>%
  melt() %>%
  mutate(id=1:nrow(.)) %>%
  left_join(asfr1 %>%
              group_by(period, agegr, id.agegr.period) %>%
              summarise() %>%
              ungroup(), by = c("id" = "id.agegr.period")) %>%
  ggplot(aes(x=period, y=value, group=agegr)) +
  geom_line(aes(color=agegr))

#### Period random effect
  
rw2_mod <- mod$marginals.random$id.period %>%
  lapply(t_emarg) %>%
  melt() %>%
  mutate(period = 1995:2016) %>%
  left_join(mod$marginals.random$id.period %>%
              lapply(t_qmarg) %>%
              lapply(exp) %>%
              bind_rows() %>%
              t() %>%
              `colnames<-`(value=c("lower", "upper")) %>%
              data.frame() %>%
              mutate(period = 1995:2016), by = "period") %>%
  ggplot(aes(x=period, y=value, ymin=lower, ymax=upper)) +
  geom_line() +
  geom_ribbon(alpha=0.5) +
  xlab("") +
  ylab("")


rw1_mod <- mod$marginals.random$id.period %>%
  lapply(t_emarg) %>%
  melt() %>%
  mutate(period = 1995:2016) %>%
  left_join(mod$marginals.random$id.period %>%
              lapply(t_qmarg) %>%
              lapply(exp) %>%
              bind_rows() %>%
              t() %>%
              `colnames<-`(value=c("lower", "upper")) %>%
              data.frame() %>%
              mutate(period = 1995:2016), by = "period") %>%
  ggplot(aes(x=period, y=value, ymin=lower, ymax=upper)) +
  geom_line() +
  geom_ribbon(alpha=0.5) +
  xlab("") +
  ylab("")

gridExtra::grid.arrange(rw1_mod, rw2_mod)

asfr1 %>%
  group_by(period, agegr, id.agegr.period) %>%
  summarise()

  mutate(Var2 = as.character(Var2)) %>%
  pivot_wider(names_from=Var2) %>%
  mutate(period = as.numeric(factor(L1))+1994) %>%
  ggplot(aes(x=x, y=y, group=agegr)) +
    geom_line(aes(color=agegr)) +
    xlab("ASFR") +
    ylab("")+
    xlim(0,3)



lapply(tmp, t_tmarg)

View(mod$marginals.random$id.agegr)

lapply(mod$marginals.random$id.agegr, inla.tmarginal(exp)) %>%
       melt() %>%
       pivot_wider(names_from=Var2) %>%
       mutate(agegr = factor(as.numeric(factor(L1)), labels=unique(asfr1$agegr))) %>%
  ggplot(aes(x=x, y=y, group=agegr)) +
  geom_line(aes(color=agegr)) +
  xlab("ASFR") +
  ylab("")+
  xlim(0,3)

lapply(x, function(y){exp(inla.hpdmarginal(0.95, y))}) %>%
  melt() %>%
  pivot_wider(names_from=Var2) %>%
  mutate(L1 = factor(as.numeric(factor(L1)), labels=unique(asfr1$agegr)))

exp.b0.mean <- inla.emarginal(exp, mod$marginals.fixed[[1]])
exp.b0.CI <- inla.qmarginal(c(0.025, 0.975), inla.tmarginal(exp, mod$marginals.fixed[[1]]))
exp.b0.CI

 mod$summary.fitted.values

n.samples = 1000
samples = inla.posterior.sample(n.samples, result = mod)

mod$misc$configs$contents


round(mod$summary.hyperpar, 3)

  mod2 <- glmer(births ~ agegr + period + v025 + (1 |agegr) , family = poisson, data = asfr2, offset = log(pys))
  mod2.int <- glmer(births ~ agegr*period + v025 + (1 + period | agegr), family = poisson, data = asfr2, offset = log(pys))
  janova(mod2.int, mod2)
  
  View(coef(mod2.inter))
  
  df_pred2 <- crossing(period = 1995:2017,
                       agegr = levels(asfr$agegr),
                       pys = 1,
                       v025 = levels(asfr$v025))
  
  pred2 <- predict(mod2.int, newdata = df_pred2)
  
  df_pred2 <- data.frame(df_pred2, pred2) %>%
    mutate(asfr = exp(pred2))
  
  ggplot(asfr2, aes(period, asfr, color = surveyid)) +
    geom_point() +
    geom_line() +
    #geom_ribbon(aes(ymin = lower, ymax = upper), data = df_pred, color = NA, fill = "darkred", alpha = 0.3) +
    geom_line(data = df_pred2, aes(y=asfr), color = "darkred") +
    # facet_wrap(~v025)
    facet_grid(v025~agegr)
  
