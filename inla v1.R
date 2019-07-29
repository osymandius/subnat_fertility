library(tidyverse)
library(rdhs)
library(demogsurv)
library(INLA)
library(reshape2)

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

tfr <- Map(calc_tfr, ir,
           by = list(~surveyid + country + survyear + v025),
           tips = tips_surv,
           period = list(1995:2017))

tfr <- tfr %>%
  bind_rows %>%
  type.convert %>% 
  mutate(v025 = factor(v025, levels= c(1,2), labels = c("Urban", "Rural"))) %>%
  filter(period <= survyear) %>%
  mutate(lower = tfr - qnorm(0.975) * se_tfr,
         upper = tfr + qnorm(0.975) * se_tfr)

asfr <- Map(calc_asfr, ir,
              by = list(~surveyid + country + survyear+ v025),
              tips = tips_surv,
              period = list(1995:2017),
              counts = TRUE)

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
  


asfr2 <- asfr %>%
    bind_rows %>%
    type.convert %>%
    mutate(v025 = factor(v025, levels= c(1,2), labels = c("Urban", "Rural"))) %>%
    filter(period <= survyear) %>%
    mutate(births = round(births),
           pys = round(pys),
           sum_var = ifelse(agegr=="40-44" | agegr=="45-49", TRUE, FALSE)) %>%
    group_by(surveyid, country,  survyear,tips,  v025, period, sum_var) %>%
    mutate(births = ifelse(sum_var, sum(births), births),
           pys = ifelse(sum_var, sum(pys), pys),
           asfr = births/pys,
           agegr = plyr::revalue(agegr, c("40-44" = "40-49")),
           agegr = factor(agegr, levels=c("15-19", "20-24", "25-29", "30-34", "35-39", "40-49"))
    ) %>%
    droplevels() %>%
    ungroup() %>%
    select(-sum_var) %>%
    distinct(surveyid, country, survyear, tips,v025, period, births, pys, asfr, .keep_all=TRUE)

asfr1$agegr2 <- asfr1$agegr
asfr1$agegr3 <- asfr1$agegr

lme_mod <- glmer(births ~ agegr + period + v025 + (1 | agegr), family = poisson, data = asfr2, offset = log(pys))
summary(lme_mod)

vc0<-as.data.frame(VarCorr(lme_mod))
vc0$init<-log(1/vc0$vcov); vc0

 f(id.agegr, model="rw1") + f(id.period, model="rw1") + f(id.agegr.period, model="iid")
formula.2 <- births ~ v025 + f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr.period, model="iid")

formula.5 <- births ~ v025 + f(id.period, model="rw2") + f(id.agegr.period, model="iid") + f(id.agegr, model="rw1")

mod <- inla(formula.5, family="poisson", data=asfr_pred, offset=log(pys),
            control.family=list(link='log'),
            control.predictor=list(compute=TRUE, link=1),
            control.compute=list(config = TRUE, dic= TRUE, cpo=TRUE))
            # control.fixed=list(mean=0, prec=0.00001,
                              # mean.intercept=0, prec.intercept=0.00001))

mod_rw2 <- inla(formula.2, family="poisson", data=asfr1, offset=log(pys),
                control.family=list(link='log'),
                control.predictor=list(compute=TRUE, link=1),
                control.compute=list(config = TRUE, dic= TRUE, cpo=TRUE))
# control.fixed=list(mean=0, prec=0.00001,
# mean.intercept=0, prec.intercept=0.00001))

summary(mod)


asfr_pred %>%
  filter(id<309) %>%
  select("v025", "agegr", "period","pys", "id") %>%
  left_join(exp(mod$summary.fitted.values[1:308,]) %>%
              mutate(id=1:nrow(.)), by="id") %>%
  arrange(period, v025, agegr) %>%
  ggplot(aes(x=period, y=mean, ymin=`0.025quant`, max=`0.975quant`, group=agegr))+
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
  
