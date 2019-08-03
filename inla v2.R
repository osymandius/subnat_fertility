library(tidyverse)
library(rdhs)
library(demogsurv)
library(INLA)
library(reshape2)
library(survival)
library(sf)
library(spdep)

setwd("~/Documents/GitHub/subnat_fertility")
# load("~/Documents/GitHub/subnat_fertility/asfr_singleage.Rda")
load("~/Documents/GitHub/naomi-dev/data/shapefile.rda")

calc_asfr1 <- function(data,
                      by = NULL,
                      agegr = NULL,
                      period = NULL,
                      cohort = NULL,
                      tips = c(0, 3),
                      clusters=~v021,
                      strata=~v024+v025,
                      id="caseid",
                      dob="v011",
                      intv = "v008",
                      weight= "v005",
                      varmethod = "none",
                      bvars = grep("^b3\\_[0-9]*", names(data), value=TRUE),
                      birth_displace = 1e-6,
                      origin=1900,
                      scale=12,
                      bhdata = NULL,
                      counts=FALSE,
                      clustcounts = FALSE){
  
  data$id <- data[[id]]
  data$dob <- data[[dob]]
  data$intv <- data[[intv]]
  data$weights <- data[[weight]] / mean(data[[weight]])
  
  if(is.null(by))
    by <- ~1
  
  vars <- unique(unlist(lapply(c(by, strata, clusters), all.vars)))
  f <- formula(paste("~", paste(vars, collapse = "+")))
  mf <- model.frame(formula = f, data = data, na.action = na.pass,
                    id = id, weights = weights, dob = dob, intv = intv)
  
  if(is.null(bhdata)) {
    births <- reshape(model.frame(paste("~", paste(bvars, collapse="+")),
                                  data, na.action=na.pass, id=id),
                      idvar="(id)", timevar="bidx",
                      varying=bvars, v.names="bcmc", direction="long")
  } else {
    if(length(bvars) > 1)
      stop("If `bhdata' is provided, bvars must provide a single variable name (length(bvars) = 1)")
    
    bhdata$id <- bhdata[[id]]
    bhdata$bcmc <- bhdata[[bvars]]
    births <- model.frame(~bcmc, data = bhdata, id = id)
    births$bidx <- ave(births$bcmc, births$`(id)`, FUN = seq_along)
  }
  births <- births[!is.na(births$bcmc), ]
  births$bcmc <- births$bcmc + births$bidx * birth_displace
  
  epis <- tmerge(mf, mf, id=`(id)`, tstart=`(dob)`, tstop=`(intv)`)
  epis <- tmerge(epis, births, id=`(id)`, birth = event(bcmc))
  
  aggr <- demog_pyears(f, epis, period=period, agegr=agegr, cohort=cohort, tips=tips,
                       event="birth", weights="(weights)", origin=origin, scale=scale)$data
  
  ## construct interaction of all factor levels that appear
  byvar <- intersect(c(all.vars(by), "agegr", "period", "cohort", "tips"),
                     names(aggr))
  aggr$byf <- interaction(aggr[byvar], drop=TRUE)
  
  ## prediction for all factor levels that appear
  pred <- data.frame(aggr[c(byvar, "byf")])[!duplicated(aggr$byf),]
  pred <- pred[order(pred$byf), ]
  
  if(counts || varmethod == "none"){
    mc <- model.matrix(~-1+byf, aggr)
    clong <- aggr[c("event", "pyears")]
    pred[c("births", "pys")] <- t(mc) %*% as.matrix(clong)
  }
  
  if(varmethod == "none") {
    
    pred$asfr <- pred$births / pred$pys
    pred$byf <- NULL
    if(!counts)
      pred[c("births", "pys")] <- NULL
    
  } else if(varmethod == "lin") {
    
    des <- survey::svydesign(ids=clusters, strata=strata, data=aggr, weights=~1)
    class(des) <- c("svypyears", class(des))
    
    ## fit model
    f <- if(length(levels(aggr$byf)) == 1)
      event ~ offset(log(pyears))
    else
      event ~ -1 + byf + offset(log(pyears))
    
    mod <- survey::svyglm(f, des, family=quasipoisson)
    
    ## prediction for all factor levels that appear
    pred$pyears <- 1
    
    asfr <- predict(mod, pred, type="response", vcov=TRUE)
    v <- vcov(asfr)
    dimnames(v) <- list(pred$byf, pred$byf)
    
    pred$asfr <- as.numeric(asfr)
    pred$se_asfr <- sqrt(diag(v))
    pred[c("byf", "pyears")] <- NULL
    attr(pred, "var") <- v
  } else if(varmethod %in% c("jkn", "jk1")) {
    
    ## Convert to array with events and PYs for each cluster
    ## reshape2::acast is MUCH faster than stats::reshape
    events_clust <- reshape2::acast(aggr, update(clusters, byf ~ .), value.var="event")
    pyears_clust <- reshape2::acast(aggr, update(clusters, byf ~ .), value.var="pyears")
    
    if(varmethod == "jkn"){
      aggr$strataid <- as.integer(interaction(aggr[all.vars(strata)], drop=TRUE))
      strataid <- drop(reshape2::acast(unique(aggr[c(all.vars(clusters), "strataid")]),
                                       update(clusters,  1 ~ .), value.var="strataid"))
    } else
      strataid <- NULL
    
    estdf <- jackknife(events_clust, pyears_clust, strataid)
    
    pred$asfr <- estdf$est
    pred$se_asfr <- estdf$se
    attr(pred, "var") <- vcov(estdf)
    pred$byf <- NULL
    if(clustcounts){
      attr(pred, "events_clust") <- events_clust
      attr(pred, "pyears_clust") <- pyears_clust
      attr(pred, "strataid") <- strataid
    }
  } else
    stop(paste0("varmethod = \"", varmethod, "\" is not recognized."))
  
  
  rownames(pred) <- NULL
  
  return(pred)
}

set_rdhs_config(email="o.stevens@imperial.ac.uk", project="Subnational fertility", config_path = "~/.rdhs.json")

##+ datasets
surveys <- dhs_surveys(countryIds = c("MW", "ZW", "KE"), surveyYearStart=1995)
ird <- dhs_datasets(fileType = "IR", fileFormat = "flat", surveyIds = surveys$SurveyId)
ird$path <- unlist(get_datasets(ird))
# 
# 
ir <- lapply(ird$path, readRDS) %>%
  Map(data.frame,
      surveyid = surveys$SurveyId,
      country = surveys$CountryName,
      survyear = surveys$SurveyYear,
      survtype = surveys$SurveyType,
      .,
      stringsAsFactors = FALSE)


geo_2000_dhs <- rdhs::get_datasets("MWGE43FL.zip")[[1]] %>% readRDS %>% st_as_sf
geo_2000_dhs <- st_join(geo_2000_dhs, sh32, join = st_intersects, suffix = c("", "_sh32")) %>%
  mutate(district_32 = ifelse(district %in% c("Blantyre", "Lilongwe", "Mzimba", "Zomba"), 
                              ifelse(URBAN_RURA == "U", paste(district, "City"), as.character(district)), as.character(district)))

ir[["MW2000DHS"]] <- ir[["MW2000DHS"]] %>%
  left_join(geo_2000_dhs, by = c("v001" = "DHSCLUST"))

ir[["MW2004DHS"]] <- ir[["MW2004DHS"]] %>%
  mutate(district = as_factor(sdist2) %>%
           sub("nkhota kota", "nkhotakota", .) %>%
           str_to_title)


ir[["MW2010DHS"]] <- ir[["MW2010DHS"]] %>%
  left_join(
    data.frame("distcode" = attr(ir[["MW2010DHS"]]$sdistrict, "labels"), "district" = names(attr(ir[["MW2010DHS"]]$sdistrict, "labels")) %>% 
                 sub("nkhota kota", "nkhotakota", .) %>% 
                 sub("nkhatabay", "nkhata bay", .) %>% 
                 str_to_title()
    ),
    by = c("sdistrict" = "distcode")
  ) %>%
  mutate(district_32 = ifelse(district %in% c("Blantyre", "Lilongwe", "Mzimba", "Zomba"), 
                              ifelse(v025 == 1, paste(district, "City"), as.character(district)), as.character(district)))

ir[["MW2015DHS"]] <- ir[["MW2015DHS"]] %>%
  mutate(district = sub("(.*) - (.*)", "\\1", as_factor(v022)),
         district = fct_recode(district,
                               "nkhotakota" = "nkhota kota",
                               "nkhata bay" = "nkhatabay") %>%
           str_to_title)

geo_2012_mis<- rdhs::get_datasets("MWGE6AFL.ZIP")[[1]] %>% readRDS %>% st_as_sf
geo_2012_mis <- st_join(geo_2012_mis, sh32, join = st_intersects, suffix = c("", "_sh32")) %>%
  type.convert() ## Check you haven't fucked this factor conversion up like last time.

geo_2012_mis %>%
  mutate(district_32 = ifelse(district %in% c("Blantyre", "Lilongwe", "Mzimba", "Zomba"), 
                              ifelse(URBAN_RURA == "U", paste(district, "City"), as.character(district)), as.character(district))) %>%
  dplyr::select(district, district32)

ir[["MW2012MIS"]] <- ir[["MW2012MIS"]] %>%
  left_join(geo_2012_mis, by = c("v001" = "DHSCLUST"))


geo_2014_mis<- rdhs::get_datasets("MWGE71FL.ZIP")[[1]] %>% readRDS %>% st_as_sf
geo_2014_mis <- st_join(geo_2014_mis, sh32, join = st_intersects, suffix = c("", "_sh32"))

ir[["MW2014MIS"]] <- ir[["MW2014MIS"]] %>%
  left_join(geo_2014_mis, by = c("v001" = "DHSCLUST"))


geo_2017_mis<- rdhs::get_datasets("MWGE7IFL.ZIP")[[1]] %>% readRDS %>% st_as_sf
geo_2017_mis <- st_join(geo_2017_mis, sh32, join = st_intersects, suffix = c("", "_sh32"))

ir[["MW2017MIS"]] <- ir[["MW2017MIS"]] %>%
  left_join(geo_2017_mis, by = c("v001" = "DHSCLUST"))

# foo <- ir %>%
#   bind_rows() %>%
#   mutate(district_fac = factor(district)) %>%
#   group_split(surveyid) %>%
#   Map(calc_asfr1, .,
#       by = list(~surveyid + country + survyear + district_fac),
#       tips = tips_surv,
#       agegr= list(15:50),
#       period = list(1995:2017),
#       counts = TRUE)
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
tips_surv <- list("DHS" = c(0, 7), "MIS" = c(0, 5))[surveys$SurveyType]

asfr <- Map(calc_asfr1, ir,
            by = list(~surveyid + country + survyear + v024),
            tips = tips_surv,
            agegr= list(15:50),
            period = list(1995:2017),
            counts = TRUE)

region_names <- lapply(ir, function(x) {x %>%
    dplyr::select(surveyid, v024) %>%
    unique() %>%
    arrange(v024) %>%
    mutate(region_name = names(attr(x$v024, "labels")) %>% str_to_title())
}) %>%
  bind_rows()

asfr_int <- lapply(asfr, function(x) {
  x %>%
    type.convert %>%
    left_join(region_names, by=c("v024", "surveyid")) %>%
    mutate(region_name = fct_recode(region_name, 
                                    "Central" = "Central Region",
                                    "South" = "Southern",
                                    "South" = "Southern Region",
                                    "North" = "Northern",
                                    "North" = "Northern Region",
                                    "Matabeleland North" = "Matebeleland North",
                                    "Matabeleland South" = "Matebeleland South",
                                    "North Eastern" = "Northeastern"
          )
    )
})

asfr1 <- asfr_int %>%
  bind_rows() %>%
  group_by(country) %>%
  group_split() %>%
  lapply(function(x) {
    x %>%
      mutate(id.region = group_indices(., region_name))
  }) %>%
  bind_rows() %>%
  filter(period <= survyear) %>%
  mutate(id.period = group_indices(., period),
         id.period2 = id.period,
         id.agegr = group_indices(., agegr),
         id.agegr2 = id.agegr,
         id.agegr.period = group_indices(., period, agegr),
         # id.district = group_indices(., district),
         id = 1:nrow(.))


#   
# View(lapply(foo, function(x){x %>%
#   mutate(new_v024_fac = as.numeric(new_v024))}
# ) %>%
#   bind_rows() %>%
#   select(country, v024, new_v024, new_v024_fac) %>%
#   arrange(country, new_v024_fac) %>%
#   unique()
# )
# 
# foo[[1]] %>%
#   mutate(new_v024_fac = as.numeric(region_name)) %>%
#   dplyr::select(v024, region_name, new_v024_fac)

  



asfr_pred <- crossing(period = asfr1$period, agegr = asfr1$agegr, v024 = asfr1$v024, pys=1) %>%
  bind_rows(asfr %>%
              bind_rows %>%
              type.convert %>%
              mutate(v025 = factor(v025, levels= c(1,2), labels = c("Urban", "Rural"))) %>%
              mutate(v024 = factor(v024, levels= c(1,2,3), labels = c("North", "Central", "South"))) %>%
              filter(period <= survyear)) %>%
  mutate(id.period = group_indices(., period),
         id.period2 = id.period,
         id.agegr = group_indices(., agegr),
         id.agegr2 = id.agegr,
         id.agegr.period = group_indices(., period, agegr),
         id.district = group_indices(., district),
         id.agegr.period.district = group_indices(., agegr, period, district),
         id = 1:nrow(.))

# asfr_pred <- asfr %>%
#   bind_rows %>%
#   type.convert %>%
#   mutate(v025 = factor(v025, levels= c(1,2), labels = c("Urban", "Rural"))) %>%
#   filter(period <= survyear) %>%
#   mutate(births = NA_integer_,
#          pys = 1,
#          asfr = NA_integer_,
#          se_asfr = NA_integer_) %>%
#   select("v025", "agegr", "period", "pys") %>%
#   group_by_all() %>%
#   summarise() %>%
#   ungroup() %>%
#   bind_rows(asfr %>%
#               bind_rows %>%
#               type.convert %>%
#               mutate(v025 = factor(v025, levels= c(1,2), labels = c("Urban", "Rural"))) %>%
#               filter(period <= survyear)) %>%
#   mutate(id.period = group_indices(., period),
#          id.period2 = id.period,
#          id.agegr = group_indices(., agegr),
#          id.agegr2 = id.agegr,
#          id.agegr.period = group_indices(., period, agegr),
#          id = 1:nrow(.))

# next things to do:
#   - Make the `agegr x period` term a RW2 on period (instead of iid). I think the way to do this with `f(id.agegr.period, model = "rw2", group = agegr)` or something like that.
# - Stratify the data by single-year of age and make a RW2 over age
# - Fit an interaction of RW2 on age and RW2 on period
# - Get the district-level data sorted out and start doing space

test <- poly2nb(sh32 %>% filter(!area_id %in% c(7, 17, 31, 32)), row.names = sh32$area_id, queen = TRUE)
nb2INLA("mw.adj", test)

formula.0 <- births ~ v025 + f(id.agegr.period, model="rw2") + f(id.agegr, model="rw2", group=id.period, control.group=list(model="rw2"))
formula.1 <- births ~ v025 + f(id.agegr.period, model="rw2", group=id.agegr, control.group=list(model="rw2")) + f(id.agegr, model="rw2", group=id.period, control.group=list(model="rw2"))

formula.2 <- births ~ f(id.district, model="bym2", graph="mw.adj") + f(id.agegr.period, model="rw2") + f(id.agegr, model="rw2", group=id.period, control.group=list(model="rw2"))
## 21561

formula.3 <- births ~ f(id.district, model="bym2", graph="mw.adj") + f(id.agegr.period, model="rw2", group=id.agegr, control.group=list(model="rw2")) + f(id.agegr, model="rw2", group=id.period, control.group=list(model="rw2"))
## 21561

formula.4 <- births ~ f(id.period, model="rw2") + f(id.agegr, model="rw2") + f(id.agegr2, model="rw2", group=id.period, control.group=list(model="rw2"))

formula.5 <- births ~ f(id.period, model="rw2") + f(id.agegr, model="rw2")

mod <- inla(formula.4, family="poisson", data=asfr_pred, E=pys,
            control.family=list(link='log'),
            control.predictor=list(compute=TRUE, link=1),
            control.inla = list(strategy = "gaussian", int.strategy = "eb"),
            control.compute=list(config = TRUE, dic= TRUE, cpo=TRUE))
            # control.fixed=list(mean=0, prec=0.00001,
                              # mean.intercept=0, prec.intercept=0.00001))

summary(mod5)

pred_size <- nrow(crossing(period = asfr1$period, agegr = asfr1$agegr, region = asfr1$v024, pys=1))

pred5 <- asfr_pred %>%
  filter(id<pred_size+1) %>%
  dplyr::select("agegr", "period", "pys", "v024","id") %>%
  left_join(mod$summary.fitted.values[1:pred_size, ] %>%
              mutate(id = 1:pred_size), by="id") %>%
  arrange(period, agegr) %>%
  mutate(agegroup = rep(1:7, each=5, times=pred_size/35),
         agegroup = factor(agegroup, levels=1:7, labels=c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49")))
  ggplot(aes(x=period, y=`0.5quant`, ymin=`0.025quant`, ymax=`0.975quant`, group=agegr))+
  geom_line(aes(color=agegroup)) +
  labs(title="Without period/age interaction")
  ylim(0,0.4) 
  facet_wrap(~district)
  
pred4 <- pred4 %>%
  mutate(source="With interaction")
  
pred5 <- pred5 %>%
  mutate(source="Without interaction")

pred4 %>%
  bind_rows(pred5) %>%
  filter(period %in% c(2000, 2005, 2010, 2015)) %>%
  ggplot(aes(x=agegr, y=log(`0.5quant`), group=source))+
    geom_line(aes(color=source)) +
    facet_grid(v024~period)
    
  
gridExtra::grid.arrange(p4, p5)

foo <- mod$summary.random$id.agegr

mod$marginals.random$id.agegr %>%
  lapply(t_emarg) %>%
  melt() %>%
  mutate(age = rep(15:49, 22),
         agegroup = rep(1:7, each=5, times=22),
         agegroup = factor(agegr, labels=c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49")),
         period = rep(1995:2016, each=35)) %>%
  ggplot(aes(x=period, y=value, group=age)) +
    geom_line(aes(color=agegr)) +
    xlab("ASFR") +
    ylab("")





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
  mutate(Var2 = as.character(Var2))
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

