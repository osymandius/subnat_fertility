library(tidyverse)
library(magrittr)
library(rdhs)
library(demogsurv)
library(INLA)
library(reshape2)
library(survival)
library(sf)
library(spdep)

calc_asfr1 <- function(data,
                       by = NULL,
                       agegr = NULL,
                       period = NULL,
                       cohort = NULL,
                       tips = NULL,
                       clusters=~v021,
                       strata=NULL,
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

ir_to_area <- function(ir, areas2) {
  ir %>%
    left_join(areas2, by=c("v001" = "cluster_id")) %>%
    filter(!is.na(area_id))
}

setwd("~/Documents/GitHub/subnat_fertility")
set_rdhs_config(email="o.stevens@imperial.ac.uk", project="Subnational fertility", config_path = "~/.rdhs.json")

oli_surveys <- readRDS("~/Documents/GitHub/naomi-data-edit/oli_cluster.rds") %>%
  dplyr::select(survey_id) %>%
  mutate(iso3 = survey_id) %>%
  separate(col="iso3", into="iso3", sep=3) %>%
  distinct

#' Choose a country
iso2_code <- "KE"
iso3_code <- "KEN"

clusters <- readRDS("~/Documents/GitHub/naomi-data-edit/oli_cluster.rds") %>%
  mutate(iso3 = survey_id) %>%
  separate(col="iso3", into="iso3", sep=3) %>%
  filter(iso3 == iso3_code)

admin1_areas <- readRDS("~/Documents/GitHub/naomi-data/data/areas/areas_wide.rds") %>%
  filter(iso3 == iso3_code) %>%
  left_join(clusters, by=c("id5" = "geoloc_area_id", "iso3")) %>%
  dplyr::select(iso3, id1, name1, survey_id, cluster_id) %>%
  distinct

# areas <- readRDS("~/Documents/GitHub/naomi-data/data/areas/areas_long.rds") %>%
#   filter(iso3 ==iso3_code) %>%
#   #filter(parent_area_id %in% c("ZWE")) %>%
#   left_join(clusters, by=c("area_id" = "geoloc_area_id", "iso3"))
# 

# areas <- readRDS("~/Documents/GitHub/naomi-data/data/areas/areas_long.rds") %>%
#   filter(iso3 ==iso3_code) %>%
#   filter(area_level==1) %>%
#   #filter(area_level == ifelse(iso3_code=="MWI", 4, max(area_level))) %>%
#   left_join(clusters, by=c("area_id" = "geoloc_area_id", "iso3"))

boundaries <- readRDS("~/Documents/GitHub/naomi-data/data/areas/boundaries.rds")

##+ datasets
surveys <- dhs_surveys(countryIds = c(iso2_code), surveyYearStart=1995) %>%
  filter(SurveyType == "DHS")

foo <- surveys %>%
  separate(SurveyId, into=c("cc", "year", "type"), sep = c(2,6)) %>%
  type.convert %>%
  dplyr::select(year)

admin1_areas <- admin1_areas %>%
  separate(survey_id, into=c("cc", "year", "type"), sep = c(3,7), remove=FALSE) %>%
  mutate(year = as.numeric(year)) %>%
  filter(year %in% foo$year) %>%
  dplyr::select(-c(cc, year, type)) %>%
  rename(area_id = id1, area_name = name1)
  

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

head(ir[["ET2005DHS"]]$v008)

test_ir <- ir[["ET2005DHS"]]
test_ir$v008 <- ir[["ET2005DHS"]]$v008+92
test_ir$v011 <- ir[["ET2005DHS"]]$v011+92

if (iso3_code == "ETH") {
  
  calc_asfr
  
  ## 2005 DHS
  ## 2011 DHS
  ## 2016 DHS
  ir[[3]]$v007 <- floor(ir[[3]]$v007+92/12) %>% as.integer()
  ir[[3]]$v008 <- ir[[3]]$v008+92 %>% as.integer()
  ir[[3]]$v010 <- floor(ir[[3]]$v007+92/12) %>% as.integer()
  ir[[3]]$v011 <- ir[[3]]$v011+92 %>% as.integer()
  
  ir[[2]]$v008 <- ir[[2]]$v008+92 %>% as.integer()
  ir[[3]]$v008 <- ir[[3]]$v008+92 %>% as.integer()
  ir[[4]]$v008 <- ir[[4]]$v008+92 %>% as.integer()
}

# Tanzania 2017 MIS- shock because TIPS going back too far


areas2 <- admin1_areas %>%
  group_by(survey_id) %>%
  group_split(keep=TRUE)

# areas2 <- areas %>%
#   group_by(survey_id) %>%
#   group_split(keep=TRUE)

ir_area <- Map(ir_to_area, ir, areas2) %>% 
  bind_rows() %>%
  group_by(survey_id, area_id) %>%
  group_split(keep = TRUE)

names(ir_area) <- sapply(ir_area,  function(x) {
    paste(x[["area_id"]][1], x[["survey_id"]][1])
})

tips_surv <- list("DHS" = c(0:10), "MIS" = c(0:5), "AIS" = c(0:5))[surveys$SurveyType]

debugonce(demog_pyears)

## ETH testing - why are there 0 births?
# asfr <- calc_asfr1(test_ir, by = ~surveyid + survyear,
#                    tips = c(0,7),
#                    agegr= 3:10*5,
#                    #period = list(seq(1995, 2017, by=0.5)),
#                    period = c(1995:2017),
#                    counts = TRUE)

asfr <- Map(calc_asfr1, ir,
            by = list(~surveyid + survyear),
            tips = tips_surv,
            agegr= list(3:10*5),
            #period = list(seq(1995, 2017, by=0.5)),
            period = list((1995*1):(2017*1)*(1/1)),
            counts = TRUE)

# asfr_15to19 <- Map(calc_asfr1, ir_area,
#                by = list(~surveyid + survyear),
#                tips = tips_surv,
#                agegr= list(15:20),
#                period = list(1995:2017),
#                counts = TRUE)
# 
# asfr_20to49 <- Map(calc_asfr1, ir_area,
#             by = list(~surveyid + survyear),
#             tips = tips_surv,
#             agegr= list(4:10*5),
#             period = list(1995:2017),
#             counts = TRUE)

### Plot ASFR data at national level as check:
# asfr %>%
#   bind_rows() %>%
#   type.convert() %>%
#   ggplot(aes(period, asfr, group=surveyid)) +
#     geom_point(aes(color=surveyid)) +
#     geom_line(aes(color=surveyid)) +
#     #ylim(0,0.5) +
#     facet_wrap(~agegr)

asfr1 <- asfr %>%
  bind_rows %>%
  separate(period, into=c("period", "stuff"), sep="-") %>%
  dplyr::select(-stuff) %>%
  type.convert

asfr1 <- asfr %>%
  bind_rows %>%
  type.convert

# asfr1 <- asfr_15to19 %>%
#   bind_rows %>%
#   type.convert %>%
#   bind_rows(asfr_20to49 %>%
#               bind_rows %>%
#               type.convert %>%
#               mutate(agegr = as.character(agegr)) %>%
#               separate(agegr, into = c("agegr", "rest"), sep="-") %>%
#               dplyr::select(-rest) %>%
#               type.convert %>%
#               mutate(agegr = agegr + 2.5))

#area_name = asfr1$area_name,
asfr_pred <- crossing(surveyid = asfr1$surveyid, period = asfr1$period, agegr = asfr1$agegr,  pys=1) %>%
  bind_rows(asfr1) %>%
  mutate(id.period = group_indices(., period),
         id.period2 = id.period,
         id.agegr = group_indices(., agegr),
         id.agegr2 = id.agegr,
         id.agegr.period = group_indices(., period, agegr),
         id.tips = (group_indices(., tips)),
         dummy_step = ifelse(tips>5, 1, 0),
         tips.iid = ifelse(tips>5, 2, 1),
         # id.district = group_indices(., area_name),
         # id.district2 = id.district,
         # id.district3 = id.district,
         # id.survey = group_indices(., surveyid),
         # id.agegr.period.district = group_indices(., agegr, period, district),
         # id.region = group_indices(., region_name),
         id = 1:nrow(.),
  )

# sh <- areas %>%
#   # filter(parent_area_id %in% c("ZWE")) %>%
#   #filter(iso3 == iso3_code, area_level == ifelse(iso3_code=="MWI", 3, max(area_level))) %>%
#   # select(-parent_area_id) %>%
#   mutate(area_idx = row_number())

sh <- readRDS("~/Documents/GitHub/naomi-data/data/areas/areas_long.rds") %>%
  filter(parent_area_id == iso3_code) %>%
  mutate(area_idx = row_number())

#' Neighbor list
nb <- sh %>%
  left_join(boundaries) %>%
  st_as_sf %>%
  as("Spatial") %>%
  spdep::poly2nb() %>%
  `names<-`(sh$area_idx)

# adj <- poly2nb(adj_data, row.names = adj_data$area_id, queen = TRUE)
nb2INLA(paste0(iso3_code, ".adj"), nb)

formulae <- list()
formulae[[1]] <- births ~  f(id.district, model="bym2", graph=paste0(iso3_code, ".adj")) + f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2")) + f(id.agegr.period, model="iid")

formulae[[2]] <- formula.2 <- births ~  f(id.district, model="bym2", graph=paste0(iso3_code, ".adj")) + f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2")) + f(id.agegr.period, model="iid") + f(id.district2, model="bym2", graph=paste0(iso3_code, ".adj"), group=id.period, control.group=list(model="rw1"))

formulae[[3]] <- formula.3 <- births ~ f(id.district, model="bym2", graph=paste0(iso3_code, ".adj")) + f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2")) + f(id.agegr.period, model="iid") + f(id.district2, model="bym2", graph=paste0(iso3_code, ".adj"), group=id.period, control.group=list(model="rw1")) +f(id.district3, model="bym2", graph=paste0(iso3_code, ".adj"), group=id.agegr, control.group=list(model="rw1"))

## Irregular spacing RW2 testing
formulae[[4]] <- births ~ f(agegr, model="rw2", scale.model=TRUE, values=c(15:19, 4.5:9.5*5)) + f(id.period, model="rw2")
formulae[[5]] <- births ~ f(agegr, model="rw2", scale.model=TRUE, values=c(15:19, 4.5:9.5*5)) + f(id.period, model="rw2") + f(id.survey, model="iid") ## Something terrible happens here
formulae[[6]] <- births ~ f(agegr, model="rw2", scale.model=TRUE, values=c(30:100*0.5)) + f(id.period, model="rw2")

## Smaller period testing
formulae[[7]] <- births ~ f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2")) + f(id.agegr.period, model="iid")
## Quarter fixed effects
formulae[[8]] <- births ~ id.quarter + f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2"))

formulae[[9]] <- births ~ f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2"))

formulae[[1]] <- births ~ f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2")) + f(id.tips, model="rw1")

formulae[[2]] <- births ~ f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2")) + f(id.tips, dummy_step, model="rw1")

formulae[[3]] <- births ~ dummy_step + f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2")) + f(id.tips, model="rw1")

formulae[[4]] <- births ~ f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2")) + f(id.tips, model="rw1", group=dummy_step, control.group=list(model="iid"))

# 
# formulae[2] <- births ~  
#   f(id.district, model="bym2", graph=paste0(iso3_code, ".adj")) + 
#   f(id.agegr, model="rw1") +
#   f(id.period, model="rw2") +
#   f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2")) +
#   f(id.agegr.period, model="iid") +
#   f(id.district2, model="bym2", graph=paste0(iso3_code, ".adj"), group=id.period, control.group=list(model="rw1")) 
# 
# formula[3] <- births ~  
#   f(id.district, model="bym2", graph=paste0(iso3_code, ".adj")) + 
#   f(id.agegr, model="rw1") +
#   f(id.period, model="rw2") +
#   f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2")) +
#   f(id.agegr.period, model="iid") +
#   f(id.district2, model="bym2", graph=paste0(iso3_code, ".adj"), group=id.period, control.group=list(model="rw1")) +
#   f(id.district3, model="bym2", graph=paste0(iso3_code, ".adj"), group=id.agegr, control.group=list(model="rw1"))


  # f(id.agegr2, id.period2, model="rw2", constr = TRUE) +
  # f(id.district2, model="bym2", graph=paste0(iso3_code, ".adj"), group=id.period, control.group=list(model="rw1")) +
  # f(id.district3, model="bym2", graph=paste0(iso3_code, ".adj"), group=id.agegr, control.group=list(model="rw1")) +

plots <- lapply(formulae, function(formulae) {
  
  mod1 <- inla(formulae, family="poisson", data=asfr_pred, E=pys,
               control.family=list(link='log'),
               control.predictor=list(compute=TRUE, link=1),
               control.inla = list(strategy = "gaussian", int.strategy = "eb"),
               control.compute=list(config = TRUE, dic= TRUE, cpo=TRUE),
               verbose=TRUE)
  
  pred_size <- nrow(asfr_pred) - nrow(asfr1)
  
  pred <- asfr_pred %>%
    filter(id<pred_size+1) %>%
    dplyr::select("agegr", "period", "pys","id") %>%
    left_join(mod1$summary.fitted.values[1:pred_size, ] %>%
                mutate(id = 1:pred_size), by="id") %>%
    arrange(period, agegr) %>%
    #filter(agegr <20) %>%
    mutate(agegr = factor(agegr))
  # mutate(agegroup = rep(1:7, each=5, times=pred_size/35),
  #        agegroup = factor(agegroup, levels=1:7, labels=c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49")))
  
  plot <- ggplot()+
    geom_line(data = pred, aes(x=period, y=`0.5quant`, group=agegr, color=agegr)) +
    geom_vline(aes(xintercept=1995:2016), linetype=3) +
    geom_rect(data=foo, aes(xmin = year-1, xmax = year, ymin = -Inf, ymax=Inf), fill="blue", alpha=0.2) +
    labs(y="ASFR", x=element_blank(), title="") +
    ylim(0, 0.34)
  
  return(plot)
  
})

  
gridExtra::grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]])

  ggplot()+
    geom_line(data = pred, aes(x=period, y=`0.5quant`, group=agegr, color=agegr)) +
    geom_vline(aes(xintercept=1995:2016), linetype=3) +
    geom_rect(data=foo, aes(xmin = year-1, xmax = year, ymin = -Inf, ymax=Inf), fill="blue", alpha=0.2) +
    labs(y="ASFR", x=element_blank(), title=paste(formula))
    ylim(0, 0.34)
    
    
  form12 <- exp(mod1$summary.random$id.tips) %>% 
      mutate(id = 1:11) %>%
      ggplot(aes(x=id)) +
      geom_point(aes(y=`0.5quant`)) +
      geom_errorbar(aes(ymin = `0.025quant`, ymax=`0.975quant`), width=0.2)
    
  quarter_summary <- exp(mod1$summary.fixed) %>%
    mutate(country = "ZW") %>%
    bind_rows(quarter_summary2)
  
  quarter_summary %>%
    mutate(id = rownames(.)) %>%
    ggplot(aes(x=id)) +
    geom_point(aes(y=`0.5quant`)) +
    geom_errorbar(aes(ymin = `0.025quant`, ymax=`0.975quant`), width=0.2)
    
    int_mnth <- ir %>%
      bind_rows
    
    svy_mnth <- int_mnth %>% 
      count(surveyid, v008) %>%
      mutate(month = 1+(12*(v008/12 - floor(v008/12)))) %>%
      type.convert %>%
      group_by(surveyid) %>%
      mutate(prop = round(100*(n/sum(n)))) %>%
      ungroup %>%
      mutate(month = factor(month, labels=c("Jan", "Feb", "Mar", "Apr" ,"May", "Jun", "Jul" ,"Aug", "Sep", "Oct", "Nov","Dec"))) %>%
      filter(prop>0)
    
