library(tidyverse)
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

setwd("~/Documents/GitHub/subnat_fertility")

clusters <- readRDS("~/Documents/GitHub/naomi-data/data/survey/survey_clusters.rds")
areas <- readRDS("~/Documents/GitHub/naomi-data/data/areas/areas_long.rds")
boundaries <- readRDS("~/Documents/GitHub/naomi-data/data/areas/boundaries.rds")

clusters %>%
  filter(iso3 =="MWI") %>%
  select(survey_id) %>%
  distinct

#' Choose a country
iso2_code <- "ET"
iso3_code <- "ETH"

clusters <- clusters %>%
  mutate(iso3 = survey_id) %>%
  separate(col="iso3", into="iso3", sep=3)

areas <- areas %>%
  filter(area_level == 2) %>% 
  filter(iso3 %in% c(iso3_code)) %>%
  left_join(clusters, by=c("area_id" = "geoloc_area_id", "iso3"))



set_rdhs_config(email="o.stevens@imperial.ac.uk", project="Subnational fertility", config_path = "~/.rdhs.json")

##+ datasets
surveys <- dhs_surveys(countryIds = c(iso2_code), surveyYearStart=2000)
ird <- dhs_datasets(fileType = "IR", fileFormat = "flat", surveyIds = surveys$SurveyId)
#ged <- dhs_datasets(fileType = "GE", fileFormat = "flat", surveyIds = surveys$SurveyId)

ird$path <- unlist(get_datasets(ird))
#ged$path <- unlist(get_datasets(ged))

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

areas2 <- areas %>%
  # mutate(survey_id = recode(survey_id, 
  #                           "ZWE2005DHS" = "ZW2005DHS",
  #                           "ZWE2010DHS" = "ZW2010DHS",
  #                           "ZWE2015DHS" = "ZW2015DHS",)
  # ) %>% 
  group_by(survey_id) %>%
  group_split(keep=TRUE)

 # names(areas2) <-sapply(areas2,  function(x) {
 #    paste(x[1,2], x[1,6])
 #  })


ir_to_area <- function(ir, areas2) {
  ir %>%
    left_join(areas2, by=c("v001" = "cluster_id")) %>%
    filter(!is.na(area_id))
}

ir <- Map(ir_to_area, ir, areas2) %>% 
  bind_rows() %>%
  group_by(survey_id, area_id) %>%
  group_split(keep = TRUE)

  names(ir) <- sapply(ir,  function(x) {
    paste(x[["area_id"]][1], x[["survey_id"]][1])
  })

tips_surv <- list("DHS" = c(0, 7), "MIS" = c(0, 5))[surveys$SurveyType]

asfr <- Map(calc_asfr1, ir[1:250],
               by = list(~surveyid + survyear + area_id + area_name),
               tips = tips_surv,
               agegr= list(15:50),
               period = list(1995:2017),
               counts = TRUE)

asfr1 <- asfr %>%
  bind_rows() %>%
  type.convert() %>%
  mutate(id.period = group_indices(., period),
         id.period2 = id.period,
         id.agegr = group_indices(., agegr),
         id.agegr2 = id.agegr,
         id.agegr.period = group_indices(., period, agegr),
         id.district = group_indices(., area_name),
         # id.agegr.period.district = group_indices(., agegr, period, district),
         # id.region = group_indices(., region_name),
         id = 1:nrow(.))

asfr_pred <- crossing(period = asfr1$period, agegr = asfr1$agegr, area_name = asfr1$area_name, pys=1) %>%
  bind_rows(asfr1)%>%
  mutate(id.period = group_indices(., period),
         id.period2 = id.period,
         id.agegr = group_indices(., agegr),
         id.agegr2 = id.agegr,
         id.agegr.period = group_indices(., period, agegr),
         id.district = group_indices(., area_name),
         # id.agegr.period.district = group_indices(., agegr, period, district),
         # id.region = group_indices(., region_name),
         id = 1:nrow(.))

# adj_data <- areas %>%
#   filter(iso3 == "MWI", area_level == 2) %>%
#   select(area_id) %>%
#   distinct() %>%
#   left_join(boundaries, by="area_id") %>%
#   st_as_sf

sh <- areas %>%
  filter(iso3 == iso3_code, area_level == 2) %>%
  select(-parent_area_id) %>%
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

formula.1 <- births ~ f(id.district, model="bym2", graph=paste0(iso3_code, ".adj")) + f(id.agegr.period, model="iid") + f(id.agegr, model="rw2", group=id.period, control.group=list(model="rw2"))

mod1 <- inla(formula.1, family="poisson", data=asfr_pred, E=pys,
             control.family=list(link='log'),
             control.predictor=list(compute=TRUE, link=1),
             control.inla = list(strategy = "gaussian", int.strategy = "eb"),
             control.compute=list(config = TRUE, dic= TRUE, cpo=TRUE))

pred_size <- nrow(asfr_pred) - nrow(asfr1)

pred <- asfr_pred %>%
  filter(id<pred_size+1) %>%
  dplyr::select("agegr", "period", "pys", "area_name","id") %>%
  left_join(mod1$summary.fitted.values[1:pred_size, ] %>%
              mutate(id = 1:pred_size), by="id") %>%
  arrange(period, agegr) %>%
  mutate(agegroup = rep(1:7, each=5, times=pred_size/35),
         agegroup = factor(agegroup, levels=1:7, labels=c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49")))
head(
pred %>%
  filter(agegr==22, period==2016)
)

pred %>%
  ggplot(aes(x=period, y=`0.5quant`, group=agegr))+
  geom_line(aes(color=agegr)) +
  facet_wrap(~area_name)
