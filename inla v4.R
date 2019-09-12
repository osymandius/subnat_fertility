library(countrycode)
library(tidyverse)
library(magrittr)
library(rdhs)
library(demogsurv)
library(INLA)
library(reshape2)
library(survival)
library(sf)
library(spdep)

setwd("~/Documents/GitHub/subnat_fertility")
set_rdhs_config(email="o.stevens@imperial.ac.uk", project="Subnational fertility", config_path = "~/.rdhs.json")

source("fertility_funs.R")

## Set flags. Multicountry == TRUE for running several countries in a single model. Subnational currently only set up for running 1 country at a time, setting the iso codes as strings manually {see lines 131 onwards for ZWE admin 1 example}
subnational <- TRUE
multicountry <- FALSE

## If Namibia is in this list, ensure it's first. Crap but it'll do for now.
iso2 <- c("NM", "RW", "UG", "MZ", "MW", "KE", "ZM", "ZW", "TZ")

iso2_list <- as.list(iso2)

iso3_list <- lapply(iso2_list, function(x) {x %>% countrycode("iso2c", "iso3c")})

if (iso2[1]=="NM") {
  iso3_list[1] <- "NAM"
}


formulae <- list()

formulae[[1]] <- births ~ f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2")) + f(id.tips, model="rw1")

formulae[[2]] <- births ~ f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2")) + f(id.tips, dummy_step, model="rw1")

formulae[[3]] <- births ~ dummy_step + f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2")) + f(id.tips, model="rw1")

formulae[[4]] <- births ~ f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2")) + f(id.tips, model="rw1", group=dummy_step, control.group=list(model="iid"))

formulae_number <- length(formulae)

## Get cluster coordinates
clusters <- readRDS("~/Documents/GitHub/naomi-data-edit/oli_cluster.rds") %>%
  mutate(iso3 = survey_id) %>%
  separate(col="iso3", into="iso3", sep=3) %>%
  filter(iso3 %in% unlist(iso3_list))


## Get surveys and filter w.r.t clusters. Split into country list.
surveys <- dhs_surveys(countryIds = c(iso2), surveyYearStart=1995) %>%
  separate(SurveyId, into=c("code", "rest"), sep=2, remove=FALSE) %>%
  mutate(iso3 = countrycode(CountryName, "country.name", "iso3c"),
         survey_id = paste0(iso3, rest)) %>%
  select(-c(code, rest)) %>%
  filter(survey_id %in% clusters$survey_id,
         SurveyType != "AIS") %>%
  group_split(CountryName)

## For testing take ZW at admin 2 level
surveys <- surveys[8]

areas <- readRDS("~/Documents/GitHub/naomi-data/data/areas/areas_wide.rds") %>%
  filter(iso3 == "ZWE") %>%
  left_join(clusters, by=c("id2" = "geoloc_area_id", "iso3")) %>%
  dplyr::select(iso3, id1, name1, survey_id, cluster_id) %>%
  rename(area_id = id1, area_name = name1)

area_list <- areas %>%
  group_by(survey_id) %>%
  group_split(keep=TRUE)

# area_list <- areas %>%
#   group_by(survey_id) %>%
#   group_split(keep=TRUE)

boundaries <- readRDS("~/Documents/GitHub/naomi-data/data/areas/boundaries.rds")

ird <- lapply(surveys, function(surveys) {
  dhs_datasets(fileType = "IR", fileFormat = "flat", surveyIds = surveys$SurveyId)
})

ird <- lapply(ird, function(x) {
  x %>%
    mutate(path = unlist(get_datasets(x))) %>%
    bind_rows
})

ir <- unlist(lapply(ird, "[", "path")) %>%
  lapply(readRDS) %>%
  Map(function(ir, surveys) {
    mutate(ir,
           surveyid = surveys$SurveyId,
           country = surveys$CountryName,
           survyear = surveys$SurveyYear,
           survtype = surveys$SurveyType)
  }, ., surveys %>% 
    bind_rows %>%
    group_split(SurveyId))

ir_area <- Map(ir_by_area, ir, area_list) %>% 
  bind_rows() %>%
  group_by(survey_id, area_id) %>%
  group_split(keep = TRUE)

names(ir_area) <- sapply(ir_area,  function(x) {
  paste(x[["area_id"]][1], x[["survey_id"]][1])
})

test <- lapply(ir_area, function(x){
  x %>% select(survtype) %>% distinct
}) %>%
  bind_rows

tips_surv <- list("DHS" = c(0:10), "MIS" = c(0:5), "AIS" = c(0:5))[test$survtype]

# tips_surv <- list("DHS" = c(0:10), "MIS" = c(0:5), "AIS" = c(0:5))[surveys %>%
#                                                                      bind_rows %>%
#                                                                      .$SurveyType]

asfr <- Map(calc_asfr1, ir_area,
            by = list(~country + surveyid + survyear + area_name),
            tips = tips_surv,
            agegr= list(3:10*5),
            #period = list(seq(1995, 2017, by=0.5)),
            period = list((1995*1):(2017*1)*(1/1)),
            counts = TRUE)

#### ZWE ADMIN 1 TEST #####

asfr1 <- asfr %>%
  bind_rows %>%
  type.convert

asfr_pred <- crossing(area_name = asfr1$area_name, period = asfr1$period, agegr = asfr1$agegr,  pys=1) %>%
  bind_rows(asfr1) %>%
  mutate(id.period = group_indices(., period),
         id.period2 = id.period,
         id.agegr = group_indices(., agegr),
         id.agegr2 = id.agegr,
         id.agegr.period = group_indices(., period, agegr),
         id.tips = (group_indices(., tips)),
         dummy_step = ifelse(tips>5, 1, 0),
         tips.iid = ifelse(tips>5, 2, 1), ## do na_for_cats or whatever
         id.district = group_indices(., area_name),
         # id.district2 = id.district,
         # id.district3 = id.district,
         # id.region = group_indices(., region_name),
         id = 1:nrow(.),
  )

# sh <- readRDS("~/Documents/GitHub/naomi-data/data/areas/areas_long.rds") %>%
#   #filter(parent_area_id %in% c("MWI")) %>%
#   filter(iso3 == "ZWE", area_level == 1) %>%
#   select(-parent_area_id) %>%
#   mutate(area_idx = row_number())

## FOR ADMIN 1
sh <- readRDS("~/Documents/GitHub/naomi-data/data/areas/areas_long.rds") %>%
  filter(parent_area_id == "ZWE") %>%
  mutate(area_idx = row_number())

#' Neighbor list
nb <- sh %>%
  left_join(boundaries) %>%
  st_as_sf %>%
  as("Spatial") %>%
  spdep::poly2nb() %>%
  `names<-`(sh$area_idx)

nb2INLA(paste0("ZWE", ".adj"), nb)

formulae <- lapply(formulae, function(x) {
  update(x,    ~ . + f(id.district, model="bym2", graph="ZWE.adj"))
})

mod_list <- Map(run_mod, formulae, list(asfr_pred))

pred_list <- Map(get_pred, mod_list, list(asfr_pred), list(asfr1))

pred_split <- Map(function(x, y) x %>% mutate(source=paste("formula", y)), x=pred_list, y=(1:formulae_number))

plots <- lapply(pred_split, function(pred_split) {
  pred_split %>%
    ggplot()+
    geom_line(aes(x=period, y=`0.5quant`, group=agegr, color=agegr)) +
    labs(y="ASFR", x=element_blank(), title=unique(pred_split$source)) +
    facet_wrap(~area_name)
  
})

ggsave("ZWE admin 1, formula 1.png", plot=plots[[1]], device="png")
ggsave("ZWE admin 1, formula 2.png", plot=plots[[2]], device="png")
ggsave("ZWE admin 1, formula 3.png", plot=plots[[3]], device="png")
ggsave("ZWE admin 1, formula 4.png", plot=plots[[4]], device="png")

#### COUNTRIES RUN TOGETHER #####

    asfr1 <- asfr %>%
      bind_rows %>%
      type.convert
    
    #area_name = asfr1$area_name,
    asfr_pred <- crossing(country = asfr1$country, period = asfr1$period, agegr = asfr1$agegr,  pys=1) %>%
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
    
    formulae <- lapply(formulae, function(x) {
      update(x,    ~ . + country)
    })
    
    mod_list <- Map(run_mod, formulae, list(asfr_pred))
    
    pred_list <- Map(get_pred, mod_list, list(asfr_pred), list(asfr1))
    
    pred_split <- Map(function(x, y) x %>% mutate(source=paste("formula", y)), x=pred_list, y=(1:formulae_number))
    
    plots <- lapply(pred_split, function(pred_split) {
      pred_split %>%
        ggplot()+
        geom_line(aes(x=period, y=`0.5quant`, group=agegr, color=agegr)) +
        labs(y="ASFR", x=element_blank(), title=unique(pred_split$source)) +
        facet_wrap(~country)
      
    })
    
    ggsave("Together, formula 1.png", plot=plots[[1]], device="png")
    ggsave("Together, formula 2.png", plot=plots[[2]], device="png")
    ggsave("Together, formula 3.png", plot=plots[[3]], device="png")
    ggsave("Together, formula 4.png", plot=plots[[4]], device="png")
    
##### COUNTRIES RUN SEPARATELY #####
  
  asfr1_country <- asfr %>%
    bind_rows %>%
    type.convert %>%
    group_split(country, keep=TRUE) %>%
    lapply(droplevels)
  
  asfr_pred_country <- lapply(asfr1_country, function(asfr1_country) {
    crossing(country = asfr1_country$country, period = asfr1_country$period, agegr = asfr1_country$agegr,  pys=1) %>%
      bind_rows(asfr1_country) %>%
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
  })
  
  formulae <- repeat_formulae(formulae, length(iso2))
  
  mod_list <- Map(run_mod, formulae, asfr_pred = rep(asfr_pred_country, formulae_number))
  
  mod_split <- Map(function(x, y) {
    mod_list[x:(x+length(iso2)-1)]
  }, x=seq(1, formulae_number*length(iso2), by=length(iso2)), y=1:formulae_number)
  
  pred_list <- Map(get_pred, mod_list, rep(asfr_pred_country, formulae_number), rep(asfr1_country, formulae_number))
  
  pred_split <- Map(function(x, y) {
    pred_list[x:(x+length(iso2)-1)] %>%
      bind_rows %>%
      mutate(source=paste("formula", y))
  }, x=seq(1, formulae_number*length(iso2), by=length(iso2)), y=1:formulae_number)
  
  plots <- lapply(pred_split, function(pred_split) {
    pred_split %>%
      ggplot()+
      geom_line(aes(x=period, y=`0.5quant`, group=agegr, color=agegr)) +
      labs(y="ASFR", x=element_blank(), title=unique(pred_split$source)) +
      facet_wrap(~country)
    
  })
  
  ggsave("Separate, formula 1.png", plot=plots[[1]], device="png")
  ggsave("Separate, formula 2.png", plot=plots[[2]], device="png")
  ggsave("Separate, formula 3.png", plot=plots[[3]], device="png")
  ggsave("Separate, formula 4.png", plot=plots[[4]], device="png")
  
  
  


