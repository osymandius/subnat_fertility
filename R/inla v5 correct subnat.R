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
library(parallel)
library(viridis)
library(tidyr)
library(parallel)

setwd("~/Documents/GitHub/subnat_fertility")
source("fertility_funs.R")
boundaries <- readRDS("~/Documents/GitHub/naomi-data/data/area_boundaries.rds")
areas_long <- readRDS("~/Documents/GitHub/naomi-data/data/area_hierarchy.rds")
areas_wide <- readRDS("~/Documents/GitHub/naomi-data/data/areas_wide.rds")

dhs_iso3 <- dhs_countries(returnFields=c("CountryName", "DHS_CountryCode")) %>%
  mutate(iso3 = countrycode(CountryName, "country.name", "iso3c"),
         iso3 = ifelse(CountryName == "Eswatini", "SWZ", iso3))

## Get cluster coordinates
clusters <- readRDS("~/Documents/GitHub/naomi-data-edit/oli_cluster2.rds") %>%
  mutate(iso3 = survey_id) %>%
  separate(col="iso3", into="iso3", sep=3) %>%
  left_join(dhs_iso3 %>% select(-CountryName), by="iso3") %>%
  separate(survey_id, into=c(NA, "surv"), sep=3, remove=FALSE) %>%
  mutate(DHS_survey_id = paste0(DHS_CountryCode, surv)) %>%
  select(-surv) %>%
  filter(survey_id != "UGA2011AIS")

iso3 <- as.list(clusters %>% .$iso3 %>% unique)

## Get surveys for which we have clusters. Split into country list.
surveys <- dhs_surveys(surveyIds = unique(clusters$DHS_survey_id)) %>%
  left_join(clusters %>% select(c(DHS_survey_id, survey_id)) %>% distinct, by=c("SurveyId" = "DHS_survey_id")) %>%
  filter(CountryName != "Ethiopia", SurveyType != "AIS") %>%
  group_split(CountryName)

areas <- readRDS("~/Documents/GitHub/naomi-data/data/areas/areas_long.rds") %>%
  inner_join(clusters, by=c("area_id" = "geoloc_area_id", "iso3")) %>%
  filter(!iso3 %in% c("MWI", "ETH"))

areas_mwi <- readRDS("~/Documents/GitHub/naomi-data/data/areas/areas_long.rds") %>%
  inner_join(clusters, by=c("area_id" = "geoloc_area_id", "iso3")) %>%
  filter(iso3 == "MWI") %>%
  left_join(areas_wide %>% select(area_id, name4, id4), by=c("area_id")) %>%
  select(-c(area_id, area_name, parent_area_id, area_level)) %>%
  rename(area_name = name4, area_id = id4)

areas <- areas %>%
  bind_rows(areas_mwi) %>%
  separate(survey_id, into=c(NA, "SurveyType"), sep=-3, remove=FALSE) %>%
  filter(SurveyType != "AIS") %>%
  select(-SurveyType)

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

area_list <- areas %>%
  group_by(survey_id) %>%
  group_split(keep=TRUE)

ir_area <- Map(ir_by_area2, ir, area_list, n=1:length(ir), total=length(ir)) %>% 
  unlist(recursive = FALSE)


test <- lapply(ir_area, function(x) {
  x$foo <- x %>%
    .$survtype %>%
    unique
}) %>%
  unlist

test <- data.frame("survtype" = test)

names(ir_area) <- sapply(ir_area,  function(x) {
  paste(x[["area_id"]][1], x[["survey_id"]][1])
})


tips_surv <- list("DHS" = c(0:15), "MIS" = c(0:5), "AIS" = c(0:5))[test$survtype]
tips_surv <- list("DHS" = c(0:10), "MIS" = c(0:5), "AIS" = c(0:5))[surveys %>%
                                                                     bind_rows %>%
                                                                     .$SurveyType]

asfr <- Map(calc_asfr1, ir_area,
            y=1:length(ir_area),
            by = list(~country + surveyid + survtype + survyear + area_name + area_id),
            tips = list(0:10),
            agegr= list(3:10*5),
            #period = list(seq(1995, 2017, by=0.5)),
            period = list((1995*1):(2017*1)*(1/1)),
            counts = TRUE) %>%
  bind_rows

asfr1_country <- list(asfr_mwi, asfr_lso, asfr_rwa, asfr_zwe)

asfr1_country <- lapply(asfr1_country, type.convert)

asfr_pred_country <- lapply(asfr1_country, function(asfr1_country) {
  crossing(country = asfr1_country$country, period = 1995:2020, agegr = asfr1_country$agegr,  pys=1) %>%
    bind_rows(asfr1_country) %>%
    mutate(id.period = group_indices(., period),
           id.period2 = id.period,
           id.period3 = id.period,
           id.agegr = group_indices(., agegr),
           id.agegr2 = id.agegr,
           id.agegr3 = id.agegr,
           id.agegr.period = group_indices(., period, agegr),
           # tips = factor(tips),
           id.tips = (group_indices(., tips)),
           id.tips = ifelse(is.na(tips), NA, id.tips),
           tips_dummy = ifelse(tips>5, 1, 0),
           # survey_dummy = (group_indices(., survtype)),
           # survey_dummy = ifelse(is.na(survtype), NA, survey_dummy),
           # id.district = group_indices(., area_name),
           # id.district2 = id.district,
           # id.district3 = id.district,
           # id.survey = group_indices(., surveyid),
           # id.agegr.period.district = group_indices(., agegr, period, district),
           # id.region = group_indices(., region_name),
           id = 1:nrow(.),
           
    )
})

asfr_pred_country <- lapply(asfr_pred_country, function(x) {
  x %>%
    mutate_if(is.factor, as.character)
})

formula <- births ~ tips_dummy + f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2")) + f(id.tips, model="rw1")

mod_list <- lapply(asfr_pred_country, run_mod, formulae = formula, model_family="poisson")
subnational = FALSE
multicountry = FALSE

pred_list <- Map(get_pred, mod_list, asfr_pred = asfr_pred_country, asfr1 = asfr1_country)

pred <- pred_list %>%
  bind_rows

saveRDS(asfr_pred_country_subnat, file="asfr_pred_country_subnat.rds")

ZWE_poisson <- stats[[4]]
RWA_poisson <- stats[[3]]
UGA_poisson <- stats

saveRDS(MWI_poisson, file="MWI_poisson_mod.rds")
saveRDS(RWA_poisson, file="RWA_poisson_mod.rds")
saveRDS(ZWE_poisson, file="ZWE_poisson_mod.rds")
saveRDS(ZWE_poisson, file="ZWE_poisson_mod.rds")
saveRDS(UGA_poisson, file="UGA_poisson_mod.rds")

saveRDS(mod, file="ZWE_0inflate.rds")
saveRDS(mod2, file="RWA_0inflate.rds")

rwa_mod <- readRDS("RWA_poisson_mod.rds")
mwi_mod <- readRDS("MWI_poisson_mod.rds")
zwe_mod <- readRDS("ZWE_poisson_mod.rds")
lso_mod <- readRDS("LSO_poisson_mod.rds")
uga_mod <- readRDS("UGA_poisson_mod.rds")

pred <- asfr_pred_country_subnat[[5]] %>%
  filter(id<23296+1) %>%
  dplyr::select("area_name", "agegr", "period", "pys","id") %>%
  left_join(uga_mod$summary.fitted.values[1:23296, ] %>%
              mutate(id = 1:23296), by="id") %>%
  arrange(area_name, period, agegr) %>%
  mutate(agegr = factor(agegr))

pred %>%
  ggplot()+
  geom_line(aes(x=period, y=`0.5quant`, group=agegr, color=agegr)) +
  #geom_vline(aes(xintercept=SurveyYear), linetype=3) +
  labs(y="ASFR", x=element_blank(), title="Zim admin-2") +
  facet_wrap(~area_name)

pred %>%
  filter(agegr == "20-24", as.integer(factor(area_name))<20) %>%
  ggplot() +
  geom_line(aes(x=period, y=`0.5quant`, group=agegr, color=agegr)) +
  geom_point(data=asfr_pred_country_subnat[[1]] %>% filter(agegr == "20-24", as.integer(factor(area_name))<20), aes(x=period, y=asfr, group=surveyid, color=surveyid), size=1) +
  labs(y="ASFR", x=element_blank(), title="Zim admin-2") +
  facet_wrap(~area_name)+
  ylim(0,0.35)

plots <- lapply(pred_split, function(pred_split) {
  pred_split %>%
    left_join(
      lapply(surveys, function(surveys) surveys %>% bind_rows %>% select(CountryName, SurveyYear)) %>% bind_rows,
      by = c("country" = "CountryName")
    ) %>%
    ggplot()+
    geom_line(aes(x=period, y=`0.5quant`, group=agegr, color=agegr)) +
    geom_vline(aes(xintercept=SurveyYear), linetype=3) +
    labs(y="ASFR", x=element_blank(), title="Estimated ASFR from MIS only") +
    facet_wrap(~country)
  
})
