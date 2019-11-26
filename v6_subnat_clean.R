library(countrycode)
library(tidyverse)
library(magrittr)
library(rdhs)
library(demogsurv)
library(INLA)
library(reshape2)
library(survival)
library(geojsonsf)
library(sf)
library(spdep)
library(tidyr)
library(naomi)
library(viridis)
library(parallel)

setwd("~/Documents/GitHub/subnat_fertility")
source("fertility_funs.R")

#boundaries_old <- readRDS("~/Documents/GitHub/naomi-data-old/data/area_boundaries.rds")
boundaries <- readRDS("~/Documents/GitHub/subnat_fertility/area_boundaries.rds")

areas_long <- readRDS("~/Documents/GitHub/subnat_fertility/areas_long.RDS")

area_paths <- paste0("~/Documents/GitHub/naomi-data/", iso3_code, "/data/", tolower(iso3_code), "_areas.geojson")

areas_wide <- lapply(area_paths, read_sf) %>%
  lapply(function(x) {spread_areas(as.data.frame(x))}) %>%
  lapply(function(x) {x %>% mutate(iso3 = area_id0)}) %>%
  bind_rows

set_rdhs_config(email="o.stevens@imperial.ac.uk", project="Subnational fertility", config_path = "~/.rdhs.json")

dhs_iso3 <- dhs_countries(returnFields=c("CountryName", "DHS_CountryCode")) %>%
  mutate(iso3 = countrycode(CountryName, "country.name", "iso3c"),
         iso3 = ifelse(CountryName == "Eswatini", "SWZ", iso3))


## Get cluster coordinates
clusters <- readRDS("~/Documents/GitHub/naomi-data-edit/clusters_2019_11_21.rds") %>%
  mutate(iso3 = survey_id) %>%
  separate(col="iso3", into="iso3", sep=3) %>%
  left_join(dhs_iso3 %>% select(-CountryName), by="iso3") %>%
  separate(survey_id, into=c(NA, "surv"), sep=3, remove=FALSE) %>%
  mutate(DHS_survey_id = paste0(DHS_CountryCode, surv)) %>%
  separate(surv, into=c(NA, "SurveyType"), sep=-3) %>%
  filter(iso3 == "MWI", DHS_CountryCode != "OS") %>%
  filter(!survey_id %in% c("MOZ2009AIS", "TZA2003AIS", "UGA2011AIS")) 

## Get surveys for which we have clusters. Split into country list.
surveys <- dhs_surveys(surveyIds = unique(clusters$DHS_survey_id)) %>%
  left_join(clusters %>% select(c(DHS_survey_id, survey_id)) %>% distinct, by=c("SurveyId" = "DHS_survey_id")) %>%
  filter(!SurveyId %in% c("MZ2009AIS", "TZ2003AIS", "UG2011AIS"))  %>%
  filter(CountryName == "Tanzania") %>%
  group_split(CountryName)

areas <- areas_long %>%
  inner_join(clusters, by=c("area_id" = "geoloc_area_id", "iso3"))

area_list <- areas %>%
  group_by(survey_id) %>%
  group_split(keep=TRUE)

names(area_list) <- area_list %>% lapply("[", "survey_id") %>% lapply(unique) %>% bind_rows %>% .$survey_id

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

names(ir) <- ir %>% 
  lapply("[", "surveyid") %>% 
  lapply(unique) %>% 
  bind_rows %>% 
  left_join(clusters %>% 
              select(survey_id, DHS_survey_id) %>% 
              unique, 
            by=c("surveyid" = "DHS_survey_id")) %>% 
  select(-surveyid) %>% 
  .$survey_id

# cols <- c("b3_01", "b3_02", "b3_03", "b3_04", "b3_05", "b3_06", "b3_07", "b3_08", "b3_09", "b3_10", "b3_11", "b3_12", "b3_13", "b3_14", "b3_15", "b3_16", "b3_17", "b3_18", "b3_19", "b3_20")
# 
# ir[[1]]$v008 %<>% add(92)
# ir[[1]]$v011 %<>% add(92)
# ir[[1]][, cols] %<>% add(92)
# 
# ir[[2]]$v008 %<>% add(92)
# ir[[2]]$v011 %<>% add(92)
# ir[[2]][, cols] %<>% add(92)
# 
# ir[[3]]$v008 %<>% add(92)
# ir[[3]]$v011 %<>% add(92)
# ir[[3]][, cols] %<>% add(92)
# 
# ir[[4]]$v008 %<>% add(92)
# ir[[4]]$v011 %<>% add(92)
# ir[[4]][, cols] %<>% add(92)

ir_area <- Map(ir_by_area2, ir, area_list[names(ir)], n=1:length(ir), total=length(ir)) %>%
   unlist(recursive = FALSE)

survey_type <- ir_area %>%
  lapply("[", "survtype") %>%
  lapply(unique) %>%
  bind_rows

tips_surv <- list("DHS" = c(0:15), "MIS" = c(0:5), "AIS" = c(0:5))[survey_type$survtype]

asfr <- Map(calc_asfr1, ir_area,
                y=1:length(ir_area),
                by = list(~country + surveyid + survtype + survyear + area_name + area_id),
                tips = tips_surv,
                #tips = list(c(0,10)),
                agegr= list(3:10*5),
                period = list(1995:2017),
                counts = TRUE)

asfr1_country <- asfr %>%
  bind_rows %>%
  type.convert %>%
  filter(period<=survyear) %>%
  group_split(country, keep=TRUE) %>%
  lapply(droplevels)


asfr_pred_subnat_15_2020 <- lapply(asfr1_country, function(asfr1_country) {
  crossing(country = asfr1_country$country, area_id = asfr1_country$area_id,  period = min(asfr1_country$period):2020, agegr = asfr1_country$agegr,  pys=1) %>%
    bind_rows(asfr1_country) %>%
    mutate(id.period = group_indices(., period),
           id.period2 = id.period,
           id.period3 = id.period,
           id.agegr = group_indices(., agegr),
           id.agegr2 = id.agegr,
           id.agegr3 = id.agegr,
           id.agegr.period = group_indices(., period, agegr),
           id.tips = (group_indices(., tips)),
           id.tips = ifelse(is.na(tips), NA, id.tips),
           id.tips = factor(id.tips),
           tips_dummy = ifelse(tips>5, 1, 0),
           id.district = group_indices(., area_id),
           id.district2 = id.district,
           id.district3 = id.district,
           id = 1:nrow(.),
           
    ) %>%
    mutate_if(is.factor, as.character)
})


saveRDS(asfr_pred_subnat_15_2020, file="asfr_pred_subnat_15_2020_NEW.rds")

max_level <- areas %>% filter(naomi_level == TRUE) %>% select(iso3, area_level) %>% unique

iso3_list <- as.list(max_level$iso3)

Map(function(iso3_list, max_level, areas_long, boundaries) {
  
  sh <- areas_long %>%
    filter(iso3 == iso3_list, area_level == max_level$area_level[max_level$iso3 == iso3_list]) %>%
    mutate(area_idx = row_number())
  
  #' Neighbor list
  nb <- sh %>%
    left_join(boundaries) %>%
    st_as_sf %>%
    as("Spatial") %>%
    spdep::poly2nb() %>%
    `names<-`(sh$area_idx)
  
  nb2INLA(paste0(iso3_list, ".adj"), nb)
  
}, iso3_list, list(max_level), list(areas_long), list(boundaries))


### RUN MODEL ON CLUSTER ####
#' Recommend:
#' 1) Writing single formula within the run_mod function called by the cluster. 
#' 2) Can either write the .adj files ahead of time and dropping into the directory, or write them one by one in the cluster loop
#' 3) lapply(asfr_pred_country, run_mod), rather than the expand.grid method on cluster.
#' Below code requires completed models saved as individual .rds in a list "mod_list"
#' 

mod_list <- readRDS("~/Downloads/ar1_newh_noTZA.rds")

asfr_pred_subnat_15_2020 <- readRDS("asfr_pred_subnat_15_2020_NEW.rds")

asfr_pred_subnat_15_2020 <- asfr_pred_subnat_15_2020[-8]

iso3_code <- asfr_pred_subnat_15_2020 %>%
  lapply("[", "country") %>%
  lapply(unique) %>%
  bind_rows %>%
  mutate(country = ifelse(country== "Eswatini", "SWZ", countrycode(country, "country.name", "iso3c"))) %>%
  .$country

names(mod_list) <- iso3_code
names(asfr_pred_subnat_15_2020) <- iso3_code

calendar_quarter_targets <- convert_calendar_quarter(c(2000:2020), 2)

population_age_female <- load_population_agesex("~/Documents/GitHub/subnat_fertility/population_agesex_wide_wpp2019_raked.rds", "~/Documents/GitHub/subnat_fertility/areas_long.RDS") %>%
  filter(sex == "female", iso3 != "UGA") %>%
  mutate(calendar_quarter = quarter_id_to_calendar_quarter(quarter_id)) %>%
  left_join(get_age_groups()) %>%
  interpolate_fertility_population(., calendar_quarter_targets) %>%
  left_join(areas_long) %>%
  filter(naomi_level == TRUE)



pop_areas <- population_age_female %>%
  left_join(areas_wide) %>%
  group_by(period, agegr, area_id0) %>%
  mutate(id0_agepop = sum(population))  %>%
  
  group_by(period, agegr, area_id1) %>%
  mutate(id1_agepop = sum(population)) %>%
  
  group_by(period, agegr, area_id2) %>%
  mutate(id2_agepop = sum(population)) %>%
  
  group_by(period, agegr, area_id3) %>%
  mutate(id3_agepop = sum(population)) %>%
  
  group_by(period, agegr, area_id4) %>%
  mutate(id4_agepop = sum(population)) %>%
  
  group_by(period, agegr, area_id5) %>%
  mutate(id5_agepop = sum(population)) %>%
  
  mutate(id2_agepop = ifelse(is.na(area_id2), NA, id2_agepop),
         id3_agepop = ifelse(is.na(area_id3), NA, id3_agepop),
         id4_agepop = ifelse(is.na(area_id4), NA, id4_agepop),
         id5_agepop = ifelse(is.na(area_id5), NA, id5_agepop),
  ) %>%
  filter(!is.na(area_name0))

new_iso3 <- pop_areas$iso3 %>% unique

############ Get ASFR, TFR, Births, Births by age at all administrative levels from model.


mod_results <- Map(get_mod_results, mod_list %>% .[new_iso3], asfr_pred_subnat_15_2020 %>% .[new_iso3], list(pop_areas), list(areas_long), list(population_age_female))

mod_results[[5]] <- get_mod_results(mod_list[["SWZ"]], asfr_pred_subnat_15_2020[["SWZ"]], pop_areas, areas_long, population_age_female)

debugonce(get_mod_results)

