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
devtools::load_all("~/Documents/GitHub/naomi")

setwd("~/Documents/GitHub/subnat_fertility")
source("fertility_funs.R")

#boundaries_old <- readRDS("~/Documents/GitHub/naomi-data-old/data/area_boundaries.rds")
boundaries <- readRDS("~/Documents/GitHub/subnat_fertility/area_boundaries.rds")

iso3_codes <- c("LSO", "MOZ", "MWI", "NAM", "SWZ", "TZA", "UGA", "ZMB", "ZWE")

paths <- paste0("~/Documents/GitHub/naomi-data/", iso3_codes, "/data/", tolower(iso3_codes), "_areas.geojson")
area_cols <- c("iso3", "area_id", "area_name", "area_level", "parent_area_id", "naomi_level")

areas_long <- lapply(paths, read_sf) %>% lapply(function(x) {
  iso3_code <- x %>%
    filter(area_level == 0) %>%
    select(area_id) %>%
    unique %>%
    .$area_id
  
  x <- x %>%
    mutate(iso3 = iso3_code) %>%
    st_drop_geometry() %>%
    select(area_cols)
  
  return(x)
}) %>% 
  bind_rows

#areas_long <- readRDS("~/Documents/GitHub/subnat_fertility/areas_long.RDS")

areas_wide <- lapply(paths, read_sf) %>%
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
    filter(CountryName == "Malawi") %>%
    group_split(CountryName)
  
  areas <- areas_long %>%
    inner_join(clusters, by=c("area_id" = "geoloc_area_id", "iso3")) %>%
    left_join(areas_wide %>% select(area_id, area_id2)) %>%
    select(-c(area_id, area_name, area_level, parent_area_id, naomi_level)) %>%
    rename(area_id = area_id2)
  
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
  
  tips_surv <- list("DHS" = c(0,15), "MIS" = c(0,5), "AIS" = c(0,5))[survey_type$survtype]
  tips_surv <- list("DHS" = c(0:15), "MIS" = c(0:5), "AIS" = c(0:5))[surveys[[1]]$SurveyType]

  asfr <- Map(calc_asfr1, ir,
              y=1:length(ir),
              by = list(~country + surveyid + survtype + survyear),
              tips = tips_surv,
              #tips = list(c(0,10)),
              agegr= list(3:10*5),
              period = list(1995:2017),
              counts = TRUE)
  
  asfr <- asfr %>%
    bind_rows %>%
    type.convert %>%
    filter(period<=survyear) %>%
    group_split(country, keep=TRUE) %>%
    lapply(droplevels)

  asfr_pred <- list()
  
  asfr_pred <- asfr[[1]] %>%
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
           id = 1:nrow(.)
    ) %>%
    mutate_if(is.factor, as.character)
  
moz_tips_ar1_mod <- run_mod4(asfr_pred[[1]])

asfr_tza_admin1 <- lapply(asfr1_country, function(asfr1_country) {
  
  iso3_code <- ifelse(unique(asfr1_country$country) == "Eswatini", "SWZ",
                      unique(countrycode(asfr1_country$country, "country.name", "iso3c")))
  
  area_df <- areas_long %>% filter(iso3 == iso3_code, area_level ==  2)
  
  pred_df <- crossing(country = asfr1_country$country, area_id = area_df$area_id, period = min(asfr1_country$period):2020, agegr = asfr1_country$agegr,  pys=1)
  
  pred_df <- pred_df %>%
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
  
  return(pred_df)
})


saveRDS(asfr_pred_subnat_15_2020, file= "~/Documents/GitHub/subnat_fertility/asfr_pred_subnat_15_2020_NEW.rds")

max_level <- areas_long %>% filter(naomi_level == TRUE) %>% select(iso3, area_level) %>% unique

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
 
mwi_comparison <- run_mod(formulae = births ~ f(id.period, model = "rw1") + f(id.agegr, model="rw1") + f(id.district, model="bym", graph="MWI.adj"), asfr_pred = asfr, model_family = "poisson")

asfr_pred_subnat_15_2020 <- readRDS("~/Documents/GitHub/subnat_fertility/asfr_pred_subnat_15_2020_NEW.rds")
mod_list <- readRDS("2019_12_2_mods.rds")
mod_results <- readRDS("2019_12_3_results.rds")

iso3_code <- asfr_pred_subnat_15_2020 %>%
  lapply("[", "country") %>%
  lapply(unique) %>%
  bind_rows %>%
  mutate(country = ifelse(country== "Eswatini", "SWZ", countrycode(country, "country.name", "iso3c"))) %>%
  .$country

mod_list[9] <- NULL
iso3_code <- names(mod_list)

calendar_quarter_targets <- convert_calendar_quarter(c(2000:2020), 2)

population_age_female <- load_population_agesex("~/Documents/GitHub/subnat_fertility/population_agesex_wide_wpp2019_raked.rds", "~/Documents/GitHub/subnat_fertility/areas_long.RDS") %>%
  filter(sex == "female", iso3 %in% iso3_code, iso3 != "UGA") %>%
  mutate(calendar_quarter = quarter_id_to_calendar_quarter(quarter_id)) %>%
  left_join(get_age_groups()) %>%
  interpolate_fertility_population(., calendar_quarter_targets) %>%
  left_join(areas_long) %>%
  filter(naomi_level == TRUE)

### Add Uganda back in from NSO numbers..

uga_pop <- read.csv("~/Documents/GitHub/naomi-data/UGA/data/uga_population_ubos.csv") %>%
  mutate(period = year_labels(calendar_quarter_to_quarter_id(calendar_quarter))) %>%
  select(-c(source, calendar_quarter))

uga_pop_raked <- readRDS("uga_pop_raked.rds") 

uga_scale_pop <- crossing(area_id = uga_pop_raked $area_id, sex = "female", period = 2000:2015, age_group_id = 1:17) %>%
  left_join(wpp19 %>% 
              filter(iso3 == "UGA", sex=="female") %>% 
              group_by(iso3, year, age_group_id) %>% 
              summarise(wpp_age_population = sum(population)),
            by = c("period" = "year", "age_group_id")) %>%
  left_join(uga_pop_raked %>% 
              filter(period == 2015, sex=="female") %>%
              group_by(age_group_id, sex) %>%
              mutate(nat_age_pop = sum(population),
                     area_age_ratio = population/nat_age_pop) %>%
              ungroup %>%
              select(area_id, age_group_id, area_age_ratio)) %>%
  mutate(population = wpp_age_population *  area_age_ratio)

uga_reconstructed_pop <- uga_scale_pop %>%
  bind_rows(uga_pop_raked %>% filter(period >2015, sex=="female")) %>%
  left_join(get_age_groups() %>% select(age_group_id, age_group_label)) %>%
  rename(agegr = age_group_label)

### THIS IS A BODGE WHILE UGA AREA HIERARCHY IS WRONG
uga_reconstructed_pop$area_id <- gsub("_5_", "_3_", uga_reconstructed_pop$area_id)

population_age_female <- population_age_female %>%
  bind_rows(uga_reconstructed_pop %>%
  select(intersect(colnames(uga_reconstructed_pop), colnames(population_age_female))) %>%
  left_join(areas_long) %>%
  select(-parent_area_id))

###

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

############ Get ASFR, TFR, Births, Births by age at all administrative levels from model.

asfr_pred_subnat_15_2020 <- readRDS("~/Documents/GitHub/subnat_fertility/asfr_pred_subnat_15_2020_NEW.rds")

pred_df <- asfr %>%
  filter(is.na(surveyid)) %>%
  select(area_id, period, agegr,  pys) %>%
  mutate(surveyid = NA) %>%
  bind_rows(asfr %>%
              filter(!is.na(surveyid)) %>%
              group_by(area_id, period, agegr) %>%
              summarise(births = sum(births),
                        pys = sum(births)
              ) %>%
              ungroup %>%
              mutate(surveyid = 1)
  ) %>%
  mutate(id.period = group_indices(., period),
         id.period2 = id.period,
         id.period3 = id.period,
         id.agegr = group_indices(., agegr),
         id.agegr2 = id.agegr,
         id.agegr3 = id.agegr,
         id.agegr.period = group_indices(., period, agegr),
         id.district = group_indices(., area_id),
         id.district2 = id.district,
         id.district3 = id.district,
         id = 1:nrow(.)
  )


mod_comp <- get_mod_results(mwi_test_mod, pred_df %>% mutate(country = "Malawi"), pop_areas, areas_long, population_age_female)

mod_results <- Map(get_mod_results, mod_list, asfr_pred_subnat_15_2020[names(mod_list)], list(pop_areas), list(areas_long), list(population_age_female))

nam_fe_res <- get_mod_results(int, asfr_pred_subnat_15_2020[["NAM"]], pop_areas, areas_long, population_age_female)

View(mods[["MWI"]])

mwi_cub_res[["asfr"]] %>%
  filter(naomi_level) %>%
  ggplot(aes(x=period, y=median, group = agegr, color=agegr)) +
    geom_line()+
    facet_wrap(~area_id)

mod_results[["MWI"]]$asfr %>%
  filter(naomi_level) %>%
  ggplot(aes(x=period, y=median, group = agegr, color=agegr)) +
  geom_line()+
  facet_wrap(~area_id)

saveRDS(mod_results, "2019_12_3_results.rds")

nam <- names(mod_results)

asfr1_admin1 <- readRDS("asfr_admin1_newh.rds")

tfr_maps <- lapply(nam, function(x) {
  mod_results[[x]]$tfr %>%
    filter(naomi_level, period== 2019) %>%
    left_join(boundaries) %>%
    st_as_sf %>%
    ggplot() +
      geom_sf(aes(fill=median, geometry = geometry)) +
      labs(title=paste("2019 TFR |", x)) +
      scale_fill_viridis()
  
})

tfr_trend <- lapply(nam, function(x) {
  
  naomi_level_number <- mod_results[[x]]$tfr %>%
    ungroup %>%
    select(naomi_level, area_level) %>%
    filter(naomi_level) %>%
    unique %>%
    .$area_level
  
  mod_results[[x]]$tfr %>%
    filter(area_level == naomi_level_number -1) %>%
    ggplot() +
    geom_line(aes(x=period, y=median)) +
    labs(title=paste("TFR |", x)) +
    facet_wrap(~area_id)
  
})


tfr_plot <- lapply(nam, function(x) {
  mod_results[[x]]$tfr %>%
    ungroup %>%
    mutate(area_num = as.numeric(as.factor(area_id))) %>%
    filter(naomi_level, area_num < 11) %>%
    ggplot(aes(x=period)) +
    geom_line(aes(y=median)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5)+
    #geom_point(data=asfr1_admin1[[x]] %>% bind_rows %>% filter(agegr %in% c("20-24")), aes(y=asfr)) +
    labs(title=paste("TFR |", x), y="TFR", x="Year") +
    facet_wrap(~area_id)
  
})

tfr_plot[[4]]

nam <- names(mod_results)

mod_results[[2]]$tfr %>%
  ungroup %>%
  filter(naomi_level) %>%
  mutate(area_num = as.numeric(as.factor(area_id))) %>%
  filter(naomi_level, area_num < 11) %>%
  ggplot(aes(x=period)) +
  geom_line(aes(y=median)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5)+
  #geom_point(data=asfr1_admin1[[x]] %>% bind_rows %>% filter(agegr %in% c("20-24")), aes(y=asfr)) +
  labs(title=paste("TFR |", x), y="TFR", x="Year") +
  facet_wrap(~area_id)



asfr_pl <- mod_results[["TZA"]]$asfr %>%
  filter(area_level == 2, agegr %in% c("20-24")) %>%
  ggplot(aes(x=period)) +
  geom_line(aes(y=median)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5)+
  geom_point(data=asfr1_admin1[["TZA"]] %>% bind_rows %>% filter(agegr %in% c("20-24")), aes(y=asfr)) +
  labs(title=paste("ASFR in 20-24 |", "TZA")) +
  facet_wrap(~area_id) +
  ylim(0,0.5)

pwalk(list(paste0(nam, "_tfr_map.png"), tfr_maps), ggsave, path = getwd())
pwalk(list(paste0(nam, "_asfr_overlay.png"), asfr_data_overlay), ggsave, path=getwd())

plot5 <- lapply(nam, function(x) {
  mod_results[[x]]$tfr %>%
    filter(naomi_level) %>%
    ggplot(aes(x=period, y=median)) +
      geom_line() +
      facet_wrap(~area_id) +
      labs(title=paste("TFR trends |", x))
})

plot6 <- lapply(nam, function(x) {
  
  naomi_level_number <- mod_results[[x]]$tfr %>%
    ungroup %>%
    select(naomi_level, area_level) %>%
    filter(naomi_level) %>%
    unique %>%
    .$area_level
  
  mod_results[[x]]$asfr %>%
    filter(area_level == naomi_level_number - 1) %>%
    ggplot(aes(x=period, y=median, group=agegr, color=agegr)) +
    geom_line() +
    facet_wrap(~area_name) +
    labs(title=paste("ASFR trends |", x, "admin-1"), x="Year", y="ASFR")
})

plot6[[3]]

pwalk(list(paste0(nam, "_asfr_admin1.png"), plot6), ggsave, path = getwd())

