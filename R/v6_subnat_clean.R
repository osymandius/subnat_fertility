library(countrycode)
library(tidyverse)
library(rdhs)
library(demogsurv)
library(INLA)
library(reshape2)
library(survival)
library(geojsonsf)
library(sf)
library(spdep)
library(parallel)
library(naomi)

library(here)

# naomi_data_path <- "~/naomi-data"
naomi_data_path <- "~/Documents/GitHub/naomi-data"

source(here("R/inputs.R"))
source(here("R/fertility_funs.R"))

iso3_current <- "MWI"
# iso3 <- c("LSO", "MOZ", "MWI", "NAM", "TZA", "UGA", "ZMB", "ZWE")
# 
list2env(make_areas_population(iso3_current, naomi_data_path), globalenv())

## set_rdhs_config(email="o.stevens@imperial.ac.uk", project="Subnational fertility", config_path = "~/.rdhs.json")

dhs_iso3 <- dhs_countries(returnFields=c("CountryName", "DHS_CountryCode")) %>%
  mutate(iso3 = countrycode(CountryName, "country.name", "iso3c"),
         iso3 = ifelse(CountryName == "Eswatini", "SWZ", iso3))

clusters <- readRDS(here("input_data/clusters_2019_11_21.rds")) %>%
  mutate(iso3 = survey_id) %>%
  separate(col="iso3", into="iso3", sep=3) %>%
  left_join(dhs_iso3 %>% select(-CountryName), by="iso3") %>%
  separate(survey_id, into=c(NA, "surv"), sep=3, remove=FALSE) %>%
  mutate(DHS_survey_id = paste0(DHS_CountryCode, surv)) %>%
  separate(surv, into=c(NA, "SurveyType"), sep=-3) %>%
  filter(iso3 == iso3_current, DHS_CountryCode != "OS") %>%
  filter(!survey_id %in% c("MOZ2009AIS", "TZA2003AIS", "UGA2011AIS")) %>%
  filter(survey_id == "MWI2015DHS")

## Get surveys for which we have clusters. Split into country list.
surveys <- dhs_surveys(surveyIds = unique(clusters$DHS_survey_id)) %>%
  left_join(clusters %>% 
              select(c(DHS_survey_id, survey_id, iso3)) %>% 
              distinct, 
            by=c("SurveyId" = "DHS_survey_id")) %>%
  filter(!SurveyId %in% c("MZ2009AIS", "TZ2003AIS", "UG2011AIS"),
         iso3 == iso3_current
  )

## Needs check to ensure level is < max_level
cluster_areas <- assign_cluster_area(clusters, 0)

dat <- clusters_to_surveys(surveys, cluster_areas, single_tips = FALSE)

asfr <- Map(calc_asfr1, dat$ir,
              y=1:length(dat$ir),
              by = list(~country + surveyid + survtype + survyear + area_id),
              tips = dat$tips_surv,
              agegr= list(3:10*5),
              period = list(1995:2017),
              counts = TRUE) %>%
  bind_rows %>%
  type.convert %>%
  filter(period<=survyear) %>%
  rename(age_group = agegr) %>%
  mutate(iso3 = countrycode(country, "country.name", "iso3c"),
         iso3 = ifelse(country == "Eswatini", "SWZ", iso3)) %>%
  select(-country)


mics_asfr <- Map(calc_asfr_mics, mics_data$wm, y=list(1),
                 by = list(~area_id + survyear + surveyid + survtype),
                 tips = list(c(0:5)),
                 agegr= list(3:10*5),
                 period = list(1995:2017),
                 counts = TRUE,
                 bhdata = mics_data$bh_df) %>%
  bind_rows

asfr_pred <- make_asfr_pred_df(asfr)

get_neighbourhood_structure(asfr, areas_long, boundaries)

##################################

dat <- get_asfr_pred_df("ZWE", area_level = 0, project = FALSE)

dat_list <- lapply(c("LSO", "MOZ", "MWI", "NAM", "TZA", "UGA", "ZMB", "ZWE"), get_asfr_pred_df, area_level = "naomi", project=FALSE)

####################### NATIONAL MODS

dat <- lapply(c("LSO", "MOZ", "MWI", "NAM", "TZA", "UGA", "ZMB", "ZWE"), function(iso3_current) {
  areas_long <- filter(areas_long, iso3==iso3_current)
  dat <- readRDS(here("countries", paste0(iso3_current, "/data/", iso3_current, "_asfr_admin0.rds"))) %>%
    filter(period<2016)
  
  return(dat)
})

names(dat) <- c("LSO", "MOZ", "MWI", "NAM", "TZA", "UGA", "ZMB", "ZWE")

mod_list <- lapply( c("LSO", "MOZ", "MWI", "NAM", "TZA", "UGA", "ZMB", "ZWE"),
                   function(x) readRDS(here("countries", paste0(x, "/mods/", x, "_nat_mod.rds"))))

mod_list <- Map(function(x, y) {
  
  formula <- births ~ 
    f(id.age_group, model="rw1") + 
    f(id.period, model="rw2") + 
    f(id.age_group2, model = "rw1", group = id.period, control.group = list(model = "rw2")) +
    tips_dummy + 
    f(id.tips, model="rw1")
  
  mod <- run_mod_nat(formula, y)
  
  # saveRDS(mod, paste0("countries/", x, "/mods/", x, "_nat_mod.rds"))
  
  return(mod)
}, x=c("LSO", "MOZ", "MWI", "NAM", "TZA", "UGA", "ZMB", "ZWE"), y=dat)

mod <- run_mod_nat(formula, dat)


res_list <- Map(function(x, y,  z) {
  
  res <- get_mod_results_test(z, y)
  
  # saveRDS(res, paste0("countries/", x, "/mods/", x, "_nat_res.rds"))
  
  return(res)
}, x=c("LSO", "MOZ", "MWI", "NAM", "TZA", "UGA", "ZMB", "ZWE"), y=dat, z=mod_list)

res_list %>%
  bind_rows %>%
ggplot(aes(x=period, y=median, color=age_group, group=age_group)) +
  geom_line() +
  facet_wrap(~iso3)

mod_list <- list()
mod_list[[1]] <- zwe.r

r_list <- lapply(mod_list, "[[", "summary.random") %>%
  lapply(function(x){
  x %>%
    lapply(function(x) mutate(x, ID = as.integer(ID))) %>%
    lapply(arrange, ID) %>%
    lapply("[[", "mean") %>%
    lapply(exp)
})

zwe.r$summary.fixed

exp(0.05003866)

names(r_list) <- c("LSO", "MOZ", "MWI", "NAM", "TZA", "UGA", "ZMB", "ZWE")


##########################

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

population_age_female <- load_population_agesex("~/Documents/GitHub/subnat_fertility/input_data/population_agesex_wide_wpp2019_raked.rds", "~/Documents/GitHub/subnat_fertility/input_data/areas_long.RDS") %>%
  filter(sex == "female", iso3 == "ZWE", iso3 != "UGA") %>%
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

