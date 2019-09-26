library(tidyverse)

setwd("~/Documents/GitHub/subnat_fertility")

source("fertility_funs.R")

saveRDS(national_pred, file="pred_df_national.rds")

asfr_pred_country_subnat <- readRDS("asfr_pred_country_subnat.rds")

mwi_mod <- readRDS("MWI_poisson_mod.rds")
lso_mod <- readRDS("LSO_poisson_mod.rds")
rwa_mod <- readRDS("RWA_poisson_mod.rds")
zwe_mod <- readRDS("ZWE_poisson_mod.rds")

mod_list <- list(mwi_mod, lso_mod, rwa_mod, zwe_mod)

asfr1_country_subnat <- asfr_pred_country_subnat %>%
  lapply(function(x) {
    x <- x %>%
      filter(!is.na(surveyid))
  })

subnational <- TRUE
multicountry <- FALSE

pred_list <- Map(get_pred, mod_list, asfr_pred_country_subnat, asfr1_country_subnat)

pred <- pred_list %>%
  bind_rows

### WORLD POP. Filter areas on cluster area. Convert MWI from level 5 to level 4.

worldpop <- read_csv("WorldPop_agesex.csv", col_types = cols(X1 = col_skip(), X = col_skip(), sex = col_character())) %>%
  filter(sex=="f")

iso3_code <- c("MWI", "LSO", "RWA", "ZWE")

areas_long <- readRDS("~/Documents/GitHub/naomi-data/data/areas/areas_long.rds")
areas_wide <- readRDS("~/Documents/GitHub/naomi-data/data/areas/areas_wide.rds")

worldpop_rep <- worldpop %>%
  left_join(areas_wide, by=c("iso3", "id" = "area_id")) %>%
  group_split(year) %>%
  lapply(function(x) {
    n <- nrow(x)
    year_start <- unique(x$year)
    x <- x %>%
      bind_rows(., ., ., ., .) %>%
      mutate(year = rep(year_start:(year_start+4), each=n))
    return(x)
  })  %>%
  bind_rows

pop_areas_wide <- worldpop_rep %>%
  filter(id %in% unique(clusters$geoloc_area_id)) %>%
  group_by(year, age, id0) %>%
  mutate(id0_agepop = sum(population))  %>%
  group_by(year, age, id1) %>%
  mutate(id1_agepop = sum(population)) %>%
  group_by(year, age, id2) %>%
  mutate(id2_agepop = sum(population)) %>%
  group_by(year, age, id3) %>%
  mutate(id3_agepop = sum(population)) %>%
  group_by(year, age, id4) %>%
  mutate(id4_agepop = sum(population)) %>%
  mutate(id2_agepop = ifelse(is.na(id2), NA, id2_agepop),
         id3_agepop = ifelse(is.na(id3), NA, id3_agepop),
         id4_agepop = ifelse(is.na(id4), NA, id4_agepop),
  ) %>%
  select(-c(name5, id5))


# pop_areas_wide <- worldpop %>%
#   left_join(areas_wide, by=c("iso3", "id" = "area_id")) %>%
#   filter(id %in% unique(clusters$geoloc_area_id)) %>%
#   group_by(year, age, id0) %>%
#   mutate(id0_agepop = sum(population))  %>%
#   group_by(year, age, id1) %>%
#   mutate(id1_agepop = sum(population)) %>%
#   group_by(year, age, id2) %>%
#   mutate(id2_agepop = sum(population)) %>%
#   group_by(year, age, id3) %>%
#   mutate(id3_agepop = sum(population)) %>%
#   group_by(year, age, id4) %>%
#   mutate(id4_agepop = sum(population)) %>%
#   group_by(year, id0) %>%
#   mutate(id0_totpop = sum(population))  %>%
#   group_by(year, id1) %>%
#   mutate(id1_totpop = sum(population)) %>%
#   group_by(year, id2) %>%
#   mutate(id2_totpop = sum(population)) %>%
#   group_by(year, id3) %>%
#   mutate(id3_totpop = sum(population)) %>%
#   group_by(year, id4) %>%
#   mutate(id4_totpop = sum(population)) %>%
#   mutate(id2_agepop = ifelse(is.na(id2), NA, id2_agepop),
#          id3_agepop = ifelse(is.na(id3), NA, id3_agepop),
#          id4_agepop = ifelse(is.na(id4), NA, id4_agepop),
#          id2_totpop = ifelse(is.na(id2), NA, id2_totpop),
#          id3_totpop = ifelse(is.na(id3), NA, id3_totpop),
#          id4_totpop = ifelse(is.na(id4), NA, id4_totpop)
#   ) %>%
#   select(-c(name5, id5))

pop_areas_wide_mwi <- pop_areas_wide %>%
  filter(iso3 == "MWI") %>%
  mutate(id = id4,
         name = name4,
         level = 4
  ) %>%
  group_by(year, age, id4) %>%
  mutate(population = sum(population)) %>%
  distinct()

pop_areas <- pop_areas_wide %>%
  filter(iso3 != "MWI") %>%
  bind_rows(pop_areas_wide_mwi)

############ Calculating aggregated ASFR and TFR against worldpop age/sex data

asfr_aggr <- pred %>%
  select(-id) %>%
  left_join(pop_areas %>% filter(startage>=15, startage<50, agespan==5), by=c("area_name" = "name", "period" = "year", "agegr" = "age")) %>%
  mutate(asfr_ratio0 = mean*(population/id0_agepop), asfr_ratio0_u = `0.975quant`*(population/id0_agepop), asfr_ratio0_l = `0.025quant`*(population/id0_agepop),
         asfr_ratio1 = mean*(population/id1_agepop), asfr_ratio1_u = `0.975quant`*(population/id1_agepop), asfr_ratio1_l = `0.025quant`*(population/id1_agepop),
         asfr_ratio2 = mean*(population/id2_agepop), asfr_ratio2_u = `0.975quant`*(population/id2_agepop), asfr_ratio2_l = `0.025quant`*(population/id2_agepop),
         asfr_ratio3 = mean*(population/id3_agepop), asfr_ratio3_u = `0.975quant`*(population/id3_agepop), asfr_ratio3_l = `0.025quant`*(population/id3_agepop),
         asfr_ratio4 = mean*(population/id4_agepop), asfr_ratio4_u = `0.975quant`*(population/id4_agepop), asfr_ratio4_l = `0.025quant`*(population/id4_agepop)
  ) %>%
  filter(period>1999)

asfr_admin0 <- asfr_aggr %>%
  group_by(period, agegr, id0) %>%
  summarise(val = sum(asfr_ratio0),
            lower = sum(asfr_ratio0_l),
            upper = sum(asfr_ratio0_u)
  ) %>%
  mutate(variable = "asfr") %>%
  rename(area_id = id0)

asfr_admin1 <- asfr_aggr %>%
  group_by(period, agegr, id1) %>%
  summarise(val = sum(asfr_ratio1),
            lower = sum(asfr_ratio1_l),
            upper = sum(asfr_ratio1_u)
  ) %>%
  mutate(variable = "asfr") %>%
  rename(area_id = id1)

asfr_admin2 <- asfr_aggr %>%
  group_by(period, agegr, id2) %>%
  summarise(val = sum(asfr_ratio2),
            lower = sum(asfr_ratio2_l),
            upper = sum(asfr_ratio2_u)
  ) %>%
  mutate(variable = "asfr") %>%
  rename(area_id = id2)

asfr_admin3 <- asfr_aggr %>%
  group_by(period, agegr, id3) %>%
  summarise(val = sum(asfr_ratio3),
            lower = sum(asfr_ratio3_l),
            upper = sum(asfr_ratio3_u)
  ) %>%
  mutate(variable = "asfr") %>%
  rename(area_id = id3)

asfr_admin4 <- asfr_aggr %>%
  group_by(period, agegr, id4) %>%
  summarise(val = sum(asfr_ratio4),
            lower = sum(asfr_ratio4_l),
            upper = sum(asfr_ratio4_u)
  ) %>%
  mutate(variable = "asfr") %>%
  rename(area_id = id4)

tfr_admin0 <- asfr_aggr %>%
  group_by(period, id0) %>%
  summarise(val = 5*sum(asfr_ratio0),
            lower = 5*sum(asfr_ratio0_l),
            upper = 5*sum(asfr_ratio0_u)
  ) %>%
  mutate(variable="tfr") %>%
  rename(area_id = id0)

tfr_admin1 <- asfr_aggr %>%
  group_by(period, id1) %>%
  summarise(val = 5*sum(asfr_ratio1),
            lower = 5*sum(asfr_ratio1_l),
            upper = 5*sum(asfr_ratio1_u)
  ) %>%
  mutate(variable="tfr") %>%
  rename(area_id = id1)

tfr_admin2 <- asfr_aggr %>%
  group_by(period, id2) %>%
  summarise(val = 5*sum(asfr_ratio2),
            lower = 5*sum(asfr_ratio2_l),
            upper = 5*sum(asfr_ratio2_u)
  ) %>%
  mutate(variable="tfr") %>%
  rename(area_id = id2)

tfr_admin3 <- asfr_aggr %>%
  group_by(period, id3) %>%
  summarise(val = 5*sum(asfr_ratio3),
            lower = 5*sum(asfr_ratio3_l),
            upper = 5*sum(asfr_ratio3_u)
  ) %>%
  mutate(variable="tfr") %>%
  rename(area_id = id3)

tfr_admin4 <- asfr_aggr %>%
  group_by(period, id4) %>%
  summarise(val = 5*sum(asfr_ratio4),
            lower = 5*sum(asfr_ratio4_l),
            upper = 5*sum(asfr_ratio4_u)
  ) %>%
  mutate(variable="tfr") %>%
  rename(area_id = id4)

model_asfr <- areas_long %>%
  select(-parent_area_id) %>%
  filter(iso3 %in% iso3_code) %>%
  left_join(asfr_admin0 %>% bind_rows(asfr_admin1, asfr_admin2, asfr_admin3, asfr_admin4), by="area_id") %>%
  mutate(source = "Model")

model_births_age <- model_asfr %>%
  select(-variable) %>%
  left_join(worldpop_rep %>% select(c(id, age, population, year)), by=c("agegr" = "age", "period" = "year", "area_id" = "id")) %>%
  mutate(val = val*population,
         lower = lower*population,
         upper = upper*population,
         variable = "births",
         source = "Model"
  ) %>%
  select(-population) 


model_tfr <- areas_long %>%
  select(-parent_area_id) %>%
  filter(iso3 %in% iso3_code) %>%
  left_join(tfr_admin0 %>% bind_rows(tfr_admin1, tfr_admin2, tfr_admin3, tfr_admin4), by="area_id") %>%
  mutate(source = "Model")


model_births <- model_births_age %>%
  group_by(iso3, area_id, area_name, area_level, period, source, variable) %>%
  summarise(val = sum(val),
            lower = sum(lower),
            upper = sum(upper))



worldpop_U1 <- read_csv("WorldPop_agesex.csv", col_types = cols(X1 = col_skip(), X = col_skip(), sex = col_character())) %>%
  filter(startage==0, iso3 %in% c("MWI", "LSO", "RWA", "ZWE")) %>%
  group_by(iso3, id, level, name, source, year, age, startage, agespan) %>%
  summarise(population = sum(population)) %>%
  ungroup %>%
  mutate(variable = "U1_pop") %>%
  rename(val = population, area_id = id, area_level = level, area_name = name, period = year) %>%
  select(iso3, area_id, area_level, area_name, source, period, source, variable, val)


######### National comparison of ASFR and TFR vs WPP 2019

wpp_asfr <- read_excel("WPP2019_FERT_F07_AGE_SPECIFIC_FERTILITY.xlsx")

wpp_asfr <- wpp_asfr %>%
  separate(period, into=c("period", NA), sep="-") %>%
  type.convert %>%
  mutate(iso3 = countrycode(country, "country.name", "iso3c"),
         area_name = countrycode(iso3, "iso3c", "country.name"),
         area_id = iso3
  ) %>%
  select(-c(countrycode, country)) %>%
  melt(id=c("iso3", "area_name", "area_id", "period"), variable.name="agegr", value.name = "val") %>%
  filter(iso3 %in% iso3_code) %>%
  group_split(period) %>%
  lapply(function(x) {
    n <- nrow(x)
    year_start <- unique(x$period)
    x <- x %>%
      bind_rows(., ., ., ., .) %>%
      mutate(period = rep(year_start:(year_start+4), each=n))
    return(x)
  })  %>%
  bind_rows %>%
  mutate(val = as.numeric(val)/1000,
         source = "WPP2019",
         variable = "asfr",
         area_level = 0)

wpp_tfr <- read_excel("WPP2019_FERT_F04_TOTAL_FERTILITY(1).xlsx")

wpp_tfr <- wpp_tfr %>%
  melt(id="area_name", variable.name="period", value.name = "val") %>%
  separate(period, into=c("period", NA), sep="-") %>%
  type.convert() %>%
  mutate(iso3 = countrycode(area_name, "country.name", "iso3c")) %>%
  filter(iso3 %in% iso3_code) %>%
  group_split(period) %>%
  lapply(function(x) {
    n <- nrow(x)
    year_start <- unique(x$period)
    x <- x %>%
      bind_rows(., ., ., ., .) %>%
      mutate(period = rep(year_start:(year_start+4), each=n))
    return(x)
  })  %>%
  bind_rows %>%
  mutate(source = "WPP2019",
         area_level = 0,
         area_id = iso3,
         variable = "tfr")


######### GBD 2017 national ASFR and TFR

gbd_asfr <- read_csv("IHME_GBD_2017_FERT_ESTIMATES_1950_2017_Y2018M11D08.CSV", col_types = cols(age_group_id = col_skip(), location_id = col_skip(), measure_id = col_skip(), measure_name = col_skip(), metric_name = col_skip(), sex_id = col_skip(), sex_name = col_skip()))

gbd_asfr <- gbd_asfr %>%
  mutate(iso3 = countrycode(location_name, "country.name", "iso3c")) %>%
  filter(iso3 %in% iso3_code) %>%
  separate(age_group_name, into=c("agegr", "up"), sep=" to ") %>%
  mutate(agegr = paste0(agegr, "-", up),
         area_id = iso3,
         area_level = 0,
         source = "GBD2017",
         variable = "asfr") %>%
  select(-up) %>%
  rename(area_name = location_name, period = year_id) 


gbd_tfr <- gbd_asfr %>%
  group_by(iso3, area_name, area_id, area_level,  period, source) %>%
  summarise(val = 5*sum(val),
            upper = 5*sum(upper),
            lower = 5*sum(lower),
  ) %>%
  mutate(variable = "tfr")


#### GBD Under 1 population

gbd_2000_pop <- read_csv("IHME_GBD_2017_POP_2000_2004_Y2018M11D08.CSV", col_types = cols(location_id = col_skip(), measure_id = col_skip(), measure_name = col_skip(), metric_name = col_skip()))

gbd_2005_pop <- read_csv("IHME_GBD_2017_POP_2005_2009_Y2018M11D08.CSV", col_types = cols(location_id = col_skip(), measure_id = col_skip(), measure_name = col_skip(), metric_name = col_skip()))

gbd_2010_pop <- read_csv("IHME_GBD_2017_POP_2010_2014_Y2018M11D08.CSV", col_types = cols(location_id = col_skip(), measure_id = col_skip(), measure_name = col_skip(), metric_name = col_skip()))

gbd_2015_pop <- read_csv("IHME_GBD_2017_POP_2015_2017_Y2018M11D08.CSV", col_types = cols(location_id = col_skip(), measure_id = col_skip(), measure_name = col_skip(), metric_name = col_skip()))

gbd_2000_pop <- gbd_2000_pop %>%
  filter(age_group_id == 28, sex_id ==3)

gbd_2005_pop <- gbd_2005_pop %>%
  filter(age_group_id == 28, sex_id ==3)

gbd_2010_pop <- gbd_2010_pop %>%
  filter(age_group_id == 28, sex_id ==3)

gbd_2015_pop <- gbd_2015_pop %>%
  filter(age_group_id == 28, sex_id ==3)

gbd_U1 <- gbd_2000_pop %>%
  bind_rows(gbd_2005_pop, gbd_2010_pop, gbd_2015_pop) %>%
  select(-c(sex_id, sex_name, age_group_id, age_group_name)) %>%
  rename(area_name = location_name, period = year_id) %>%
  mutate(iso3 = countrycode(area_name, "country.name", "iso3c"),
         area_level = 0,
         area_id = iso3,
         source = "GBD2017",
         variable = "U1_pop") %>%
  filter(iso3 %in% iso3_code)

## Spectrum TFR and ASFR

# spec_tfr <- data.frame(
#   "iso3" = iso3_code,
#   "area_name" = countrycode(iso3_code, "iso3c", "country.name"),
#   "area_id" = iso3_code,
#   "area_level" = 0,
#   "period" = 1970:2023,
#   "source" = "Spectrum18",
#   "variable" = "tfr",
#   "val" = c(5.8,5.8,5.8,5.79817,5.78494,5.76197,5.4,5.39,5.38,5.37,5.36,5.35,5.34,5.33,5.32,5.31,5.3,5.18,5.06,4.94,4.82,4.7,4.58,4.46,4.34,4.22,4.1,4.04,3.98,3.92,3.86,3.8,3.74,3.68,3.62,3.56,3.5,3.47,3.44,3.41,3.38,3.35,3.32,3.29,3.26,3.23,3.2,3.03699,2.99126,2.94807,2.90641,2.86637,2.82802,2.79136)
# )
# 
# spec_asfr <- data.frame(
#   "iso3" = iso3_code,
#   "area_name" = countrycode(iso3_code, "iso3c", "country.name"),
#   "area_id" = iso3_code,
#   "area_level" = 0,
#   "period" = 1970:2023,
#   "source" = "Spectrum18",
#   "variable" = "asfr_prop",
#   "1" = c(8.1479,8.1479,8.1479,8.1479,8.1479,8.1478,8.1478,8.1478,8.1461,8.1357,8.1239,8.1197,8.1324,8.1687,8.2226,8.296,8.3933,8.5189,8.6851,8.9356,9.2529,9.6057,9.9553,10.2659,10.5748,10.8841,11.1851,11.469,11.7211,11.9185,12.0825,12.2432,12.4341,12.6835,12.9568,13.2475,13.5302,13.7456,13.8346,13.8479,13.8246,13.795,13.7895,13.8424,13.9857,14.1927,14.43,14.6494,14.8215,14.9777,15.1261,15.2652,15.3944,15.5128),
#   
#   "2" = c(22.8469,22.8469,22.8469,22.8469,22.8469,22.8469,22.8469,22.8469,22.8462,22.8426,22.8386,22.8371,22.8416,22.855,22.8811,22.9218,22.9792,23.0554,23.1577,23.3076,23.4871,23.6751,23.8565,24.023,24.1825,24.3364,24.49,24.6505,24.8232,24.9948,25.1557,25.3028,25.4421,25.608,25.8925,26.2196,26.5019,26.6716,26.6962,26.6953,26.6937,26.6928,26.6941,26.7058,26.762,26.8478,26.9415,27.024,27.0829,27.1341,27.1803,27.2211,27.2557,27.2841),
#   
#   "3" = c(22.8568,22.8568,22.8568,22.8568,22.8568,22.8568,22.8568,22.8568,22.8575,22.8617,22.8666,22.8683,22.8631,22.8486,22.8285,22.8016,22.765,22.7159,22.6499,22.5554,22.4368,22.3052,22.1766,22.0665,21.9549,21.8426,21.7398,21.6588,21.6306,21.747,21.9488,22.1568,22.2972,22.3265,22.3453,22.3511,22.3338,22.3023,22.2754,22.2694,22.2769,22.2875,22.2909,22.2787,22.2577,22.2314,22.2028,22.1769,22.1602,22.1473,22.1378,22.1319,22.1301,22.1327),
#   
#   "4" = c(19.1198,19.1198,19.1198,19.1198,19.1198,19.1198,19.1198,19.1199,19.1199,19.1207,19.1216,19.122,19.121,19.119,19.1216,19.1277,19.1354,19.1428,19.1507,19.1715,19.1951,19.211,19.214,19.2072,19.2064,19.1991,19.1719,19.1131,18.9922,18.7128,18.3337,17.936,17.6005,17.3955,17.2805,17.2114,17.1636,17.1277,17.1027,17.0976,17.1052,17.1153,17.118,17.1028,17.065,17.0127,16.9558,16.9055,16.8721,16.8452,16.8229,16.8063,16.7963,16.7932),
#   
#   "5" = c(15.5213,15.5213,15.5213,15.5213,15.5213,15.5213,15.5212,15.5212,15.522,15.5262,15.531,15.5327,15.5274,15.5124,15.488,15.4531,15.4062,15.345,15.2644,15.1455,14.9964,14.8323,14.6712,14.5291,14.3892,14.2495,14.115,13.9914,13.9009,13.9159,13.9817,14.0253,13.9788,13.7653,13.3217,12.7764,12.2736,11.9457,11.8741,11.8701,11.8817,11.8948,11.8953,11.8573,11.719,11.511,11.2769,11.0591,10.8878,10.7309,10.581,10.438,10.302,10.1732),
#   
#   "6" = c(8.5121,8.5121,8.5121,8.5121,8.5121,8.5121,8.5121,8.5121,8.5125,8.515,8.5179,8.5191,8.5158,8.5065,8.4908,8.468,8.4379,8.399,8.349,8.2774,8.1879,8.0876,7.9858,7.8918,7.804,7.7153,7.6202,7.515,7.3918,7.2377,7.0761,6.9325,6.8301,6.7865,6.7685,6.7606,6.7568,6.7556,6.7565,6.7577,6.7586,6.7586,6.7572,6.7531,6.7426,6.7264,6.7071,6.687,6.6687,6.6502,6.6312,6.6112,6.5902,6.5678),
#   
#   "7" = c(2.9955,2.9955,2.9955,2.9955,2.9955,2.9955,2.9955,2.9954,2.9958,2.9979,3.0004,3.0014,2.9987,2.9898,2.9676,2.9315,2.883,2.8227,2.7432,2.6069,2.4437,2.283,2.1405,2.0163,1.8883,1.7731,1.6782,1.6024,1.5402,1.4733,1.4215,1.4036,1.4172,1.4347,1.4348,1.4335,1.4401,1.4515,1.4606,1.4622,1.4594,1.4558,1.455,1.46,1.4681,1.4778,1.4882,1.4983,1.5071,1.5144,1.5209,1.5263,1.5314,1.5362)
# 
# ) %>%
#   melt(id=c("iso3", "area_id", "area_level", "area_name", "period", "variable", "source"), variable.name="agegr", value.name="val") %>%
#   mutate(agegr = rep(c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49"), each=54))
# 
# 
# 
# asfr_prop <- model_asfr %>%
#   filter(area_level ==0) %>%
#   group_by(period) %>%
#   mutate(val = 100*(val/sum(val)),
#          variable = "asfr_prop") %>%
#   bind_rows(spec_asfr)
# 
# asfr_prop %>%
#   ggplot(aes(x=period, y=val, group=agegr, color=agegr, linetype=source)) +
#   geom_line(data=asfr_prop %>% filter(source =="Model", period>1999)) +
#   geom_line(data=asfr_prop %>% filter(source =="Spectrum18", period>1999)) +
#   labs(y="Age distribution of fertility (%)", title="Age distribution of fertility | Malawi")


## Combine

all_age_outputs <- model_births %>%
  bind_rows(model_tfr, worldpop_U1, wpp_tfr, gbd_U1, gbd_tfr, spec_tfr) %>%
  filter(!is.na(val))

age_specific_outputs <- model_asfr %>%
  bind_rows(model_births_age, wpp_asfr, gbd_asfr, asfr_prop) %>%
  filter(!is.na(val))

saveRDS(age_specific_outputs, file="age_specific_outputs.rds")
saveRDS(all_age_outputs, file="allage_outputs.rds")

age_specific_outputs <- readRDS(file="age_specific_outputs.rds")
all_age_outputs <- readRDS(file="allage_outputs.rds")



####### MAPS

all_age_outputs %>%
  filter(iso3 =="LSO", area_level %in% c(0, 1, 2), variable=="tfr", period==2015, source =="Model") %>%
  left_join(boundaries) %>%
  ggplot() +
    geom_sf(aes(geometry = geometry, fill=val)) +
    facet_wrap(~area_level) +
    labs(title="TFR by admin level in 2015 | Lesotho")


ggplot(data=data, aes(x=period, y=val, group=agegr, color=agegr, linetype=source)) +
  geom_line(data=age_specific_outputs %>% filter(variable == "asfr_prop", source =="Model", period>1999, area_level==0)) +
  geom_line(data=age_specific_outputs %>% filter(variable == "asfr_prop", source =="Spectrum18", period>1999, area_level==0)) +
  labs(y="Age distribution of fertility (%)", title="Age distribution of fertility") +
  facet_wrap(~area_name)
  

all_age_outputs %>%
  filter((variable=="U1_pop" | variable =="births"), area_level==0, period>1999) %>%
  ggplot(aes(x=period, y=val, group=source, fill=source)) +
    geom_line(aes(color=source)) +
    labs(y="U1 pop or births", title="U1 pop [WorldPop & GBD] compared to model births") +
    facet_wrap(~area_name, scales="free") 

age_specific_outputs %>%
  filter(variable=="asfr", period %in% c(2000, 2005, 2010, 2015), area_level==0, agegr != "10-14", agegr != "50-54") %>%
  ggplot(aes(x=agegr, y=val, group = source, color=source)) +
    geom_line()+
    facet_grid(period ~ area_name, scales="free") +
  labs(y="ASFR")

age_specific_outputs %>%
  filter(variable == "asfr", area_level == 0, source == "Model", agegr =="20-24") %>%
  ggplot(aes(x=period, y=val, group=agegr, color=agegr)) +
    geom_line() +
    geom_point(data=asfr1_country %>% bind_rows %>% rename(area_name = country) %>% filter(agegr=="20-24") %>% type.convert, aes(y=asfr, group=surveyid, color=surveyid))+
    facet_wrap(~area_name) +
  labs(y="ASFR")

age_specific_outputs %>%
  filter(variable == "asfr", area_level == 0, source == "Model") %>%
  ggplot(aes(x=period, y=val, group=agegr, color=agegr)) +
  geom_line() +
  geom_line(data=pred_df_national %>% bind_rows %>% rename(area_name = country) %>% filter(period>1999), aes(y=`0.5quant`), linetype=2) +
  facet_wrap(~area_name) +
  labs(y="ASFR")

all_age_outputs %>%
  filter((variable=="U1_pop" | variable == "births"), ((area_level==2 & iso3 != "MWI") | (area_level==4 & iso3 == "MWI")), period==2015) %>%
  group_by(iso3, area_id, area_name) %>%
  summarise(births_to_U1 = val[variable=="births"]/val[variable=="U1_pop"]) %>%
  group_by(iso3) %>%
  mutate(mean_diff = mean(births_to_U1)) %>%
  ggplot() +
    geom_point(aes(x=area_name, y=births_to_U1)) +
    geom_hline(aes(yintercept = 1), linetype=2) +
    geom_hline(aes(yintercept = mean_diff), color="red") +
    labs(x="", y="Births:U1 population", title="Ratio of model births to WorldPop U1 pop") +
    coord_flip() +
    facet_wrap(~iso3, ncol=2)

all_age_outputs %>%
  filter((variable=="U1_pop" | variable == "births"), ((area_level==2 & iso3 != "MWI") | (area_level==4 & iso3 == "MWI")), period==2015) %>%
  group_by(iso3, area_id, area_name) %>%
  summarise(births_to_U1 = val[variable=="births"]/val[variable=="U1_pop"]) %>%
  group_by(iso3) %>%
  mutate(mean_diff = mean(births_to_U1)) %>%
  distinct(mean_diff)

  ggplot() +
    geom_col(aes(reorder(area_name, val, sum), val, group=source, fill=source), position = position_dodge()) +
    coord_flip()
    
View(all_age_outputs %>%
       filter(variable != "U1_pop") %>%
       bind_rows(worldpop_U1) %>%
       filter((variable=="U1_pop" | variable == "births"), area_level==2, iso3 == "RWA", period==2015)
     )

