library(tidyverse)

setwd("~/Documents/GitHub/subnat_fertility")

source("fertility_funs.R")

saveRDS(national_pred, file="pred_df_national.rds")

asfr_pred_country_subnat <- readRDS("asfr_pred_country_subnat.rds")

mod_list <- list(readRDS("MWI_poisson_mod.rds"), readRDS("LSO_poisson_mod.rds"), readRDS("RWA_poisson_mod.rds"), readRDS("ZWE_poisson_mod.rds"), readRDS("UGA_poisson_mod.rds"))

### WORLD POP. Filter areas on cluster area. Convert MWI from level 5 to level 4.

worldpop <- read_csv("WorldPop_agesex.csv", col_types = cols(X1 = col_skip(), X = col_skip(), sex = col_character())) %>%
  filter(sex=="f")

iso3_code <- c("MWI", "RWA", "ZWE", "LSO", "UGA")

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

############ Get ASFR, TFR, Births, Births by age at all administrative levels from model.

debugonce(get_mod_results)
mod_results <- Map(get_mod_results, mod_list, asfr_pred_country_subnat, list(pop_areas), list(areas_long), list(worldpop_rep))

########

worldpop_U1 <- read_csv("WorldPop_agesex.csv", col_types = cols(X1 = col_skip(), X = col_skip(), sex = col_character())) %>%
  filter(startage==0, iso3 %in% iso3_code) %>%
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

spec_tfr_uga <- data.frame(
  "iso3" = "UGA",
  "area_name" = countrycode("UGA", "iso3c", "country.name"),
  "area_id" = "UGA",
  "area_level" = 0,
  "period" = 1970:2022,
  "source" = "Spectrum18",
  "variable" = "tfr",
  "val" = c(7.11296,7.10704,7.10208,7.1,7.1,7.1,7.1,7.1,7.1,7.1,7.1,7.1,7.1,7.1,7.1,7.1,7.1,7.1,7.1,7.09795,7.09237,7.08412,7.07404,7.063,7.04837,7.02792,7.00308,6.9753,6.946,6.91446,6.87884,6.83888,6.79435,6.745,6.68697,6.61854,6.54262,6.46213,6.38,6.29371,6.2008,6.1043,6.00723,5.9126,5.8202,5.72801,5.63635,5.54554,5.4559,5.36733,5.27958,5.19279,5.10705)
)

spec_asfr_uga <- data.frame(
  "iso3" = "UGA",
  "area_name" = countrycode("UGA", "iso3c", "country.name"),
  "area_id" = "UGA",
  "area_level" = 0,
  "period" = 1970:2022,
  "source" = "Spectrum18",
  "variable" = "asfr_prop",
  "1" = c(12.8001,12.8001,12.8001,12.8001,12.8001,12.8001,12.8001,12.8001,12.8001,12.8001,12.8001,12.8001,12.8001,12.8001,12.8133,12.8505,12.9068,12.9781,13.0601,13.2304,13.5174,13.8405,14.1143,14.2499,14.2431,14.1658,14.0438,13.8997,13.7499,13.5892,13.4062,13.206,12.9969,12.79,12.544,12.2877,12.073,11.8972,11.7501,11.6315,11.5492,11.5223,11.5444,11.5107,11.3738,11.2132,11.0432,10.876,10.7188,10.5706,10.4268,10.2865,10.1489),

  "2" = c(22.2601,22.2601,22.2601,22.2601,22.2601,22.2601,22.2601,22.2601,22.2601,22.2601,22.2601,22.2601,22.2601,22.2601,22.2651,22.2795,22.3012,22.3286,22.36,22.4245,22.5278,22.6375,22.7343,22.81,22.9315,23.1361,23.3736,23.5981,23.77,23.9296,24.1298,24.3526,24.5793,24.79,24.9986,25.2114,25.4154,25.613,25.81,25.9993,26.1703,26.3188,26.4476,26.5895,26.7592,26.9319,27.1029,27.2681,27.4244,27.5723,27.7129,27.847,27.9747),

  "3" = c(21.8301,21.83,21.83,21.83,21.83,21.83,21.83,21.83,21.83,21.83,21.83,21.83,21.83,21.83,21.8342,21.846,21.8636,21.8855,21.9101,21.969,22.0653,22.1592,22.223,22.24,22.1548,21.9714,21.7639,21.6092,21.59,21.7102,21.9004,22.1362,22.3917,22.64,22.9097,23.2173,23.5332,23.8422,24.1299,24.3987,24.6588,24.904,25.1348,25.3804,25.653,25.9298,26.2069,26.4808,26.7496,27.0122,27.2693,27.5216,27.77),

  "4" = c(18.71,18.71,18.71,18.71,18.71,18.71,18.71,18.71,18.71,18.71,18.71,18.71,18.71,18.71,18.7065,18.6966,18.6817,18.6623,18.64,18.6007,18.5349,18.4541,18.3795,18.3401,18.359,18.4222,18.5045,18.5852,18.65,18.6865,18.6985,18.697,18.6939,18.7,18.7179,18.7273,18.721,18.7032,18.6801,18.6543,18.623,18.5806,18.5271,18.4828,18.458,18.438,18.4201,18.4026,18.3842,18.366,18.3491,18.3337,18.3192),

  "5" = c(15.36,15.36,15.36,15.36,15.36,15.36,15.36,15.36,15.36,15.36,15.36,15.36,15.36,15.36,15.3398,15.2839,15.1992,15.0922,14.97,14.7025,14.2514,13.7583,13.3562,13.17,13.1636,13.2022,13.2601,13.3153,13.3499,13.3303,13.249,13.1309,13.0002,12.8799,12.7698,12.6461,12.5069,12.3591,12.2101,12.0619,11.9119,11.7569,11.5973,11.4464,11.3103,11.1792,11.052,10.9275,10.8055,10.6862,10.5706,10.4583,10.3489),

  "6" = c(6.55,6.55,6.55,6.55,6.55,6.55,6.55,6.55,6.55,6.55,6.55,6.55,6.55,6.55,6.551,6.5536,6.5577,6.5632,6.5699,6.596,6.645,6.6945,6.7237,6.71,6.6427,6.5362,6.4072,6.2727,6.1501,6.0409,5.9357,5.8331,5.7316,5.6301,5.5327,5.438,5.3417,5.2423,5.1399,5.0346,4.9276,4.8191,4.7105,4.6083,4.5154,4.4268,4.3413,4.2578,4.1762,4.0965,4.0189,3.9437,3.8705),

  "7" = c(2.49,2.49,2.49,2.49,2.49,2.49,2.49,2.49,2.49,2.49,2.49,2.49,2.49,2.49,2.4901,2.4901,2.4901,2.49,2.4899,2.4768,2.4584,2.4561,2.4691,2.48,2.5053,2.566,2.6468,2.7197,2.7401,2.713,2.6803,2.6442,2.6067,2.5699,2.5273,2.4722,2.4087,2.3429,2.2799,2.2199,2.1592,2.0986,2.0385,1.982,1.9302,1.881,1.8335,1.7871,1.7413,1.6962,1.6522,1.6092,1.5674)

) %>%
  melt(id=c("iso3", "area_id", "area_level", "area_name", "period", "variable", "source"), variable.name="agegr", value.name="val") %>%
  mutate(agegr = rep(c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49"), each=53))



asfr_prop_uga <- model_asfr %>%
  filter(area_level ==0, iso3 == "UGA") %>%
  group_by(period) %>%
  mutate(val = 100*(val/sum(val)),
         variable = "asfr_prop") %>%
  bind_rows(spec_asfr_uga)

asfr_prop %>%
  ggplot(aes(x=period, y=val, group=agegr, color=agegr, linetype=source)) +
  geom_line(data=asfr_prop %>% filter(source =="Model", period>1999)) +
  geom_line(data=asfr_prop %>% filter(source =="Spectrum18", period>1999)) +
  labs(y="Age distribution of fertility (%)", title="Age distribution of fertility | Malawi")


## Combine

asfr_prop <- age_specific_outputs %>% 
  filter(variable=="asfr_prop") %>%
  bind_rows(asfr_prop_uga)

spec_tfr <- all_age_outputs %>%
  filter(variable=="tfr", source =="Spectrum18") %>%
  bind_rows(spec_tfr_uga)

all_age_outputs <- model_births %>%
  bind_rows(model_tfr, worldpop_U1, wpp_tfr, gbd_U1, gbd_tfr) %>%
  filter(!is.na(val))

age_specific_outputs <- model_asfr %>%
  bind_rows(model_births_age, wpp_asfr, gbd_asfr) %>%
  filter(!is.na(val))

saveRDS(age_specific_outputs, file="age_specific_outputs.rds")
saveRDS(all_age_outputs, file="allage_outputs.rds")

age_specific_outputs <- readRDS(file="age_specific_outputs.rds")
all_age_outputs <- readRDS(file="allage_outputs.rds")

####### MAPS

age_specific_outputs %>%
  filter(iso3 =="UGA", area_level %in% c(0, 1, 2), variable=="asfr", period==2015, source =="Model") %>%
  left_join(boundaries) %>%
  ggplot() +
    geom_sf(aes(geometry = geometry, fill=val)) +
    facet_grid(area_level~agegr) +
    labs(title="ASFR by admin level in 2015 | Uganda")


ggplot(data=data, aes(x=period, y=val, group=agegr, color=agegr, linetype=source)) +
  geom_line(data=age_specific_outputs %>% filter(variable == "asfr_prop", source =="Model", period>1999, area_level==0)) +
  geom_line(data=age_specific_outputs %>% filter(variable == "asfr_prop", source =="Spectrum18", period>1999, area_level==0)) +
  labs(y="Age distribution of fertility (%)", title="Age distribution of fertility") +
  facet_wrap(~area_name)
  

all_age_outputs %>%
  filter(iso3 =="ZWE", area_level %in% c(0, 1, 2), variable=="tfr", period==2015, source =="Model") %>%
  left_join(boundaries) %>%
  ggplot() +
    geom_sf(aes(geometry = geometry, fill=val)) +
    facet_wrap(~area_level) +
    labs(title="TFR by admin level in 2015 | Zimbabwe")

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

all_age_outputs %>%
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

