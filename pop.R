library(tidyverse)

setwd("~/Documents/GitHub/subnat_fertility")

### WORLD POP. Filter areas on cluster area. Convert MWI from level 5 to level 4.

worldpop <- read_csv("WorldPop_agesex.csv", col_types = cols(X1 = col_skip(), X = col_skip(), sex = col_character())) %>%
  filter(sex=="f")

areas_long <- readRDS("~/Documents/GitHub/naomi-data/data/areas/areas_long.rds")
areas_wide <- readRDS("~/Documents/GitHub/naomi-data/data/areas/areas_wide.rds")

worldpop_rep <- worldpop %>%
  left_join(areas_wide, by=c("iso3", "id" = "area_id")) %>%
  filter(year != 2020) %>%
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
  )

asfr_admin0 <- asfr_aggr %>%
  group_by(period, agegr, id0) %>%
  summarise(asfr = sum(asfr_ratio0),
            asfr_l = sum(asfr_ratio0_l),
            asfr_u = sum(asfr_ratio0_u)
  ) %>%
  rename(area_id = id0)

asfr_admin1 <- asfr_aggr %>%
  group_by(period, agegr, id1) %>%
  summarise(asfr = sum(asfr_ratio1),
            asfr_l = sum(asfr_ratio1_l),
            asfr_u = sum(asfr_ratio1_u)
  ) %>%
  rename(area_id = id1)

asfr_admin2 <- asfr_aggr %>%
  group_by(period, agegr, id2) %>%
  summarise(asfr = sum(asfr_ratio2),
            asfr_l = sum(asfr_ratio2_l),
            asfr_u = sum(asfr_ratio2_u)
  ) %>%
  rename(area_id = id2)

asfr_admin3 <- asfr_aggr %>%
  group_by(period, agegr, id3) %>%
  summarise(asfr = sum(asfr_ratio3),
            asfr_l = sum(asfr_ratio3_l),
            asfr_u = sum(asfr_ratio3_u)
  ) %>%
  rename(area_id = id3)

asfr_admin4 <- asfr_aggr %>%
  group_by(period, agegr, id4) %>%
  summarise(asfr = sum(asfr_ratio4),
            asfr_l = sum(asfr_ratio4_l),
            asfr_u = sum(asfr_ratio4_u)
  ) %>%
  rename(area_id = id4)

tfr_admin0 <- asfr_aggr %>%
  group_by(period, id0) %>%
  summarise(tfr = 5*sum(asfr_ratio0),
            tfr_l = 5*sum(asfr_ratio0_l),
            tfr_u = 5*sum(asfr_ratio0_u)
  ) %>%
  rename(area_id = id0)

tfr_admin1 <- asfr_aggr %>%
  group_by(period, id1) %>%
  summarise(tfr = 5*sum(asfr_ratio1),
            tfr_l = 5*sum(asfr_ratio1_l),
            tfr_u = 5*sum(asfr_ratio1_u)
  ) %>%
  rename(area_id = id1)

tfr_admin2 <- asfr_aggr %>%
  group_by(period, id2) %>%
  summarise(tfr = 5*sum(asfr_ratio2),
            tfr_l = 5*sum(asfr_ratio2_l),
            tfr_u = 5*sum(asfr_ratio2_u)
  ) %>%
  rename(area_id = id2)

tfr_admin3 <- asfr_aggr %>%
  group_by(period, id3) %>%
  summarise(tfr = 5*sum(asfr_ratio3),
            tfr_l = 5*sum(asfr_ratio3_l),
            tfr_u = 5*sum(asfr_ratio3_u)
  ) %>%
  rename(area_id = id3)

tfr_admin4 <- asfr_aggr %>%
  group_by(period, id4) %>%
  summarise(tfr = 5*sum(asfr_ratio4),
            tfr_l = 5*sum(asfr_ratio4_l),
            tfr_u = 5*sum(asfr_ratio4_u)
  ) %>%
  rename(area_id = id4)

age_disag_outputs <- areas_long %>%
  filter(iso3=="MWI") %>%
  left_join(asfr_admin0 %>% bind_rows(asfr_admin1, asfr_admin2, asfr_admin3, asfr_admin4), by="area_id") %>%
  left_join(worldpop_rep %>% select(c(id, age, population, year)), by=c("agegr" = "age", "period" = "year", "area_id" = "id")) %>%
  mutate(births = asfr*population,
         births_l = asfr_l*population,
         births_u = asfr_u*population)


nonage_outputs <- areas_long %>%
  filter(iso3 == "MWI") %>%
  left_join(tfr_admin0 %>% bind_rows(tfr_admin1, tfr_admin2, tfr_admin3, tfr_admin4), by="area_id") %>%
  left_join(age_disag_outputs %>%
              group_by(iso3, area_id, area_name, area_level, parent_area_id, period) %>%
              summarise(births = sum(births),
                        births_l = sum(births_l),
                        births_u = sum(births_u)) %>%
              ungroup %>%
              select(area_id, period, births, births_l, births_u), 
            by=c("area_id", "period")) %>%
  left_join(worldpop_rep %>%
              filter(startage==0) %>%
              select(id, year, population) %>%
              rename(worldpop_U1 = population),
            by=c("area_id" = "id", "period" = "year"))


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
  melt(id=c("iso3", "area_name", "area_id", "period"), variable.name="agegr", value.name = "wpp19_asfr") %>%
  filter(iso3 == "MWI") %>%
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
  mutate(wpp19_asfr = as.numeric(wpp19_asfr)/1000)

age_disag_outputs <- age_disag_outputs %>%
  left_join(wpp_asfr)

wpp_tfr <- read_excel("WPP2019_FERT_F04_TOTAL_FERTILITY(1).xlsx")

wpp_tfr <- wpp_tfr %>%
  melt(id="area_name", variable.name="period", value.name = "wpp19_tfr") %>%
  separate(period, into=c("period", NA), sep="-") %>%
  type.convert() %>%
  mutate(iso3 = countrycode(area_name, "country.name", "iso3c")) %>%
  filter(iso3 == "MWI") %>%
  group_split(period) %>%
  lapply(function(x) {
    n <- nrow(x)
    year_start <- unique(x$period)
    x <- x %>%
      bind_rows(., ., ., ., .) %>%
      mutate(period = rep(year_start:(year_start+4), each=n))
    return(x)
  })  %>%
  bind_rows

nonage_outputs <- nonage_outputs %>%
  left_join(wpp_tfr, by=c("area_name", "period", "iso3"))

######### GBD 2017 national ASFR and TFR

gbd_asfr <- read_csv("IHME_GBD_2017_FERT_ESTIMATES_1950_2017_Y2018M11D08.CSV", col_types = cols(age_group_id = col_skip(), location_id = col_skip(), measure_id = col_skip(), measure_name = col_skip(), metric_name = col_skip(), sex_id = col_skip(), sex_name = col_skip()))

gbd_asfr <- gbd_asfr %>%
  mutate(iso3 = countrycode(location_name, "country.name", "iso3c")) %>%
  filter(iso3 == "MWI") %>%
  separate(age_group_name, into=c("agegr", "up"), sep=" to ") %>%
  mutate(agegr = paste0(agegr, "-", up)) %>%
  select(-up) %>%
  rename(area_name = location_name, period = year_id, gbd_asfr = val, gbd_asfr_u = upper, gbd_asfr_l = lower)

age_disag_outputs <- age_disag_outputs %>%
  left_join(gbd_asfr)

gbd_tfr <- gbd_asfr %>%
  group_by(iso3, area_name, period) %>%
  summarise(gbd_tfr = 5*sum(gbd_asfr),
            gbd_tfr_u = 5*sum(gbd_asfr_u),
            gbd_tfr_l = 5*sum(gbd_asfr_l),
  )

nonage_outputs <- nonage_outputs %>%
  left_join(gbd_tfr)

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
  rename(area_name = location_name, period = year_id, gbd_U1 = val, gbd_U1_u = upper, gbd_U1_l = lower) %>%
  mutate(iso3 = countrycode(area_name, "country.name", "iso3c"))

nonage_outputs <- nonage_outputs %>%
  left_join(gbd_U1)
