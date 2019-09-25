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
  mutate(asfr_ratio0 = mean*(population/id0_agepop),
         asfr_ratio1 = mean*(population/id1_agepop),
         asfr_ratio2 = mean*(population/id2_agepop),
         asfr_ratio3 = mean*(population/id3_agepop),
         asfr_ratio4 = mean*(population/id4_agepop)
  )

asfr_admin0 <- asfr_aggr %>%
  group_by(period, agegr, id0) %>%
  summarise(asfr = sum(asfr_ratio0)) %>%
  rename(area_id = id0)

asfr_admin1 <- asfr_aggr %>%
  group_by(period, agegr, id1) %>%
  summarise(asfr = sum(asfr_ratio1)) %>%
  rename(area_id = id1)

asfr_admin2 <- asfr_aggr %>%
  group_by(period, agegr, id2) %>%
  summarise(asfr = sum(asfr_ratio2)) %>%
  rename(area_id = id2)

asfr_admin3 <- asfr_aggr %>%
  group_by(period, agegr, id3) %>%
  summarise(asfr = sum(asfr_ratio3)) %>%
  rename(area_id = id3)

asfr_admin4 <- asfr_aggr %>%
  group_by(period, agegr, id4) %>%
  summarise(asfr = sum(asfr_ratio4)) %>%
  rename(area_id = id4)

tfr_admin0 <- asfr_aggr %>%
  group_by(period, id0) %>%
  summarise(tfr = 5*sum(asfr_ratio0)) %>%
  rename(area_id = id0)

tfr_admin1 <- asfr_aggr %>%
  group_by(period, id1) %>%
  summarise(tfr = 5*sum(asfr_ratio1)) %>%
  rename(area_id = id1)

tfr_admin2 <- asfr_aggr %>%
  group_by(period, id2) %>%
  summarise(tfr = 5*sum(asfr_ratio2)) %>%
  rename(area_id = id2)

tfr_admin3 <- asfr_aggr %>%
  group_by(period, id3) %>%
  summarise(tfr = 5*sum(asfr_ratio3)) %>%
  rename(area_id = id3)

tfr_admin4 <- asfr_aggr %>%
  group_by(period, id4) %>%
  summarise(tfr = 5*sum(asfr_ratio4)) %>%
  rename(area_id = id4)

age_disag_outputs <- areas_long %>%
  filter(iso3=="MWI") %>%
  left_join(asfr_admin0 %>% bind_rows(asfr_admin1, asfr_admin2, asfr_admin3, asfr_admin4), by="area_id") %>%
  left_join(worldpop_rep %>% select(c(id, age, population, year)), by=c("agegr" = "age", "period" = "year", "area_id" = "id")) %>%
  mutate(births = asfr*population)





######### National comparison of ASFR and TFR vs WPP 2019

wpp_asfr <- read_excel("WPP2019_FERT_F07_AGE_SPECIFIC_FERTILITY.xlsx")

wpp_asfr <- wpp_asfr %>%
  separate(period, into=c("period", NA), sep="-") %>%
  type.convert %>%
  mutate(period = period + 2.5,
         iso3 = countrycode(country, "country.name", "iso3c"),
         area_name = countrycode(iso3, "iso3c", "country.name"),
         area_id = iso3
  ) %>%
  filter(iso3 %in% unique(clusters$iso3)) %>%
  select(-c(countrycode, country)) %>%
  melt(value.name="asfr", id=c("iso3", "area_name", "area_id", "period"), variable.name="agegr") %>%
  mutate(asfr = as.numeric(asfr)/1000)

asfr_area %>%
  filter(area_level==0) %>%
  mutate(source = "mod") %>%
  select(-c(area_level, parent_area_id, population, births)) %>%
  rbind(wpp_asfr %>% filter(period == 2012.5) %>% mutate(source="WPP2019 2012.5", period=2015)) %>%
  rbind(wpp_asfr %>% filter(period == 2017.5) %>% mutate(source="WPP2019 2017.5", period=2015)) %>%
  filter(iso3 == "MWI") %>%
  ggplot(aes(group=source)) +
    geom_line(aes(x=agegr, y=asfr, color=source))

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

nonage_ouputs <- areas_long %>%
  filter(iso3 == "MWI") %>%
  left_join(tfr_admin0 %>% bind_rows(tfr_admin1, tfr_admin2, tfr_admin3, tfr_admin4), by="area_id") %>%
  left_join(asfr_area %>%
              group_by(iso3, area_id, area_name, area_level, parent_area_id, period) %>%
              summarise(est_births = sum(births)) %>%
              ungroup %>%
              select(area_id, period, est_births), 
            by=c("area_id", "period")) %>%
  left_join(worldpop_rep %>%
              filter(startage==0) %>%
              select(id, year, population) %>%
              rename(worldpop_U1 = population),
            by=c("area_id" = "id", "period" = "year")) %>%
  left_join(wpp_tfr, by=c("area_name", "period", "iso3"))

