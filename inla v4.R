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
library(extrafont)

devtools::install_github("tidyverse/tidyr")
library(tidyr)

setwd("~/Documents/GitHub/subnat_fertility")
load("asfr.rds")
load("asfr_dhs.rds")

set_rdhs_config(email="o.stevens@imperial.ac.uk", project="Subnational fertility", config_path = "~/.rdhs.json")

source("fertility_funs.R")

## Set flags. Multicountry == TRUE for running several countries in a single model. Subnational currently only set up for running 1 country at a time, setting the iso codes as strings manually {see lines 62-69 and 131-198 for ZWE admin 1 example}
subnational <- TRUE
multicountry <- FALSE

# ## If Namibia is in this list, ensure it's first. Crap but it'll do for now.
# iso2 <- c("NM", "RW", "UG", "MZ", "MW", "KE", "ZM", "ZW", "TZ", "LS")
# 
# iso2_list <- as.list(iso2)
# 
# iso3_list <- lapply(iso2_list, function(x) {x %>% countrycode("iso2c", "iso3c")})
# 
# if (iso2[1]=="NM") {
#   iso3_list[1] <- "NAM"
# }


formulae <- list()

formulae[[1]] <- births ~ tips_dummy + f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2")) + f(id.tips, model="rw1")

formulae_number <- length(formulae)

dhs_iso3 <- dhs_countries(returnFields=c("CountryName", "DHS_CountryCode")) %>%
  mutate(iso3 = countrycode(CountryName, "country.name", "iso3c"),
         iso3 = ifelse(CountryName == "Eswatini", "SWZ", iso3))

## Get cluster coordinates
clusters <- readRDS("~/Documents/GitHub/naomi-data-edit/oli_cluster.rds") %>%
  mutate(iso3 = survey_id) %>%
  separate(col="iso3", into="iso3", sep=3) %>%
  left_join(dhs_iso3 %>% select(-CountryName), by="iso3") %>%
  separate(survey_id, into=c(NA, "surv"), sep=3, remove=FALSE) %>%
  mutate(DHS_survey_id = paste0(DHS_CountryCode, surv)) %>%
  select(-surv) %>%
  filter(survey_id != "CIV2005AIS")

clusters <- readRDS("oli_clusters.rds")

iso3 <- as.list(clusters %>% .$iso3 %>% unique)

## Get surveys for which we have clusters. Split into country list.
surveys <- dhs_surveys(surveyIds = unique(clusters$DHS_survey_id)) %>%
  left_join(clusters %>% select(c(DHS_survey_id, survey_id)) %>% distinct, by=c("SurveyId" = "DHS_survey_id")) %>%
  filter(CountryName == "Ethiopia") %>%
  group_split(CountryName)

## NATIONAL
surveys <- dhs_surveys(surveyIds = unique(clusters$DHS_survey_id)) %>%
  filter(!SurveyId %in% c("MZ2015AIS", "TZ2007AIS", "UG2014MIS"), CountryName != "Malawi")
  group_split(CountryName)
  
surveys <- dhs_surveys(countryIds = c("MW", "KE", "ZW", "RW", "LS", "ZM",  "UG"), surveyYearStart = 1995) %>%
  filter(SurveyType == "DHS") %>%
  group_split(CountryName)


# surveys <- dhs_surveys(countryIds = iso3_to_dhs(unique(clusters$iso3)), surveyYearStart=1995) %>%
#   separate(SurveyId, into=c(NA, "surv"), sep =2, remove=FALSE) %>%
#   mutate(survey_id = paste0(iso3, surv)) %>%
#   select(-surv) %>%
#   filter(SurveyType != "AIS") %>%
#   filter(survey_id %in% unique(clusters$survey_id))
#   group_split(CountryName)

# 
# areas <- readRDS("~/Documents/GitHub/naomi-data/data/areas/areas_wide.rds") %>%
#   filter(iso3 == "ZWE") %>%
#   left_join(clusters, by=c("id2" = "geoloc_area_id", "iso3")) %>%
#   dplyr::select(iso3, id1, name1, survey_id, cluster_id) %>%
#   rename(area_id = id1, area_name = name1)

areas <- readRDS("~/Documents/GitHub/naomi-data/data/areas/areas_long.rds") %>%
  inner_join(clusters, by=c("area_id" = "geoloc_area_id", "iso3")) %>%
  filter(iso3 == "ETH")

area_list <- areas %>%
  group_by(survey_id) %>%
  group_split(keep=TRUE)

boundaries <- readRDS("~/Documents/GitHub/naomi-data/data/areas/boundaries.rds")

ird <- lapply(surveys, function(surveys) {
  dhs_datasets(fileType = "IR", fileFormat = "flat", surveyIds = surveys$SurveyId)
})

### FOR WHEN API IS DOWN ####

dhs_datasets_cached <- readRDS("~/Downloads/dhs_datasets_cached.rds")

dhs_datasets_cached <- dhs_datasets_cached %>%
  filter(SurveyId %in% surveys$SurveyId, FileType == "Individual Recode", FileFormat == "Flat ASCII data (.dat)") %>%
  separate(FileName, into=c("FileName", NA), sep=-4, remove=TRUE) %>%
  mutate(path = paste0("/Users/os210/Library/Caches/rdhs/datasets/", FileName, ".rds"))

ird <- dhs_datasets_cached %>%
  group_split(DHS_CountryCode)

###########

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

ir_by_area2 <- function(ir, area_list) {
  
  print("run done")
  
  ir_int <- ir %>%
    left_join(area_list, by=c("v001" = "cluster_id")) %>%
    filter(!is.na(area_id)) %>%
    group_split(area_id)
  
  return(ir_int)
  
}


ir_area <- Map(ir_by_area2, ir, area_list) %>%
  unlist(recursive = FALSE)

save(ir_area, file="ir_area.RData")
load("ir_area.RData")

test <- lapply(ir_area, function(x) {
  x$foo <- x %>%
    .$survtype %>%
    unique
}) %>%
  unlist

names(ir_area) <- sapply(ir_area,  function(x) {
  paste(x[["area_id"]][1], x[["survey_id"]][1])
})

test <- data.frame("survtype" = test)

tips_surv <- list("DHS" = c(0:10), "MIS" = c(0:5), "AIS" = c(0:5))[test$survtype]

tips_surv <- list("DHS" = c(0,10), "MIS" = c(0,5), "AIS" = c(0,5))[test$survtype]

tips_surv <- list("DHS" = c(0:10), "MIS" = c(0:5), "AIS" = c(0:5))[surveys %>%
                                                                     bind_rows %>%
                                                                     .$SurveyType]

tips_surv <- list("DHS" = c(0,10), "MIS" = c(0,5), "AIS" = c(0,5))[surveys %>%
                                                                     bind_rows %>%
                                                                     .$SurveyType]

asfr_eth <- Map(calc_asfr1, ir_area,
            y=1:length(ir_area),
            by = list(~country + surveyid + survtype + survyear + area_name + area_id),
            tips = tips_surv,
            agegr= list(3:10*5),
            #period = list(seq(1995, 2017, by=0.5)),
            period = list(1985:2009),
            counts = TRUE)

asfr_eth <- asfr_eth %>%
  bind_rows() %>%
  type.convert() %>%
  mutate(period = period +8)

asfr <- asfr_eth

asfr <- Map(calc_asfr1, ir,
                     y=1:length(ir),
                     by = list(~country + surveyid + survtype + survyear),
                     tips = tips_surv,
                     agegr= list(3:10*5),
                     #period = list(seq(1995, 2017, by=0.5)),
                     period = list((1995*1):(2017*1)*(1/1)),
                     counts = TRUE)

tfr <- Map(calc_tfr, ir, by=list(~country + surveyid + survyear), tips = tips_surv, period=list(1995:2017))

x %>%
  group_by(country, surveyid, period) %>%
  summarise(tfr = mean(tfr)) %>%
  mutate(id = 1:nrow(.),
         tfr_stan = tfr/tfr[id==1]
  )

tfr %>%
  bind_rows %>%
  type.convert %>%
  filter(period<=survyear, tips == "0-9") %>%
  ggplot(aes(x=period, y=tfr)) +
  geom_point()+
  #scale_x_continuous(breaks=0:9, labels=as.character(0:9))+
  facet_wrap(~surveyid)

tfr_ZM2013DHS <- tfr %>%
  bind_rows %>%
  type.convert %>%
  filter(period<=survyear) %>%
  filter(surveyid == "ZM2013DHS", period>2003)

saveRDS(tfr_ZM2013DHS, file="tfr_ZM2013DHS.rds")

tfr_ZM2013DHS %>%
  ggplot(aes(x=period, y=tfr)) +
  scale_x_continuous(breaks=seq(2004,2013, 1), labels=as.character(seq(2004,2013, 1)))+
    geom_point(size=2) +
    geom_line(size=0.75) +
    annotate(geom = "text", x = 2004:2013, y = 4.75, label = as.character(9:0)) + 
    geom_vline(aes(xintercept = 2007.5), color="red", linetype=2)+
    coord_cartesian(ylim = c(5, 7), expand = FALSE, clip = "off") +
    labs(y="TFR") +
    theme(
      axis.title.x = element_blank(),
      plot.margin = unit(c(1, 1, 2, 1), "lines")
    )

tfr[[3]] %>%
  type.convert %>%
  filter(period<=survyear) %>%
  group_by(country, surveyid, tips) %>%
  summarise(tfr = mean(tfr)) %>%
  ungroup %>%
  mutate(tfr_stan = tfr/tfr[tips==0])

tfr_plot %>%
  bind_rows %>%
  mutate(id = as.factor(id)) %>%
  filter(is.na(id))
  mutate(id = factor(id, levels=0:9, labels=as.character(0:9))) %>%
  ggplot(aes(x=id, y=tfr_stan)) +
    geom_line() +
    facet_wrap(~surveyid)

tfr_eth <- tfr_eth %>%
  bind_rows %>%
  type.convert() %>%
  filter(period<survyear)

tfr_eth %>% 
  mutate(lower = tfr+(qnorm(0.95)*se_tfr),
         upper = tfr-(qnorm(0.95)*se_tfr)
  ) %>%
  ggplot(aes(x=period, y=tfr, group=surveyid, color=surveyid)) +
    geom_point() +
    geom_line()
    geom_ribbon(aes(ymin=lower, ymax=upper, fill=surveyid), alpha=0.2)

asfr_zwe <- asfr_zwe %>%
  bind_rows %>%
  type.convert

saveRDS(asfr_zwe_subnat, file="asfr_zwe_subnat.rds")

asfr_zwe_subnat <- asfr_zwe %>%
  bind_rows() %>%
  type.convert()



calc_asfr1(ir[[28]], y=836, by = ~country + surveyid + survtype + survyear,
           tips = tips_surv[[28]],
           agegr= 3:10*5,
           #period = list(seq(1995, 2017, by=0.5)),
           period = 1995:2017,
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

# summary(mod_list[[3]])
# 
# exp(mod_list[[3]]$summary.fixed)
# 
# exp(mod_list[[3]]$summary.random$id.tips) %>%
#   mutate(id = 1:11) %>%
#   ggplot(aes(x=id)) +
#   geom_point(aes(y=`0.5quant`)) +
#   geom_errorbar(aes(ymin = `0.025quant`, ymax=`0.975quant`), width=0.2) +
#   labs(title="TIPS RW1, FE Dummy | Zimbabwe admin1")

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
             id.period3 = id.period,
             id.agegr = group_indices(., agegr),
             id.agegr2 = id.agegr,
             id.agegr3 = id.agegr,
             id.agegr.period = group_indices(., period, agegr),
             id.tips = (group_indices(., tips)),
             id.tips2 = id.tips,
             tips_dummy = ifelse(tips>5, 1, 0),
             id.country = group_indices(., country),
             #survey_dummy = (group_indices(., survtype)),
             # id.district = group_indices(., area_name),
             # id.district2 = id.district,
             # id.district3 = id.district,
             # id.survey = group_indices(., surveyid),
             # id.agegr.period.district = group_indices(., agegr, period, district),
             # id.region = group_indices(., region_name),
             id = 1:nrow(.),
      )
    
    formulae <- lapply(formulae, function(x) {
      update(x,    ~ . + country + f(id.agegr3, model = "rw1", group = id.country, control.group = list(model = "iid") + f(id.period2, model = "rw2", group = id.country, control.group = list(model = "iid"))))
    })
    
    mod_list <- Map(run_mod, formulae, list(asfr_pred), model_family = "poisson")
    
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
    
  load("asfr_subnational.RData")
  
  asfr1_country <- asfr %>%
    bind_rows %>%
    type.convert %>%
    filter(period<=survyear) %>%
    # separate(surveyid, into=c(NA, "type"), sep=6, remove=FALSE) %>%
    # filter(type == "MIS") %>%
    # select(-type) %>%
    group_split(country, keep=TRUE) %>%
    lapply(droplevels)
  
  
  
  asfr_pred_country <- lapply(asfr1_country, function(asfr1_country) {
    crossing(country = asfr1_country$country, area_name = asfr1_country$area_name, period = asfr1_country$period, agegr = asfr1_country$agegr,  pys=1) %>%
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
             id.tips = factor(id.tips),
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
  asfr_pred_country_eth <- asfr_pred_country
  asfr_pred_country <- readRDS("asfr_pred_country_subnat.rds")
  
  asfr_pred_country[6]  <- asfr_pred_country_eth
  
  saveRDS(asfr_pred_country, file="asfr_pred_country_subnat.rds")
  
  # formulae <- repeat_formulae(formulae, length(iso2))
  #  formulae <- rep(formulae, length(iso2))
  #formulae <- formulae[-10]
  
  max_level <- areas %>% group_by(iso3) %>% summarise(max_level = max(area_level)) %>%
    mutate(max_level = ifelse(iso3 == "MWI", 4, max_level))
  
  saveRDS(max_level, file="max_level.rds")
  
  iso3_list <- as.list(max_area %>% .$iso3)
  
  areas_long <- readRDS("~/Documents/GitHub/naomi-data/data/areas/areas_long.rds")
  
  Map(function(iso3_list, max_level) {
    
    sh <- areas_long %>%
      filter(iso3 == "ETH", area_level == 2) %>%
      mutate(area_idx = row_number())
    
    #' Neighbor list
    nb <- sh %>%
      left_join(boundaries) %>%
      st_as_sf %>%
      as("Spatial") %>%
      spdep::poly2nb() %>%
      `names<-`(sh$area_idx)
    
    nb2INLA(paste0("ETH.adj"), nb)
  }, iso3_list, max_level, list(areas_long), list(boundaries))
  
  sh <- readRDS("~/Documents/GitHub/naomi-data/data/areas/areas_long.rds") %>%
    filter(iso3 == "LSO", area_level == 2) %>%
    mutate(area_idx = row_number())
  
  #' Neighbor list
  nb <- sh %>%
    left_join(boundaries) %>%
    st_as_sf %>%
    as("Spatial") %>%
    spdep::poly2nb() %>%
    `names<-`(sh$area_idx)
  
  nb2INLA(paste0("LSO", ".adj"), nb)
  
  formula <- births ~ tips_dummy + f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2")) + f(id.tips, model="rw1")
  
  mod <- run_mod(formula, asfr_pred_country[[1]], model_family = "poisson")
  
tips <- Map(function(mod_list, asfr_pred_country) {
  
  mod_list$summary.fixed %>%
    exp %>%
    filter(mean>0.2) %>%
    mutate(tips = 1:nrow(.)) %>%
    bind_rows(data.frame("mean" = 1, tips = 0)) %>%
    mutate(country = unique(asfr_pred_country$country))
}, mod_list, asfr_pred_country) %>%
  bind_rows

saveRDS(tips, file="tips.rds")

tips %>%
  filter(country %in% c("Malawi", "Zambia", "Zimbabwe")) %>%
  ggplot(aes(x=tips, y=mean, ymin=`0.025quant`, ymax=`0.975quant`)) +
    geom_point() +
    geom_errorbar(width=0.5) +
    scale_x_continuous(breaks=seq(0,10, 1), labels = as.character(0:10)) +
    labs(x="Time Preceding Survey (TIPS)", y="Mean") +
    facet_wrap(~country)
    
  mod_list <- Map(run_mod, list(formula), asfr_pred_country, model_family="poisson")
  
tips_rw1 <- mod_list[[1]]$summary.random$id.tips %>%
    type.convert %>%
    mutate(tips = 0:9,
           mean = ifelse(tips>5, exp(mean+mod_list[[1]]$summary.fixed$mean[2]), exp(mean)),
           `0.025quant` = ifelse(tips>5, exp(`0.025quant`+mod_list[[1]]$summary.fixed$mean[2]), exp(`0.025quant`)),
           `0.975quant` = ifelse(tips>5, exp(`0.975quant`+mod_list[[1]]$summary.fixed$mean[2]), exp(`0.975quant`))
    )

saveRDS(tips_rw1, file="tips_rw1.rds")
  
data.frame(4) %>%
    ggplot(aes(x=tips, y=mean, ymin=`0.025quant`, ymax=`0.975quant`)) +
    geom_line(data = tips_rw1  %>% filter(tips<6)) +
    geom_line(data = tips_rw1  %>% filter(tips>5)) +
    geom_ribbon(data = tips_rw1 %>% filter(tips<6), alpha=0.2) +
    geom_ribbon(data = tips_rw1  %>% filter(tips>5), alpha=0.2) +
    scale_x_continuous(breaks=seq(0,10, 1), labels = as.character(0:10)) +
    labs(x="Time Preceding Survey (TIPS)", y="Mean")
  
  
  tfr_mwi %>%
    bind_rows %>%
    type.convert %>%
    ggplot(aes(x=period, y=tfr, group = surveyid, color=surveyid)) +
      geom_point() +
      geom_line()
  
  mod_list_all_survey <- mod_list
  mod_list_all_survey[[2]]$summary.fixed %>%
    exp() %>%
    filter(mean>0.2) %>%
    ggplot(aes(x=1:9, y=mean)) +
      geom_line()
  
  mod_list <- Map(run_mod, formulae, asfr_pred = rep(asfr_pred_country, formulae_number), model_family="poisson")
  
  poisson_summary <- lapply(mod_list, summary)
  

  plot_df_fac <- Map(function(mod_list, iso3_list) {
    exp(mod_list$summary.fixed) %>%
      filter(mean>0.3) %>%
      mutate(tips = 1:9,
             country = countrycode(iso3_list, "iso3c", "country.name")) %>%
      bind_rows(data.frame("mean" = 1, "tips" = 0, "country" = countrycode(iso3_list, "iso3c", "country.name")))
  }, mod_list, iso3_list) %>%
      bind_rows
  
  plot_df <- Map(function(mod_list, iso3_list) {
    mod_list$summary.random$id.tips %>%
      mutate(tips = (0:(nrow(.)-1))) %>%
      mutate(mean = ifelse(tips>5, exp(mean+mod_list$summary.fixed[2,1]), exp(mean)),
             `0.025quant` = ifelse(tips>5, exp(`0.025quant`+mod_list$summary.fixed[2,1]), exp(`0.025quant`)),
             `0.975quant` = ifelse(tips>5, exp(`0.975quant`+mod_list$summary.fixed[2,1]), exp(`0.975quant`)),
             country = countrycode(iso3_list, "iso3c", "country.name"))
    }, mod_list, iso3_list) %>%
      bind_rows
  
  plot_df <- Map(function(mod_list, iso3_list) {
    exp(mod_list$summary.random$id.tips) %>%
      mutate(tips = (0:(nrow(.)-1))) %>%
      mutate(country = iso3_list)
  }, mod_list, as.list(unique(asfr1_country %>% bind_rows %>% .$country))) %>%
    bind_rows
    
      
  plot_df %>%
      ggplot(aes(x=tips, y=mean, ymin=`0.025quant`, ymax = `0.975quant`)) +
        geom_line(data=plot_df %>% filter(tips<6)) +
        geom_ribbon(data=plot_df %>% filter(tips<6), alpha=0.5, fill="blue") +
        geom_line(data=plot_df %>% filter(tips>5)) +
        geom_ribbon(data=plot_df %>% filter(tips>5), alpha=0.5, fill="blue")+
        geom_point(data=plot_df_fac, size=0.3) +
        geom_errorbar(data = plot_df_fac, aes(width=0.2)) +
        geom_hline(aes(yintercept=1), linetype=3)+
        scale_x_continuous(breaks=0:9) +
        labs(title="TIPS in DHS surveys: Smooth vs factor")+
        facet_wrap(~country)
  
  if(formulae_number>1) {
  
    mod_split <- Map(function(x, y) {
      mod_list[x:(x+length(iso2)-1)]
    }, x=seq(1, formulae_number*length(iso2), by=length(iso2)), y=1:formulae_number)
    
  }
  
  pred_list <- Map(get_pred, mod_list, rep(asfr_pred_country, formulae_number), rep(asfr1_country, formulae_number))
  
  
    pred_split <- Map(function(x, y) {
      pred_list[x:(x+length(iso2)-1)] %>%
        bind_rows %>%
        mutate(source=paste("formula", y))
    }, x=seq(1, formulae_number*length(iso2), by=length(iso2)), y=1:formulae_number)

  
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
  
  plots[[1]]
  
  ggsave("Separate, formula 1.png", plot=plots[[1]], device="png")
  ggsave("Separate, formula 2.png", plot=plots[[2]], device="png")
  ggsave("Separate, formula 3.png", plot=plots[[3]], device="png")
  ggsave("Separate, formula 4.png", plot=plots[[4]], device="png")
  
  
  


