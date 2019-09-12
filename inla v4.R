library(countrycode)
library(tidyverse)
library(magrittr)
library(rdhs)
library(demogsurv)
library(INLA)
library(reshape2)
library(survival)
library(sf)

setwd("~/Documents/GitHub/subnat_fertility")


source("fertility_funs.R")

multicountry <- FALSE

## If Namibia is in this list, ensure it's first. Crap but it'll do for now.
iso2 <- c("NM", "RW", "UG", "MZ", "MW", "KE", "ZM", "ZW", "TZ")
iso2_list <- as.list(sapply(iso2, rep, formulae_number))

# names(iso2_list) <- sapply(iso2_list, countrycode, "iso2c", "country.name")

iso3_list <- lapply(iso2_list, function(x) {x %>% countrycode("iso2c", "iso3c")})

if (iso2[1]=="NM") {
  iso3_list[1:formulae_number] <- "NAM"
}

names(iso3_list) <- sapply(iso3_list, countrycode, "iso3c", "country.name")

formulae <- list()

formulae[[1]] <- births ~ f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2")) + f(id.tips, model="rw1")

formulae[[2]] <- births ~ f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2")) + f(id.tips, dummy_step, model="rw1")

formulae[[3]] <- births ~ dummy_step + f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2")) + f(id.tips, model="rw1")

formulae[[4]] <- births ~ f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2")) + f(id.tips, model="rw1", group=dummy_step, control.group=list(model="iid"))

formulae_number <- length(formulae)

formulae <- repeat_formulae(formulae, formulae_number)

surveys <- dhs_surveys(countryIds = c(iso2), surveyYearStart=1995) %>%
  filter(SurveyType == "DHS") %>%
  group_split(CountryName)

names(surveys) <- surveys %>%
  bind_rows %>%
  group_keys(CountryName) %>%
  pull(1)


ird <- lapply(surveys, function(surveys) {
  dhs_datasets(fileType = "IR", fileFormat = "flat", surveyIds = surveys$SurveyId)
})

ird <- lapply(ird, function(x) {
  x %>%
    mutate(path = unlist(get_datasets(x))) %>%
    bind_rows
})

# 
# 

ir <- unlist(lapply(ird, "[", "path")) %>%
  lapply(readRDS) %>%
  Map(function(test, surveys) {
    mutate(test,
           surveyid = surveys$SurveyId,
           country = surveys$CountryName,
           survyear = surveys$SurveyYear,
           survtype = surveys$SurveyType)
  }, ., surveys %>% 
    bind_rows %>%
    group_split(SurveyId))

tips_surv <- list("DHS" = c(0:10), "MIS" = c(0:5), "AIS" = c(0:5))[surveys %>%
                                                                     bind_rows %>%
                                                                     .$SurveyType]

asfr <- Map(calc_asfr1, ir,
            by = list(~country + surveyid + survyear),
            tips = tips_surv,
            agegr= list(3:10*5),
            #period = list(seq(1995, 2017, by=0.5)),
            period = list((1995*1):(2017*1)*(1/1)),
            counts = TRUE)



if (multicountry == FALSE) {

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

} else {
  
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
  
  mod_list <- Map(run_mod, formulae = rep(formulae, length(iso2)), asfr_pred = rep(asfr_pred_country, formulae_number))
  
  pred_list <- Map(get_pred, mod_list, asfr_pred_country, rep(asfr1_country, 4))
    
  
  
}
  
  # ggplot()+
  #   geom_line(data = foo, aes(x=period, y=`0.5quant`, group=agegr, color=agegr)) +
  #   geom_vline(aes(xintercept=1995:2016), linetype=3) +
  #   geom_rect(data=foo, aes(xmin = year-1, xmax = year, ymin = -Inf, ymax=Inf), fill="blue", alpha=0.2) +
  #   labs(y="ASFR", x=element_blank(), title=paste(formula)) +
  #   ylim(0, 0.34)

test <- lapply(formulae, rep, 4)
