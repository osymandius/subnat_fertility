

#' Choose a country
iso2_code <- "ZM"
iso3_code <- "ZMB"

surveys <- dhs_surveys(countryIds = c(iso2_code), surveyYearStart=1995) %>%
  filter(SurveyType != "AIS")
  
survey_text <- surveys$SurveyYear %>%
  lapply(print)
  
ird <- dhs_datasets(fileType = "IR", fileFormat = "flat", surveyIds = surveys$SurveyId)
ird$path <- unlist(get_datasets(ird))

ir <- lapply(ird$path, readRDS) %>%
  Map(data.frame,
      surveyid = surveys$SurveyId,
      country = surveys$CountryName,
      survyear = surveys$SurveyYear,
      survtype = surveys$SurveyType,
      .,
      stringsAsFactors = FALSE) 

tips_surv <- list("DHS" = c(0, 7), "MIS" = c(0, 5), "AIS" = c(0, 5))[surveys$SurveyType]

ir_loo <- lapply(1:nrow(surveys), function(x) {
  ir[-x]
})

tips_loo <- lapply(1:nrow(surveys), function(x) {
  tips_surv[-x]
})

loo_fun <- function(ir_loo, tips_loo, survey_text) {
  
  asfr <- Map(calc_asfr1, ir_loo[[1]],
              by = list(~surveyid + survyear),
              tips = tips_loo[[1]],
              agegr= list(3:10*5),
              period = list((1995*1):(2017*1)*(1/1)),
              counts = TRUE)
  
  asfr1 <- asfr %>%
    bind_rows %>%
    type.convert %>%
    filter(period<=survyear)
  
  asfr_pred <- crossing(period = asfr1$period, agegr = asfr1$agegr,  pys=1) %>%
    bind_rows(asfr1)%>%
    mutate(id.period = group_indices(., period),
           id.period2 = id.period,
           id.agegr = group_indices(., agegr),
           id.agegr2 = id.agegr,
           agegr2 = agegr,
           id.agegr.period = group_indices(., period, agegr),
           # id.district = group_indices(., area_name),
           # id.district2 = id.district,
           # id.district3 = id.district,
           # id.survey = group_indices(., surveyid),
           # id.agegr.period.district = group_indices(., agegr, period, district),
           # id.region = group_indices(., region_name),
           id = 1:nrow(.))
  
  formulae <- list()
  formulae[[7]] <- births ~ + f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2")) + f(id.agegr.period, model="iid")
  
  mod1 <- inla(formulae[[7]], family="poisson", data=asfr_pred, E=pys,
               control.family=list(link='log'),
               control.predictor=list(compute=TRUE, link=1),
               control.inla = list(strategy = "gaussian", int.strategy = "eb"),
               control.compute=list(config = TRUE, dic= TRUE, cpo=TRUE),
               verbose=TRUE)
  
  pred_size <- nrow(asfr_pred) - nrow(asfr1)
  
  pred <- asfr_pred %>%
    filter(id<pred_size+1) %>%
    dplyr::select("agegr", "period", "pys","id") %>%
    left_join(mod1$summary.fitted.values[1:pred_size, ] %>%
                mutate(id = 1:pred_size), by="id") %>%
    arrange(period, agegr) %>%
    mutate(agegr = factor(agegr),
           source = paste("No", survey_text))
  
  print(survey_text)
  
  pred
}

pred_loo <- Map(loo_fun, ir_loo, tips_loo, survey_text)

pred_loo %>%
  bind_rows %>%
  filter(agegr == "20-24") %>%
  ggplot()+
    geom_line(aes(x=period, y=`0.5quant`, group=source, color=source)) +
    geom_vline(data=surveys, aes(xintercept=SurveyYear), linetype=3) +
    geom_rect(data=surveys, aes(xmin = SurveyYear-1, xmax = SurveyYear, ymin = -Inf, ymax=Inf), fill="blue", alpha=0.2) +
    labs(y="ASFR", x=element_blank(), title="Zambia") +
    ylim(0, 0.4)
  

ggplot()+
  geom_line(data = pred, aes(x=period, y=`0.5quant`, group=agegr, color=agegr)) +
  geom_vline(data=foo, aes(xintercept=year), linetype=3) +
  geom_rect(data=foo, aes(xmin = year-1, xmax = year, ymin = -Inf, ymax=Inf), fill="blue", alpha=0.2) +
  geom_label() +
  labs(y="ASFR", x=element_blank(), title="Malawi") +
  ylim(0, 0.34)