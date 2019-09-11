#' Choose a country
iso2_code <- "NM"
iso3_code <- "NAM"

iso2 <- c("MW", "RW", "UG", "MZ", "NM", "KE", "ZM", "ZW", "TZ")
iso2_list <- lapply(iso2, print)
iso3_list <- lapply(iso2_list, function(x) {x %>% countrycode("iso2c", "iso3c")})
iso3_list[5] <- "NAM"


fixed_eff_fun <- function(iso2_code, iso3_code, survey_type) {

    surveys <- dhs_surveys(countryIds = c(iso2_code), surveyYearStart=1995) %>%
      filter(SurveyType == survey_type)
    
    # foo <- surveys %>%
    #   separate(SurveyId, into=c("cc", "year", "type"), sep = c(2,6)) %>%
    #   type.convert %>%
    #   dplyr::select(year)
    
    ird <- dhs_datasets(fileType = "IR", fileFormat = "flat", surveyIds = surveys$SurveyId)
    ird$path <- unlist(get_datasets(ird))
    
    # 
    # 
    ir <- lapply(ird$path, readRDS) %>%
      Map(data.frame,
          surveyid = surveys$SurveyId,
          country = surveys$CountryName,
          survyear = surveys$SurveyYear,
          survtype = surveys$SurveyType,
          .,
          stringsAsFactors = FALSE) 
    
    tips_surv <- list("DHS" = c(0:10), "MIS" = c(0:5), "AIS" = c(0:5))[surveys$SurveyType]
    
    asfr <- Map(calc_asfr1, ir,
                by = list(~surveyid + survyear),
                tips = tips_surv,
                agegr= list(3:10*5),
                #period = list(seq(1995, 2017, by=0.5)),
                period = list((1995*1):(2017*1)*(1/1)),
                counts = TRUE)
    
    asfr1 <- asfr %>%
      bind_rows %>%
      separate(period, into=c("period", "stuff"), sep="-") %>%
      dplyr::select(-stuff) %>%
      type.convert %>%
      mutate(id.period = group_indices(., period),
             id.period2 = id.period,
             id.agegr = group_indices(., agegr),
             id.agegr2 = id.agegr,
             agegr2 = agegr,
             id.agegr.period = group_indices(., period, agegr),
             id.tips = group_indices(., tips),
             # id.district = group_indices(., area_name),
             # id.district2 = id.district,
             # id.district3 = id.district,
             # id.survey = group_indices(., surveyid),
             # id.agegr.period.district = group_indices(., agegr, period, district),
             # id.region = group_indices(., region_name),
             id = 1:nrow(.),
             split_period = as.character(period)) %>%
      separate(split_period, into=c("quot", "remainder"), sep=5, remove=TRUE) %>% 
      mutate(id.quarter = ifelse(remainder==25, 2, 
                                 ifelse(remainder==5, 3, 
                                        ifelse(remainder==75, 4, 1
                                               
                                        )
                                 )),
             id.quarter = factor(id.quarter),
             id.tips = factor(id.tips)
      ) %>%
      dplyr::select(-c(quot, remainder))
    
    # asfr_pred <- crossing(surveyid = asfr1$surveyid, period = asfr1$period, agegr = asfr1$agegr,  pys=1) %>%
    #   bind_rows(asfr1) %>%
    #   mutate(id.period = group_indices(., period),
    #          id.period2 = id.period,
    #          id.agegr = group_indices(., agegr),
    #          id.agegr2 = id.agegr,
    #          agegr2 = agegr,
    #          id.agegr.period = group_indices(., period, agegr),
    #          id.tips = group_indices(., tips),
    #          # id.district = group_indices(., area_name),
    #          # id.district2 = id.district,
    #          # id.district3 = id.district,
    #          # id.survey = group_indices(., surveyid),
    #          # id.agegr.period.district = group_indices(., agegr, period, district),
    #          # id.region = group_indices(., region_name),
    #          id = 1:nrow(.),
    #          split_period = as.character(period)) %>%
    #   separate(split_period, into=c("quot", "remainder"), sep=5, remove=TRUE) %>% 
    #   mutate(id.quarter = ifelse(remainder==25, 2, 
    #                              ifelse(remainder==5, 3, 
    #                                     ifelse(remainder==75, 4, 1
    #                                            
    #                                     )
    #                              )),
    #          id.quarter = factor(id.quarter)
    #   ) %>%
    #   dplyr::select(-c(quot, remainder))
    
    
    formulae <- list()
    formulae[[8]] <- births ~ id.quarter + f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2"))
    formulae[[9]] <- births ~ id.tips + f(id.agegr, model="rw1") + f(id.period, model="rw2") + f(id.agegr2, model="rw1", group=id.period, control.group=list(model="rw2"))
    
    mod1 <- inla(formulae[[9]], family="poisson", data=asfr1, E=pys,
                 control.family=list(link='log'),
                 control.predictor=list(compute=TRUE, link=1),
                 control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                 control.compute=list(config = TRUE, dic= TRUE, cpo=TRUE),
                 verbose=TRUE)
    
    # pred_size <- nrow(asfr_pred) - nrow(asfr1)
    # 
    # pred <- asfr_pred %>%
    #   filter(id<pred_size+1) %>%
    #   dplyr::select("agegr", "period", "pys","id") %>%
    #   left_join(mod1$summary.fitted.values[1:pred_size, ] %>%
    #               mutate(id = 1:pred_size), by="id") %>%
    #   arrange(period, agegr) %>%
    #   #filter(agegr <20) %>%
    #   mutate(agegr = factor(agegr))
    # # mutate(agegroup = rep(1:7, each=5, times=pred_size/35),
    # #        agegroup = factor(agegroup, levels=1:7, labels=c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49")))
    
    # ggplot()+
    #   geom_line(data = pred, aes(x=period, y=`0.5quant`, group=agegr, color=agegr)) +
    #   geom_vline(aes(xintercept=1995:2016), linetype=3) +
    #   geom_rect(data=foo, aes(xmin = year-1, xmax = year, ymin = -Inf, ymax=Inf), fill="blue", alpha=0.2) +
    #   labs(y="ASFR", x=element_blank(), title=paste(iso3_code, "0.25 year periods, without iid"))
    
    
    fixed_eff <- exp(mod1$summary.fixed) %>%
      mutate(country = iso2_code)
    
    fixed_eff 
}

fixed_eff <- Map(fixed_eff_fun, iso2_list, iso3_list, list("DHS"))

fixed_eff %>%
  bind_rows %>%
  filter(mean>0.2) %>%
  mutate(tips = rep(1:9, times=length(iso2_list))) %>%
  bind_rows(data.frame("mean" = 1, "tips" = 0, "country"=unlist(iso2_list))) %>%
  mutate(tips = factor(tips)) %>%
  ggplot(aes(x=tips)) +
  geom_point(aes(y=mean)) +
  geom_errorbar(aes(ymin = `0.025quant`, ymax=`0.975quant`), width=0.2) +
  geom_hline(aes(yintercept=1), linetype=3)+
  facet_wrap(~country)

debugonce(fixed_eff_fun)
fixed_eff_fun(iso2_list[[3]], iso3_list[[3]], "DHS")

quarter_summary %>%
  filter(mean>0.2) %>%
  mutate(tips = rep(1:6, times=length(unique(quarter_summary$country)))) %>%
  bind_rows(data.frame("mean" = 1, "tips" = 0, "country"=unique(quarter_summary$country))) %>%
  mutate(tips = factor(tips)) %>%
  ggplot(aes(x=tips)) +
  geom_point(aes(y=mean)) +
  geom_errorbar(aes(ymin = `0.025quant`, ymax=`0.975quant`), width=0.2) +
  geom_hline(aes(yintercept=1), linetype=3)+
  facet_wrap(~country)
