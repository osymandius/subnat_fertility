#' ---
#' title: "Estimation of subnational fertility in SSA"
#' author: "Oli Stevens"
#' date: "`r Sys.Date()`"
#' output: 
#'   rmarkdown::html_vignette:
#'   keep_md: true
#' ---
  
##+ setup, include = FALSE
knitr::opts_chunk$set(
  collapse = TRUE,
  cache = TRUE,
  comment = "#>"
)

#' This is some intro text about how I went and did the thing.

#" ## Load packages

##+ load packages, results=FALSE, message=FALSE, warning=FALSE
library(tidyverse)
library(rdhs)
library(demogsurv)
library(lme4)

#' ## Identify and load all surveys from Malawi since 1995

#' Use rdhs to identify and retrieve individual recode (IR) datasets for all DHS surveys conducted in Malawi since 1995.
set_rdhs_config(email="o.stevens@imperial.ac.uk", project="Subnational fertility", config_path = "~/.rdhs.json")

##+ datasets
surveys <- dhs_surveys(countryIds = "MW", surveyYearStart=1995)
ird <- dhs_datasets(fileType = "IR", fileFormat = "flat", surveyIds = surveys$SurveyId)
ird$path <- unlist(get_datasets(ird))

#' Load the datasets into R as a list and asssign survey-level variables

ir <- lapply(ird$path, readRDS) %>%
  Map(data.frame,
      surveyid = surveys$SurveyId,
      country = surveys$CountryName,
      survyear = surveys$SurveyYear,
      survtype = surveys$SurveyType,
      .,
      stringsAsFactors = FALSE)
      

#' ## Calculate TFR 
#'
#' Calculate annual (calendar periods) TFR estimates for 0-10 years preceding
#' each DHS survey or 0 to 5 years preceding each MIS.

##+ warning = FALSE
## Shorten the TIPS for 2015 DHS to 0-8 due to outliers in urban for 10 year recall.
tips_surv <- list("DHS" = c(0, 7), "MIS" = c(0, 5))[surveys$SurveyType]

# names(tips_surv) <- unique(surveys$SurveyId)

# tips_surv2 <- list(c(0,5))


tfr <- Map(calc_tfr, ir,
           by = list(~surveyid + country + survyear + v025),
           tips = tips_surv,
           period = list(1995:2017))

tfr <- tfr %>%
  bind_rows %>%
  type.convert %>% 
  mutate(v025 = factor(v025, levels= c(1,2), labels = c("Urban", "Rural"))) %>%
  filter(period <= survyear) %>%
  mutate(lower = tfr - qnorm(0.975) * se_tfr,
         upper = tfr + qnorm(0.975) * se_tfr)

#' Plot annual TFR estimates

##+ fig.show='hold', fig.height=3.8, fig.width=5.0, fig.align = "center"
ggplot(tfr, aes(period, tfr, ymin = lower, ymax = upper,
                color=surveyid, fill=surveyid)) +
  geom_point() +
  geom_line() +
  geom_ribbon(alpha=0.2, linetype="blank") +
  ggtitle("TFR") +
  theme(legend.position="bottom") +
  ylim(0, 12.5) +
  facet_grid(~v025)


#' ## Calculate ASFR and births / PYs
#'
#' Calculate annual (calendar periods) TFR estimates for 0-10 years preceding
#' each DHS survey or 0 to 5 years preceding each MIS. Argument `counts = TRUE`
#' indicates to return the counts of weighted births and person-years.

### FIXED EFFECTS ONLY
{
##+ warning = FALSE
asfr <- Map(calc_asfr, ir,
            by = list(~surveyid + country + survyear+ v025),
            tips = tips_surv,
            period = list(1995:2017),
            counts = TRUE)
  
asfr <- asfr %>%
  bind_rows %>%
  type.convert %>%
  mutate(v025 = factor(v025, levels= c(1,2), labels = c("Urban", "Rural"))) %>%
  filter(period <= survyear)

#' ## Fit model
#'
#' Fit Poisson GLM with log-linear time trend for each age group

##+ warning = FALSE
mod <- glm(births ~ agegr +  period + v025, family = poisson, data = asfr, offset = log(pys))
mod.inter <- glm(births ~ agegr*period + v025, family = poisson, data = asfr, offset = log(pys))

anova(mod, mod.inter, test="LRT")

summary(mod)
summary(mod2)


#' ## Generate predictions for ASFR

#' Create a data frame for values to predict---all age groups in each year from 1995 through 2017.
df_pred <- crossing(period = 1995:2017,
                    agegr = levels(asfr$agegr),
                    pys = 1,
                    v025 = unique(asfr$v025))

#' Generate predictions log ASFR and return standard error
pred<- predict(mod, newdata = df_pred, se.fit = TRUE)


df_pred <- data.frame(df_pred, pred) %>%
  mutate(asfr = exp(fit),
         lower = exp(fit - qnorm(0.975) * se.fit),
         upper = exp(fit + qnorm(0.975) * se.fit))


#' Plot predicted vs. observed ASFR

##+ fig.show='hold', fig.height=5.0, fig.width=7, fig.align = "center"
plot <- ggplot(asfr, aes(period, asfr, color = surveyid)) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), data = df_pred, color = NA, fill = "darkred", alpha = 0.3) +
  geom_line(data = df_pred, color = "darkred") +
  #ylim(0, 0.5) +
  facet_grid(v025~agegr)
  }

gridExtra::grid.arrange(plot, plot.inter, ncol=2)

### MIXED EFFECTS
{
  
  foo2 <- Map(calc_asfr, ir,
              by = list(~surveyid + country + survyear+ v025),
              tips = list("DHS" = c(0, 7), "MIS" = c(0, 5))[surveys$SurveyType],
              period = list(1995:2017),
              counts = TRUE)
  
  ## Combine births and pys from 40-44 and 45-49 age groups to reduce 0 counts?
  asfr2 <- foo2 %>%
    bind_rows %>%
    type.convert %>%
    mutate(v025 = factor(v025, levels= c(1,2), labels = c("Urban", "Rural"))) %>%
    filter(period <= survyear) %>%
    mutate(births = round(births),
           pys = round(pys),
           sum_var = ifelse(agegr=="40-44" | agegr=="45-49", TRUE, FALSE)) %>%
    group_by(surveyid, country,  survyear,tips,  v025, period, sum_var) %>%
    mutate(births = ifelse(sum_var, sum(births), births), 
           pys = ifelse(sum_var, sum(pys), pys),
           asfr = births/pys,
           agegr = plyr::revalue(agegr, c("40-44" = "40-49")),
           agegr = factor(agegr, levels=c("15-19", "20-24", "25-29", "30-34", "35-39", "40-49"))
    ) %>%
    droplevels() %>%
    ungroup() %>%
    select(-sum_var) %>%
    distinct(surveyid, country, survyear, tips,v025, period, births, pys, asfr, .keep_all=TRUE)
  
  asfr %>%
    filter(surveyid == "MW2010DHS") %>%
    arrange(survyear, v025)
  
  unique(asfr2$agegr)
  
  mod2 <- glmer(births ~ agegr + period + v025 + (1 |agegr) , family = poisson, data = asfr2, offset = log(pys))
  mod2.int <- glmer(births ~ agegr*period + v025 + (1 + period | agegr), family = poisson, data = asfr2, offset = log(pys))
  janova(mod2.int, mod2)
  
  View(coef(mod2.int))

  df_pred2 <- crossing(period = 1995:2017,
                       agegr = levels(asfr2$agegr),
                       pys = 1,
                       v025 = levels(asfr2$v025))
  
  pred2 <- predict(mod2.int, newdata = df_pred2)
  
  df_pred2 <- data.frame(df_pred2, pred2) %>%
    mutate(asfr = exp(pred2))
  
  ggplot(asfr2, aes(period, asfr, color = surveyid)) +
    geom_point() +
    geom_line() +
    #geom_ribbon(aes(ymin = lower, ymax = upper), data = df_pred, color = NA, fill = "darkred", alpha = 0.3) +
    geom_line(data = df_pred2, aes(y=asfr), color = "darkred") +
    # facet_wrap(~v025)
    facet_grid(v025~agegr)
  
  
}

#' ## Calculate predicted TFR from predicted ASFR
#'
#' TFR is a linear transformation of ASFR given by:
#' $$ TFR = \sum_{ao \in \{15-19, \ldots, 45-49\}} 5\cdot ASFR_a.$$
#'
#' The linear transformation can be defined as a matrix $A$ such that $TFR = A \cdot ASFR.$
#' The function `.mm_aggr()` is an internal function in `demogsurv` created for this purpose.
#' (As it has broader usage, we should develop it more and make it an exported function.)
#' `.mm_aggr()` returns a list of two objects:
#' * A data frame `df` defining the stratification of the resulting aggregated transformation.
#' * A design matrix `mm` defining the transformation.
#'
#' Convert columns of `df_pred` to factors for internal workings of `.mm_aggr().`

df_predf <- df_pred %>%
  transmute(period = as.factor(period),
            agegr = as.factor(agegr))
mm <- demogsurv:::.mm_aggr(df_predf, seq(15, 50, 5))

tfr_pred <- data.frame(mm$df,
                       tfr = c(df_pred$asfr %*% mm$mm))

#' Observe that this gives the same results as using `group_by()` and `summarise()`.
tfr_pred %>%
  head

df_pred %>%
  group_by(period) %>%
  summarise(tfr = sum(5 * asfr)) %>%
  head


#' #### Variance for TFR
#'
#' If $y = Ax$ is a linear tranformation, then $$var(y) = var(Ax) = A\cdot var(x) \cdot A^T.$$
#' Thus, to estimate the variance of TFR, we need to:
#' 
#' 1. Obtain the covariance matrix $V_{log(asfr)}$ for the log ASFR values generated by `predict.glm()`.
#' 1. Use the delta method to approximate the covariance matrix $V_{asfr}$ for $asfr = exp(log(asfr))$.
#' 1. Use the formula above to obtain the covariance matrix for the predicted TFR $V_{tfr} = A\cdot V_{asfr} \cdot A^T.$
#'
#' First, create a design matrix for the predictions. This is the first internal step of `predict.lm()`.

Xpred <- model.matrix(delete.response(terms(mod)), data = df_pred)

#' `Xpred %*% coef(mod)` returns the same values as the log ASFR returned by `predict.glm`.

log_asfr_pred <- c(Xpred %*% coef(mod))

all.equal(log_asfr_pred,
          unname(predict(mod, df_pred, type = "link")))

#' Calculate the covariance matrix of the predicted log ASFR values
V_log_asfr <- Xpred %*% vcov(mod) %*% t(Xpred)

#' This covariance matrix returns the same standard errors as `predict(..., se = TRUE)`.
all.equal(sqrt(diag(V_log_asfr)),
          predict(mod, df_pred, type = "link", se = TRUE)$se.fit)

#' Since $\frac{d}{dx} exp(x) = exp(x)$, the delta method approximation for $V_{asfr}$ is
asfr_pred <- exp(log_asfr_pred)
V_asfr <- diag(asfr_pred) %*% V_log_asfr %*% diag(asfr_pred)

#' Finally, use the linear transformation $A$ above to calculate TFR from ASFR

V_tfr <- t(mm$mm) %*% V_asfr %*% mm$mm

tfr_pred <- tfr_pred %>%
  type.convert %>%
  mutate(se = sqrt(diag(V_tfr)),
         lower = tfr - qnorm(0.975) * se,
         upper = tfr + qnorm(0.975) * se)
         
#' Plot modelled TFR compared to direct survey estimates

##+ fig.show='hold', fig.height=3.8, fig.width=5.0, fig.align = "center"
ggplot(tfr, aes(period, tfr, ymin = lower, ymax = upper,
                color=surveyid, fill=surveyid)) +
  geom_point() +
  geom_line() +
  geom_ribbon(alpha=0.1, linetype="blank") +
  geom_line(aes(fill = NA), data = tfr_pred, color = "darkred", size = 1.2) +
  geom_ribbon(data = tfr_pred, color = NA, fill = "darkred", alpha = 0.2) +
  ggtitle("TFR") +
  theme(legend.position="bottom")
