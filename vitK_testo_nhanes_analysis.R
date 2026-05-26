library(broom)
library(dplyr)
library(dietaryindex)
library(emmeans)
library(ggplot2)
library(gtsummary)
library(haven)
library(here)
library(janitor)
library(labelled)
library(marginaleffects)
library(nhanesA)
library(openxlsx)
library(purrr)
library(splines)
library(srvyr)
library(survey)
library(tableone)
library(tidyverse)

#--------------------------------------------------
# Analysis settings
#--------------------------------------------------
options(survey.lonely.psu = "adjust")

set.seed(2026)

dir.create(here("results"), showWarnings = FALSE)

dir.create(here("figures"), showWarnings = FALSE)

read_cycle <- function(file, cycle_label){
  read_xpt(here("data_raw", file)) %>%
    clean_names() %>%
    mutate(cycle = cycle_label)
}

#------------------------data import--------------------------

demo_1112 <- read_cycle("DEMO_G.xpt", "2011-2012")
demo_1314 <- read_cycle("DEMO_H.xpt", "2013-2014")
demo_1516 <- read_cycle("DEMO_I.xpt", "2015-2016")

demo_allagesex <- bind_rows(demo_1112, demo_1314, demo_1516)

demo_all <- demo_allagesex %>%
  filter(
    riagendr == 1,
    ridageyr >= 18
  ) %>%
  select(
    seqn,
    cycle,
    ridageyr,
    ridreth1,
    indfmpir,
    sdmvpsu,
    sdmvstra,
    wtmec2yr
  )

diet_1112 <- read_cycle("DR1TOT_G.xpt", "2011-2012")
diet_1314 <- read_cycle("DR1TOT_H.xpt", "2013-2014")
diet_1516 <- read_cycle("DR1TOT_I.xpt", "2015-2016")

diet_all <- bind_rows(diet_1112, diet_1314, diet_1516) %>%
  select(seqn, cycle, dr1tvk, dr1tkcal)


tst_1112 <- read_cycle("TST_G.xpt", "2011-2012")
tst_1314 <- read_cycle("TST_H.xpt", "2013-2014")
tst_1516 <- read_cycle("TST_I.xpt", "2015-2016")

tst_all <- bind_rows(tst_1112, tst_1314, tst_1516) %>%
  select(seqn, cycle, lbxtst, lbxshbg)


bmx_1112 <- read_cycle("BMX_G.xpt", "2011-2012")
bmx_1314 <- read_cycle("BMX_H.xpt", "2013-2014")
bmx_1516 <- read_cycle("BMX_I.xpt", "2015-2016")

bmx_all <- bind_rows(bmx_1112, bmx_1314, bmx_1516) %>%
  select(seqn, cycle, bmxbmi, bmxwaist)


glu_1112 <- read_cycle("GLU_G.xpt", "2011-2012") %>%
  select(seqn, cycle, lbxglu, lbxin)

glu_1314 <- read_cycle("GLU_H.xpt", "2013-2014") %>%
  select(seqn, cycle, lbxglu)

glu_1516 <- read_cycle("GLU_I.xpt", "2015-2016") %>%
  select(seqn, cycle, lbxglu)

in_1314 <- read_cycle("INS_H.xpt", "2013-2014") %>%
  select(seqn, cycle, lbxin)

in_1516 <- read_cycle("INS_I.xpt", "2015-2016") %>%
  select(seqn, cycle, lbxin)

glu_in_1314 <- left_join(glu_1314, in_1314, by = c("seqn","cycle"))
glu_in_1516 <- left_join(glu_1516, in_1516, by = c("seqn","cycle"))

glu_all <- bind_rows(glu_1112,
                     glu_in_1314,
                     glu_in_1516)


smq_1112 <- read_cycle("SMQ_G.xpt", "2011-2012")
smq_1314 <- read_cycle("SMQ_H.xpt", "2013-2014")
smq_1516 <- read_cycle("SMQ_I.xpt", "2015-2016")

smq_all <- bind_rows(smq_1112, smq_1314, smq_1516) %>%
  transmute(seqn, cycle,
            smoking =
              case_when(
                smq020 == 2 ~ "Never",
                smq020 == 1 & smq040 == 3 ~ "Former",
                smq020 == 1 & smq040 %in% c(1,2) ~ "Current",
                TRUE ~ NA_character_
              ))


alq_1112 <- read_cycle("ALQ_G.xpt", "2011-2012")
alq_1314 <- read_cycle("ALQ_H.xpt", "2013-2014")
alq_1516 <- read_cycle("ALQ_I.xpt", "2015-2016")

alq_all <- bind_rows(alq_1112, alq_1314, alq_1516) %>%
  transmute(
    seqn,
    cycle,
    alcohol = case_when(
      alq110 == 2 ~ "Never",
      alq101 == 2 & alq110 == 1 ~ "Former",
      alq101 == 1 ~ "Current",
      TRUE ~ NA_character_
    )
  )

paq_1112 <- read_cycle("PAQ_G.xpt", "2011-2012")
paq_1314 <- read_cycle("PAQ_H.xpt", "2013-2014")
paq_1516 <- read_cycle("PAQ_I.xpt", "2015-2016")

paq_all <- bind_rows(paq_1112, paq_1314, paq_1516) %>%
  transmute(
    seqn,
    cycle,
    phys_active = case_when(
      paq650 == 1 | paq665 == 1 ~ "Active",
      paq650 == 2 & paq665 == 2 ~ "Inactive",
      TRUE ~ NA_character_
    ),
    phys_active = factor(
      phys_active,
      levels = c("Inactive", "Active")
    )
  )


fast_1112 <- read_cycle("FASTQX_G.xpt", "2011-2012")
fast_1314 <- read_cycle("FASTQX_H.xpt", "2013-2014")
fast_1516 <- read_cycle("FASTQX_I.xpt", "2015-2016")

fast_all <- bind_rows(fast_1112, fast_1314, fast_1516) %>%
  select(seqn, cycle, phdsesn)


mcq_1112 <- read_cycle("MCQ_G.xpt", "2011-2012")
mcq_1314 <- read_cycle("MCQ_H.xpt", "2013-2014")
mcq_1516 <- read_cycle("MCQ_I.xpt", "2015-2016")

mcq_all <- bind_rows(mcq_1112, mcq_1314, mcq_1516) %>%
  select(seqn, cycle, mcq230a)


rx_1112 <- read_cycle("RXQ_RX_G.xpt", "2011-2012")
rx_1314 <- read_cycle("RXQ_RX_H.xpt", "2013-2014")
rx_1516 <- read_cycle("RXQ_RX_I.xpt", "2015-2016")


rx_all <- bind_rows(rx_1112, rx_1314, rx_1516) %>%
  group_by(seqn, cycle) %>%
  summarise(
    testosterone_med = any(
      str_detect(
        str_to_upper(rxddrug),
        "TESTOSTERONE"
      ),
      na.rm = TRUE
    ),
    .groups = "drop"
  )


analysis_raw <- demo_all %>%
  left_join(diet_all, by=c("seqn","cycle")) %>%
  left_join(tst_all, by=c("seqn","cycle")) %>%
  left_join(bmx_all, by=c("seqn","cycle")) %>%
  left_join(glu_all, by=c("seqn","cycle")) %>%
  left_join(smq_all, by=c("seqn","cycle")) %>%
  left_join(alq_all, by=c("seqn","cycle")) %>%
  left_join(paq_all, by=c("seqn","cycle")) %>%
  left_join(fast_all, by=c("seqn","cycle")) %>%
  left_join(mcq_all, by=c("seqn","cycle")) %>%
  left_join(rx_all, by=c("seqn","cycle"))



#------------------------participant flowchart--------------------------

# All NHANES participants
flow_0 <- demo_allagesex

n0 <- nrow(flow_0)

# Adult men only
flow_1 <- flow_0 %>%
  filter(
    riagendr == 1,
    ridageyr >= 18
  )

excluded_agesex <- n0 - nrow(flow_1)
n1 <- nrow(flow_1)

# Merge analytic datasets
flow_2 <- flow_1 %>%
  left_join(diet_all, by = c("seqn","cycle")) %>%
  left_join(tst_all, by = c("seqn","cycle")) %>%
  left_join(bmx_all, by = c("seqn","cycle")) %>%
  left_join(glu_all, by = c("seqn","cycle")) %>%
  left_join(smq_all, by = c("seqn","cycle")) %>%
  left_join(alq_all, by = c("seqn","cycle")) %>%
  left_join(paq_all, by = c("seqn","cycle")) %>%
  left_join(fast_all, by = c("seqn","cycle")) %>%
  left_join(mcq_all, by = c("seqn","cycle")) %>%
  left_join(rx_all, by = c("seqn","cycle"))

# Exclude missing testosterone
flow_3 <- flow_2 %>%
  filter(!is.na(lbxtst))

excluded_t <- nrow(flow_2) - nrow(flow_3)
n3 <- nrow(flow_3)

# Exclude missing vitamin K intake
flow_4 <- flow_3 %>%
  filter(
    !is.na(dr1tvk),
    !is.na(dr1tkcal)
  )

excluded_vk <- nrow(flow_3) - nrow(flow_4)
n4 <- nrow(flow_4)

# Exclude missing covariates
flow_5 <- flow_4 %>%
  filter(
    !is.na(indfmpir),
    !is.na(smoking),
    !is.na(alcohol),
    !is.na(phys_active),
    !is.na(phdsesn)
  )

excluded_covars <- nrow(flow_4) - nrow(flow_5)
n5 <- nrow(flow_5)

# Exclude testosterone therapy
flow_6 <- flow_5 %>%
  filter(
    testosterone_med == FALSE |
      is.na(testosterone_med)
  )

excluded_trt <- nrow(flow_5) - nrow(flow_6)
n_final <- nrow(flow_6)

# Flowchart summary table
flowchart_table <- tibble(
  Step = c(
    "NHANES 2011-2016 participants",
    "Excluded women and age <18 years",
    "Adult men aged ≥18 years",
    "Excluded missing testosterone",
    "Participants with testosterone data",
    "Excluded missing vitamin K intake",
    "Participants with vitamin K data",
    "Excluded missing covariates",
    "Participants with complete covariate data",
    "Excluded testosterone therapy",
    "Final analytic sample"
  ),
  
  N = c(
    n0,
    excluded_agesex,
    n1,
    excluded_t,
    n3,
    excluded_vk,
    n4,
    excluded_covars,
    n5,
    excluded_trt,
    n_final
  )
)

write.xlsx(
  flowchart_table,
  here("results", "Figure1_Analysis_Population.xlsx"),
  rowNames = FALSE
)


#-----------------------------cleanup--------------------------------
analysis_clean <- analysis_raw %>%
  mutate(
    wtmec6yr = wtmec2yr / 3,
    wtmec4yr = wtmec2yr / 2,
    vk_density = dr1tvk / (dr1tkcal / 1000),
    homa_ir = (lbxin * lbxglu) / 405,
    lowt = if_else(lbxtst < 300, 1, 0),
    
    race = factor(
      ridreth1,
      levels = c(1,2,3,4,5),
      labels = c(
        "Mexican American",
        "Other Hispanic",
        "Non-Hispanic White",
        "Non-Hispanic Black",
        "Other race"
      )
    ),
    
    smoking = factor(
      smoking,
      levels = c("Never", "Former", "Current")
    ),
    
    alcohol = factor(
      alcohol,
      levels = c("Never", "Former", "Current")
    ),
    
    exam_session = factor(
      phdsesn,
      levels = c(0,1,2),
      labels = c(
        "Morning",
        "Afternoon",
        "Evening"
      )
    )
  ) %>%
  
  filter(
    !is.na(lbxtst),
    !is.na(vk_density),
    !is.na(ridageyr),
    !is.na(sdmvpsu),
    !is.na(sdmvstra),
    !is.na(wtmec6yr)
  ) %>%
  
  filter(
    testosterone_med == FALSE |
      is.na(testosterone_med)
  )


design_1116 <- svydesign(
  ids = ~sdmvpsu,
  strata = ~sdmvstra,
  weights = ~wtmec6yr,
  data = analysis_clean,
  nest = TRUE
)


#===========================================================
# MAIN ANALYSIS
#===========================================================
vk_sd <- sqrt(
  as.numeric(
    svyvar(~vk_density, design_1116)
  )
)

vk_mean <- as.numeric(
  svymean(~vk_density, design_1116, na.rm = TRUE)
)

analysis_clean <- analysis_clean %>%
  mutate(
    vk_density_sd =
      (vk_density - vk_mean) / vk_sd
  )

design_1116 <- update(
  design_1116,
  vk_density_sd = analysis_clean$vk_density_sd
)

m1 <- svyglm(
  lbxtst ~ vk_density_sd + ridageyr,
  design = design_1116
)

m1_tab <- tidy(m1, conf.int = TRUE) %>%
  filter(term == "vk_density_sd") %>%
  transmute(
    Model = "Model 1 (Age-adjusted)",
    Exposure = "Vitamin K density (per 1 SD)",
    Beta = estimate,
    CI_low = conf.low,
    CI_high = conf.high,
    P_value = p.value
  )

m2 <- svyglm(
  lbxtst ~ vk_density_sd +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session,
  design = design_1116
)

m2_tab <- tidy(m2, conf.int = TRUE) %>%
  filter(term == "vk_density_sd") %>%
  transmute(
    Model = "Model 2 (Fully adjusted)",
    Exposure = "Vitamin K density (per 1 SD)",
    Beta = estimate,
    CI_low = conf.low,
    CI_high = conf.high,
    P_value = p.value
  )

m3 <- svyglm(
  lbxtst ~ vk_density_sd +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session +
    bmxbmi +
    bmxwaist,
  design = design_1116
)

m3_tab <- tidy(m3, conf.int = TRUE) %>%
  filter(term == "vk_density_sd") %>%
  transmute(
    Model = "Model 3 (+ Adiposity)",
    Exposure = "Vitamin K density (per 1 SD)",
    Beta = estimate,
    CI_low = conf.low,
    CI_high = conf.high,
    P_value = p.value
  )

main_results <- bind_rows(m1_tab, m2_tab, m3_tab)

write.xlsx(
  main_results,
  here("results", "M1to3_Continuous_VK.xlsx")
)


#------------------------quantify adiposity attenuation--------------------------
beta_m2 <- tidy(m2) %>%
  filter(term == "vk_density_sd") %>%
  pull(estimate)

beta_m3 <- tidy(m3) %>%
  filter(term == "vk_density_sd") %>%
  pull(estimate)

attenuation_pct <- 100 * (beta_m2 - beta_m3) / beta_m2

attenuation_table <- tibble(
  Beta_Model2 = beta_m2,
  Beta_Model3 = beta_m3,
  Percent_attenuation = attenuation_pct
) %>%
  mutate(
    Beta_Model2 = round(Beta_Model2, 3),
    Beta_Model3 = round(Beta_Model3, 3),
    Percent_attenuation = round(Percent_attenuation, 1)
  )

write.xlsx(
  attenuation_table,
  here("results", "Adiposity_Attenuation.xlsx")
)


#------------------------quartile/table 1 creation--------------------------
vk_quants <- svyquantile(
  ~vk_density,
  design_1116,
  quantiles = c(0.25, 0.5, 0.75),
  na.rm = TRUE
)

vk_cut <- as.numeric(coef(vk_quants))

analysis_clean <- analysis_clean %>%
  mutate(
    vk_q = cut(
      vk_density,
      breaks = c(-Inf,
                 vk_cut[1],
                 vk_cut[2],
                 vk_cut[3],
                 Inf),
      labels = c("Q1","Q2","Q3","Q4")
    )
  )

analysis_clean$vk_q <- factor(analysis_clean$vk_q,
                              levels = c("Q1","Q2","Q3","Q4"))

design_1116 <- update(design_1116,
                      vk_q = analysis_clean$vk_q)

# Unweighted N
n_overall <- nrow(analysis_clean)

n_by_q <- analysis_clean %>%
  count(vk_q)

# Table 1
table_1 <- tbl_svysummary(
  design_1116,
  by = vk_q,
  include = c(
    vk_density,
    ridageyr,
    race,
    indfmpir,
    dr1tkcal,
    bmxbmi,
    bmxwaist,
    lbxtst,
    smoking,
    alcohol,
    phys_active,
    testosterone_med
  ),
  statistic = list(
    vk_density ~ "{median} ({p25}, {p75})",
    ridageyr ~ "{mean} ({sd})",
    indfmpir ~ "{mean} ({sd})",
    dr1tkcal ~ "{mean} ({sd})",
    bmxbmi ~ "{mean} ({sd})",
    bmxwaist ~ "{mean} ({sd})",
    lbxtst ~ "{mean} ({sd})",
    all_categorical() ~ "{p}%"
  ),
  digits = list(
    vk_density ~ 1,
    all_continuous() ~ 1
  ),
  missing = "no"
) %>%
  add_overall() %>%
  bold_labels()

# Add N to headers
header_map <- c(
  stat_0 = paste0("Overall\nN = ", n_overall),
  setNames(
    paste0(n_by_q$vk_q, "\nN = ", n_by_q$n),
    paste0("stat_", seq_len(nrow(n_by_q)))
  )
)

table_1 <- table_1 %>%
  modify_header(!!!header_map)

# Export
table_1_df <- as_tibble(table_1, col_labels = TRUE)

write.xlsx(
  table_1_df,
  here("results", "Table1_StudyPopulation_by_VK_Quartiles.xlsx"),
  overwrite = TRUE
)


#------------------------quartile analyses--------------------------
m1_q <- svyglm(
  lbxtst ~ vk_q + ridageyr,
  design = design_1116
)

m2_q <- svyglm(
  lbxtst ~ vk_q +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session,
  design = design_1116
)

extract_q <- function(model, model_name) {
  tidy(model, conf.int = TRUE) %>%
    filter(str_detect(term, "vk_q")) %>%
    mutate(
      Quartile = case_when(
        term == "vk_qQ2" ~ "Q2 vs Q1",
        term == "vk_qQ3" ~ "Q3 vs Q1",
        term == "vk_qQ4" ~ "Q4 vs Q1"
      ),
      Model = model_name
    ) %>%
    select(Model,
           Quartile,
           Beta = estimate,
           CI_low = conf.low,
           CI_high = conf.high,
           P_value = p.value)
}

tab_m1_q <- extract_q(m1_q, "Model 1 (Age-adjusted)")
tab_m2_q <- extract_q(m2_q, "Model 2 (Fully adjusted)")

quartile_results <- bind_rows(tab_m1_q, tab_m2_q)

vk_medians <- analysis_clean %>%
  group_by(vk_q) %>%
  summarise(
    vk_q_median = median(vk_density, na.rm = TRUE),
    .groups = "drop"
  )

analysis_clean <- analysis_clean %>%
  left_join(vk_medians, by = "vk_q")

design_1116 <- update(
  design_1116,
  vk_q_median = analysis_clean$vk_q_median
)

trend_m2 <- svyglm(
  lbxtst ~ vk_q_median +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session,
  design = design_1116
)

p_trend <- tidy(trend_m2) %>%
  filter(term == "vk_q_median") %>%
  pull(p.value)

trend_row <- tibble(
  Model = "Model 2 (Fully adjusted)",
  Quartile = "P for trend",
  Beta = NA,
  CI_low = NA,
  CI_high = NA,
  P_value = p_trend
)

quartile_results <- bind_rows(quartile_results, trend_row)

quartile_results_clean <- quartile_results %>%
  mutate(
    `β (95% CI)` =
      ifelse(is.na(Beta),
             "",
             sprintf("%.2f (%.2f, %.2f)",
                     Beta, CI_low, CI_high)),
    P_value = ifelse(is.na(P_value),
                     "",
                     sprintf("%.3f", P_value))
  ) %>%
  select(Model, Quartile, `β (95% CI)`, P_value)

write.xlsx(
  quartile_results_clean,
  here("results", "Quartile_Model1and2.xlsx")
)

#===========================================================
# SENSITIVITY ANALYSES
#===========================================================
#------------------------sensitivity, vk log--------------------------

analysis_clean <- analysis_clean %>%
  mutate(
    log_vk = log(vk_density + 1)
  )

design_1116 <- update(
  design_1116,
  log_vk = analysis_clean$log_vk
)

log_vk_sd <- sqrt(
  as.numeric(
    svyvar(~log_vk, design_1116)
  )
)

log_vk_mean <- as.numeric(
  svymean(~log_vk, design_1116, na.rm = TRUE)
)

analysis_clean <- analysis_clean %>%
  mutate(
    log_vk_sd_var =
      (log_vk - log_vk_mean) / log_vk_sd
  )

design_1116 <- update(
  design_1116,
  log_vk_sd_var = analysis_clean$log_vk_sd_var
)

m1_logvk <- svyglm(
  lbxtst ~ log_vk_sd_var + ridageyr,
  design = design_1116
)

m2_logvk <- svyglm(
  lbxtst ~ log_vk_sd_var +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session,
  design = design_1116
)

# Extract results
extract_log <- function(model, model_name) {
  tidy(model, conf.int = TRUE) %>%
    filter(term == "log_vk_sd_var") %>%
    transmute(
      Model = model_name,
      Exposure = "Log vitamin K density (per 1 SD)",
      Beta = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      P_value = p.value
    )
}

tab_m1_logvk <- extract_log(m1_logvk, "Model 1 (Age-adjusted)")
tab_m2_logvk <- extract_log(m2_logvk, "Model 2 (Fully adjusted)")

log_results <- bind_rows(tab_m1_logvk, tab_m2_logvk)

log_results_clean <- log_results %>%
  mutate(
    `β (95% CI)` =
      sprintf("%.2f (%.2f, %.2f)",
              Beta, CI_low, CI_high),
    P_value = sprintf("%.3f", P_value)
  ) %>%
  select(Model, Exposure, `β (95% CI)`, P_value)

write.xlsx(
  log_results_clean,
  here("results", "Log_VK_Model_1and2.xlsx")
)

#------------------------sensitivity, T log--------------------------

analysis_clean <- analysis_clean %>%
  mutate(
    log_t = log(lbxtst)
  )

design_1116 <- update(design_1116,
                      log_t = analysis_clean$log_t)

m1_logT <- svyglm(
  log_t ~ vk_density_sd +
    ridageyr,
  design = design_1116
)

m2_logT <- svyglm(
  log_t ~ vk_density_sd +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session,
  design = design_1116
)


extract_logT <- function(model, model_name) {
  tidy(model, conf.int = TRUE) %>%
    filter(term == "vk_density_sd") %>%
    mutate(
      Model = model_name,
      Beta = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      Ratio = exp(estimate),
      Ratio_low = exp(conf.low),
      Ratio_high = exp(conf.high)
    ) %>%
    select(Model,
           Beta,
           CI_low,
           CI_high,
           Ratio,
           Ratio_low,
           Ratio_high,
           p.value)
}

tab_m1_logT <- extract_logT(m1_logT, "Model 1 (Age-adjusted)")
tab_m2_logT <- extract_logT(m2_logT, "Model 2 (Fully adjusted)")

logT_results <- bind_rows(tab_m1_logT, tab_m2_logT)

logT_results_clean <- logT_results %>%
  mutate(
    `β (log scale)` = sprintf("%.4f (%.4f, %.4f)",
                              Beta, CI_low, CI_high),
    `Ratio (95% CI)` = sprintf("%.3f (%.3f, %.3f)",
                               Ratio, Ratio_low, Ratio_high),
    P_value = sprintf("%.3f", p.value)
  ) %>%
  select(Model, `β (log scale)`, `Ratio (95% CI)`, P_value)

write.xlsx(
  logT_results_clean,
  here("results", "Log_Testosterone_Model1_and_Model2.xlsx")
)

#-----------------------sensitivity, low T only--------------------------

m1_lowT <- svyglm(
  lowt ~ vk_density_sd + ridageyr,
  design = design_1116,
  family = quasibinomial()
)

m2_lowT <- svyglm(
  lowt ~ vk_density_sd +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session,
  design = design_1116,
  family = quasibinomial()
)

extract_or <- function(model, model_name) {
  tidy(model, conf.int = TRUE) %>%
    filter(term == "vk_density_sd") %>%
    mutate(
      OR = exp(estimate),
      CI_low = exp(conf.low),
      CI_high = exp(conf.high),
      Model = model_name
    ) %>%
    select(Model, OR, CI_low, CI_high, p.value)
}

tab_m1_lowT <- extract_or(m1_lowT, "Model 1 (Age-adjusted)")
tab_m2_lowT <- extract_or(m2_lowT, "Model 2 (Fully adjusted)")

lowT_results <- bind_rows(tab_m1_lowT, tab_m2_lowT) %>%
  mutate(
    `OR (95% CI)` = sprintf("%.2f (%.2f, %.2f)",
                            OR, CI_low, CI_high),
    P_value = sprintf("%.3f", p.value)
  ) %>%
  select(Model, `OR (95% CI)`, P_value)

write.xlsx(
  lowT_results,
  here("results", "LowT_Logistic_Models.xlsx")
)

#-----------------------sensitivity, hei--------------------------

cycles_hei <- tibble(
  cycle = c("2011-2012","2013-2014","2015-2016"),
  suffix = c("G","H","I")
)

hei_scores <- pmap_dfr(cycles_hei, function(cycle, suffix) {
  
  fped_file <- case_when(
    suffix == "G" ~ "fped_dr1tot_1112.sas7bdat",
    suffix == "H" ~ "fped_dr1tot_1314.sas7bdat",
    suffix == "I" ~ "fped_dr1tot_1516.sas7bdat"
  )
  
  hei_cycle <- HEI2020_NHANES_FPED(
    FPED_PATH     = here("data_raw", fped_file),
    NUTRIENT_PATH = here("data_raw", paste0("DR1TOT_", suffix, ".xpt")),
    DEMO_PATH     = here("data_raw", paste0("DEMO_", suffix, ".xpt"))
  )
  
  hei_cycle$cycle <- cycle
  
  return(hei_cycle)
})

hei_final <- hei_scores %>%
  select(seqn = SEQN,
         cycle,
         hei_total = HEI2020_ALL)

analysis_clean_hei <- analysis_clean %>%
  left_join(hei_final, by = c("seqn","cycle"))

design_1116_hei <- svydesign(
  ids = ~sdmvpsu,
  strata = ~sdmvstra,
  weights = ~wtmec6yr,
  data = analysis_clean_hei,
  nest = TRUE
)

m2_hei <- svyglm(
  lbxtst ~ vk_density_sd +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session +
    hei_total,
  design = design_1116_hei
)

tab_m2 <- tidy(m2, conf.int = TRUE) %>%
  filter(term == "vk_density_sd") %>%
  mutate(
    Model = "Model 2 (Fully adjusted)"
  )

tab_m2_hei <- tidy(m2_hei, conf.int = TRUE) %>%
  filter(term == "vk_density_sd") %>%
  mutate(
    Model = "Model 2 + HEI"
  )

beta_m2 <- tab_m2$estimate
beta_m2_hei <- tab_m2_hei$estimate

attenuation_pct <- 100 * (beta_m2 - beta_m2_hei) / beta_m2

comparison_table <- bind_rows(tab_m2, tab_m2_hei) %>%
  mutate(
    `β (95% CI)` = sprintf("%.2f (%.2f, %.2f)",
                           estimate, conf.low, conf.high),
    P_value = sprintf("%.3f", p.value)
  ) %>%
  select(Model, `β (95% CI)`, P_value)

# Add attenuation row
comparison_table <- comparison_table %>%
  bind_rows(
    tibble(
      Model = "% Attenuation after HEI",
      `β (95% CI)` = paste0(round(attenuation_pct,1), "%"),
      P_value = ""
    )
  )

write.xlsx(
  comparison_table,
  here("results", "Model2_vs_Model2plusHEI.xlsx"),
  rowNames = FALSE
)

#---------------------sensitivity, exclude implausible kcal------------------------

analysis_energy <- analysis_clean %>%
  filter(dr1tkcal > 800,
         dr1tkcal < 5000)

design_energy <- svydesign(
  ids = ~sdmvpsu,
  strata = ~sdmvstra,
  weights = ~wtmec6yr,
  data = analysis_energy,
  nest = TRUE
)

design_energy <- update(design_energy,
                        vk_density_sd = analysis_energy$vk_density_sd)

m2_energy <- svyglm(
  lbxtst ~ vk_density_sd +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session,
  design = design_energy
)

tab_m2_energy <- tidy(m2_energy, conf.int = TRUE) %>%
  filter(term == "vk_density_sd")

beta_energy <- tab_m2_energy$estimate
beta_original <- tab_m2$estimate

attenuation_energy <- 100 * (beta_original - beta_energy) / beta_original

energy_table <- tibble(
  Model = c("Model 2 (Fully adjusted)",
            "Model 2 – Excluding implausible energy",
            "Percent change"),
  `β (95% CI)` = c(
    sprintf("%.2f (%.2f, %.2f)",
            tab_m2$estimate,
            tab_m2$conf.low,
            tab_m2$conf.high),
    sprintf("%.2f (%.2f, %.2f)",
            tab_m2_energy$estimate,
            tab_m2_energy$conf.low,
            tab_m2_energy$conf.high),
    paste0(round(attenuation_energy,1), "%")
  ),
  P_value = c(
    sprintf("%.3f", tab_m2$p.value),
    sprintf("%.3f", tab_m2_energy$p.value),
    ""
  )
)

write.xlsx(
  energy_table,
  here("results", "Robustness_Excluding_Implausible_Energy.xlsx"),
  rowNames = FALSE
)

#---------------------sensitivity, exclude extreme vk------------------------

vk_p1  <- quantile(analysis_clean$vk_density, 0.01, na.rm = TRUE)
vk_p99 <- quantile(analysis_clean$vk_density, 0.99, na.rm = TRUE)

vk_p1
vk_p99

analysis_trim_vk <- analysis_clean %>%
  filter(
    vk_density >= vk_p1,
    vk_density <= vk_p99
  )

design_trim_vk <- svydesign(
  ids = ~sdmvpsu,
  strata = ~sdmvstra,
  weights = ~wtmec6yr,
  data = analysis_trim_vk,
  nest = TRUE
)

vk_sd_trim <- sqrt(
  as.numeric(
    svyvar(~vk_density, design_trim_vk)
  )
)

vk_mean_trim <- as.numeric(
  svymean(~vk_density, design_trim_vk, na.rm = TRUE)
)

analysis_trim_vk <- analysis_trim_vk %>%
  mutate(
    vk_density_sd_trim =
      (vk_density - vk_mean_trim) / vk_sd_trim
  )

design_trim_vk <- update(
  design_trim_vk,
  vk_density_sd_trim = analysis_trim_vk$vk_density_sd_trim
)

m2_trim_vk <- svyglm(
  lbxtst ~ vk_density_sd_trim +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session,
  design = design_trim_vk
)

# Extract trimmed model
tab_m2_trim_vk <- tidy(m2_trim_vk, conf.int = TRUE) %>%
  filter(term == "vk_density_sd_trim")

# Compute percent change vs original Model 2
beta_trim <- tab_m2_trim_vk$estimate
beta_original <- tab_m2$estimate

attenuation_vk <- 100 * (beta_original - beta_trim) / beta_original

# Build table (analogous to energy table)
vk_trim_table <- tibble(
  Model = c(
    "Model 2 (Fully adjusted)",
    "Model 2 – Excluding extreme VK (1–99%)",
    "Percent change"
  ),
  `β (95% CI)` = c(
    sprintf("%.2f (%.2f, %.2f)",
            tab_m2$estimate,
            tab_m2$conf.low,
            tab_m2$conf.high),
    sprintf("%.2f (%.2f, %.2f)",
            tab_m2_trim_vk$estimate,
            tab_m2_trim_vk$conf.low,
            tab_m2_trim_vk$conf.high),
    paste0(round(attenuation_vk, 1), "%")
  ),
  P_value = c(
    sprintf("%.3f", tab_m2$p.value),
    sprintf("%.3f", tab_m2_trim_vk$p.value),
    ""
  )
)

write.xlsx(
  vk_trim_table,
  here("results", "Robustness_Excluding_Extreme_VK_1to99.xlsx"),
  rowNames = FALSE
)

#-----------------------sensitivity, extreme T --------------------------

t_p1  <- quantile(analysis_clean$lbxtst, 0.01, na.rm = TRUE)
t_p99 <- quantile(analysis_clean$lbxtst, 0.99, na.rm = TRUE)

t_p1
t_p99

analysis_trim_t <- analysis_clean %>%
  filter(
    lbxtst >= t_p1,
    lbxtst <= t_p99
  )

design_trim_t <- svydesign(
  ids = ~sdmvpsu,
  strata = ~sdmvstra,
  weights = ~wtmec6yr,
  data = analysis_trim_t,
  nest = TRUE
)

vk_sd_trim_t <- sqrt(
  as.numeric(
    svyvar(~vk_density, design_trim_t)
  )
)

vk_mean_trim_t <- as.numeric(
  svymean(~vk_density, design_trim_t, na.rm = TRUE)
)

analysis_trim_t <- analysis_trim_t %>%
  mutate(
    vk_density_sd_trim_t =
      (vk_density - vk_mean_trim_t) / vk_sd_trim_t
  )

design_trim_t <- update(
  design_trim_t,
  vk_density_sd_trim_t =
    analysis_trim_t$vk_density_sd_trim_t
)

m2_trim_t <- svyglm(
  lbxtst ~ vk_density_sd_trim_t +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session,
  design = design_trim_t
)

tab_m2_trim_t <- tidy(m2_trim_t, conf.int = TRUE) %>%
  filter(term == "vk_density_sd_trim_t")

beta_trim_t <- tab_m2_trim_t$estimate
beta_original <- tab_m2$estimate

attenuation_t <- 100 * (beta_original - beta_trim_t) / beta_original

t_trim_table <- tibble(
  Model = c(
    "Model 2 (Fully adjusted)",
    "Model 2 – Excluding extreme testosterone (1–99%)",
    "Percent change"
  ),
  `β (95% CI)` = c(
    sprintf("%.2f (%.2f, %.2f)",
            tab_m2$estimate,
            tab_m2$conf.low,
            tab_m2$conf.high),
    sprintf("%.2f (%.2f, %.2f)",
            tab_m2_trim_t$estimate,
            tab_m2_trim_t$conf.low,
            tab_m2_trim_t$conf.high),
    paste0(round(attenuation_t, 1), "%")
  ),
  P_value = c(
    sprintf("%.3f", tab_m2$p.value),
    sprintf("%.3f", tab_m2_trim_t$p.value),
    ""
  )
)

write.xlsx(
  t_trim_table,
  here("results", "Robustness_Excluding_Extreme_Testosterone_1to99.xlsx"),
  rowNames = FALSE
)

#-----------------------sensitivity, SHBG --------------------------

analysis_shbg <- analysis_clean %>%
  filter(
    cycle %in% c("2013-2014","2015-2016"),
    !is.na(lbxshbg)
  )

design_shbg <- svydesign(
  ids = ~sdmvpsu,
  strata = ~sdmvstra,
  weights = ~wtmec4yr,
  data = analysis_shbg,
  nest = TRUE
)

m2_restricted <- svyglm(
  lbxtst ~ vk_density_sd +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session,
  design = design_shbg
)

tab_m2_restricted <- tidy(m2_restricted, conf.int = TRUE) %>%
  filter(term == "vk_density_sd")

m2_shbg <- svyglm(
  lbxtst ~ vk_density_sd +
    lbxshbg +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session,
  design = design_shbg
)

tab_m2_shbg <- tidy(m2_shbg, conf.int = TRUE) %>%
  filter(term == "vk_density_sd")

beta_restricted <- tab_m2_restricted$estimate
beta_shbg <- tab_m2_shbg$estimate

attenuation_shbg <- 100 * (beta_restricted - beta_shbg) / beta_restricted

shbg_table <- tibble(
  Model = c(
    "Model 2 (Restricted 2013–2016)",
    "Model 2 + SHBG",
    "Percent change"
  ),
  `β (95% CI)` = c(
    sprintf("%.2f (%.2f, %.2f)",
            tab_m2_restricted$estimate,
            tab_m2_restricted$conf.low,
            tab_m2_restricted$conf.high),
    sprintf("%.2f (%.2f, %.2f)",
            tab_m2_shbg$estimate,
            tab_m2_shbg$conf.low,
            tab_m2_shbg$conf.high),
    paste0(round(attenuation_shbg,1), "%")
  ),
  P_value = c(
    sprintf("%.3f", tab_m2_restricted$p.value),
    sprintf("%.3f", tab_m2_shbg$p.value),
    ""
  )
)

write.xlsx(
  shbg_table,
  here("results", "Robustness_Adjusted_for_SHBG_2013to2016.xlsx"),
  rowNames = FALSE
)

#-----------------------sensitivity, exclude P cancer --------------------------

analysis_no_pc <- analysis_clean %>%
  filter(mcq230a != 30 | is.na(mcq230a))

design_no_pc <- svydesign(
  ids = ~sdmvpsu,
  strata = ~sdmvstra,
  weights = ~wtmec6yr,
  data = analysis_no_pc,
  nest = TRUE
)

m2_no_pc <- svyglm(
  lbxtst ~ vk_density_sd +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session,
  design = design_no_pc
)

tab_m2_no_pc <- tidy(m2_no_pc, conf.int = TRUE) %>%
  filter(term == "vk_density_sd")

# Compute percent change like before
beta_no_pc <- tab_m2_no_pc$estimate
beta_original <- tab_m2$estimate

attenuation_no_pc <- 100 * (beta_original - beta_no_pc) / beta_original

pc_table <- tibble(
  Model = c(
    "Model 2 (Fully adjusted)",
    "Model 2 – Excluding prostate cancer",
    "Percent change"
  ),
  `β (95% CI)` = c(
    sprintf("%.2f (%.2f, %.2f)",
            tab_m2$estimate,
            tab_m2$conf.low,
            tab_m2$conf.high),
    sprintf("%.2f (%.2f, %.2f)",
            tab_m2_no_pc$estimate,
            tab_m2_no_pc$conf.low,
            tab_m2_no_pc$conf.high),
    paste0(round(attenuation_no_pc,1), "%")
  ),
  P_value = c(
    sprintf("%.3f", tab_m2$p.value),
    sprintf("%.3f", tab_m2_no_pc$p.value),
    ""
  )
)

write.xlsx(
  pc_table,
  here("results", "Robustness_Excluding_Prostate_Cancer_Specific.xlsx"),
  rowNames = FALSE
)

#-----------------------sensitivity, test nonlinearity--------------------------

m2_linear_raw <- svyglm(
  lbxtst ~ vk_density +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session,
  design = design_1116
)

m2_spline <- svyglm(
  lbxtst ~ ns(vk_density, df = 3) +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session,
  design = design_1116
)

extract_spline_test <- function(linear_model, spline_model) {
  
  overall_test <- regTermTest(
    spline_model,
    ~ ns(vk_density, df = 3)
  )
  
  nonlinear_test <- anova(
    linear_model,
    spline_model
  )
  
  tibble(
    N = nobs(spline_model),
    `P overall association` = signif(overall_test$p, 3),
    `P for non-linearity` = signif(nonlinear_test$p, 3)
  )
}

spline_results_table <- extract_spline_test(
  m2_linear_raw,
  m2_spline
)

write.xlsx(
  spline_results_table,
  file = here("results", "Sensitivity_RestrictedCubicSplines_Model2.xlsx"),
  overwrite = TRUE
)


#--------------------nonlinearity, generate predicted curve-----------------------

m2_spline <- svyglm(
  lbxtst ~ ns(vk_density, df = 3) +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session,
  design = design_1116
)

regTermTest(m2_spline, ~ ns(vk_density, df = 3))

library(marginaleffects)
library(ggplot2)

vk_grid <- datagrid(
  model = m2_spline,
  vk_density = seq(
    quantile(analysis_clean$vk_density, 0.01, na.rm = TRUE),
    quantile(analysis_clean$vk_density, 0.99, na.rm = TRUE),
    length.out = 200
  )
)

pred <- predictions(
  m2_spline,
  newdata = vk_grid
)

spline_plot <- ggplot(pred, aes(x = vk_density, y = estimate)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.2
  ) +
  labs(
    x = "Vitamin K density (µg/1000 kcal)",
    y = "Predicted serum testosterone (ng/dL)"
  )

ggsave(
  filename = here(
    "figures",
    "Spline_VitaminK_Testosterone.png"
  ),
  plot = spline_plot,
  width = 7,
  height = 5,
  dpi = 300
)

#-----------------------effect modification, obesity--------------------------
analysis_clean <- analysis_clean %>%
  mutate(
    obesity = if_else(bmxbmi >= 30, "Obese", "Non-obese"),
    obesity = factor(obesity,
                     levels = c("Non-obese", "Obese"))
  )

design_1116 <- update(design_1116,
                      obesity = analysis_clean$obesity)

m2_interaction_obesity <- svyglm(
  lbxtst ~ vk_density_sd * obesity +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session,
  design = design_1116
)

interaction_p_obesity <- tidy(m2_interaction_obesity) %>%
  filter(grepl("vk_density_sd:obesity", term)) %>%
  pull(p.value)

m_nonobese <- svyglm(
  lbxtst ~ vk_density_sd +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session,
  design = subset(design_1116, obesity == "Non-obese")
)

m_obese <- svyglm(
  lbxtst ~ vk_density_sd +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session,
  design = subset(design_1116, obesity == "Obese")
)

extract_strata <- function(model, group_label, interaction_p) {
  tidy(model, conf.int = TRUE) %>%
    filter(term == "vk_density_sd") %>%
    mutate(
      Group = group_label,
      N = nobs(model),
      `P for interaction` = interaction_p
    ) %>%
    select(Group, estimate, conf.low, conf.high, p.value, N, `P for interaction`)
}

em_obesity <- bind_rows(
  extract_strata(m_nonobese, "Non-obese", interaction_p_obesity),
  extract_strata(m_obese, "Obese", interaction_p_obesity)
)

#-----------------------effect modification, age--------------------------
analysis_clean <- analysis_clean %>%
  mutate(
    age_group = if_else(ridageyr < 50, "<50", "≥50"),
    age_group = factor(age_group,
                       levels = c("<50", "≥50"))
  )

design_1116 <- update(design_1116,
                      age_group = analysis_clean$age_group)

m2_interaction_age <- svyglm(
  lbxtst ~ vk_density_sd * age_group +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session,
  design = design_1116
)

interaction_p_age <- tidy(m2_interaction_age) %>%
  filter(grepl("vk_density_sd:age_group", term)) %>%
  pull(p.value)

m_age_under50 <- svyglm(
  lbxtst ~ vk_density_sd +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session,
  design = subset(design_1116, age_group == "<50")
)

m_age_50plus <- svyglm(
  lbxtst ~ vk_density_sd +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session,
  design = subset(design_1116, age_group == "≥50")
)

em_age <- bind_rows(
  extract_strata(m_age_under50, "<50", interaction_p_age),
  extract_strata(m_age_50plus, "≥50", interaction_p_age)
)


#-----------------------effect modification, HOMA-IR--------------------------
analysis_clean_ir <- analysis_clean %>%
  filter(!is.na(homa_ir))

design_ir <- svydesign(
  ids = ~sdmvpsu,
  strata = ~sdmvstra,
  weights = ~wtmec6yr,
  data = analysis_clean_ir,
  nest = TRUE
)

homa_cut <- as.numeric(
  coef(svyquantile(~homa_ir,
                   design_ir,
                   quantiles = 0.75,
                   na.rm = TRUE))
)

analysis_clean_ir <- analysis_clean_ir %>%
  mutate(
    ir_group = if_else(homa_ir >= homa_cut,
                       "Top quartile",
                       "Lower 3 quartiles"),
    ir_group = factor(ir_group,
                      levels = c("Lower 3 quartiles",
                                 "Top quartile"))
  )

design_ir <- update(design_ir,
                    ir_group = analysis_clean_ir$ir_group)

m2_interaction_ir <- svyglm(
  lbxtst ~ vk_density_sd * ir_group +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal,
  design = design_ir
)

interaction_p_ir <- tidy(m2_interaction_ir) %>%
  filter(grepl("vk_density_sd:ir_group", term)) %>%
  pull(p.value)

m_ir_low <- svyglm(
  lbxtst ~ vk_density_sd +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal,
  design = subset(design_ir,
                  ir_group == "Lower 3 quartiles")
)

m_ir_high <- svyglm(
  lbxtst ~ vk_density_sd +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal,
  design = subset(design_ir,
                  ir_group == "Top quartile")
)

em_ir <- bind_rows(
  extract_strata(m_ir_low,
                 "Lower 3 quartiles",
                 interaction_p_ir),
  extract_strata(m_ir_high,
                 "Top quartile",
                 interaction_p_ir)
)

effect_mod_all <- bind_rows(
  em_obesity %>% mutate(Modifier = "Obesity"),
  em_age %>% mutate(Modifier = "Age"),
  em_ir %>% mutate(Modifier = "Insulin resistance")
)

effect_mod_final <- effect_mod_all %>%
  mutate(
    `β (95% CI)` = sprintf(
      "%.2f (%.2f, %.2f)",
      estimate, conf.low, conf.high
    ),
    `P value` = sprintf("%.3f", p.value),
    `P for interaction` = sprintf("%.3f", `P for interaction`)
  ) %>%
  select(
    Modifier,
    Group,
    `β (95% CI)`,
    `P value`,
    N,
    `P for interaction`
  )

write.xlsx(
  effect_mod_final,
  here("results", "Effect_Modification_Model2_VK_Testosterone.xlsx"),
  overwrite = TRUE
)

#------------------effect modification, age, further investigation-----------------
m3_interaction_age <- svyglm(
  lbxtst ~ vk_density_sd * age_group +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session +
    bmxbmi +
    bmxwaist,
  design = design_1116
)

interaction_p_age_m3 <- tidy(m3_interaction_age) %>%
  filter(grepl("vk_density_sd:age_group", term)) %>%
  pull(p.value)

interaction_p_age_m3



m2_interaction_age_cont <- svyglm(
  lbxtst ~ vk_density * age_group +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session,
  design = design_1116
)

interaction_p_age_cont <- tidy(m2_interaction_age_cont) %>%
  filter(grepl("vk_density:age_group", term)) %>%
  pull(p.value)

interaction_p_age_cont



m2_lowT_interaction_age <- svyglm(
  lowt ~ vk_density_sd * age_group +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session,
  design = design_1116,
  family = quasibinomial()
)

interaction_p_lowT_age <- tidy(m2_lowT_interaction_age) %>%
  filter(grepl("vk_density_sd:age_group", term)) %>%
  pull(p.value)

interaction_p_lowT_age

analysis_trim_age <- analysis_trim_vk %>%
  mutate(
    age_group = if_else(ridageyr < 50, "<50", "≥50"),
    age_group = factor(age_group,
                       levels = c("<50", "≥50"))
  )

design_trim_age <- svydesign(
  ids = ~sdmvpsu,
  strata = ~sdmvstra,
  weights = ~wtmec6yr,
  data = analysis_trim_age,
  nest = TRUE
)

m2_trim_interaction_age <- svyglm(
  lbxtst ~ vk_density_sd_trim * age_group +
    ridageyr +
    race +
    indfmpir +
    smoking +
    alcohol +
    phys_active +
    dr1tkcal +
    exam_session,
  design = design_trim_age
)

interaction_p_trim_age <- tidy(m2_trim_interaction_age) %>%
  filter(grepl("vk_density_sd_trim:age_group", term)) %>%
  pull(p.value)

interaction_p_trim_age


#------------------sensitivity analyses N summary table-----------------

N_main_m2        <- nobs(m2)
N_energy_m2      <- nobs(m2_energy)
N_trim_vk_m2     <- nobs(m2_trim_vk)
N_trim_t_m2      <- nobs(m2_trim_t)
N_shbg_m2        <- nobs(m2_restricted)
N_shbg_adj_m2    <- nobs(m2_shbg)
N_no_pc_m2       <- nobs(m2_no_pc)
N_hei_m2         <- nobs(m2_hei)

N_lowT_m2 <- nobs(m2_lowT)

N_interaction_age <- nobs(m2_interaction_age)
N_interaction_obesity <- nobs(m2_interaction_obesity)
N_interaction_ir <- nobs(m2_interaction_ir)

sensitivity_N_table <- tibble(
  Specification = c(
    "Model 2 (Fully adjusted)",
    "Model 2 + HEI",
    "Excluding implausible energy",
    "Excluding extreme VK (1–99%)",
    "Excluding extreme testosterone (1–99%)",
    "Restricted 2013–2016",
    "Restricted + SHBG adjustment",
    "Excluding prostate cancer"
  ),
  Analyzed_N = c(
    N_main_m2,
    N_hei_m2,
    N_energy_m2,
    N_trim_vk_m2,
    N_trim_t_m2,
    N_shbg_m2,
    N_shbg_adj_m2,
    N_no_pc_m2
  )
)

write.xlsx(
  sensitivity_N_table,
  here("results", "Supplementary_Sensitivity_Unweighted_N.xlsx"),
  overwrite = TRUE
)


svyvar(~vk_density, design_1116)
sqrt(as.numeric(svyvar(~vk_density, design_1116)))

svyvar(~vk_density, design_trim_vk)
sqrt(as.numeric(svyvar(~vk_density, design_trim_vk)))


writeLines(
  capture.output(sessionInfo()),
  here("results", "sessionInfo.txt")
)
