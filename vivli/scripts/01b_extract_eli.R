
# Load packages and functions
source("Survival_attrition/Scripts/00a_packages_ipd_functions.R")

### 6484 NCT00051558 -----------------------------------------------------------

ipdname <- "V:/From 8697/6484/Analysis Datasets_GHBZ/"
eli_6484 <- list.files("V:/From 8697/6484/Analysis Datasets_GHBZ")

eli_6484_attr <- read_sas(paste0(ipdname, eli_6484[14])) %>%
  rename_with(., tolower) %>%
  filter(!is.na(randsdte) & !is.na(thstrsdt)) %>%
  select(patient, protocol, therapy, randsdte) %>%
  inner_join(
    read_sas(paste0(ipdname, eli_6484[18])) %>%
      rename_with(., tolower) %>%
      select(patient, protocol, thstpsdt, endsdt, reasont)
  ) %>%
  mutate(
    trtgrp = if_else(therapy == "ALN10", "Alendronate", "Teriparatide"),
    trttype = if_else(grepl("PTH", trtgrp), "Exp", "AC"),
    randt_unit = "date",
    disct_unit = "date",
    reason = if_else(reasont == "Protocol completed", "Completed", reasont),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = protocol, usubjid = patient, 
    trtgrp, trttype, 
    randt = randsdte,  randt_unit,
    disct = endsdt, disct_unit,
    compl,
    reason_orig = reasont,
    reason
  )

table(eli_6484_attr$reason, useNA = "ifany")
table(eli_6484_attr$trtgrp, eli_6484_attr$reason, useNA = "ifany")

### 9797 NCT00384930 -----------------------------------------------------------

ipdname <- "V:/From 8697/9797/Analysis Datasets_LVHG/acute/"
eli_9797 <- list.files("V:/From 8697/9797/Analysis Datasets_LVHG/acute")

eli_9797_attr <- read_sas(paste0(ipdname, eli_9797[17])) %>%
  rename_with(., tolower) %>%
  filter(subjitt == 1) %>%
  select(usubjid, sdyid, trt, rndmdt) %>%
  inner_join(
    read_sas(paste0(ipdname, eli_9797[4])) %>%
      rename_with(., tolower) %>%
      select(usubjid, dsdt, dssnm)
  ) %>%
  mutate(
    trttype = if_else(grepl("Tada", trt), "Exp", "AC"),
    randt_unit = "date",
    disct_unit = "date",
    compl = if_else(dssnm == "Completed", 1, 0)
  ) %>%
  select(
    studyid = sdyid, usubjid, 
    trtgrp = trt, trttype, 
    randt = rndmdt,  randt_unit,
    disct = dsdt, disct_unit,
    compl,
    reason_orig = dssnm,
    reason = dssnm
  )

table(eli_9797_attr$reason, useNA = "ifany")
table(eli_9797_attr$trtgrp, eli_9797_attr$reason, useNA = "ifany")

### 10303 NCT00125918 ----------------------------------------------------------

list.files("V:/From 8697/10303")
ipdname <- c("V:/From 8697/10303/Analysis Datasets_LVGY.zip")
eli_10303 <- ListZip(ipdname, separate = FALSE)

eli_10303_attr <- read_sas(unz(ipdname, eli_10303$file[17])) %>%
  rename_with(., tolower) %>%
  filter(!is.na(rndmdt)) %>%
  select(usubjid, sdyid, trt, rndmdt) %>%
  inner_join(
    read_sas(unz(ipdname, "disposit.sas7bdat")) %>%
      rename_with(., tolower) %>%
      filter(subjitt == 1) %>%
      select(usubjid, sdyid, discdt, dslnm)
  ) %>%
  mutate(
    trttype = if_else(trt == "Placebo", "Plc", "Exp"),
    randt_unit = "date",
    disct_unit = "date",
    compl = if_else(dslnm == "Completed", 1, 0)
  ) %>%
  select(
    studyid = sdyid, usubjid, 
    trtgrp = trt, trttype, 
    randt = rndmdt,  randt_unit,
    disct = discdt, disct_unit,
    compl,
    reason_orig = dslnm,
    reason = dslnm
  )

table(eli_10303_attr$reason, useNA = "ifany")
table(eli_10303_attr$trtgrp, eli_10303_attr$reason, useNA = "ifany")

### 10487 NCT00861757 ----------------------------------------------------------

list.files("V:/From 8697/10487")
ipdname <- c("V:/From 8697/10487/Analysis Datasets_LVHB.zip")
eli_10487 <- ListZip(ipdname, separate = FALSE)

eli_10487_attr <- read_sas(unz(ipdname, eli_10487$file[20])) %>%
  rename_with(., tolower) %>%
  filter(!is.na(rndmdt) & subjfastr == 1) %>%
  select(usubjid, sdyid, trt, rndmdt) %>%
  inner_join(
    read_sas(unz(ipdname, eli_10487$file[5])) %>%
      rename_with(., tolower) %>%
      select(usubjid, sdyid, dsdt, dslnm)
  ) %>%
  mutate(
    trttype = 
      case_when(
        grepl("placebo", trt) ~ "Plc",
        grepl("tada", trt) ~ "Exp",
        TRUE ~ "AC"
      ),
    randt_unit = "date",
    disct_unit = "date",
    compl = if_else(dslnm == "Completed", 1, 0)
  ) %>%
  select(
    studyid = sdyid, usubjid, 
    trtgrp = trt, trttype, 
    randt = rndmdt,  randt_unit,
    disct = dsdt, disct_unit,
    compl,
    reason_orig = dslnm,
    reason = dslnm
  )

table(eli_10487_attr$reason, useNA = "ifany")
table(eli_10487_attr$trtgrp, eli_10487_attr$reason, useNA = "ifany")


### 10893 NCT00827242 ----------------------------------------------------------

list.files("V:/From 8697/10893")
ipdname <- c("V:/From 8697/10893/Analysis Datasets_LVHJ.zip")
eli_10893 <- ListZip(ipdname, separate = FALSE)

eli_10893_attr <- read_sas(unz(ipdname, eli_10893$file[23])) %>%
  rename_with(., tolower) %>%
  filter(!is.na(rndmdt)) %>%
  select(usubjid, sdyid, trt, rndmdt) %>%
  inner_join(
    read_sas(unz(ipdname, eli_10893$file[6])) %>%
      rename_with(., tolower) %>%
      select(usubjid, sdyid, dsdt, dslnm)
  ) %>%
  mutate(
    trttype = if_else(trt == "Placebo", "Plc", "Exp"),
    randt_unit = "date",
    disct_unit = "date",
    compl = if_else(dslnm == "Completed", 1, 0)
  ) %>% 
  select(
    studyid = sdyid, usubjid, 
    trtgrp = trt, trttype, 
    randt = rndmdt,  randt_unit,
    disct = dsdt, disct_unit,
    compl,
    reason_orig = dslnm,
    reason = dslnm
  )

table(eli_10893_attr$reason, useNA = "ifany")
table(eli_10893_attr$trtgrp, eli_10893_attr$reason, useNA = "ifany")

### 11373 NCT01064687 ----------------------------------------------------------

ipdname <- "V:/From 8697/11373/Analysis Datasets_GBDA/"
eli_11373 <- list.files("V:/From 8697/11373/Analysis Datasets_GBDA")

eli_11373_attr <- read_sas(paste0(ipdname, eli_11373[46])) %>%
  rename_with(., tolower) %>%
  filter(!is.na(rndmdt)  & subjitttr == 1) %>%
  select(usubjid, sdyid, trt, rndmdt) %>%
  left_join(
    read_sas(paste0(ipdname, eli_11373[15])) %>%
      rename_with(., tolower) %>% 
      select(usubjid, sdyid, dsdt, dslnm, visid)
  ) %>% 
  mutate(
    trttype = case_when(
      grepl("Placebo", trt) ~ "Plc",
      trt == "Exenatide" ~ "AC",
      TRUE ~ "Exp"
    ),
    randt_unit = "date",
    disct_unit = "date",
    compl = if_else(dslnm == "Completed", 1, 0)
  ) %>% 
  select(
    studyid = sdyid, usubjid, 
    trtgrp = trt, trttype, 
    randt = rndmdt,  randt_unit,
    disct = dsdt, disct_unit,
    compl,
    reason_orig = dslnm,
    reason = dslnm
  )

table(eli_11373_attr$reason, useNA = "ifany")
table(eli_11373_attr$trtgrp, eli_11373_attr$reason, useNA = "ifany")

### 11374 NCT01075282 ----------------------------------------------------------

ipdname <- "V:/From 8697/11374/Analysis Datasets_GBDB/"
eli_11374 <- list.files("V:/From 8697/11374/Analysis Datasets_GBDB")

eli_11374_attr <- read_sas(paste0(ipdname, eli_11374[45])) %>%
  rename_with(., tolower) %>%
  filter(!is.na(rndmdt)  & subjitttr == 1) %>%
  select(usubjid, sdyid, trt, rndmdt) %>%
  inner_join(
    read_sas(paste0(ipdname, eli_11374[17])) %>%
      rename_with(., tolower) %>% 
      select(usubjid, sdyid, dsdt, dslnm)
  ) %>% 
  mutate(
    trttype = 
      case_when(
        grepl("Dula", trt) ~ "Exp",
        TRUE ~ "AC"
      ),
    randt_unit = "date",
    disct_unit = "date",
    compl = if_else(dslnm == "Completed", 1, 0)
  ) %>% 
  select(
    studyid = sdyid, usubjid, 
    trtgrp = trt, trttype, 
    randt = rndmdt,  randt_unit,
    disct = dsdt, disct_unit,
    compl,
    reason_orig = dslnm,
    reason = dslnm
  )

table(eli_11374_attr$reason, useNA = "ifany")
table(eli_11374_attr$trtgrp, eli_11374_attr$reason, useNA = "ifany")

### 11375 NCT01126580 ----------------------------------------------------------

ipdname <- "V:/From 8697/11375/Analysis Datasets_GBDC/"
eli_11375 <- list.files("V:/From 8697/11375/Analysis Datasets_GBDC")

eli_11375_attr <- read_sas(paste0(ipdname, eli_11375[45])) %>%
  rename_with(., tolower) %>%
  filter(!is.na(rndmdt)  & subjitttr == 1) %>%
  select(usubjid, sdyid, trt, rndmdt) %>%
  inner_join(
    read_sas(paste0(ipdname, eli_11375[15])) %>%
      rename_with(., tolower) %>% 
      select(usubjid, sdyid, dsdt, dslnm)
  ) %>% 
  mutate(
    trttype = 
      case_when(
        grepl("Dula", trt) ~ "Exp",
        TRUE ~ "AC"
      ),
    randt_unit = "date",
    disct_unit = "date",
    compl = if_else(dslnm == "Completed", 1, 0)
  ) %>% 
  select(
    studyid = sdyid, usubjid, 
    trtgrp = trt, trttype, 
    randt = rndmdt,  randt_unit,
    disct = dsdt, disct_unit,
    compl,
    reason_orig = dslnm,
    reason = dslnm
  )

table(eli_11375_attr$reason, useNA = "ifany")
table(eli_11375_attr$trtgrp, eli_11375_attr$reason, useNA = "ifany")

### 11376 NCT01191268 ----------------------------------------------------------

ipdname <- "V:/From 8697/11376/Analysis Datasets_GBDD/"
eli_11376 <- list.files("V:/From 8697/11376/Analysis Datasets_GBDD")

eli_11376_attr <- read_sas(paste0(ipdname, eli_11376[45])) %>%
  rename_with(., tolower) %>%
  filter(!is.na(rndmdt)  & subjitt == 1) %>%
  select(usubjid, sdyid, trt, rndmdt) %>%
  left_join(
    read_sas(paste0(ipdname, eli_11376[19])) %>%
      rename_with(., tolower) %>% 
      select(usubjid, sdyid, dsdt, dslnm)
  ) %>% 
  mutate(
    trttype = 
      case_when(
        grepl("Dula", trt) ~ "Exp",
        TRUE ~ "AC"
      ),
    randt_unit = "date",
    disct_unit = "date",
    compl = if_else(dslnm == "Completed", 1, 0)
  ) %>% 
  filter(dslnm != "") %>% 
  select(
    studyid = sdyid, usubjid, 
    trtgrp = trt, trttype, 
    randt = rndmdt,  randt_unit,
    disct = dsdt, disct_unit,
    compl,
    reason_orig = dslnm,
    reason = dslnm
  )

table(eli_11376_attr$reason, useNA = "ifany")
table(eli_11376_attr$trtgrp, eli_11376_attr$reason, useNA = "ifany")

### 11377 NCT01624259 ----------------------------------------------------------

ipdname <- "V:/From 8697/11377/Analysis Datasets_GBDE/"
eli_11377 <- list.files("V:/From 8697/11377/Analysis Datasets_GBDE")

eli_11377_attr <- read_sas(paste0(ipdname, eli_11377[47])) %>%
  rename_with(., tolower) %>%
  filter(!is.na(rndmdt) & subjitttr == 1) %>%
  select(usubjid, sdyid, trt, rndmdt) %>%
  left_join(
    read_sas(paste0(ipdname, eli_11377[17])) %>%
      rename_with(., tolower) %>% 
      filter(dstp == "Final") %>% 
      select(usubjid, sdyid, dsdt, dslnm)
  ) %>% 
  mutate(
    trttype = 
      case_when(
        grepl("Dula", trt) ~ "Exp",
        TRUE ~ "AC"
      ),
    randt_unit = "date",
    disct_unit = "date",
    compl = if_else(dslnm == "Completed", 1, 0)
  ) %>% 
  filter(dslnm != "") %>% 
  select(
    studyid = sdyid, usubjid, 
    trtgrp = trt, trttype, 
    randt = rndmdt,  randt_unit,
    disct = dsdt, disct_unit,
    compl,
    reason_orig = dslnm,
    reason = dslnm
  )

table(eli_11377_attr$reason, useNA = "ifany")
table(eli_11377_attr$trtgrp, eli_11377_attr$reason, useNA = "ifany")

### 11422 NCT00734474 ----------------------------------------------------------

ipdname <- "V:/From 8697/11422/Analysis Datasets_GBCF/"
eli_11422 <- list.files("V:/From 8697/11422/Analysis Datasets_GBCF")

eli_11422_attr <- read_sas(paste0(ipdname, eli_11422[41])) %>%
  rename_with(., tolower) %>%
  filter(!is.na(rndmdt) & subjitttr == 1) %>%
  select(usubjid, sdyid, trt, rndmdt) %>%
  left_join(
    read_sas(paste0(ipdname, eli_11422[12])) %>%
      rename_with(., tolower) %>% 
      filter(dstp == "Final") %>% 
      select(usubjid, sdyid, dsdt, dslnm)
  ) %>% 
  mutate(
    trttype = 
      case_when(
        grepl("Dula", trt) ~ "Exp",
        TRUE ~ "AC"
      ),
    randt_unit = "date",
    disct_unit = "date",
    compl = if_else(dslnm == "Completed", 1, 0)
  ) %>% 
  filter(
    dslnm != ""
    & trt %in% c("LY 0.75mg", "LY 1.5mg", "Placebo/Sitagliptin", "Sitagliptin")
  ) %>% 
  select(
    studyid = sdyid, usubjid, 
    trtgrp = trt, trttype, 
    randt = rndmdt,  randt_unit,
    disct = dsdt, disct_unit,
    compl,
    reason_orig = dslnm,
    reason = dslnm
  )

table(eli_11422_attr$reason, useNA = "ifany")
table(eli_11422_attr$trtgrp, eli_11422_attr$reason, useNA = "ifany")

### 11667 NCT00855582 ----------------------------------------------------------

list.files("V:/From 8697/11667")
ipdname <- c("V:/From 8697/11667/Analysis Datasets_LVHR.zip")
eli_11667 <- ListZip(ipdname, separate = FALSE)

eli_11667_attr <- read_sas(unz(ipdname, eli_11667$file[27])) %>%
  rename_with(., tolower) %>%
  filter(!is.na(rndmdt)) %>%
  select(usubjid, sdyid, trt, rndmdt) %>%
  inner_join(
    read_sas(unz(ipdname, eli_11667$file[6])) %>%
      rename_with(., tolower) %>%
      select(usubjid, sdyid, dsdt, dslnm)
  ) %>%
  mutate(
    trttype = 
      case_when(
        grepl("Tada", trt) ~ "Exp",
        TRUE ~ "Plc"
      ),
    randt_unit = "date",
    disct_unit = "date",
    compl = if_else(dslnm == "Completed", 1, 0)
  ) %>%
  select(
    studyid = sdyid, usubjid, 
    trtgrp = trt, trttype, 
    randt = rndmdt,  randt_unit,
    disct = dsdt, disct_unit,
    compl,
    reason_orig = dslnm,
    reason = dslnm
  )

table(eli_11667_attr$reason, useNA = "ifany")
table(eli_11667_attr$trtgrp, eli_11667_attr$reason, useNA = "ifany")


### 11668 NCT00848081 ----------------------------------------------------------

list.files("V:/From 8697/11668")
ipdname <- c("V:/From 8697/11668/Analysis Datasets_LVHS.zip")
eli_11668 <- ListZip(ipdname, separate = FALSE)

eli_11668_attr <- read_sas(unz(ipdname, eli_11668$file[20])) %>%
  rename_with(., tolower) %>%
  filter(!is.na(rndmdt)) %>%
  select(usubjid, sdyid, trt, rndmdt) %>%
  inner_join(
    read_sas(unz(ipdname, eli_11668$file[4])) %>%
      rename_with(., tolower) %>%
      select(usubjid, sdyid, dsdt, dslnm)
  ) %>%
  mutate(
    trttype = 
      case_when(
        grepl("Tada", trt) ~ "Exp",
        TRUE ~ "Plc"
      ),
    randt_unit = "date",
    disct_unit = "date",
    compl = if_else(dslnm == "Completed", 1, 0)
  ) %>%
  select(
    studyid = sdyid, usubjid, 
    trtgrp = trt, trttype, 
    randt = rndmdt,  randt_unit,
    disct = dsdt, disct_unit,
    compl,
    reason_orig = dslnm,
    reason = dslnm
  )

table(eli_11668_attr$reason, useNA = "ifany")
table(eli_11668_attr$trtgrp, eli_11668_attr$reason, useNA = "ifany")

### 12932 NCT00970632 ----------------------------------------------------------

list.files("V:/From 8697/12932")
ipdname <- c("V:/From 8697/12932/Analysis Datasets_LVID.zip")
eli_12932 <- ListZip(ipdname, separate = FALSE)

eli_12932_attr <- read_sas(unz(ipdname, eli_12932$file[23])) %>%
  rename_with(., tolower) %>%
  filter(!is.na(rndmdt)) %>%
  select(usubjid, sdyid, trt, rndmdt) %>%
  inner_join(
    read_sas(unz(ipdname, eli_12932$file[6])) %>%
      rename_with(., tolower) %>%
      select(usubjid, sdyid, dsdt, dslnm)
  ) %>%
  mutate(
    trttype = 
      case_when(
        grepl("Tada", trt) ~ "Exp",
        grepl("Tamsu", trt) ~ "AC",
        TRUE ~ "Plc"
      ),
    randt_unit = "date",
    disct_unit = "date",
    compl = if_else(dslnm == "Completed", 1, 0)
  ) %>%
  select(
    studyid = sdyid, usubjid, 
    trtgrp = trt, trttype, 
    randt = rndmdt,  randt_unit,
    disct = dsdt, disct_unit,
    compl,
    reason_orig = dslnm,
    reason = dslnm
  )

table(eli_12932_attr$reason, useNA = "ifany")
table(eli_12932_attr$trtgrp, eli_12932_attr$reason, useNA = "ifany")

### 13193 NCT01769378 ----------------------------------------------------------

list.files("V:/From 8697/13193")
ipdname <- c("V:/From 8697/13193/Analysis Datasets_GBDG.zip")
eli_13193 <- ListZip(ipdname, separate = FALSE)

eli_13193_attr <- read_sas(unz(ipdname, eli_13193$file[44])) %>%
  rename_with(., tolower) %>%
  filter(!is.na(rndmdt) & subjitttr == 1) %>%
  select(usubjid, sdyid, trt, rndmdt) %>%
  left_join(
    read_sas(unz(ipdname, eli_13193$file[16])) %>%
      rename_with(., tolower) %>%
      select(usubjid, sdyid, dsdt, dslnm)
  ) %>%
  mutate(
    trttype = 
      case_when(
        grepl("Tada", trt) ~ "Exp",
        grepl("Tamsu", trt) ~ "AC",
        TRUE ~ "Plc"
      ),
    randt_unit = "date",
    disct_unit = "date",
    compl = if_else(dslnm == "Completed", 1, 0)
  ) %>%
  filter(dslnm != "") %>% 
  select(
    studyid = sdyid, usubjid, 
    trtgrp = trt, trttype, 
    randt = rndmdt,  randt_unit,
    disct = dsdt, disct_unit,
    compl,
    reason_orig = dslnm,
    reason = dslnm
  )

table(eli_13193_attr$reason, useNA = "ifany")
table(eli_13193_attr$trtgrp, eli_13193_attr$reason, useNA = "ifany")

### 13195 NCT02152371 ----------------------------------------------------------

list.files("V:/13195")
ipdname <- c("V:/13195/ADaM_GBDI.zip")
eli_13195 <- ListZip(ipdname, separate = FALSE)

eli_13195_attr <- read_sas(unz(ipdname, eli_13195$file[22])) %>%
  rename_with(., tolower) %>%
  filter(randfl == "Y" & ittfl == "Y") %>%
  select(usubjid, studyid, trt01a, randdt) %>%
  left_join(
    read_sas(unz(ipdname, eli_13195$file[9])) %>%
      rename_with(., tolower) %>%
      filter(dsscat == "STUDY DISPOSITION") %>% 
      select(usubjid, studyid, adt, dsdecod)
  ) %>%
  distinct() %>% 
  mutate(
    trttype = if_else(grepl("Dula", trt01a), "Exp", "Plc"),
    randt_unit = "date",
    disct_unit = "date",
    compl = if_else(dsdecod == "COMPLETED", 1, 0)
  ) %>%
  select(
    studyid, usubjid, 
    trtgrp = trt01a, trttype, 
    randt = randdt,  randt_unit,
    disct = adt, disct_unit,
    compl,
    reason_orig = dsdecod,
    reason = dsdecod
  )

table(eli_13195_attr$reason, useNA = "ifany")
table(eli_13195_attr$trtgrp, eli_13195_attr$reason, useNA = "ifany")

### 13798 NCT01621178 ----------------------------------------------------------

ipdname <- "V:/From 8697/13798/ADaM_GBDX/"
eli_13798 <- list.files("V:/From 8697/13798/ADaM_GBDX")

eli_13798_attr <- read_sas(paste0(ipdname, eli_13798[21])) %>%
  rename_with(., tolower) %>%
  filter(randfl == "Y") %>%
  select(usubjid, studyid, trt01a, randdt) %>%
  filter(trt01a != "") %>% 
  left_join(
    read_sas(paste0(ipdname, eli_13798[10])) %>%
      rename_with(., tolower) %>%
      filter(dsscat == "STUDY DISPOSITION") %>%
      select(usubjid, studyid, adt, dsdecod)
  ) %>% 
  distinct() %>% 
  mutate(
    trttype = if_else(grepl("Dula", trt01a), "Exp", "AC"),
    randt_unit = "date",
    disct_unit = "date",
    compl = if_else(dsdecod == "COMPLETED", 1, 0)
  ) %>% 
  select(
    studyid, usubjid, 
    trtgrp = trt01a, trttype, 
    randt = randdt,  randt_unit,
    disct = adt, disct_unit,
    compl,
    reason_orig = dsdecod,
    reason = dsdecod
  )

table(eli_13798_attr$reason, useNA = "ifany")
table(eli_13798_attr$trtgrp, eli_13798_attr$reason, useNA = "ifany")

### 13990 NCT01558271 ----------------------------------------------------------

ipdname <- "V:/From 8697/13990/Analysis Datasets_GBDP/"
eli_13990 <- list.files("V:/From 8697/13990/Analysis Datasets_GBDP")

eli_13990_attr <- read_sas(paste0(ipdname, eli_13990[33])) %>%
  rename_with(., tolower) %>%
  filter(!is.na(rndmdt)) %>%
  select(usubjid, sdyid, trt, rndmdt) %>%
  left_join(
    read_sas(paste0(ipdname, eli_13990[11])) %>%
      rename_with(., tolower) %>% 
      filter(dstp == "Study") %>% 
      select(usubjid, sdyid, dsdt, dslnm)
  ) %>% 
  mutate(
    trttype = 
      case_when(
        grepl("75", trt) ~ "Exp",
        grepl("Lira", trt) ~ "AC",
        TRUE ~ "Placebo"
      ),
    randt_unit = "date",
    disct_unit = "date",
    compl = if_else(dslnm == "Completed", 1, 0)
  ) %>% 
  select(
    studyid = sdyid, usubjid, 
    trtgrp = trt, trttype, 
    randt = rndmdt,  randt_unit,
    disct = dsdt, disct_unit,
    compl,
    reason_orig = dslnm,
    reason = dslnm
  )

table(eli_13990_attr$reason, useNA = "ifany")
table(eli_13990_attr$trtgrp, eli_13990_attr$reason, useNA = "ifany")

### 15361 NCT02597049 ----------------------------------------------------------

ipdname <- "V:/From 8697/15361/ADaM_GBGE/"
eli_15361 <- list.files("V:/From 8697/15361/ADaM_GBGE")

eli_15361_attr <- read_sas(paste0(ipdname, eli_15361[17])) %>%
  rename_with(., tolower) %>%
  filter(randfl == "Y") %>%
  select(usubjid, studyid, trt01a, randdt) %>%
  filter(trt01a != "") %>% 
  left_join(
    read_sas(paste0(ipdname, eli_15361[7])) %>%
      rename_with(., tolower) %>%
      filter(grepl("Study Disposition", param)) %>%
      select(usubjid, studyid, adt, dsdecod)
  ) %>% 
  distinct() %>% 
  mutate(
    trttype = if_else(grepl("Dula", trt01a), "Exp", "Plc"),
    randt_unit = "date",
    disct_unit = "date",
    compl = if_else(dsdecod == "COMPLETED", 1, 0)
  ) %>% 
  select(
    studyid, usubjid, 
    trtgrp = trt01a, trttype, 
    randt = randdt,  randt_unit,
    disct = adt, disct_unit,
    compl,
    reason_orig = dsdecod,
    reason = dsdecod
  )

table(eli_15361_attr$reason, useNA = "ifany")
table(eli_15361_attr$trtgrp, eli_15361_attr$reason, useNA = "ifany")

### 15803 NCT02750410 ----------------------------------------------------------

ipdname <- "V:/From 8697/15803/ADaM_GBGF/"
eli_15803 <- list.files("V:/From 8697/15803/ADaM_GBGF")

eli_15803_attr <- read_sas(paste0(ipdname, eli_15803[18])) %>%
  rename_with(., tolower) %>%
  filter(randfl == "Y") %>%
  select(usubjid, studyid, trt01a, randdt) %>%
  filter(trt01a != "") %>% 
  left_join(
    read_sas(paste0(ipdname, eli_15803[7])) %>%
      rename_with(., tolower) %>%
      filter(
        dsscat == "STUDY DISPOSITION"
        & epoch != "SCREENING"
        & !is.na(aperiod)
      ) %>%
      select(usubjid, studyid, adt, dsdecod)
  ) %>% 
  distinct() %>% 
  mutate(
    trttype = if_else(grepl("Dula", trt01a), "Exp", "Plc"),
    randt_unit = "date",
    disct_unit = "date",
    compl = if_else(dsdecod == "COMPLETED", 1, 0)
  ) %>% 
  select(
    studyid, usubjid, 
    trtgrp = trt01a, trttype, 
    randt = randdt,  randt_unit,
    disct = adt, disct_unit,
    compl,
    reason_orig = dsdecod,
    reason = dsdecod
  )

table(eli_15803_attr$reason, useNA = "ifany")
table(eli_15803_attr$trtgrp, eli_15803_attr$reason, useNA = "ifany")

# Save output ------------------------------------------------------------------

# Get dataframe names and add to list
get_df_names <- Filter(function(x) grepl("_attr$", x), ls())
eli_lst <- lapply(get_df_names, function(i) get(i))

# Format all variables in each dataframe as characters
eli_lst_format <- lapply(eli_lst, function(i) i %>% mutate(across(everything(), as.character)))

# Bind together
eli_ready <- bind_rows(eli_lst_format)

# Final checks
length(unique(eli_ready$studyid)) # 86
sum(duplicated(eli_ready$usubjid)) # 566 (same IDs used between two studies)
glimpse(eli_ready)

# Save
write_rds(eli_ready, "Survival_attrition/Data/eli.rds")
