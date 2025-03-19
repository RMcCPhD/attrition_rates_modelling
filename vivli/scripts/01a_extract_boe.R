
# Load functions
source("Survival_attrition/Scripts/00a_packages_ipd_functions.R")

### 107.210 NCT02183064 --------------------------------------------------------

list.files("V:/From 8697/107.210")
ipdname <- c("V:/From 8697/107.210/107.210_VIVLI_v01.zip")
boe_107.210 <- ListZip(ipdname, separate = FALSE)

boe_107.210_attr <- read_sas(unz(ipdname, boe_107.210$file[36])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_107.210$file[25])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_107.210$file[21])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcd) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "A", "Exp", "AC"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcd == 1 ~ "Completed",
      termcd %in% c(2, 3, 4) ~ "Adverse Event",
      termcd == 5 ~ "Lack of Efficacy",
      termcd == 6 ~ "Protocol Deviation",
      termcd == 7 ~ "Lost to Follow-up",
      termcd == 8 ~ "Voluntary Withdrawal",
      termcd == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcd,
    reason
  )

table(boe_107.210_attr$reason, useNA = "ifany")
table(boe_107.210_attr$trtgrp, boe_107.210_attr$reason, useNA = "ifany")

### 205.130 NCT02172287 --------------------------------------------------------

list.files("V:/From 8697/205.130")
ipdname <- c("V:/From 8697/205.130/205.130_VIVLI_v01.zip")
boe_205.130 <- ListZip(ipdname, separate = FALSE)

boe_205.130_attr <- read_sas(unz(ipdname, boe_205.130$file[41])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.130$file[24])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.130$file[19])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdc) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A" ~ "Exp",
      tpatt == "B" ~ "AC",
      TRUE ~ "Plc"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdc == 1 ~ "Completed",
      termcdc %in% c(2, 3, 4) ~ "Adverse Event",
      termcdc == 5 ~ "Protocol Deviation",
      termcdc == 6 ~ "Lost to Follow-up",
      termcdc == 7 ~ "Voluntary Withdrawal",
      termcdc == 8 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdc,
    reason
  )

table(boe_205.130_attr$reason, useNA = "ifany")
table(boe_205.130_attr$trtgrp, boe_205.130_attr$reason, useNA = "ifany")

### 205.214 NCT00274014 --------------------------------------------------------

list.files("V:/From 8697/205.214")
ipdname <- c("V:/From 8697/205.214/205.214_VIVLI_v01.zip")
boe_205.214 <- ListZip(ipdname, separate = FALSE)

boe_205.214_attr <- read_sas(unz(ipdname, boe_205.214$file[37])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.214$file[28])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.214$file[42])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdc) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A" ~ "Exp",
      TRUE ~ "Plc"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdc == 1 ~ "Completed",
      termcdc %in% c(2, 3, 4) ~ "Adverse Event",
      termcdc == 5 ~ "Protocol Deviation",
      termcdc == 6 ~ "Lost to Follow-up",
      termcdc == 7 ~ "Voluntary Withdrawal",
      termcdc == 8 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdc,
    reason
  )

table(boe_205.214_attr$reason, useNA = "ifany")
table(boe_205.214_attr$trtgrp, boe_205.214_attr$reason, useNA = "ifany")

### 205.254 NCT00168844---------------------------------------------------------

list.files("V:/From 8697/205.254")
ipdname <- c("V:/From 8697/205.254/205.254_VIVLI_v01.zip")
boe_205.254 <- ListZip(ipdname, separate = FALSE)

boe_205.254_attr <- read_sas(unz(ipdname, boe_205.254$file[76])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.254$file[54])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.254$file[82])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdc) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      grepl("D|H", tpatt) ~ "Exp",
      TRUE ~ "Plc"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdc == 1 ~ "Completed",
      termcdc %in% c(2, 3, 4) ~ "Adverse Event",
      termcdc == 5 ~ "Protocol Deviation",
      termcdc == 6 ~ "Lost to Follow-up",
      termcdc == 7 ~ "Voluntary Withdrawal",
      termcdc == 8 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdc,
    reason
  )

table(boe_205.254_attr$reason, useNA = "ifany")
table(boe_205.254_attr$trtgrp, boe_205.254_attr$reason, useNA = "ifany")

### 205.255 NCT00168831 --------------------------------------------------------

list.files("V:/From 8697/205.255")
ipdname <- c("V:/From 8697/205.255/205.255_VIVLI_v01.zip")
boe_205.255 <- ListZip(ipdname, separate = FALSE)

boe_205.255_attr <- read_sas(unz(ipdname, boe_205.255$file[76])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.255$file[54])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.255$file[82])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdc) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      grepl("D|H", tpatt) ~ "Exp",
      TRUE ~ "Plc"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdc == 1 ~ "Completed",
      termcdc %in% c(2, 3, 4) ~ "Adverse Event",
      termcdc == 5 ~ "Protocol Deviation",
      termcdc == 6 ~ "Lost to Follow-up",
      termcdc == 7 ~ "Voluntary Withdrawal",
      termcdc == 8 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdc,
    reason
  )

table(boe_205.255_attr$reason, useNA = "ifany")
table(boe_205.255_attr$trtgrp, boe_205.255_attr$reason, useNA = "ifany")

### 205.256 NCT00274053 --------------------------------------------------------

list.files("V:/From 8697/205.256")
ipdname <- c("V:/From 8697/205.256/205.256_VIVLI_v01_nowho.zip")
boe_205.256 <- ListZip(ipdname, separate = FALSE)

boe_205.256_attr <- read_sas(unz(ipdname, boe_205.256$file[48])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.256$file[35])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.256$file[53])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdc) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      grepl("D|H", tpatt) ~ "Exp",
      TRUE ~ "Plc"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdc == 1 ~ "Completed",
      termcdc %in% c(2, 3, 4) ~ "Adverse Event",
      termcdc == 5 ~ "Protocol Deviation",
      termcdc == 6 ~ "Lost to Follow-up",
      termcdc == 7 ~ "Voluntary Withdrawal",
      termcdc == 8 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdc,
    reason
  )

table(boe_205.256_attr$reason, useNA = "ifany")
table(boe_205.256_attr$trtgrp, boe_205.256_attr$reason, useNA = "ifany")

### 205.257 NCT00274573 --------------------------------------------------------

list.files("V:/From 8697/205.257")
ipdname <- c("V:/From 8697/205.257/205.257_VIVLI_v01.zip")
boe_205.257 <- ListZip(ipdname, separate = FALSE)

boe_205.257_attr <- read_sas(unz(ipdname, boe_205.257$file[33])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.257$file[23])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.257$file[36])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdc) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A" ~ "Exp",
      TRUE ~ "Plc"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdc == 1 ~ "Completed",
      termcdc %in% c(2, 3, 4) ~ "Adverse Event",
      termcdc == 5 ~ "Protocol Deviation",
      termcdc == 6 ~ "Lost to Follow-up",
      termcdc == 7 ~ "Voluntary Withdrawal",
      termcdc == 8 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdc,
    reason
  )

table(boe_205.257_attr$reason, useNA = "ifany")
table(boe_205.257_attr$trtgrp, boe_205.257_attr$reason, useNA = "ifany")

### 205.259 NCT00277264 --------------------------------------------------------

list.files("V:/From 8697/205.259")
ipdname <- c("V:/From 8697/205.259/205.259_VIVLI_v01.zip")
boe_205.259 <- ListZip(ipdname, separate = FALSE)

boe_205.259_attr <- read_sas(unz(ipdname, boe_205.259$file[47])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.259$file[32])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.259$file[52])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdc) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A" ~ "Exp",
      TRUE ~ "Plc"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdc == 1 ~ "Completed",
      termcdc %in% c(2, 3, 4) ~ "Adverse Event",
      termcdc == 5 ~ "Protocol Deviation",
      termcdc == 6 ~ "Lost to Follow-up",
      termcdc == 7 ~ "Voluntary Withdrawal",
      termcdc == 8 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdc,
    reason
  )

table(boe_205.259_attr$reason, useNA = "ifany")
table(boe_205.259_attr$trtgrp, boe_205.259_attr$reason, useNA = "ifany")

### 205.264 NCT00274560 --------------------------------------------------------

list.files("V:/From 8697/205.264")
ipdname <- c("V:/From 8697/205.264/205.264_VIVLI_v01.zip")
boe_205.264 <- ListZip(ipdname, separate = FALSE)

boe_205.264_attr <- read_sas(unz(ipdname, boe_205.264$file[49])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.264$file[34])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.264$file[55])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdc) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A" ~ "Exp",
      TRUE ~ "Plc"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdc == 1 ~ "Completed",
      termcdc %in% c(2, 3, 4) ~ "Adverse Event",
      termcdc == 5 ~ "Protocol Deviation",
      termcdc == 6 ~ "Lost to Follow-up",
      termcdc == 7 ~ "Voluntary Withdrawal",
      termcdc == 8 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdc,
    reason
  )

table(boe_205.264_attr$reason, useNA = "ifany")
table(boe_205.264_attr$trtgrp, boe_205.264_attr$reason, useNA = "ifany")

### 205.266 NCT00274547 --------------------------------------------------------

list.files("V:/From 8697/205.266")
ipdname <- c("V:/From 8697/205.266/205.266_VIVLI_v01.zip")
boe_205.266 <- ListZip(ipdname, separate = FALSE)

boe_205.266_attr <- read_sas(unz(ipdname, boe_205.266$file[39])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.266$file[25])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.266$file[42])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdc) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A" ~ "Exp",
      TRUE ~ "Plc"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdc == 1 ~ "Completed",
      termcdc %in% c(2, 3, 4) ~ "Adverse Event",
      termcdc == 5 ~ "Protocol Deviation",
      termcdc == 6 ~ "Lost to Follow-up",
      termcdc == 7 ~ "Voluntary Withdrawal",
      termcdc == 8 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdc,
    reason
  )

table(boe_205.266_attr$reason, useNA = "ifany")
table(boe_205.266_attr$trtgrp, boe_205.266_attr$reason, useNA = "ifany")

### 205.276 NCT00274079 --------------------------------------------------------

list.files("V:/From 8697/205.276")
ipdname <- c("V:/From 8697/205.276/205.276_VIVLI_v01.zip")
boe_205.276 <- ListZip(ipdname, separate = FALSE)

boe_205.276_attr <- read_sas(unz(ipdname, boe_205.276$file[55])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.276$file[39])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.276$file[57])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdc) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A" ~ "Exp",
      TRUE ~ "AC"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdc == 1 ~ "Completed",
      termcdc %in% c(2, 3, 4) ~ "Adverse Event",
      termcdc == 5 ~ "Protocol Deviation",
      termcdc == 6 ~ "Lost to Follow-up",
      termcdc == 7 ~ "Voluntary Withdrawal",
      termcdc == 8 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdc,
    reason
  )

table(boe_205.276_attr$reason, useNA = "ifany")
table(boe_205.276_attr$trtgrp, boe_205.276_attr$reason, useNA = "ifany")

### 205.282 NCT00239408 --------------------------------------------------------

list.files("V:/From 8697/205.282")
ipdname <- c("V:/From 8697/205.282/205.282_VIVLI_v01.zip")
boe_205.282 <- ListZip(ipdname, separate = FALSE)

boe_205.282_attr <- read_sas(unz(ipdname, boe_205.282$file[38])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.282$file[26])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.282$file[40])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdc) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A" ~ "Exp",
      TRUE ~ "Plc"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdc == 1 ~ "Completed",
      termcdc %in% c(2, 3, 4) ~ "Adverse Event",
      termcdc == 5 ~ "Protocol Deviation",
      termcdc == 6 ~ "Lost to Follow-up",
      termcdc == 7 ~ "Voluntary Withdrawal",
      termcdc == 8 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdc,
    reason
  )

table(boe_205.282_attr$reason, useNA = "ifany")
table(boe_205.282_attr$trtgrp, boe_205.282_attr$reason, useNA = "ifany")

### 205.287 NCT00239421 --------------------------------------------------------

list.files("V:/From 8697/205.287")
ipdname <- c("V:/From 8697/205.287/205.287_VIVLI_v01.zip")
boe_205.287 <- ListZip(ipdname, separate = FALSE)

boe_205.287_attr <- read_sas(unz(ipdname, boe_205.287$file[49])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.287$file[37])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.287$file[52])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdc) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A" ~ "Exp",
      TRUE ~ "AC"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdc == 1 ~ "Completed",
      termcdc %in% c(2, 3, 4) ~ "Adverse Event",
      termcdc == 5 ~ "Protocol Deviation",
      termcdc == 6 ~ "Lost to Follow-up",
      termcdc == 7 ~ "Voluntary Withdrawal",
      termcdc == 8 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdc,
    reason
  )

table(boe_205.287_attr$reason, useNA = "ifany")
table(boe_205.287_attr$trtgrp, boe_205.287_attr$reason, useNA = "ifany")

### 205.301 NCT00152984 --------------------------------------------------------

list.files("V:/From 8697/205.301")
ipdname <- c("V:/From 8697/205.301/205.301_VIVLI_v01.zip")
boe_205.301 <- ListZip(ipdname, separate = FALSE)

boe_205.301_attr <- read_sas(unz(ipdname, boe_205.301$file[52])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.301$file[39])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.301$file[55])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdc) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A" ~ "Exp",
      TRUE ~ "Plc"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdc == 1 ~ "Completed",
      termcdc %in% c(2, 3, 4) ~ "Adverse Event",
      termcdc == 5 ~ "Protocol Deviation",
      termcdc == 6 ~ "Lost to Follow-up",
      termcdc == 7 ~ "Voluntary Withdrawal",
      termcdc == 8 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdc,
    reason
  )

table(boe_205.301_attr$reason, useNA = "ifany")
table(boe_205.301_attr$trtgrp, boe_205.301_attr$reason, useNA = "ifany")

### 205.368 NCT00525512 --------------------------------------------------------

list.files("V:/From 8697/205.368")
ipdname <- c("V:/From 8697/205.368/205.368_VIVLI_v01_nowho.zip")
boe_205.368 <- ListZip(ipdname, separate = FALSE)

boe_205.368_attr <- read_sas(unz(ipdname, boe_205.368$file[43])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.368$file[30])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.368$file[47])) %>%
      rename_with(., tolower) %>%
      filter(cpevent == "END BLIND") %>% 
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "AA" ~ "Exp",
      TRUE ~ "Plc"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_205.368_attr$reason, useNA = "ifany")
table(boe_205.368_attr$trtgrp, boe_205.368_attr$reason, useNA = "ifany")

### 205.389 NCT00563381 --------------------------------------------------------

list.files("V:/From 8697/205.389")
ipdname <- c("V:/From 8697/205.389/205.389_VIVLI_v01.zip")
boe_205.389 <- ListZip(ipdname, separate = FALSE)

boe_205.389_attr <- read_sas(unz(ipdname, boe_205.389$file[55])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.389$file[45])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.389$file[58])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A" ~ "Exp",
      TRUE ~ "AC"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_205.389_attr$reason, useNA = "ifany")
table(boe_205.389_attr$trtgrp, boe_205.389_attr$reason, useNA = "ifany")

### 205.416 NCT00772538 --------------------------------------------------------

list.files("V:/From 8697/205.416")
ipdname <- c("V:/From 8697/205.416/205.416_VIVLI_v01.zip")
boe_205.416 <- ListZip(ipdname, separate = FALSE)

boe_205.416_attr <- read_sas(unz(ipdname, boe_205.416$file[67])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.416$file[41])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.416$file[72])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A" ~ "Exp",
      TRUE ~ "Plc"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_205.416_attr$reason, useNA = "ifany")
table(boe_205.416_attr$trtgrp, boe_205.416_attr$reason, useNA = "ifany")

### 205.417 NCT00776984 --------------------------------------------------------

list.files("V:/From 8697/205.417")
ipdname <- c("V:/From 8697/205.417/205.417_VIVLI_v01.zip")
boe_205.417 <- ListZip(ipdname, separate = FALSE)

boe_205.417_attr <- read_sas(unz(ipdname, boe_205.417$file[67])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.417$file[41])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.417$file[72])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A" ~ "Exp",
      TRUE ~ "Plc"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_205.417_attr$reason, useNA = "ifany")
table(boe_205.417_attr$trtgrp, boe_205.417_attr$reason, useNA = "ifany")

### 205.452 NCT01126437 --------------------------------------------------------

list.files("V:/From 8697/205.452")
ipdname <- c("V:/From 8697/205.452/205.452_VIVLI_v01.zip")
boe_205.452 <- ListZip(ipdname, separate = FALSE)

boe_205.452_attr <- read_sas(unz(ipdname, boe_205.452$file[45])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.452$file[26])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  inner_join(
    read_sas(unz(ipdname, boe_205.452$file[7])) %>%
      rename_with(., tolower) %>%
      filter(popudc == "Treated Set" & popuny == 1) %>%
      select(adp_pid)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_205.452$file[49])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      grepl("A|B", tpatt) ~ "Exp",
      TRUE ~ "AC"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_205.452_attr$reason, useNA = "ifany")
table(boe_205.452_attr$trtgrp, boe_205.452_attr$reason, useNA = "ifany")

### 244.2484 NCT02177344 -------------------------------------------------------

list.files("V:/From 8697/244.2484")
ipdname <- c("V:/From 8697/244.2484/244.2484_VIVLI_v01.zip")
boe_244.2484 <- ListZip(ipdname, separate = FALSE)

boe_244.2484_attr <- read_sas(unz(ipdname, boe_244.2484$file[51])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_244.2484$file[16])) %>%
      rename_with(., tolower) %>%
      select(adp_pid, study, tpatt) %>%
      distinct()
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_244.2484$file[31])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_244.2484$file[27])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdc) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      grepl("A|B", tpatt) ~ "Exp",
      tpatt == "D" ~ "AC",
      TRUE ~ "Plc"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdc == 1 ~ "Completed",
      termcdc %in% c(2, 3, 4) ~ "Adverse Event",
      termcdc == 5 ~ "Protocol Deviation",
      termcdc == 6 ~ "Lost to Follow-up",
      termcdc == 7 ~ "Voluntary Withdrawal",
      termcdc == 8 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdc,
    reason
  )

table(boe_244.2484_attr$reason, useNA = "ifany")
table(boe_244.2484_attr$trtgrp, boe_244.2484_attr$reason, useNA = "ifany")

### 248.524 NCT00479401 --------------------------------------------------------

list.files("V:/From 8697/248.524")
ipdname <- c("V:/From 8697/248.524/248.524_VIVLI_v01_nowho.zip")
boe_248.524 <- ListZip(ipdname, separate = FALSE)

boe_248.524_attr <- read_sas(unz(ipdname, boe_248.524$file[49])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_248.524$file[28])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_248.524$file[51])) %>%
      rename_with(., tolower) %>%
      filter(dcmsubnm == "TTM1") %>% 
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt %in% c("B", "C") ~ "Exp",
      TRUE ~ "Plc"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_248.524_attr$reason, useNA = "ifany")
table(boe_248.524_attr$trtgrp, boe_248.524_attr$reason, useNA = "ifany")

### 248.525 NCT00466167 --------------------------------------------------------

list.files("V:/From 8697/248.525")
ipdname <- c("V:/From 8697/248.525/248.525v01.zip")
boe_248.525 <- ListZip(ipdname, separate = FALSE)

boe_248.525_attr <- read_sas(unz(ipdname, boe_248.525$file[57])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_248.525$file[35])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_248.525$file[59])) %>%
      rename_with(., tolower) %>%
      filter(dcmsubnm == "TTM1") %>% 
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A" ~ "Plc",
      tpatt %in% c("B", "C") ~ "Exp"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_248.525_attr$reason, useNA = "ifany")
table(boe_248.525_attr$trtgrp, boe_248.525_attr$reason, useNA = "ifany")

### 248.595 NCT00321854 --------------------------------------------------------

list.files("V:/From 8697/248.595")
ipdname <- c("V:/From 8697/248.595/248.595_VIVLI_v01_nowho.zip")
boe_248.595 <- ListZip(ipdname, separate = FALSE)

boe_248.595_attr <- read_sas(unz(ipdname, boe_248.595$file[39])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_248.595$file[23])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_248.595$file[42])) %>%
      rename_with(., tolower) %>%
      filter(dcmsubnm == "TTM2") %>% 
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = "Exp",
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_248.595_attr$reason, useNA = "ifany")
table(boe_248.595_attr$trtgrp, boe_248.595_attr$reason, useNA = "ifany")

### 248.622 NCT00402233 --------------------------------------------------------

list.files("V:/From 8697/248.622")
ipdname <- c("V:/From 8697/248.622/248.622_VIVLI_v01_nowho.zip")
boe_248.622 <- ListZip(ipdname, separate = FALSE)

boe_248.622_attr <- read_sas(unz(ipdname, boe_248.622$file[35])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, randdt = adp_randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_248.622$file[39])) %>%
      rename_with(., tolower) %>%
      select(adp_pid, atrlbl)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_248.622$file[40])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt = adp_drgstpdt, termcdc) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(atrlbl == "Placebo", "Plc", "Exp"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdc == 1 ~ "Completed",
      termcdc == 2 ~ "Worsening of Disease Under Study",
      termcdc == 3 ~ "Worsening of Other Disease",
      termcdc == 4 ~ "Adverse Event",
      termcdc == 5 ~ "Protocol Deviation",
      termcdc == 7 ~ "Voluntary Withdrawal"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = atrlbl, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdc,
    reason
  )

table(boe_248.622_attr$reason, useNA = "ifany")
table(boe_248.622_attr$trtgrp, boe_248.622_attr$reason, useNA = "ifany")

### 248.629 NCT00472199 --------------------------------------------------------

list.files("V:/From 8697/248.629")
ipdname <- c("V:/From 8697/248.629/248.629_VIVLI_v01.zip")
boe_248.629 <- ListZip(ipdname, separate = FALSE)

boe_248.629_attr <- read_sas(unz(ipdname, boe_248.629$file[39])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_248.629$file[31])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_248.629$file[46])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "Placebo", "Plc", "Exp"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Violation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      TRUE ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_248.629_attr$reason, useNA = "ifany")
table(boe_248.629_attr$trtgrp, boe_248.629_attr$reason, useNA = "ifany")

### 248.671 NCT01191944 --------------------------------------------------------

list.files("V:/From 8697/248.671")
ipdname <- c("V:/From 8697/248.671/248.671_VIVLI_v01.zip")
boe_248.671 <- ListZip(ipdname, separate = FALSE)

boe_248.671_attr <- read_sas(unz(ipdname, boe_248.671$file[47])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_248.671$file[34])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_248.671$file[48])) %>%
      rename_with(., tolower) %>%
      filter(dcmsubnm == "TTM1") %>% 
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "A", "Exp", "AC"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_248.671_attr$reason, useNA = "ifany")
table(boe_248.671_attr$trtgrp, boe_248.671_attr$reason, useNA = "ifany")

### 352.2046 NCT00975195 -------------------------------------------------------

list.files("V:/From 8697/352.2046")
ipdname <- c("V:/From 8697/352.2046/352.2046.zip")
boe_352.2046 <- ListZip(ipdname, separate = FALSE)

boe_352.2046_attr <- read_sas(unz(ipdname, boe_352.2046$file[84])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_352.2046$file[51])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_352.2046$file[90])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = "Exp",
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_352.2046_attr$reason, useNA = "ifany")
table(boe_352.2046_attr$trtgrp, boe_352.2046_attr$reason, useNA = "ifany")

### 502.254 NCT02175355 --------------------------------------------------------

list.files("V:/From 8697/502.254")
ipdname <- c("V:/From 8697/502.254/502.254_VIVLI_v01.zip")
boe_502.254 <- ListZip(ipdname, separate = FALSE)

boe_502.254_attr <- read_sas(unz(ipdname, boe_502.254$file[41])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.254$file[33])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.254$file[26])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcd) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      grepl("Telmi", tpatt) ~ "Exp",
      tpatt == "HCTZ 12.5mg" ~ "AC",
      TRUE ~ "Plc"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcd == 1 ~ "Completed",
      termcd %in% c(2, 3, 4) ~ "Adverse Event",
      termcd == 5 ~ "Lack of Efficacy",
      termcd == 6 ~ "Protocol Deviation",
      termcd == 7 ~ "Lost to Follow-up",
      termcd == 8 ~ "Voluntary Withdrawal",
      termcd == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcd,
    reason
  )

table(boe_502.254_attr$reason, useNA = "ifany")
table(boe_502.254_attr$trtgrp, boe_502.254_attr$reason, useNA = "ifany")

### 502.256 NCT02177396 --------------------------------------------------------

list.files("V:/From 8697/502.256")
ipdname <- c("V:/From 8697/502.256/502.256_VIVLI_v01.zip")
boe_502.256 <- ListZip(ipdname, separate = FALSE)

boe_502.256_attr <- read_sas(unz(ipdname, boe_502.256$file[35])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.256$file[26])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.256$file[22])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcd) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      "Telmi" %in% tpatt ~ "Exp",
      TRUE ~ "AC"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcd == 1 ~ "Completed",
      termcd %in% c(2, 3, 4) ~ "Adverse Event",
      termcd == 5 ~ "Lack of Efficacy",
      termcd == 6 ~ "Protocol Deviation",
      termcd == 7 ~ "Lost to Follow-up",
      termcd == 8 ~ "Voluntary Withdrawal",
      termcd == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcd,
    reason
  )

table(boe_502.256_attr$reason, useNA = "ifany")
table(boe_502.256_attr$trtgrp, boe_502.256_attr$reason, useNA = "ifany")

### 502.316 NCT02172586 --------------------------------------------------------

list.files("V:/From 8697/502.316")
ipdname <- c("V:/From 8697/502.316/502.316_VIVLI_v01.zip")
boe_502.316 <- ListZip(ipdname, separate = FALSE)

boe_502.316_attr <- read_sas(unz(ipdname, boe_502.316$file[36])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.316$file[29])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.316$file[22])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcd) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "A", "Exp", "AC"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcd == 1 ~ "Completed",
      termcd %in% c(2, 3, 4) ~ "Adverse Event",
      termcd == 5 ~ "Lack of Efficacy",
      termcd == 6 ~ "Protocol Deviation",
      termcd == 7 ~ "Lost to Follow-up",
      termcd == 8 ~ "Voluntary Withdrawal",
      termcd == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcd,
    reason
  )

table(boe_502.316_attr$reason, useNA = "ifany")
table(boe_502.316_attr$trtgrp, boe_502.316_attr$reason, useNA = "ifany")

### 502.317 NCT02177461 --------------------------------------------------------

list.files("V:/From 8697/502.317")
ipdname <- c("V:/From 8697/502.317/502.317_VIVLI_v01.zip")
boe_502.317 <- ListZip(ipdname, separate = FALSE)

boe_502.317_attr <- read_sas(unz(ipdname, boe_502.317$file[36])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.317$file[29])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.317$file[22])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcd) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "A", "Exp", "AC"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcd == 1 ~ "Completed",
      termcd %in% c(2, 3, 4) ~ "Adverse Event",
      termcd == 5 ~ "Lack of Efficacy",
      termcd == 6 ~ "Protocol Deviation",
      termcd == 7 ~ "Lost to Follow-up",
      termcd == 8 ~ "Voluntary Withdrawal",
      termcd == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcd,
    reason
  )

table(boe_502.317_attr$reason, useNA = "ifany")
table(boe_502.317_attr$trtgrp, boe_502.317_attr$reason, useNA = "ifany")

### 502.327 NCT00034840 --------------------------------------------------------

list.files("V:/From 8697/502.327")
ipdname <- c("V:/From 8697/502.327/502.327_VIVLI_v01.zip")
boe_502.327 <- ListZip(ipdname, separate = FALSE)

boe_502.327_attr <- read_sas(unz(ipdname, boe_502.327$file[44])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.327$file[36])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.327$file[32])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcd) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(grepl("AAC|ACA", tpatt), "Exp", "AC"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcd == 1 ~ "Completed",
      termcd %in% c(2, 3, 4) ~ "Adverse Event",
      termcd == 5 ~ "Lack of Efficacy",
      termcd == 6 ~ "Protocol Deviation",
      termcd == 7 ~ "Lost to Follow-up",
      termcd == 8 ~ "Voluntary Withdrawal",
      termcd == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcd,
    reason
  )

table(boe_502.327_attr$reason, useNA = "ifany")
table(boe_502.327_attr$trtgrp, boe_502.327_attr$reason, useNA = "ifany")

### 502.343 NCT02200640 --------------------------------------------------------

list.files("V:/From 8697/502.343")
ipdname <- c("V:/From 8697/502.343/502.343_VIVLI_v01.zip")
boe_502.343 <- ListZip(ipdname, separate = FALSE)

boe_502.343_attr <- read_sas(unz(ipdname, boe_502.343$file[42])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.343$file[31])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.343$file[27])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcd) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "A", "Exp", "AC"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcd == 1 ~ "Completed",
      termcd %in% c(2, 3, 4) ~ "Adverse Event",
      termcd == 5 ~ "Lack of Efficacy",
      termcd == 6 ~ "Protocol Deviation",
      termcd == 7 ~ "Lost to Follow-up",
      termcd == 8 ~ "Voluntary Withdrawal",
      termcd == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcd,
    reason
  )

table(boe_502.343_attr$reason, useNA = "ifany")
table(boe_502.343_attr$trtgrp, boe_502.343_attr$reason, useNA = "ifany")

### 502.344 NCT02200653 --------------------------------------------------------

list.files("V:/From 8697/502.344")
ipdname <- c("V:/From 8697/502.344/502.344_VIVLI_v01.zip")
boe_502.344 <- ListZip(ipdname, separate = FALSE)

boe_502.344_attr <- read_sas(unz(ipdname, boe_502.344$file[51])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.344$file[40])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.344$file[32])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcd) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "A", "Exp", "AC"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcd == 1 ~ "Completed",
      termcd %in% c(2, 3, 4) ~ "Adverse Event",
      termcd == 5 ~ "Lack of Efficacy",
      termcd == 6 ~ "Protocol Deviation",
      termcd == 7 ~ "Lost to Follow-up",
      termcd == 8 ~ "Voluntary Withdrawal",
      termcd == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcd,
    reason
  )

table(boe_502.344_attr$reason, useNA = "ifany")
table(boe_502.344_attr$trtgrp, boe_502.344_attr$reason, useNA = "ifany")

### 502.376 NCT02242318 --------------------------------------------------------

list.files("V:/From 8697/502.376")
ipdname <- c("V:/From 8697/502.376/502.376_VIVLI_v01.zip")
boe_502.376 <- ListZip(ipdname, separate = FALSE)

boe_502.376_attr <- read_sas(unz(ipdname, boe_502.376$file[49])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_502.376$file[4])) %>%
      rename_with(., tolower) %>%
      select(adp_pid) %>%
      distinct()
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.376$file[40])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.376$file[36])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcd) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(grepl("Telm", name), "Exp", "AC"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcd == 1 ~ "Completed",
      termcd %in% c(2, 3, 4) ~ "Adverse Event",
      termcd == 5 ~ "Lack of Efficacy",
      termcd == 6 ~ "Protocol Deviation",
      termcd == 7 ~ "Lost to Follow-up",
      termcd == 8 ~ "Voluntary Withdrawal",
      termcd == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcd,
    reason
  )

table(boe_502.376_attr$reason, useNA = "ifany")
table(boe_502.376_attr$trtgrp, boe_502.376_attr$reason, useNA = "ifany")

### 502.391 NCT00274612 --------------------------------------------------------

list.files("V:/From 8697/502.391")
ipdname <- c("V:/From 8697/502.391/502.391_VIVLI_v01.zip")
boe_502.391 <- ListZip(ipdname, separate = FALSE)

boe_502.391_attr <- read_sas(unz(ipdname, boe_502.391$file[50])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_502.391$file[4])) %>%
      rename_with(., tolower) %>%
      select(adp_pid) %>%
      distinct()
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.391$file[39])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.391$file[34])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcd) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "A", "Exp", "AC"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcd == 1 ~ "Completed",
      termcd %in% c(2, 3, 4) ~ "Adverse Event",
      termcd == 5 ~ "Lack of Efficacy",
      termcd == 6 ~ "Protocol Deviation",
      termcd == 7 ~ "Lost to Follow-up",
      termcd == 8 ~ "Voluntary Withdrawal",
      termcd == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcd,
    reason
  )

table(boe_502.391_attr$reason, useNA = "ifany")
table(boe_502.391_attr$trtgrp, boe_502.391_attr$reason, useNA = "ifany")

### 502.392 NCT00274599 --------------------------------------------------------

list.files("V:/From 8697/502.392")
ipdname <- c("V:/From 8697/502.392/502.392_VIVLI_v01.zip")
boe_502.392 <- ListZip(ipdname, separate = FALSE)

boe_502.392_attr <- read_sas(unz(ipdname, boe_502.392$file[45])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.392$file[34])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.392$file[31])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcd) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "A", "Exp", "AC"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcd == 1 ~ "Completed",
      termcd %in% c(2, 3, 4) ~ "Adverse Event",
      termcd == 5 ~ "Lack of Efficacy",
      termcd == 6 ~ "Protocol Deviation",
      termcd == 7 ~ "Lost to Follow-up",
      termcd == 8 ~ "Voluntary Withdrawal",
      termcd == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcd,
    reason
  )

table(boe_502.392_attr$reason, useNA = "ifany")
table(boe_502.392_attr$trtgrp, boe_502.392_attr$reason, useNA = "ifany")

### 502.396 NCT00153023 --------------------------------------------------------

list.files("V:/From 8697/502.396")
ipdname <- c("V:/From 8697/502.396/502.396_VIVLI_v01.zip")
boe_502.396 <- ListZip(ipdname, separate = FALSE)

boe_502.396_attr <- read_sas(unz(ipdname, boe_502.396$file[58])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.396$file[41])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.396$file[36])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcd) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "A", "Exp", "AC"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcd == 1 ~ "Completed",
      termcd %in% c(2, 3, 4) ~ "Adverse Event",
      termcd == 5 ~ "Lack of Efficacy",
      termcd == 6 ~ "Protocol Deviation",
      termcd == 7 ~ "Lost to Follow-up",
      termcd == 8 ~ "Voluntary Withdrawal",
      termcd == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcd,
    reason
  )

table(boe_502.396_attr$reason, useNA = "ifany")
table(boe_502.396_attr$trtgrp, boe_502.396_attr$reason, useNA = "ifany")

### 502.397 NCT00168857 --------------------------------------------------------

list.files("V:/From 8697/502.397")
ipdname <- c("V:/From 8697/502.397/502.397_VIVLI_v01.zip")
boe_502.397 <- ListZip(ipdname, separate = FALSE)

boe_502.397_attr <- read_sas(unz(ipdname, boe_502.397$file[48])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.397$file[32])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.397$file[29])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcd) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "A", "Exp", "AC"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcd == 1 ~ "Completed",
      termcd %in% c(2, 3, 4) ~ "Adverse Event",
      termcd == 5 ~ "Lack of Efficacy",
      termcd == 6 ~ "Protocol Deviation",
      termcd == 7 ~ "Lost to Follow-up",
      termcd == 8 ~ "Voluntary Withdrawal",
      termcd == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcd,
    reason
  )

table(boe_502.397_attr$reason, useNA = "ifany")
table(boe_502.397_attr$trtgrp, boe_502.397_attr$reason, useNA = "ifany")

### 502.413 NCT00153088 --------------------------------------------------------

list.files("V:/From 8697/502.413")
ipdname <- c("V:/From 8697/502.413/502.413_VIVLI_v01.zip")
boe_502.413 <- ListZip(ipdname, separate = FALSE)

boe_502.413_attr <- read_sas(unz(ipdname, boe_502.413$file[44])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.413$file[32])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.413$file[47])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdc) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "C", "Plc", "Exp"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdc == 1 ~ "Completed",
      termcdc %in% c(2, 3, 4) ~ "Adverse Event",
      termcdc == 5 ~ "Protocol Deviation",
      termcdc == 6 ~ "Lost to Follow-up",
      termcdc == 7 ~ "Voluntary Withdrawal",
      termcdc == 8 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdc,
    reason
  )

table(boe_502.413_attr$reason, useNA = "ifany")
table(boe_502.413_attr$trtgrp, boe_502.413_attr$reason, useNA = "ifany")

### 502.550 NCT00926289 --------------------------------------------------------

list.files("V:/From 8697/502.550")
ipdname <- c("V:/From 8697/502.550/502.550_VIVLI_v01.zip")
boe_502.550 <- ListZip(ipdname, separate = FALSE)

boe_502.550_attr <- read_sas(unz(ipdname, boe_502.550$file[29])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>%
  inner_join(
    read_sas(unz(ipdname, boe_502.550$file[5])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid) %>%
      distinct()
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.550$file[20])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_502.550$file[30])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdc) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "A", "AC", "Exp"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdc == 1 ~ "Completed",
      termcdc %in% c(2, 3, 4) ~ "Adverse Event",
      termcdc == 5 ~ "Protocol Deviation",
      termcdc == 6 ~ "Lost to Follow-up",
      termcdc == 7 ~ "Voluntary Withdrawal",
      termcdc == 8 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdc,
    reason
  )

table(boe_502.550_attr$reason, useNA = "ifany")
table(boe_502.550_attr$trtgrp, boe_502.550_attr$reason, useNA = "ifany")

### 1199.32 NCT01619085 --------------------------------------------------------

list.files("V:/From 8697/1199.32")
ipdname <- c("V:/From 8697/1199.32/1199.32_VIVLI_v01.zip")
boe_1199.32 <- ListZip(ipdname, separate = FALSE)

boe_1199.32_attr <- read_sas(unz(ipdname, boe_1199.32$file[102])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1199.32$file[67])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1199.32$file[108])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdc) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "A", "Exp", "Plc"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdc == 1 ~ "Completed",
      termcdc %in% c(2, 3, 4) ~ "Adverse Event",
      termcdc == 5 ~ "Protocol Deviation",
      termcdc == 6 ~ "Lost to Follow-up",
      termcdc == 7 ~ "Voluntary Withdrawal",
      termcdc == 8 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdc,
    reason
  )

table(boe_1199.32_attr$reason, useNA = "ifany")
table(boe_1199.32_attr$trtgrp, boe_1199.32_attr$reason, useNA = "ifany")

### 1199.34 NCT01335477 --------------------------------------------------------

list.files("V:/From 8697/1199.34")
ipdname <- c("V:/From 8697/1199.34/1199.34_VIVLI_v01.zip")
boe_1199.34 <- ListZip(ipdname, separate = FALSE)

boe_1199.34_attr <- read_sas(unz(ipdname, boe_1199.34$file[103])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1199.34$file[67])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1199.34$file[109])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdc) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A" ~ "Exp",
      tpatt == "B" ~ "Plc",
      TRUE ~ tpatt
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdc == 1 ~ "Completed",
      termcdc %in% c(2, 3, 4) ~ "Adverse Event",
      termcdc == 5 ~ "Protocol Deviation",
      termcdc == 6 ~ "Lost to Follow-up",
      termcdc == 7 ~ "Voluntary Withdrawal",
      termcdc == 8 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdc,
    reason
  )

table(boe_1199.34_attr$reason, useNA = "ifany")
table(boe_1199.34_attr$trtgrp, boe_1199.34_attr$reason, useNA = "ifany")

### 1218.15 NCT00641043 --------------------------------------------------------

list.files("V:/From 8697/1218.15")
ipdname <- c("V:/From 8697/1218.15/1218.15_VIVLI_v01.zip")
boe_1218.15 <- ListZip(ipdname, separate = FALSE)

boe_1218.15_attr <- read_sas(unz(ipdname, boe_1218.15$file[45])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.15$file[26])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.15$file[48])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "A", "Exp", "Plc"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1218.15_attr$reason, useNA = "ifany")
table(boe_1218.15_attr$trtgrp, boe_1218.15_attr$reason, useNA = "ifany")

### 1218.16 NCT00621140 --------------------------------------------------------

list.files("V:/From 8697/1218.16")
ipdname <- "V:/From 8697/1218.16/1218.16_VIVLI_v01.zip"
boe_1218.16 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1218.16_attr <- read_sas(unz(ipdname, boe_1218.16$file[45])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.16$file[25])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.16$file[48])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "A", "Exp", "Plc"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1218.16_attr$reason, useNA = "ifany")
table(boe_1218.16_attr$trtgrp, boe_1218.16_attr$reason, useNA = "ifany")

### 1218.17 NCT00601250 --------------------------------------------------------

list.files("V:/From 8697/1218.17")
ipdname <- "V:/From 8697/1218.17/1218.17_VIVLI_v01.zip"
boe_1218.17 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1218.17_attr <- read_sas(unz(ipdname, boe_1218.17$file[45])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.17$file[26])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.17$file[48])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "A", "Exp", "Plc"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1218.17_attr$reason, useNA = "ifany")
table(boe_1218.17_attr$trtgrp, boe_1218.17_attr$reason, useNA = "ifany")

### 1218.18 NCT00602472 --------------------------------------------------------

list.files("V:/From 8697/1218.18")
ipdname <- "V:/From 8697/1218.18/1218.18_VIVLI_v01.zip"
boe_1218.18 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1218.18_attr <- read_sas(unz(ipdname, boe_1218.18$file[45])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.18$file[26])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.18$file[48])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "A", "Exp", "Plc"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1218.18_attr$reason, useNA = "ifany")
table(boe_1218.18_attr$trtgrp, boe_1218.18_attr$reason, useNA = "ifany")

### 1218.20 NCT00622284 --------------------------------------------------------

list.files("V:/From 8697/1218.20")
ipdname <- "V:/From 8697/1218.20/1218.20_VIVLI_v01.zip"
boe_1218.20 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1218.20_attr <- read_sas(unz(ipdname, boe_1218.20$file[57])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.20$file[34])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.20$file[60])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "A", "Exp", "AC"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1218.20_attr$reason, useNA = "ifany")
table(boe_1218.20_attr$trtgrp, boe_1218.20_attr$reason, useNA = "ifany")

### 1218.23 NCT00654381 --------------------------------------------------------

list.files("V:/From 8697/1218.23")
ipdname <- "V:/From 8697/1218.23/1218.23v02.zip"
boe_1218.23 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1218.23_attr <- read_sas(unz(ipdname, boe_1218.23$file[48])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.23$file[34])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.23$file[51])) %>%
      rename_with(., tolower) %>%
      filter(dcmsubnm == "TTM2") %>% 
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      grepl("AAA|BBB", tpatt) ~ "Exp",
      grepl("CCA|CCB", tpatt) ~ "AC",
      TRUE ~ "Plc"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other",
      TRUE ~ "Missing"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1218.23_attr$reason, useNA = "ifany")
table(boe_1218.23_attr$trtgrp, boe_1218.23_attr$reason, useNA = "ifany")

### 1218.35 NCT00819091 --------------------------------------------------------

list.files("V:/From 8697/1218.35")
ipdname <- "V:/From 8697/1218.35/1218.35v01.zip"
boe_1218.35 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1218.35_attr <- read_sas(unz(ipdname, boe_1218.35$file[47])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.35$file[30])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.35$file[50])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "A", "Exp", "Plc"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1218.35_attr$reason, useNA = "ifany")
table(boe_1218.35_attr$trtgrp, boe_1218.35_attr$reason, useNA = "ifany")

### 1218.36 NCT00954447 --------------------------------------------------------

list.files("V:/From 8697/1218.36")
ipdname <- "V:/From 8697/1218.36/1218.36v02.zip"
boe_1218.36 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1218.36_attr <- read_sas(unz(ipdname, boe_1218.36$file[49])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.36$file[29])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.36$file[52])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "A", "Exp", "Plc"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1218.36_attr$reason, useNA = "ifany")
table(boe_1218.36_attr$trtgrp, boe_1218.36_attr$reason, useNA = "ifany")

### 1218.43 NCT00800683 --------------------------------------------------------

list.files("V:/From 8697/1218.43")
ipdname <- "V:/From 8697/1218.43/1218.43v01.zip"
boe_1218.43 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1218.43_attr <- read_sas(unz(ipdname, boe_1218.43$file[48])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.43$file[30])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.43$file[51])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "A", "Exp", "Plc"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1218.43_attr$reason, useNA = "ifany")
table(boe_1218.43_attr$trtgrp, boe_1218.43_attr$reason, useNA = "ifany")

### 1218.46 NCT00800683 --------------------------------------------------------

list.files("V:/From 8697/1218.46")
ipdname <- "V:/From 8697/1218.46/1218.46v02.zip"
boe_1218.46 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1218.46_attr <- read_sas(unz(ipdname, boe_1218.46$file[59])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1218.46$file[10])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.46$file[35])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.46$file[63])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt %in% c("A", "B", "C") ~ "Exp", 
      tpatt %in% c("D", "E") ~ "AC",
      TRUE ~ "Plc"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1218.46_attr$reason, useNA = "ifany")
table(boe_1218.46_attr$trtgrp, boe_1218.46_attr$reason, useNA = "ifany")

### 1218.50 NCT00798161 --------------------------------------------------------

list.files("V:/From 8697/1218.50")
ipdname <- "V:/From 8697/1218.50/1218.50v01.zip"
boe_1218.50 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1218.50_attr <- read_sas(unz(ipdname, boe_1218.50$file[48])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1218.50$file[4])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.50$file[28])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.50$file[52])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "A", "Exp", "Plc"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1218.50_attr$reason, useNA = "ifany")
table(boe_1218.50_attr$trtgrp, boe_1218.50_attr$reason, useNA = "ifany")

### 1218.52 NCT00915772 --------------------------------------------------------

list.files("V:/From 8697/1218.52")
ipdname <- "V:/From 8697/1218.52/1218.52v02.zip"
boe_1218.52 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1218.52_attr <- read_sas(unz(ipdname, boe_1218.52$file[34])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1218.52$file[5])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.52$file[24])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.52$file[36])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt %in% c("B", "C") ~ "Exp", 
      TRUE ~ "AC"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1218.52_attr$reason, useNA = "ifany")
table(boe_1218.52_attr$trtgrp, boe_1218.52_attr$reason, useNA = "ifany")

### 1218.60 NCT01438814 --------------------------------------------------------

list.files("V:/From 8697/1218.60")
ipdname <- "V:/From 8697/1218.60/1218.60v02.zip"
boe_1218.60 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1218.60_attr <- read_sas(unz(ipdname, boe_1218.60$file[38])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1218.60$file[6])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.60$file[27])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.60$file[41])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A" ~ "Exp", 
      TRUE ~ "AC"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1218.60_attr$reason, useNA = "ifany")
table(boe_1218.60_attr$trtgrp, boe_1218.60_attr$reason, useNA = "ifany")

### 1218.61 NCT00996658 --------------------------------------------------------

list.files("V:/From 8697/1218.61")
ipdname <- "V:/From 8697/1218.61/1218.61v01.zip"
boe_1218.61 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1218.61_attr <- read_sas(unz(ipdname, boe_1218.61$file[46])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1218.61$file[7])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.61$file[29])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.61$file[49])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "A", "Exp", "Plc"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1218.61_attr$reason, useNA = "ifany")
table(boe_1218.61_attr$trtgrp, boe_1218.61_attr$reason, useNA = "ifany")

### 1218.63 NCT01084005 --------------------------------------------------------

list.files("V:/From 8697/1218.63")
ipdname <- "V:/From 8697/1218.63/1218.63v02.zip"
boe_1218.63 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1218.63_attr <- read_sas(unz(ipdname, boe_1218.63$file[40])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1218.63$file[5])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.63$file[27])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.63$file[51])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A" ~ "Exp", 
      TRUE ~ "Plc"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1218.63_attr$reason, useNA = "ifany")
table(boe_1218.63_attr$trtgrp, boe_1218.63_attr$reason, useNA = "ifany")

### 1218.64 NCT01087502 --------------------------------------------------------

list.files("V:/From 8697/1218.64")
ipdname <- "V:/From 8697/1218.64/1218.64v02.zip"
boe_1218.64 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1218.64_attr <- read_sas(unz(ipdname, boe_1218.64$file[48])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1218.64$file[7])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.64$file[31])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.64$file[59])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A" ~ "Exp", 
      TRUE ~ "Plc/AC"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1218.64_attr$reason, useNA = "ifany")
table(boe_1218.64_attr$trtgrp, boe_1218.64_attr$reason, useNA = "ifany")

### 1218.65 NCT01215097 --------------------------------------------------------

list.files("V:/From 8697/1218.65")
ipdname <- "V:/From 8697/1218.65/1218.65v02.zip"
boe_1218.65 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1218.65_attr <- read_sas(unz(ipdname, boe_1218.65$file[40])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1218.65$file[6])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.65$file[28])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.65$file[43])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A" ~ "Exp", 
      TRUE ~ "Plc"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1218.65_attr$reason, useNA = "ifany")
table(boe_1218.65_attr$trtgrp, boe_1218.65_attr$reason, useNA = "ifany")

### 1218.66 NCT01214239 --------------------------------------------------------

list.files("V:/From 8697/1218.66")
ipdname <- "V:/From 8697/1218.66/1218.66v01.zip"
boe_1218.66 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1218.66_attr <- read_sas(unz(ipdname, boe_1218.66$file[38])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1218.66$file[6])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.66$file[26])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.66$file[41])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A" ~ "Plc", 
      TRUE ~ "Exp"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1218.66_attr$reason, useNA = "ifany")
table(boe_1218.66_attr$trtgrp, boe_1218.66_attr$reason, useNA = "ifany")

### 1218.75 NCT01194830 --------------------------------------------------------

list.files("V:/From 8697/1218.75")
ipdname <- "V:/From 8697/1218.75/1218.75v01.zip"
boe_1218.75 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1218.75_attr <- read_sas(unz(ipdname, boe_1218.75$file[54])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1218.75$file[7])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.75$file[30])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.75$file[57])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "A", "Plc", "Exp"),
    trtgrp = case_when(
      name == "Treatment A" ~ "Placebo",
      TRUE ~ "Linagliptin"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1218.75_attr$reason, useNA = "ifany")
table(boe_1218.75_attr$trtgrp, boe_1218.75_attr$reason, useNA = "ifany")

### 1218.78 NCT01204294 --------------------------------------------------------

list.files("V:/From 8697/1218.78")
ipdname <- "V:/From 8697/1218.78/1218.78v02.zip"
boe_1218.78 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1218.78_attr <- read_sas(unz(ipdname, boe_1218.78$file[40])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1218.78$file[6])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.78$file[29])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.78$file[43])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt %in% c("A", "C", "E", "F", "G") ~ "Exp", 
      tpatt %in% c("B", "D") ~ "AC"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1218.78_attr$reason, useNA = "ifany")
table(boe_1218.78_attr$trtgrp, boe_1218.78_attr$reason, useNA = "ifany")

### 1218.89 NCT01792518 --------------------------------------------------------

list.files("V:/From 8697/1218.89")
ipdname <- "V:/From 8697/1218.89/1218.89v02.zip"
boe_1218.89 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1218.89_attr <- read_sas(unz(ipdname, boe_1218.89$file[73])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1218.89$file[20])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.89$file[51])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.89$file[76])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A" ~ "Exp", 
      tpatt == "B" ~ "Plc"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1218.89_attr$reason, useNA = "ifany")
table(boe_1218.89_attr$trtgrp, boe_1218.89_attr$reason, useNA = "ifany")

### 1218.149 NCT02240680 -------------------------------------------------------

list.files("V:/From 8697/1218.149")
ipdname <- "V:/From 8697/1218.149/1218.149v02.zip"
boe_1218.149 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1218.149_attr <- read_sas(unz(ipdname, boe_1218.149$file[67])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1218.149$file[11])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.149$file[46])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1218.149$file[70])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "Y", "Exp", "Plc"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1218.149_attr$reason, useNA = "ifany")
table(boe_1218.149_attr$trtgrp, boe_1218.149_attr$reason, useNA = "ifany")

### 1222.11 NCT00782210 --------------------------------------------------------

list.files("V:/From 8697/1222.11")
ipdname <- "V:/From 8697/1222.11/1222.11_VIVLI_v01.zip"
boe_1222.11 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1222.11_attr <- read_sas(unz(ipdname, boe_1222.11$file[56])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1222.11$file[30])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1222.11$file[61])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcd) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "C", "Plc", "Exp"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcd == 1 ~ "Completed",
      termcd %in% c(2, 3, 4) ~ "Adverse Event",
      termcd == 5 ~ "Lack of Efficacy",
      termcd == 6 ~ "Protocol Deviation",
      termcd == 7 ~ "Lost to Follow-up",
      termcd == 8 ~ "Voluntary Withdrawal",
      termcd == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcd,
    reason
  )

table(boe_1222.11_attr$reason, useNA = "ifany")
table(boe_1222.11_attr$trtgrp, boe_1222.11_attr$reason, useNA = "ifany")

### 1222.12 NCT00782509 --------------------------------------------------------

list.files("V:/From 8697/1222.12")
ipdname <- "V:/From 8697/1222.12/1222.12_VIVLI_v01.zip"
boe_1222.12 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1222.12_attr <- read_sas(unz(ipdname, boe_1222.12$file[56])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1222.12$file[30])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1222.12$file[61])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcd) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "C", "Plc", "Exp"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcd == 1 ~ "Completed",
      termcd %in% c(2, 3, 4) ~ "Adverse Event",
      termcd == 5 ~ "Lack of Efficacy",
      termcd == 6 ~ "Protocol Deviation",
      termcd == 7 ~ "Lost to Follow-up",
      termcd == 8 ~ "Voluntary Withdrawal",
      termcd == 9 ~ "Other",
      TRUE ~ "Missing"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcd,
    reason
  )

table(boe_1222.12_attr$reason, useNA = "ifany")
table(boe_1222.12_attr$trtgrp, boe_1222.12_attr$reason, useNA = "ifany")

### 1222.13 NCT00793624 --------------------------------------------------------

list.files("V:/From 8697/1222.13")
ipdname <- "V:/From 8697/1222.13/1222.13_VIVLI_v01.zip"
boe_1222.13 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1222.13_attr <- read_sas(unz(ipdname, boe_1222.13$file[60])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1222.13$file[34])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1222.13$file[67])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcd) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "C" ~ "AC",
      tpatt == "D" ~ "Plc",
      TRUE ~ "Exp"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcd == 1 ~ "Completed",
      termcd %in% c(2, 3, 4) ~ "Adverse Event",
      termcd == 5 ~ "Lack of Efficacy",
      termcd == 6 ~ "Protocol Deviation",
      termcd == 7 ~ "Lost to Follow-up",
      termcd == 8 ~ "Voluntary Withdrawal",
      termcd == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcd,
    reason
  )

table(boe_1222.13_attr$reason, useNA = "ifany")
table(boe_1222.13_attr$trtgrp, boe_1222.13_attr$reason, useNA = "ifany")

### 1222.51 NCT01694771 --------------------------------------------------------

list.files("V:/From 8697/1222.51")
ipdname <- "V:/From 8697/1222.51/1222.51_VIVLI_v02.zip"
boe_1222.51 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1222.51_attr <- read_sas(unz(ipdname, boe_1222.51$file[36])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1222.51$file[27])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1222.51$file[42])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcd) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A" ~ "Exp",
      TRUE ~ "Plc"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcd == 1 ~ "Completed",
      termcd %in% c(2, 3, 4) ~ "Adverse Event",
      termcd == 5 ~ "Lack of Efficacy",
      termcd == 6 ~ "Protocol Deviation",
      termcd == 7 ~ "Lost to Follow-up",
      termcd == 8 ~ "Voluntary Withdrawal",
      termcd == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcd,
    reason
  )

table(boe_1222.51_attr$reason, useNA = "ifany")
table(boe_1222.51_attr$trtgrp, boe_1222.51_attr$reason, useNA = "ifany")

### 1245.19 NCT01210001 --------------------------------------------------------

list.files("V:/From 8697/1245.19")
ipdname <- "V:/From 8697/1245.19/1245.19_VIVLI_v01.zip"
boe_1245.19 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1245.19_attr <- read_sas(unz(ipdname, boe_1245.19$file[84])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1245.19$file[61])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1245.19$file[87])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A" ~ "Plc",
      TRUE ~ "Exp"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1245.19_attr$reason, useNA = "ifany")
table(boe_1245.19_attr$trtgrp, boe_1245.19_attr$reason, useNA = "ifany")

### 1245.20 NCT01177813 --------------------------------------------------------

list.files("V:/From 8697/1245.20")
ipdname <- "V:/From 8697/1245.20/1245.20_VIVLI_v01.zip"
boe_1245.20 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1245.20_attr <- read_sas(unz(ipdname, boe_1245.20$file[85])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1245.20$file[31])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1245.20$file[61])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1245.20$file[88])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A" ~ "Plc",
      grepl("B|C|D", tpatt) ~ "Exp",
      TRUE ~ "AC"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1245.20_attr$reason, useNA = "ifany")
table(boe_1245.20_attr$trtgrp, boe_1245.20_attr$reason, useNA = "ifany")

### 1245.23 NCT01159600 --------------------------------------------------------

list.files("V:/From 8697/1245.23")
ipdname <- "V:/From 8697/1245.23/1245.23_VIVLI_v01.zip"
boe_1245.23 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1245.23_attr <- read_sas(unz(ipdname, boe_1245.23$file[88])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1245.23$file[33])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1245.23$file[61])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not Entered/Randomised") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1245.23$file[91])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "C" ~ "Plc",
      TRUE ~ "Exp"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1245.23_attr$reason, useNA = "ifany")
table(boe_1245.23_attr$trtgrp, boe_1245.23_attr$reason, useNA = "ifany")

### 1245.25 NCT01131676 --------------------------------------------------------

list.files("V:/From 8697/1245.25")
ipdname <- "V:/From 8697/1245.25/1245-0025_v01.zip"
boe_1245.25 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1245.25_attr <- read_sas(unz(ipdname, boe_1245.25$file[94])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1245.25$file[35])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1245.25$file[70])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1245.25$file[98])) %>%
      rename_with(., tolower) %>%
      filter(actevent == "970.01") %>% 
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt %in% c("A", "B") ~ "Exp", 
      tpatt == "C" ~ "Plc"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other",
      TRUE ~ "Missing"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1245.25_attr$reason, useNA = "ifany")
table(boe_1245.25_attr$trtgrp, boe_1245.25_attr$reason, useNA = "ifany")

### 1245.28 NCT01159600 --------------------------------------------------------

list.files("V:/From 8697/1245.28")
ipdname <- "V:/From 8697/1245.28/1245-0028_v01.zip"
boe_1245.28 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1245.28_attr <- read_sas(unz(ipdname, boe_1245.28$file[98])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1245.28$file[41])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1245.28$file[72])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not rand") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1245.28$file[101])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "AA025" ~ "Exp",
      TRUE ~ "AC"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1245.28_attr$reason, useNA = "ifany")
table(boe_1245.28_attr$trtgrp, boe_1245.28_attr$reason, useNA = "ifany")

### 1245.29 NCT02182830 --------------------------------------------------------

list.files("V:/From 8697/1245.29")
ipdname <- "V:/From 8697/1245.29/1245.29v01.zip"
boe_1245.29 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1245.29_attr <- read_sas(unz(ipdname, boe_1245.29$file[74])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1245.29$file[56])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not rand") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1245.29$file[77])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A00P3" ~ "Plc",
      TRUE ~ "Exp"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1245.29_attr$reason, useNA = "ifany")
table(boe_1245.29_attr$trtgrp, boe_1245.29_attr$reason, useNA = "ifany")

### 1245.36 NCT01164501 --------------------------------------------------------

list.files("V:/From 8697/1245.36")
ipdname <- "V:/From 8697/1245.36/1245.36_VIVLI_v01.zip"
boe_1245.36 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1245.36_attr <- read_sas(unz(ipdname, boe_1245.36$file[83])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1245.36$file[60])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not rand") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1245.36$file[86])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "C" ~ "Plc",
      TRUE ~ "Exp"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1245.36_attr$reason, useNA = "ifany")
table(boe_1245.36_attr$trtgrp, boe_1245.36_attr$reason, useNA = "ifany")

### 1245.48 NCT01370005 --------------------------------------------------------

list.files("V:/From 8697/1245.48")
ipdname <- "V:/From 8697/1245.48/1245.48_VIVLI_v01.zip"
boe_1245.48 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1245.48_attr <- read_sas(unz(ipdname, boe_1245.48$file[80])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1245.48$file[64])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not rand") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1245.48$file[83])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "C" ~ "Plc",
      TRUE ~ "Exp"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other",
      TRUE ~ "Missing"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1245.48_attr$reason, useNA = "ifany")
table(boe_1245.48_attr$trtgrp, boe_1245.48_attr$reason, useNA = "ifany")

### 1245.49 NCT01306214 --------------------------------------------------------

list.files("V:/From 8697/1245.49")
ipdname <- "V:/From 8697/1245.49/1245.49v02.zip"
boe_1245.49 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1245.49_attr <- read_sas(unz(ipdname, boe_1245.49$file[99])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1245.49$file[68])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not rand") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1245.49$file[101])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "A00P1" ~ "Plc",
      TRUE ~ "Exp"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1245.49_attr$reason, useNA = "ifany")
table(boe_1245.49_attr$trtgrp, boe_1245.49_attr$reason, useNA = "ifany")

### 1245.52 NCT01368081 --------------------------------------------------------

list.files("V:/From 8697/1245.52")
ipdname <- "V:/From 8697/1245.52/1245.52v02.zip"
boe_1245.52 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1245.52_attr <- read_sas(unz(ipdname, boe_1245.52$file[78])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1245.52$file[62])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not rand") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1245.52$file[81])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "Metformin" ~ "AC",
      TRUE ~ "Exp"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1245.52_attr$reason, useNA = "ifany")
table(boe_1245.52_attr$trtgrp, boe_1245.52_attr$reason, useNA = "ifany")

### 1275.1 NCT01422876 ---------------------------------------------------------

list.files("V:/From 8697/1275.1")
ipdname <- "V:/From 8697/1275.1/1275.1v03.zip"
boe_1275.1 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1275.1_attr <- read_sas(unz(ipdname, boe_1275.1$file[79])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1275.1$file[29])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid)
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1275.1$file[61])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not rand") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1275.1$file[81])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv) %>% 
      distinct()
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt %in% c("AA010", "AA025", "AL005") ~ "AC",
      TRUE ~ "Exp"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other",
      TRUE ~ "Missing"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1275.1_attr$reason, useNA = "ifany")
table(boe_1275.1_attr$trtgrp, boe_1275.1_attr$reason, useNA = "ifany")

### 1275.9 NCT01734785 ---------------------------------------------------------

list.files("V:/From 8697/1275.9")
ipdname <- "V:/From 8697/1275.9/1275.9v02.zip"
boe_1275.9 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1275.9_attr <- read_sas(unz(ipdname, boe_1275.9$file[85])) %>%
  rename_with(., tolower) %>%
  filter(
    randel == 1
    & dcmsubnm == "RAND"
  ) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1275.9$file[32])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid) %>% 
      distinct()
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1275.9$file[64])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not rand") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1275.9$file[87])) %>%
      rename_with(., tolower) %>%
      filter(!is.na(drgstpdt)) %>%
      select(study, adp_pid, drgstpdt, termcdv)
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt == "AL005" ~ "Plc",
      TRUE ~ "Exp"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other",
      TRUE ~ "Missing"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1275.9_attr$reason, useNA = "ifany")
table(boe_1275.9_attr$trtgrp, boe_1275.9_attr$reason, useNA = "ifany")

### 1275.10 NCT01778049 --------------------------------------------------------

list.files("V:/From 8697/1275.10")
ipdname <- "V:/From 8697/1275.10/1275.10v02.zip"
boe_1275.10 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1275.10_attr <- read_sas(unz(ipdname, boe_1275.10$file[88])) %>%
  rename_with(., tolower) %>%
  filter(
    randel == 1
    & dcmsubnm == "RAND"
  ) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1275.10$file[32])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid) %>% 
      distinct()
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1275.10$file[65])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not rand") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1275.10$file[90])) %>%
      rename_with(., tolower) %>%
      filter(!is.na(drgstpdt)) %>%
      select(study, adp_pid, drgstpdt, termcdv)
  ) %>% 
  mutate(
    trttype = "Exp",
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other",
      TRUE ~ "Missing"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1275.10_attr$reason, useNA = "ifany")
table(boe_1275.10_attr$trtgrp, boe_1275.10_attr$reason, useNA = "ifany")

### 1275.13 NCT02489968 --------------------------------------------------------

list.files("V:/From 8697/1275.13")
ipdname <- "V:/From 8697/1275.13/1275.13v02.zip"
boe_1275.13 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1275.13_attr <- read_sas(unz(ipdname, boe_1275.13$file[81])) %>%
  rename_with(., tolower) %>%
  filter(
    randel == 1
    & !is.na(randdt)
    & dcmsubnm == "RAND_V5"
  ) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1275.13$file[33])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid) %>% 
      distinct()
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1275.13$file[65])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not rand") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1275.13$file[84])) %>%
      rename_with(., tolower) %>%
      filter(!is.na(drgstpdt)) %>%
      select(study, adp_pid, drgstpdt, termcdv)
  ) %>% 
  mutate(
    trttype = "Exp",
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other",
      TRUE ~ "Missing"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1275.13_attr$reason, useNA = "ifany")
table(boe_1275.13_attr$trtgrp, boe_1275.13_attr$reason, useNA = "ifany")

### 1275.19 NCT02453555 --------------------------------------------------------

list.files("V:/From 8697/1275.19")
ipdname <- "V:/From 8697/1275.19/1275.19v01.zip"
boe_1275.19 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1275.19_attr <- read_sas(unz(ipdname, boe_1275.19$file[75])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1 & dcmsubnm == "RAND") %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1275.19$file[32])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid) %>% 
      distinct()
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1275.19$file[61])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not rand") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1275.19$file[78])) %>%
      rename_with(., tolower) %>%
      filter(!is.na(drgstpdt)) %>%
      select(study, adp_pid, drgstpdt, termcdv)
  ) %>% 
  mutate(
    trttype = if_else(tpatt == "LA062", "Exp", "AC"),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1275.19_attr$reason, useNA = "ifany")
table(boe_1275.19_attr$trtgrp, boe_1275.19_attr$reason, useNA = "ifany")

### 1276.1 NCT01719003 ---------------------------------------------------------

list.files("V:/From 8697/1276.1")
ipdname <- "V:/From 8697/1276.1/1276.1_VIVLI_v01.zip"
boe_1276.1 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1276.1_attr <- read_sas(unz(ipdname, boe_1276.1$file[82])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1 & dcmsubnm == "RAND") %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1276.1$file[31])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid) %>% 
      distinct()
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1276.1$file[59])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(name != "Not rand") %>% 
  left_join(
    read_sas(unz(ipdname, boe_1276.1$file[84])) %>%
      rename_with(., tolower) %>%
      filter(!is.na(drgstpdt)) %>%
      select(study, adp_pid, drgstpdt, termcdv)
  ) %>% 
  mutate(
    trttype = "Exp",
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1276.1_attr$reason, useNA = "ifany")
table(boe_1276.1_attr$trtgrp, boe_1276.1_attr$reason, useNA = "ifany")

### 1288.18 NCT01734785 --------------------------------------------------------

list.files("V:/From 8697/1288.18")
ipdname <- "V:/From 8697/1288.18/1288.18v02.zip"
boe_1288.18 <- ListZip(ipdname = ipdname, separate = FALSE)

boe_1288.18_attr <- read_sas(unz(ipdname, boe_1288.18$file[46])) %>%
  rename_with(., tolower) %>%
  filter(randel == 1) %>% 
  select(study, adp_pid, tpatt, randdt) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1288.18$file[8])) %>%
      rename_with(., tolower) %>%
      filter(popu == "TS" & popuny == 1) %>%
      select(adp_pid) %>% 
      distinct()
  ) %>% 
  left_join(
    read_sas(unz(ipdname, boe_1288.18$file[32])) %>%
      rename_with(., tolower) %>% 
      select(study, tpatt, name)
  ) %>% 
  filter(!grepl("Add", name)) %>% 
  inner_join(
    read_sas(unz(ipdname, boe_1288.18$file[49])) %>%
      rename_with(., tolower) %>%
      select(study, adp_pid, drgstpdt, termcdv)
  ) %>% 
  mutate(
    trttype = case_when(
      tpatt %in% c("D", "E") ~ "Exp",
      TRUE ~ "AC"
    ),
    randt_unit = "date",
    disct_unit = "date",
    reason = case_when(
      termcdv == 1 ~ "Completed",
      termcdv %in% c(2, 3, 4) ~ "Adverse Event",
      termcdv == 5 ~ "Lack of Efficacy",
      termcdv == 6 ~ "Protocol Deviation",
      termcdv == 7 ~ "Lost to Follow-up",
      termcdv == 8 ~ "Voluntary Withdrawal",
      termcdv == 9 ~ "Other",
      TRUE ~ "Missing"
    ),
    compl = if_else(reason == "Completed", 1, 0)
  ) %>%
  select(
    studyid = study, usubjid = adp_pid, 
    trtgrp = name, trttype, 
    randt = randdt,  randt_unit,
    disct = drgstpdt, disct_unit,
    compl,
    reason_orig = termcdv,
    reason
  )

table(boe_1288.18_attr$reason, useNA = "ifany")
table(boe_1288.18_attr$trtgrp, boe_1288.18_attr$reason, useNA = "ifany")

# Save output ------------------------------------------------------------------

# Get dataframe names and add to list
get_df_names <- Filter(function(x) grepl("_attr$", x), ls())
boe_lst <- lapply(get_df_names, function(i) get(i))

# Format all variables in each dataframe as characters
boe_lst_format <- lapply(boe_lst, function(i) i %>% mutate(across(everything(), as.character)))

# Bind together
boe_ready <- bind_rows(boe_lst_format)

# Final checks
length(unique(boe_ready$studyid)) # 86
sum(duplicated(boe_ready$usubjid)) # 566 (same IDs used between two studies)
glimpse(boe_ready)

# Save
write_rds(boe_ready, "Survival_attrition/Data/boe.rds")
