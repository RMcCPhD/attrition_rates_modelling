
# Tribble for identifying plot thresholds
library(tibble)

thresholds <- tribble(
  ~cause,                    ~plot_id,                    ~threshold,
  
  # Adverse event
  "Adverse Event",           "BPH 03, n=606",             3,
  "Adverse Event",           "COPD 19, n=472",            6,
  "Adverse Event",           "Hyp 01, n=1039",            4,
  "Adverse Event",           "Hyp 08, n=426",             5,
  "Adverse Event",           "Park 04, n=475",            7,
  "Adverse Event",           "T2DM 04, n=1360",           4,
  "Adverse Event",           "T2DM 11, n=899",            5,
  "Adverse Event",           "T2DM 22, n=574",            6,
  "Adverse Event",           "T2DM 25, n=561",            7,
  "Adverse Event",           "T2DM 28, n=492",            5,
  "Adverse Event",           "T2DM 29, n=447",            7,
  "Adverse Event",           "T2DM 30, n=389",            5,
  "Adverse Event",           "T2DM 31, n=302",            7,
  "Adverse Event",           "T2DM 34, n=299",            7,
  
  # Lack of efficacy
  "Lack of Efficacy",        "BPH 03, n=606",             6,
  "Lack of Efficacy",        "COPD 03, n=2488",           1,
  "Lack of Efficacy",        "COPD 14, n=624",            6,
  "Lack of Efficacy",        "COPD 18, n=519",            6,
  "Lack of Efficacy",        "Hyp 07, n=490",             10,
  "Lack of Efficacy",        "T2DM 23, n=566",            10,
  
  # Lost to follow-up
  "Lost to Follow-up",       "COPD 04, n=1829",           2,
  "Lost to Follow-up",       "COPD 08, n=983",            5,
  "Lost to Follow-up",       "Hyp 04, n=860",             6,
  "Lost to Follow-up",       "T2DM 05, n=1303",           4,
  "Lost to Follow-up",       "T2DM 11, n=899",            6,
  "Lost to Follow-up",       "T2DM 16, n=791",            5,
  
  # Other/Miscellaneous
  "Other/Miscellaneous",     "COPD 10, n=906",            5,
  "Other/Miscellaneous",     "COPD 13, n=644",            5,
  "Other/Miscellaneous",     "T2DM 07, n=1162",           2.5,
  "Other/Miscellaneous",     "T2DM 14, n=807",            4,
  "Other/Miscellaneous",     "T2DM 35, n=272",            10,
  
  # PI/Sponsor Decision
  "PI/Sponsor Decision",     "BPH 01, n=1056",            2.5,
  
  # Protocol violation
  "Protocol Violation",      "BPH 04, n=511",             10,
  "Protocol Violation",      "Hyp 01, n=1039",            5,
  "Protocol Violation",      "T2DM 24, n=566",            5,
  
  # Voluntary withdrawal (your code had a copy-paste error here)
  "Voluntary Withdrawal",    "Hyp 04, n=860",             5,
  "Voluntary Withdrawal",    "T2DM 04, n=1360",           5,
  "Voluntary Withdrawal",    "T2DM 14, n=807",            4,
  "Voluntary Withdrawal",    "T2DM 35, n=272",            10
)