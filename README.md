# Dietary Vitamin K Intake and Serum Testosterone in U.S. Adult Men: NHANES 2011-2016
This repository contains the full R code used to reproduce the analyses presented in:

add citation once published


# Study Overview
Using survey-weighted analyses of 6,223 adult men from NHANES 2011–2016, this project evaluated associations between energy-adjusted dietary vitamin K intake and serum total testosterone.

Analyses included:
- survey-weighted multivariable linear regression analyses
- survey-weighted logistic regression analyses
- continuous exposure analyses
- trend analyses across vitamin K quartiles
- sensitivity and robustness analyses
- restricted cubic spline analyses to assess nonlinearity
- effect modification and interaction analyses
  

# Repository Contents

```text
.
├── .Rprofile
├── .gitignore
├── LICENSE
├── README.md
├── renv.lock                # Reproducible package versions
├── vitK_testo_nhanes_analysis.R        # Main analysis script
├── vitK_testo_nhanes_analysis.Rproj
│
├── renv/
│   ├── activate.R
│   └── settings.json
│
├── results/
│   └── sessionInfo.txt
│
└── data_raw/                # User-supplied NHANES data files
```

# Data Sources
NHANES datasets are publicly available from the CDC National Center for Health Statistics (NCHS): 
https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx

Healthy Eating Index (HEI-2020) scores were calculated using the dietaryindex R package:
https://jamesjiadazhan.github.io/dietaryindex_manual/articles/dietaryindex.html

The FPED files required for HEI-2020 calculation were obtained using the instructions provided in the dietaryindex documentation.


# Required Data Files
Users must manually download all required NHANES .xpt files and FPED .sas7bdat files into: 
```text
data_raw/
```
The analysis script assumes the exact original NHANES filenames.

Required file groups include:
- NHANES demographic files:
`DEMO_*.xpt`
- Dietary intake files:
`DR1TOT_*.xpt`
- Glycemic laboratory files:
`GLU_*.xpt`,
`INS_*.xpt`,
- Anthropometric files:
`BMX_*.xpt`
- Lifestyle questionnaire files:
`SMQ_*.xpt`,
`ALQ_*.xpt`,
`PAQ_*.xpt`
- FPED files for HEI-2020 calculation:
`fped_dr1tot_1112.sas7bdat`,
`fped_dr1tot_1314.sas7bdat`,
`fped_dr1tot_1516.sas7bdat`,
- Testosterone and SHBG laboratory files:
`TST_*.xpt`
- Examination session / fasting questionnaire files:
`FASTQX_*.xpt`
- Medical conditions questionnaire files (used for prostate cancer sensitivity analysis):
`MCQ_*.xpt`
- Prescription medication files (used to identify testosterone therapy):
`RXQ_RX_*.xpt`


# Example Required File Structure
The following image shows the expected contents of the data_raw/ directory.
<img width="1732" height="1106" alt="Screenshot 2026-05-26 at 07 26 10" src="https://github.com/user-attachments/assets/91624f7e-bb11-4650-a2e2-9bc823c683c5" />


# Reproducibility
This repository uses renv for reproducible package management.

To restore the exact R package environment used in the analyses:
```text
install.packages("renv")
renv::restore()
```

# Running the Analysis
1. Clone the repository
2. Download all required NHANES and FPED files into data_raw/
3. Restore the R environment:
```text
renv::restore()
```
4. Run:
```text
source("nhanes_analysis.R")
```
All output tables will be generated automatically in:
```text
results/
```

# Statistical Notes
Analyses used:
- 6-year MEC examination weights (`WTMEC6YR`)
- pooled survey-weight adjustment according to NHANES analytic guidelines
- stratification (`SDMVSTRA`)
- primary sampling units (`SDMVPSU`)
- Taylor series linearization variance estimation
- survey-weighted regression using the `survey` R package
To account for singleton strata after subsetting, analyses used:
```text
options(survey.lonely.psu = "adjust")
```

# Software Environment
The analysis was performed in R using a reproducible renv environment.
Detailed package/session information is available in:
```text
results/sessionInfo.txt
```

# License
This project is licensed under the MIT License.




