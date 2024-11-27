# CDISC Pilot Project GitHub: https://github.com/cdisc-org/sdtm-adam-pilot-project
# Submission GitHub: 
# Reference: https://www.bioinfo-scrounger.com/archives/sdtm_adam_cdisc/
# Reference: https://cran.r-project.org/web/packages/admiral/vignettes/adsl.html

#################################################################################


# 1. Import SDTM/ADAM datasets from R package "admiral" ("r2rtf" and "clinUtils" packages also contains CDISC pilot projects)

# 1.1 SDTM  (Import SDTM datasets from "admiral.test" package)
library(admiral.test)
data(admiral_ae)
ae <- admiral_ae
head(ae[1:5,1:5])


# 1.2 ADAM (Import ADAM dataset from "admiral" package)
library(admiral)
data(admiral_adsl)
adsl <- admiral_adsl
head(adsl[1:10,1:10])

##################################################################


# 2. Create ADAM

## 2.1 Create ADSL
## Reference:https://cran.r-project.org/web/packages/admiral/vignettes/adsl.html


### Read in SDTM
library(pharmaversesdtm) #SDTM Test Data
#library(pharmaverseadam) #ADAM Test Data
library(admiral)  #ADAM Test Data (LZ: Honestly, I don't know the difference between "pharmaverseadam" and "admiral", but could look it up later)
#library(dplyr, warn.conflicts = FALSE)
#library(lubridate)
#library(stringr)

data("dm")
data("ds")
data("ex")
data("ae")
data("lb")

# R does not take blank (missing character) cell or  "." (missing numbers), change blank strings to NA when applicable.
# convert_blanks_to_na() function is part of "admiral" package
dm <- convert_blanks_to_na(dm)  
ds <- convert_blanks_to_na(ds)
ex <- convert_blanks_to_na(ex)
ae <- convert_blanks_to_na(ae)
lb <- convert_blanks_to_na(lb)


# 2.1.1 Load DM as the basis of ADSL  "%>%" means "and then"
adsl <- dm %>%
  mutate(TRT01P = ARM, TRT01A = ACTARM) #mutation(): Add new derived variable while retain old variables
adsl


# 2.1.2 Derive/Impute Numeric Treatment Date/Time and Duration (TRTSDTM, TRTEDTM, TRTDURD)

# impute start and end time of exposure to first and last respectively,
# do not impute date
ex_ext <- ex %>%
  derive_vars_dtm(
    dtc = EXSTDTC,
    new_vars_prefix = "EXST"
  ) %>%
  derive_vars_dtm(
    dtc = EXENDTC,
    new_vars_prefix = "EXEN",
    time_imputation = "last"
  )

# Merge EX data into ADSL
## derive_vars_merged() can be used to derive the treatment start and end date/times using the ex domain. 
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 |
                    (EXDOSE == 0 &
                       str_detect(EXTRT, "PLACEBO"))) & !is.na(EXSTDTM),
    new_vars = exprs(TRTSDTM = EXSTDTM, TRTSTMF = EXSTTMF),
    order = exprs(EXSTDTM, EXSEQ),
    mode = "first",
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 |
                    (EXDOSE == 0 &
                       str_detect(EXTRT, "PLACEBO"))) & !is.na(EXENDTM),
    new_vars = exprs(TRTEDTM = EXENDTM, TRTETMF = EXENTMF),
    order = exprs(EXENDTM, EXSEQ),
    mode = "last",
    by_vars = exprs(STUDYID, USUBJID)
  )





use_ad_template("ADSL")








#LIBNAME SETUP
sdtm <- "//product/study/analysis/data/sdtm" # assign dir to object named sdtm
out <- "//ADAM"


READ FILES
dm <- read_sas(file.path(sdtm,"dm.sas7bdat")) # read sas file as a data frame
ds <- read_sas(file.path(sdtm,"ds.sas7bdat"))
sv <- read_sas(file.path(sdtm,"sv.sas7bdat"))
suppsv <- read_sas(file.path(sdtm,"suppsv.sas7bdat"))
MERGING PARENT AND SUPPLEMENTAL DATA SET
suppds_ <- suppds %>% # "%>%" read as "and then"
  mutate(idvarval_ = as.numeric(idvarval)) %>%
  select(usubjid,idvarval_,qnam,qval) %>%
  spread(.,qnam,qval)
ds_all <- left_join(ds, suppds_, 
                    by = c("usubjid"="usubjid", "dsseq"="idvarval_"))
BASELINE FLAG/CHANGE FROM BASELINE DERIVATION
advs <- advs_ %>%
  group_by(subjid,paramcd) %>%
  arrange(subjid,paramcd) %>%
  mutate ( base = aval[visitnum==1],
           ablfl = ifelse(visitnum == 1,"y",na),
           chg = ifelse ( !is.na(aval) & !is.na(base),aval-base,na),
           pchg = ifelse ( !is.na(aval) & !is.na(base),((aval- base)/base)*100,na)
           PhUSE US Connect 2018
           3
           SUBJECT LEVEL DERIVATION
           adsl <- dm %>% 
             select(studyid, subjid, age, sex, height, weight, race, scrfl) %>%
             mutate(bmi = (weight*703)/height^2 ) %>%
             filter(scrfl == “Y”) %>%
             select(-scrfl) %>%
             arrange(studyid, subjid)