# Data Loading Script 
# Load Libraries 
library(data.table)

Sys.setlocale("LC_ALL", "Hebrew")
# Define hebrew and open access connections
df <- fread("input/Corona_vaccin4_2022-04-12_Covid.csv", encoding = "UTF-8")
df_redcap <- fread("input/data_redcap.csv", encoding = "UTF-8")
df_cci <- fread("input/CCI_covid_2nd_booster_national_0404.csv", encoding = "UTF-8")
df_serology <- as.data.table(readxl::read_excel("input/dbo_Req830_LabResultYonat_29_03.xlsx"))
df_index <- as.data.table(readxl::read_excel("input/dbo_Req830_Indx.xlsx"))
df_percent_positive_israel <- as.data.table(readxl::read_xlsx("input/percent_positive.xlsx"))

df_all_markers_vacc <- as.data.table(readxl::read_xlsx("input/group_ptid_single_antigens_02_06_22.xlsx", sheet = "vaccinated", na = c("mid")))
df_all_markers_unvacc <- as.data.table(readxl::read_xlsx("input/group_ptid_single_antigens_02_06_22.xlsx", sheet = "unvaccinated", na = c("mid")))
df_all_markers_vacc_with_mid <- as.data.table(readxl::read_xlsx("input/group_ptid_single_antigens_02_06_22.xlsx", sheet = "vaccinated"))
df_groups_shosh <- as.data.table(readxl::read_xlsx("input/75_samples_low_high_Baseline_and_M1.xlsx"))

df_mag_marker_vacc <- as.data.table(readxl::read_xlsx("input/group_ptid_02_06_22.xlsx", sheet = "vaccinated", na = c("mid")))
df_mag_marker_unvacc <- as.data.table(readxl::read_xlsx("input/group_ptid_02_06_22.xlsx", sheet = "unvaccinated", na = c("mid")))


