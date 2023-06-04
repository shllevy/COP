# Load and prepare redcap data
# add hyphen to 0303 center
fix_num_dash <-
  df_redcap$`Record ID`[!df_redcap$`Record ID` %like% "-"]
df_redcap[`Record ID` %in% fix_num_dash, `:=`(`Record ID` = paste(substr(fix_num_dash, 1, 4), substr(fix_num_dash, 5, 7), sep = '-'))]
# add zero to hospital code 0101 center
fix_num_101 <-
  df_redcap$`Record ID`[df_redcap$`Record ID` %like% "101-"]
df_redcap[`Record ID` %in% fix_num_101, `:=`(`Record ID` = paste0("0", `Record ID`))]
# add zero before subject's number
fix_num_0202 <-
  df_redcap$`Record ID`[df_redcap$`Record ID` %like% "0202"]
df_redcap[`Record ID` %in% fix_num_0202, `:=`(`Record ID` = paste(substr(fix_num_0202, 1, 5), substr(fix_num_0202, 6, 7), sep = '0'))]
# ----------------------
# Save medication analysis data frame under df_meds_history
df_meds_cols <- df_redcap[, `1. Name of medication`:`20. Stop date`]
df_meds_history <-
  df_redcap[`Event Name` == "Concomitant Medication", .SD, .SDcols = c("Record ID", colnames(df_meds_cols))]
# Remove medication columns for later analysis
df_redcap[, (colnames(df_meds_cols)) := NULL]
# ----------------------
# Fix column names in redcap data frame
col_name_dupl <- colnames(df_redcap)
col_name_dupl <- col_name_dupl[duplicated(col_name_dupl)]
col_name_dupl <- unique(col_name_dupl)
for (i in col_name_dupl) {
  name = i
  cols <- colnames(df_redcap) == name
  names(df_redcap)[cols] <- paste0(name, "_" , seq.int(sum(cols)))
}
col_name_fix <- colnames(df_redcap)
col_name_fix <- fix_col_names(col_name_fix)
setnames(df_redcap, colnames(df_redcap), col_name_fix)
setnames(df_redcap, "record_id", "num")
# ----------------------
# add number to duplicated columns after corrections
col_name_dupl <- colnames(df_redcap)
col_name_dupl <- col_name_dupl[duplicated(col_name_dupl)]
for (i in col_name_dupl) {
  name = i
  cols <- colnames(df_redcap) == name
  names(df_redcap)[cols] <- paste0(name, "_" , seq.int(sum(cols)))
}
# ----------------------
# Create data frame for each redcap questionnaire and choose relevant coulmns 
df_redcap_enroll <-
  df_redcap[event_name == "Vaccination-Enrollment" &
              repeat_instrument != "ICF"]
df_redcap_first_month <- df_redcap[event_name == "1 Month"]
df_redcap_second_month <- df_redcap[event_name == "2 Month"]
# ----------------------
colnames(df_redcap_first_month) <-
  paste(colnames(df_redcap_first_month), "first_month", sep = "_")
setnames(df_redcap_first_month, "num_first_month", "num")
df_redcap_first_month[, `:=`(fill_form_first_month = 1)]
df_redcap_first_month <- df_redcap_first_month[, c(
  "num",
  "fill_form_first_month",
  "have_you_been_diagnosed_positive_to_covid_19_since_last_visit_first_month",
  "diagnosis_date_first_month",
  "did_you_receive_the_second_booster_4th_dose_since_last_visit_first_month",
  "vaccination_date_first_month",
  "did_you_start_treatment_with_immunosuppressant_therapy_nsaid_biological_treatment_chemotherapy_or_other_first_month",
  "visit_date_3_first_month"
)]
colnames(df_redcap_second_month) <-
  paste(colnames(df_redcap_second_month), "second_month", sep = "_")
setnames(df_redcap_second_month, "num_second_month", "num")
df_redcap_second_month[, `:=`(fill_form_second_month = 1)]
df_redcap_second_month <- df_redcap_second_month[, c(
  "num",
  "fill_form_second_month",
  "have_you_been_diagnosed_positive_to_covid_19_since_last_visit_second_month",
  "diagnosis_date_second_month",
  "did_you_receive_the_second_booster_4th_dose_since_last_visit_second_month",
  "vaccination_date_second_month",
  "did_you_start_treatment_with_immunosuppressant_therapy_nsaid_biological_treatment_chemotherapy_or_other_second_month",
  "visit_date_3_second_month"
)]
# ----------------------


# Fix columns' names for main data frame 
col_name_fix <- colnames(df)
col_name_fix <- fix_col_names(col_name_fix)
setnames(df, colnames(df), col_name_fix)
# ----------------------

# Set names for standardization and merging 
setnames(
  df,
  c(
    "reference_event_column1",
    "birth_date",
    "reference_event_age_at_event",
    "deathdate_deceased_date",
    "bmi_bmi",
    "bmi_measurement_date_copy_days_from_reference",
    "covid_positive_after_enrollment_first_positive_result_date_days_from_reference",
    "covid_positive_first_ever_first_positive_result_date_days_from_reference",
    "membership_membership_type"
  ),
  c(
    "num",
    "birth_year",
    "age",
    "death_date",
    "bmi",
    "bmi_days_from_measurment",
    "covid_positive_after_enroll_days",
    "covid_positive_first_ever_days",
    "clalit_member_1"
  )
)
# ----------------------
# Correct manualy subject's num
df[num == "0303-0080", `:=` (num = "0303-080")]
# Delete all white spaces in num for correct merging 
df[, `:=` (num = str_replace_all(num, " ", ""))]
df_redcap_enroll[, `:=` (num = str_replace_all(num, " ", ""))]
df_redcap_first_month[, `:=` (num = str_replace_all(num, " ", ""))]
df_redcap_second_month[, `:=` (num = str_replace_all(num, " ", ""))]
# Define working data frame
df_work <-
  df[, c(
    "num",
    "birth_year",
    "age",
    "gender",
    "socioeconomic_score_five_level_scale",
    "country_of_birth",
    "bmi",
    "bmi_days_from_measurment",
    "covid_positive_after_enroll_days",
    "covid_positive_first_ever_days",
    "covid_pcrlab_count_count",
    "clalit_member_1",
    "admission_after_enrollment_admission_start_date"
  )]
# ---------------
# Add two missing individuals manually they will receive the data from redcap 
df_temp_3 <- df_work[nrow(df_work) + 1L]
df_temp_4 <- df_work[nrow(df_work) + 1L]
l = list(df_temp_3[, `:=`(num = "0404-274")], df_temp_4[, `:=`(num = "0404-343")], df_work)
df_work <- rbindlist(l)
# ---------------
# merge redcap enrollment data and monthes followup to main data frame
df_work <- df_redcap_enroll[, c(
  "num",
  "ethnicity",
  "sector_of_occuppation",
  "daily_interaction_with_corona_subjects",
  "smoking",
  "have_you_been_diagnosed_positive_to_covid_19_since_last_visit",
  "dose_number_1",
  "dose_number_2",
  "dose_number_3",
  "dose_number_4",
  "vaccination_date"
)][df_work, on = "num"]
df_work <- df_redcap_first_month[df_work, on = "num"]
df_work <- df_redcap_second_month[df_work, on = "num"]
# ---------------
# Complete dose number 4 from visits 1 or 2
df_work[, `:=`(dose_number_4 = fifelse(
  !is.na(vaccination_date_first_month),
  vaccination_date_first_month,
  dose_number_4
))]
df_work[, `:=`(dose_number_4 = fifelse(
  !is.na(vaccination_date_second_month),
  vaccination_date_second_month,
  dose_number_4
))]
# ---------------
# Prepare and merge serology data
setnames(
  df_serology,
  colnames(df_serology),
  c(
    "num",
    "test_num",
    "test_date" ,
    "test_name",
    "result",
    "orderd_by"
  )
)
df_serology[, `:=`(test_name = fifelse(test_name == "COV-2IgGII", "Covid19IgG ArchiteS1", test_name))]
df_serology <- df_serology[!is.na(result) & !is.na(num)]
df_serology[!test_name == "CoV-2 IgG BioPlex", `:=`(result = gsub("[^0-9.]", "", result))]
setorder(df_serology, num, test_date)
df_serology$run_no <-
  data.table::rowid(df_serology$num, df_serology$test_name)
df_serology_n <-
  dcast(df_serology,
        num ~ run_no + test_name,
        value.var = c("test_date", "result"))
col_name_fix <- colnames(df_serology_n)
col_name_fix <- fix_col_names(col_name_fix)
setnames(df_serology_n, colnames(df_serology_n), col_name_fix)
df_work <- df_serology_n[df_work, on = "num"]
# ---------------
# Add eventdate and remarks columns 
col_name_fix_df_index <- tolower(colnames(df_index))
setnames(df_index, colnames(df_index), col_name_fix_df_index)
df_index[num == "0303-0080", `:=` (num = "0303-080")]
df_work <- df_index[df_work, on = "num"]
# ---------------
setnames(
  df_work,
  c("visit_date_3_first_month",
    "visit_date_3_second_month"),
  c("visit_date_first_month",
    "visit_date_second_month")
)
# Complete missing data on age and gender from redcap
df_redcap_enroll[, `:=`(age_redcap = age, gender_redcap = sex)]
df_work <-
  df_redcap_enroll[, .(num, age_redcap, gender_redcap)][df_work, on = "num"]
df_work[is.na(age), `:=`(age = age_redcap)]
df_work[is.na(gender), `:=`(gender = gender_redcap)]

# Redefine correct class for bioplex results
serology_col_names_bioplex <- colnames(df_work)[colnames(df_work) %like% "bioplx"]
serology_col_names_alinity <- colnames(df_work)[colnames(df_work) %like% "archites"]
column_convert_numeric_loop(serology_col_names_bioplex, df_work)
column_convert_numeric_loop(serology_col_names_alinity, df_work)
# ----------------
# Add days to covid diagnosis based on subjects reports 
df_work[is.na(covid_positive_after_enroll_days), `:=`(
  covid_positive_after_enroll_days = time_length(eventdate %--% diagnosis_date_first_month, "day")
)]
df_work[is.na(covid_positive_after_enroll_days), `:=`(
  covid_positive_after_enroll_days = time_length(eventdate %--% diagnosis_date_second_month, "day")
)]
df_work[num == "0404-150", `:=`(
  covid_positive_after_enroll_days = time_length(eventdate %--% diagnosis_date_first_month, "day")
)]
df_work[num == "0404-357", `:=`(covid_positive_after_enroll_days = time_length(eventdate %--% ymd("2022-01-27"), "day"))]
# ----------------
# Follow up time from event date
df_work[!is.na(dose_number_4), `:=`(vacc_four_after_enroll_days = time_length(eventdate %--% dose_number_4, "day"))]

###### 
# Define vaccinated after
# 30 days post enrollment
# as unvaccinated
df_work[vacc_four_after_enroll_days > 30, `:=`(dose_number_4 = NA)]
######
df_work[!is.na(dose_number_4),
        `:=`(covid_positive_after_enroll_days = covid_positive_after_enroll_days - vacc_four_after_enroll_days - 7)]
df_work[, `:=`(time_to_event = covid_positive_after_enroll_days)]
df_work[, `:=`(event = fifelse(is.na(covid_positive_after_enroll_days), 0, 1))]
df_work[, `:=`(vaccinated = fifelse(is.na(dose_number_4), "No", "Yes"))]
df_work[is.na(time_to_event) &
          clalit_member_1 == 1 &
          vaccinated == "No",
        `:=`(time_to_event = time_length(eventdate %--% date("2022-04-03"), "day"))]
df_work[is.na(time_to_event) &
          clalit_member_1 == 1  &
          vaccinated == "Yes",
        `:=`(time_to_event = time_length(dose_number_4 %--% date("2022-04-03"), "day") - 7)]
df_work[is.na(time_to_event) &
          is.na(clalit_member_1) &
          vaccinated == "No",
        `:=`(time_to_event = time_length(eventdate %--% visit_date_second_month, "day"))]
df_work[is.na(time_to_event) &
          is.na(clalit_member_1) &
          vaccinated == "Yes",
        `:=`(time_to_event = time_length(dose_number_4 %--% visit_date_second_month, "day") - 7)]
df_work[is.na(time_to_event) &
          is.na(clalit_member_1) &
          vaccinated == "No",
        `:=`(time_to_event = time_length(eventdate %--% visit_date_first_month, "day"))]
df_work[is.na(time_to_event) &
          is.na(clalit_member_1) &
          vaccinated == "Yes",
        `:=`(time_to_event = time_length(dose_number_4 %--% visit_date_first_month, "day") - 7)]
# ----------------------
##### 
# Right censor at
# vaccination date for
# individuals who were
# vaccinated after 30 FU days
df_work[vacc_four_after_enroll_days > 30 , `:=`(time_to_event = vacc_four_after_enroll_days)]
##### 
# Calendar scale all period
df_work[vaccinated == "No", `:=`(time_since_first_enroll = time_length(date("2022-01-06") %--% eventdate, "day"))]
df_work[vaccinated == "Yes", `:=`(time_since_first_enroll = time_length(date("2022-01-06") %--% dose_number_4, "day") + 7)]
df_work[, `:=`(time_event_calendar = time_since_first_enroll + time_to_event)]
# Calendar scale first 30 days
df_work[, `:=`(time_to_event_30 = fifelse(time_to_event < 31 &
                                            event == 1, time_to_event, 30))]
df_work[, `:=`(event_30 = fifelse(time_to_event < 31  &
                                    event == 1, 1, 0))]
df_work[, `:=`(
  time_event_calendar_30 = fifelse(
    time_to_event < 31 & event == 1,
    time_since_first_enroll + time_to_event,
    time_since_first_enroll + 30
  )
)]
df_work[, `:=`(low_high_S2 = fifelse(result_1_cov_2_s2_igg_bioplx >= 10, "high", "low"))]
df_work[, `:=`(low_high_alinity = fifelse(result_1_covid19igg_archites1 >= 5000, "high", "low"))]
df_work[, `:=`(Positive_Covid = fifelse(covid_positive_after_enroll_days >= 0, TRUE, NA))]
df_work[, `:=`(Positive_Covid_name = fifelse(!is.na(Positive_Covid), "Yes", "No"))]
df_work[, `:=`(positive_covid_30_days = fifelse(
  covid_positive_after_enroll_days > 30 |
    is.na(covid_positive_after_enroll_days),
  0,
  1
))]
df_work[, `:=`(positive_covid_60_days = fifelse(
  covid_positive_after_enroll_days > 60 |
    is.na(covid_positive_after_enroll_days),
  0,
  1
))]
df_work[, `:=`(Positive_Covid_name_30_days = fifelse(positive_covid_30_days == 1, "Yes", "No"))]
df_work[, `:=`(Positive_Covid_name_30_days = factor(
  Positive_Covid_name_30_days,
  levels = c("No", "Yes"),
  labels = c("No", "Yes")
))]
df_work[, `:=`(Positive_Covid_name = factor(
  Positive_Covid_name,
  levels = c("No", "Yes"),
  labels = c("No", "Yes")
))]
df_work[, `:=`(Positive_Covid_num = fifelse(Positive_Covid_name == "Yes", 1, 0))]
df_work[, `:=`(time_from_third_vaccine = time_length(dose_number_3 %--% eventdate, "day"))]
df_work[, `:=`(
  ethnicity = fifelse(ethnicity == "", "Other or Not Mentioned", ethnicity),
  sector_of_occuppation = fifelse(
    sector_of_occuppation  == "" | is.na(sector_of_occuppation),
    "Other or Not Mentioned",
    sector_of_occuppation
  ),
  daily_interaction_with_corona_subjects = fifelse(
    daily_interaction_with_corona_subjects  == "" |
      is.na(daily_interaction_with_corona_subjects),
    "Not Mentioned",
    daily_interaction_with_corona_subjects
  )
)]
df_work[, `:=`(
  ethnicity = factor(
    ethnicity,
    labels = c("Jewish", "Bedouin", "Other or Not Mentioned"),
    levels = c("Jewish", "Bedouin", "Other or Not Mentioned")
  ),
  socioeconomic_score_five_level_scale = factor(
    socioeconomic_score_five_level_scale ,
    labels = c("Very High", "High", "Medium", "Low", "Very Low", "no data"),
    levels = c("Very High", "High", "Medium", "Low", "Very Low", "no data")
  ),
  sector_of_occuppation_grp = fifelse(
    sector_of_occuppation == "Nurses" |
      sector_of_occuppation == "Doctors",
    "Physicians or Nurses",
    "Administrative or Support staff"
  ),
  daily_interaction_with_corona_subjects   = factor(
    daily_interaction_with_corona_subjects ,
    labels = c("Yes", "No", "Not Mentioned"),
    levels = c("Yes", "No", "Not Mentioned"),
  ),
  vaccinated   = factor(
    vaccinated ,
    labels = c("No", "Yes"),
    levels = c("No", "Yes")
  )
)]
df_work[, `:=`(
  daily_interaction_with_corona_subjects_grp = fifelse(
    daily_interaction_with_corona_subjects == "Yes",
    "Yes",
    "Other"
  )
)]
df_work[, `:=`(age_grp = cut(
  age,
  breaks = c(-Inf, 34, 49, 64, Inf),
  labels = c("18-34", "35-49", "50-64", "65+")
))]
df_work[, `:=`(sex = gender)]
df_work[, `:=`(hospital = str_split_fixed(num, "-", 2)[, 1])]
df_work[hospital == "0101", `:=`(hospital = "Emek")]
df_work[hospital == "0202", `:=`(hospital = "Meir")]
df_work[hospital == "0303", `:=`(hospital = "Carmel")]
df_work[hospital == "0404", `:=`(hospital = "Soroka")]
df_work[, `:=`(
  sector_of_occuppation_n = fifelse(
    !sector_of_occuppation == "Doctors" &
      !sector_of_occuppation == "Nurses",
    "Administration and support staff",
    sector_of_occuppation
  )
)]
df_work[sector_of_occuppation_n == "Doctors", `:=`(sector_of_occuppation_n = "Physician")]
df_work[sector_of_occuppation_n == "Nurses", `:=`(sector_of_occuppation_n = "Nurse")]
df_work[, `:=`(sector_of_occuppation_n = factor(
  sector_of_occuppation_n ,
  labels = c("Physician", "Nurse", "Administration and support staff"),
  levels = c("Physician", "Nurse", "Administration and support staff")
))]
# Time in the study calculation:
df_work[clalit_member_1 == 1 &
          vaccinated == "No", `:=`(time_FU = time_length(eventdate %--% date("2022-04-03"), "day"))]
df_work[is.na(time_FU) &
          clalit_member_1 == 1 &
          vaccinated == "Yes", `:=`(time_FU = time_length(dose_number_4 %--% date("2022-04-03"), "day") - 7)]
df_work[is.na(time_FU) &
          clalit_member_1 == 0 &
          vaccinated == "No", `:=`(time_FU = time_length(dose_number_4 %--% visit_date_second_month, "day"))]
df_work[is.na(time_FU) &
          clalit_member_1 == 0 &
          vaccinated == "Yes", `:=`(time_FU = time_length(dose_number_4 %--% visit_date_second_month, "day") - 7)]

df_work[, `:=`(time_from_third_vaccine_month = round(time_from_third_vaccine /
                                                       30.417, digit = 0))]

setcolorder(
  df_work,
  c(
    "num",
    "age",
    "gender",
    "ethnicity",
    "sector_of_occuppation",
    "socioeconomic_score_five_level_scale",
    "country_of_birth",
    "dose_number_1",
    "dose_number_2",
    "dose_number_3",
    "dose_number_4",
    "covid_positive_after_enroll_days",
    "visit_date_first_month",
    "visit_date_second_month",
    "clalit_member_1"
  )
)
df_work[, `:=`(
  eventdate = ymd(eventdate),
  dose_number_1 = ymd(dose_number_1),
  dose_number_2 = ymd(dose_number_2),
  dose_number_3 = ymd(dose_number_3),
  dose_number_4 = ymd(dose_number_4),
  visit_date_first_month = ymd(visit_date_first_month),
  visit_date_second_month = ymd(visit_date_second_month)
)]
df_work[, `:=`(sector_of_occuppation_grp = as.character(sector_of_occuppation_grp))]
delete_col_sample = colnames(df_work)[colnames(df_work) %like% "sample"]
df_work[, (delete_col_sample) := NULL]


## Last Filters
df_work <-
  df_work[!is.na(dose_number_1) &
            !is.na(dose_number_2) &
            !is.na(dose_number_3)]
df_work <- df_work[!is.na(time_to_event)]

df_work <- df_work[time_to_event > 0]
df_work <- df_work[gender == "", `:=`(gender = NA)]



df_work <- df_mag_marker_vacc[df_work, on = "num"]
df_work <- df_mag_marker_unvacc[df_work, on = "num"]

# Prepare groups without mid
df_work[,`:=`(s2_group_unvaccinated_no_mid = s2_group_unvaccinated)]
df_work[, `:=`(
  s2_group_unvaccinated_no_mid = factor(
    s2_group_unvaccinated_no_mid,
    levels = c("high", "low"),
    labels = c("IgG BioPlex S2 High", "IgG BioPlex S2 Low")
  )
)]
df_work[,`:=`(alinity_group_unvaccinated_no_mid = alinity_group_unvaccinated)]
df_work[, `:=`(
  alinity_group_unvaccinated_no_mid = factor(
    alinity_group_unvaccinated_no_mid,
    levels = c("high", "low"),
    labels = c("IgG Alinity RBD High", "IgG Alinity RBD Low")
  )
)]
df_work[, `:=`(Variants_magnitude_response_group_unvaccinated_no_mid = Variants_magnitude_response_group_unvaccinated)]
df_work[Variants_magnitude_response_group_unvaccinated == "mid", `:=`(Variants_magnitude_response_group_unvaccinated_no_mid = NA)]
df_work[, `:=`(
  Variants_magnitude_response_group_unvaccinated_no_mid = factor(
    Variants_magnitude_response_group_unvaccinated_no_mid,
    levels = c("high", "low"),
    labels = c("IgA Variants High", "IgA Variants Low")
  )
)]
df_work[, `:=`(SARS_Cov_2_magnitude_response_group_unvaccinated_no_mid = SARS_Cov_2_magnitude_response_group_unvaccinated)]
df_work[SARS_Cov_2_magnitude_response_group_unvaccinated == "mid", `:=`(SARS_Cov_2_magnitude_response_group_unvaccinated_no_mid = NA)]
df_work[, `:=`(
  SARS_Cov_2_magnitude_response_group_unvaccinated_no_mid = factor(
    SARS_Cov_2_magnitude_response_group_unvaccinated_no_mid,
    levels = c("high", "low"),
    labels = c("IgA Wuhan High", "IgA Wuhan Low")
  )
)]
df_work[, `:=`(Mutants_magnitude_response_group_unvaccinated_no_mid = Mutants_magnitude_response_group_unvaccinated)]
df_work[Mutants_magnitude_response_group_unvaccinated == "mid", `:=`(Mutants_magnitude_response_group_unvaccinated_no_mid = NA)]
df_work[, `:=`(
  Mutants_magnitude_response_group_unvaccinated_no_mid = factor(
    Mutants_magnitude_response_group_unvaccinated_no_mid,
    levels = c("high", "low"),
    labels = c("IgG Mutants RBD High", "IgG Mutants RBD Low")
  )
)]


# Vaccinated section ------------------------------------------------
df_work[, `:=`(s2_group_vaccinated_no_mid = factor(
  s2_group_vaccinated,
  levels = c("high", "low"),
  labels = c("IgG BioPlex S2 High", "IgG BioPlex S2 Low")
))]

df_work[, `:=`(alinity_group_vaccinated_no_mid = factor(
  alinity_group_vaccinated,
  levels = c("high", "low"),
  labels = c("IgG Alinity RBD High", "IgG Alinity RBD Low")
))]

df_work[, `:=`(Variants_magnitude_response_group_vaccinated_no_mid = Variants_magnitude_response_group_vaccinated)]
df_work[Variants_magnitude_response_group_vaccinated_no_mid == "mid", `:=`(Variants_magnitude_response_group_vaccinated_no_mid = NA)]
df_work[, `:=`(
  Variants_magnitude_response_group_vaccinated_no_mid = factor(
    Variants_magnitude_response_group_vaccinated_no_mid,
    levels = c("high", "low"),
    labels = c("IgA Variants High", "IgA Variants Low")
  )
)]
df_work[, `:=`(SARS_Cov_2_magnitude_response_group_vaccinated_no_mid = SARS_Cov_2_magnitude_response_group_vaccinated)]
df_work[SARS_Cov_2_magnitude_response_group_vaccinated == "mid", `:=`(SARS_Cov_2_magnitude_response_group_vaccinated_no_mid = NA)]
df_work[, `:=`(
  SARS_Cov_2_magnitude_response_group_vaccinated_no_mid = factor(
    SARS_Cov_2_magnitude_response_group_vaccinated_no_mid,
    levels = c("high", "low"),
    labels = c("IgA Wuhan High", "IgA Wuhan Low")
  )
)]
df_work[, `:=`(Mutants_magnitude_response_group_vaccinated_no_mid = Mutants_magnitude_response_group_vaccinated)]
df_work[Mutants_magnitude_response_group_vaccinated == "mid", `:=`(Mutants_magnitude_response_group_vaccinated_no_mid = NA)]
df_work[, `:=`(
  Mutants_magnitude_response_group_vaccinated_no_mid = factor(
    Mutants_magnitude_response_group_vaccinated_no_mid,
    levels = c("high", "low"),
    labels = c("IgG Mutants RBD High", "IgG Mutants RBD Low")
  )
)]





# Create combinations

# High Low alinity_variants_iga_unvacc 90 days
df_work[, `:=`(
  alinity_variants_iga_comb_unvacc = paste0(
    Variants_magnitude_response_group_unvaccinated,
    alinity_group_unvaccinated
  )
)]
df_work <-
  df_work[alinity_variants_iga_comb_unvacc == "lowhigh" |
            alinity_variants_iga_comb_unvacc == "highlow", `:=`(alinity_variants_iga_comb_unvacc = NA)]
df_work[, `:=`(
  alinity_variants_iga_comb_unvacc = factor(
    alinity_variants_iga_comb_unvacc,
    levels = c("highhigh", "lowlow"),
    labels = c("IgG Alinity RBD & IgA Variants High", "IgG Alinity RBD & IgA Variants Low")
  )
)]
# High Low s2_variants_iga_unvaccinated 90 days

df_work[, `:=`(
  s2_variants_iga_comb_unvacc = paste0(
    Variants_magnitude_response_group_unvaccinated,
    s2_group_unvaccinated
  )
)]
df_work <-
  df_work[s2_variants_iga_comb_unvacc == "lowhigh" |
            s2_variants_iga_comb_unvacc == "highlow", `:=`(s2_variants_iga_comb_unvacc = NA)]
df_work[, `:=`(
  s2_variants_iga_comb_unvacc = factor(
    s2_variants_iga_comb_unvacc,
    levels = c("highhigh", "lowlow"),
    labels = c("BioPlex IgG S2 & IgA Variants High", "BioPlex IgG S2 & IgA Variants Low")
  )
)]

# High Low sars_cov_2_variants_iga_unvaccinated 90 days

df_work[, `:=`(
  sars_cov_2_variants_iga_comb_unvacc = paste0(
    Variants_magnitude_response_group_unvaccinated,
    SARS_Cov_2_magnitude_response_group_unvaccinated
  )
)]
df_work <-
  df_work[sars_cov_2_variants_iga_comb_unvacc == "lowhigh" |
            sars_cov_2_variants_iga_comb_unvacc == "highlow", `:=`(sars_cov_2_variants_iga_comb_unvacc = NA)]
df_work[, `:=`(
  sars_cov_2_variants_iga_comb_unvacc = factor(
    sars_cov_2_variants_iga_comb_unvacc,
    levels = c("highhigh", "lowlow"),
    labels = c("IgA Wuhan & IgA Variants High", "IgA Wuhan & IgA Variants Low")
  )
)]

# High Low mutants_variants_iga_unvaccinated 90 days

df_work[, `:=`(
  Mutants_variants_iga_comb_unvacc = paste0(
    Variants_magnitude_response_group_unvaccinated,
    Mutants_magnitude_response_group_unvaccinated
  )
)]
df_work <-
  df_work[Mutants_variants_iga_comb_unvacc == "lowhigh" |
            Mutants_variants_iga_comb_unvacc == "highlow", `:=`(Mutants_variants_iga_comb_unvacc = NA)]
df_work[, `:=`(
  Mutants_variants_iga_comb_unvacc = factor(
    Mutants_variants_iga_comb_unvacc,
    levels = c("highhigh", "lowlow"),
    labels = c("IgG Mutants RBD & IgA Variants High", "IgG Mutants RBD & IgA Variants Low")
  )
)]
# High Low s2_alinity_unvaccinated 90 days

df_work[, `:=`(
  s2_alinity_comb_unvacc = paste0(
    alinity_group_unvaccinated,
    s2_group_unvaccinated
  )
)]
df_work[, `:=`(
  s2_alinity_comb_unvacc = factor(
    s2_alinity_comb_unvacc,
    levels = c("highhigh", "lowlow"),
    labels = c("BioPlex IgG S2 & IgG Alinity RBD High", "BioPlex IgG S2 & IgG Alinity RBD Low")
  )
)]
# High Low s2_wuhan_unvaccinated 90 days

df_work[, `:=`(
  s2_wuhan_comb_unvacc = paste0(
    SARS_Cov_2_magnitude_response_group_unvaccinated,
    s2_group_unvaccinated
  )
)]
df_work[, `:=`(
  s2_wuhan_comb_unvacc = factor(
    s2_wuhan_comb_unvacc,
    levels = c("highhigh", "lowlow"),
    labels = c("BioPlex IgG S2 & IgA Wuhan High", "BioPlex IgG S2 & IgA Wuhan Low")
  )
)]

# High Low s2_mutants_unvaccinated 90 days

df_work[, `:=`(
  s2_mutants_comb_unvacc = paste0(
    Mutants_magnitude_response_group_unvaccinated,
    s2_group_unvaccinated
  )
)]
df_work[, `:=`(
  s2_mutants_comb_unvacc = factor(
    s2_mutants_comb_unvacc,
    levels = c("highhigh", "lowlow"),
    labels = c("BioPlex IgG S2 & IgG Mutants RBD High", "BioPlex IgG S2 & IgG Mutants RBD Low")
  )
)]
# High Low alinity Wuhan_unvaccinated 90 days

df_work[, `:=`(
  wuhan_alinity_comb_unvacc = paste0(
    alinity_group_unvaccinated,
    SARS_Cov_2_magnitude_response_group_unvaccinated
  )
)]
df_work[, `:=`(
  wuhan_alinity_comb_unvacc = factor(
    wuhan_alinity_comb_unvacc,
    levels = c("highhigh", "lowlow"),
    labels = c("IgG Alinity RBD & IgA Wuhan High", "IgG Alinity RBD & IgA Wuhan Low")
  )
)]
# High Low alinity mutants unvaccinated 90 days

df_work[, `:=`(
  mutants_alinity_comb_unvacc = paste0(
    alinity_group_unvaccinated,
    Mutants_magnitude_response_group_unvaccinated
  )
)]
df_work[, `:=`(
  mutants_alinity_comb_unvacc = factor(
    mutants_alinity_comb_unvacc,
    levels = c("highhigh", "lowlow"),
    labels = c("IgG Alinity RBD & IgG Mutants RBD High", "IgG Alinity RBD & IgG Mutants RBD Low")
  )
)]

# High Low wuhan mutants_unvaccinated 90 days

df_work[, `:=`(
  mutants_wuhan_comb_unvacc = paste0(
    SARS_Cov_2_magnitude_response_group_unvaccinated,
    Mutants_magnitude_response_group_unvaccinated
  )
)]
df_work[, `:=`(
  mutants_wuhan_comb_unvacc = factor(
    mutants_wuhan_comb_unvacc,
    levels = c("highhigh", "lowlow"),
    labels = c("IgG Mutants RBD & IgA Wuhan High", "IgG Mutants RBD & IgA Wuhan Low")
  )
)]

df_work[, `:=`(
  Mutants_variants_iga_comb_unvacc = paste0(
    Variants_magnitude_response_group_unvaccinated,
    Mutants_magnitude_response_group_unvaccinated
  )
)]
df_work <-
  df_work[Mutants_variants_iga_comb_unvacc == "lowhigh" |
            Mutants_variants_iga_comb_unvacc == "highlow", `:=`(Mutants_variants_iga_comb_unvacc = NA)]
df_work[, `:=`(
  Mutants_variants_iga_comb_unvacc = factor(
    Mutants_variants_iga_comb_unvacc,
    levels = c("highhigh", "lowlow"),
    labels = c("IgG Mutants RBD & IgA Variants High", "IgG Mutants RBD & IgA Variants Low")
  )
)]
# Vaccinated section 
# High Low s2_variants_iga_unvaccinated 90 days

df_work[, `:=`(
  s2_variants_iga_comb_vacc = paste0(
    Variants_magnitude_response_group_vaccinated,
    s2_group_vaccinated
  )
)]
df_work <-
  df_work[s2_variants_iga_comb_vacc == "lowhigh" |
            s2_variants_iga_comb_vacc == "highlow", `:=`(s2_variants_iga_comb_vacc = NA)]
df_work[, `:=`(
  s2_variants_iga_comb_vacc = factor(
    s2_variants_iga_comb_vacc,
    levels = c("highhigh", "lowlow"),
    labels = c("BioPlex IgG S2 & IgA Variants High", "BioPlex IgG S2 & IgA Variants Low")
  )
)]

df_work[, `:=`(
  Mutants_variants_iga_comb_vacc = paste0(
    Variants_magnitude_response_group_vaccinated,
    Mutants_magnitude_response_group_vaccinated
  )
)]
df_work <-
  df_work[Mutants_variants_iga_comb_vacc == "lowhigh" |
            Mutants_variants_iga_comb_vacc == "highlow", `:=`(Mutants_variants_iga_comb_vacc = NA)]
df_work[, `:=`(
  Mutants_variants_iga_comb_vacc = factor(
    Mutants_variants_iga_comb_vacc,
    levels = c("highhigh", "lowlow"),
    labels = c("IgG Mutants RBD & IgA Variants High", "IgG Mutants RBD & IgA Variants Low")
  )
)]

# High Low alinity_variants_iga_unvacc 90 days
df_work[, `:=`(
  alinity_variants_iga_comb_vacc = paste0(
    Variants_magnitude_response_group_vaccinated,
    alinity_group_vaccinated
  )
)]
df_work <-
  df_work[alinity_variants_iga_comb_vacc == "lowhigh" |
            alinity_variants_iga_comb_vacc == "highlow", `:=`(alinity_variants_iga_comb_vacc = NA)]
df_work[, `:=`(
  alinity_variants_iga_comb_vacc = factor(
    alinity_variants_iga_comb_vacc,
    levels = c("highhigh", "lowlow"),
    labels = c("IgG Alinity RBD & IgA Variants High", "IgG Alinity RBD & IgA Variants Low")
  )
)]

# High Low s2_variants_iga_vaccinated 90 days

df_work[, `:=`(
  s2_variants_iga_comb_vacc = paste0(
    Variants_magnitude_response_group_vaccinated,
    s2_group_vaccinated
  )
)]
df_work <-
  df_work[s2_variants_iga_comb_vacc == "lowhigh" |
            s2_variants_iga_comb_vacc == "highlow", `:=`(s2_variants_iga_comb_vacc = NA)]
df_work[, `:=`(
  s2_variants_iga_comb_vacc = factor(
    s2_variants_iga_comb_vacc,
    levels = c("highhigh", "lowlow"),
    labels = c("BioPlex IgG S2 & IgA Variants High", "BioPlex IgG S2 & IgA Variants Low")
  )
)]

# High Low sars_cov_2_variants_iga_vaccinated 90 days

df_work[, `:=`(
  sars_cov_2_variants_iga_comb_vacc = paste0(
    Variants_magnitude_response_group_vaccinated,
    SARS_Cov_2_magnitude_response_group_vaccinated
  )
)]
df_work <-
  df_work[sars_cov_2_variants_iga_comb_vacc == "lowhigh" |
            sars_cov_2_variants_iga_comb_vacc == "highlow", `:=`(sars_cov_2_variants_iga_comb_vacc = NA)]
df_work[, `:=`(
  sars_cov_2_variants_iga_comb_vacc = factor(
    sars_cov_2_variants_iga_comb_vacc,
    levels = c("highhigh", "lowlow"),
    labels = c("IgA Wuhan & IgA Variants High", "IgA Wuhan & IgA Variants Low")
  )
)]

# High Low mutants_variants_iga_vaccinated 90 days

df_work[, `:=`(
  Mutants_variants_iga_comb_vacc = paste0(
    Variants_magnitude_response_group_vaccinated,
    Mutants_magnitude_response_group_vaccinated
  )
)]
df_work <-
  df_work[Mutants_variants_iga_comb_vacc == "lowhigh" |
            Mutants_variants_iga_comb_vacc == "highlow", `:=`(Mutants_variants_iga_comb_vacc = NA)]
df_work[, `:=`(
  Mutants_variants_iga_comb_vacc = factor(
    Mutants_variants_iga_comb_vacc,
    levels = c("highhigh", "lowlow"),
    labels = c("IgG Mutants RBD & IgA Variants High", "IgG Mutants RBD & IgA Variants Low")
  )
)]

# High Low s2_alinity_unvaccinated 90 days


df_work[, `:=`(
  s2_alinity_comb_vacc = paste0(
    alinity_group_vaccinated,
    s2_group_vaccinated
  )
)]
df_work[, `:=`(
  s2_alinity_comb_vacc = factor(
    s2_alinity_comb_vacc,
    levels = c("highhigh", "lowlow"),
    labels = c("BioPlex IgG S2 & IgG Alinity RBD High", "BioPlex IgG S2 & IgG Alinity RBD Low")
  )
)]

# High Low s2_wuhan_vaccinated 90 days

df_work[, `:=`(
  s2_wuhan_comb_vacc = paste0(
    SARS_Cov_2_magnitude_response_group_vaccinated,
    s2_group_vaccinated
  )
)]
df_work[, `:=`(
  s2_wuhan_comb_vacc = factor(
    s2_wuhan_comb_vacc,
    levels = c("highhigh", "lowlow"),
    labels = c("BioPlex IgG S2 & IgA Wuhan High", "BioPlex IgG S2 & IgA Wuhan Low")
  )
)]

# High Low s2_mutants_vaccinated 90 days

df_work[, `:=`(
  s2_mutants_comb_vacc = paste0(
    Mutants_magnitude_response_group_vaccinated,
    s2_group_vaccinated
  )
)]
df_work[, `:=`(
  s2_mutants_comb_vacc = factor(
    s2_mutants_comb_vacc,
    levels = c("highhigh", "lowlow"),
    labels = c("BioPlex IgG S2 & IgG Mutants RBD High", "BioPlex IgG S2 & IgG Mutants RBD Low")
  )
)]

# High Low alinity Wuhan_vaccinated 90 days

df_work[, `:=`(
  wuhan_alinity_comb_vacc = paste0(
    alinity_group_vaccinated,
    SARS_Cov_2_magnitude_response_group_vaccinated
  )
)]
df_work[, `:=`(
  wuhan_alinity_comb_vacc = factor(
    wuhan_alinity_comb_vacc,
    levels = c("highhigh", "lowlow"),
    labels = c("IgG Alinity RBD & IgA Wuhan High", "IgG Alinity RBD & IgA Wuhan Low")
  )
)]

# High Low alinity mutants vaccinated 90 days

df_work[, `:=`(
  mutants_alinity_comb_vacc = paste0(
    alinity_group_vaccinated,
    Mutants_magnitude_response_group_vaccinated
  )
)]
df_work[, `:=`(
  mutants_alinity_comb_vacc = factor(
    mutants_alinity_comb_vacc,
    levels = c("highhigh", "lowlow"),
    labels = c("IgG Alinity RBD & IgG Mutants RBD High", "IgG Alinity RBD & IgG Mutants RBD Low")
  )
)]

# High Low wuhan mutants_vaccinated 90 days

df_work[, `:=`(
  mutants_wuhan_comb_vacc = paste0(
    SARS_Cov_2_magnitude_response_group_vaccinated,
    Mutants_magnitude_response_group_vaccinated
  )
)]
df_work[, `:=`(
  mutants_wuhan_comb_vacc = factor(
    mutants_wuhan_comb_vacc,
    levels = c("highhigh", "lowlow"),
    labels = c("IgG Mutants RBD & IgA Wuhan High", "IgG Mutants RBD & IgA Wuhan Low")
  )
)]

