# This script clears clinical TCGA data

  insert_head()
  
# reading the relevant clinical data ------
  
  insert_msg('Rading the clinical text data')
  
  tcga_clinics <- c(drugs = './input data/TCGA/nationwidechildrens.org_clinical_drug_paad.txt', 
                    fup1 = './input data/TCGA/nationwidechildrens.org_clinical_follow_up_v4.4_nte_paad.txt', 
                    fup2 = './input data/TCGA/nationwidechildrens.org_clinical_nte_paad.txt', 
                    clinical = './input data/TCGA/nationwidechildrens.org_clinical_patient_paad.txt', 
                    radiation = './input data/TCGA/nationwidechildrens.org_clinical_radiation_paad.txt') %>% 
    map(read_tsv)
  
# Clearing the general clinical data -----
  
  insert_msg('Clearing the general clinical data')
  
  tcga_clinics$clinical <- tcga_clinics$clinical[c(-1, -2), ] %>% 
    mutate(patient_id = bcr_patient_barcode, 
           histology = car::recode(histologic_diagnosis, 
                                   "'Pancreas-Adenocarcinoma-Other Subtype' = 'Adeno other'; 
                                   'Pancreas-Adenocarcinoma Ductal Type' = 'Adeno ductal'; 
                                   'Pancreas-Undifferentiated Carcinoma' = 'Undifferentiated'; 
                                   'Pancreas-Colloid (mucinous non-cystic) Carcinoma' = 'Colloid'; 
                                   '[Discrepancy]' = NA") %>% 
             factor(c('Adeno ductal', 'Adeno other', 'Undifferentiated', 'Colloid')),
           histologic_diagnosis_other = ifelse(histologic_diagnosis_other %in% c('[Not Applicable]', 
                                                                                 '[Not Available]', 
                                                                                 'not specified'), 
                                               NA, 
                                               histologic_diagnosis_other), 
           collection_type = ifelse(prospective_collection == 'NO', 
                                    'retrospective', 
                                    'prospective') %>% 
             factor, 
           sex = car::recode(gender, 
                             "'MALE' = 'Male'; 
                             'FEMALE' = 'Female'") %>% 
             factor, 
           race = car::recode(race, 
                              "'WHITE' = 'White'; 
                              'BLACK OR AFRICAN AMERICAN' = 'Black';
                              'ASIAN' = 'Asian'") %>% 
             factor('White', 'Black', 'Asian'), 
           other_malignancy = car::recode(history_other_malignancy, 
                                          "'Yes, History of Prior Malignancy' = 'Yes'") %>% 
             factor, 
           neoadjuvant = car::recode(history_neoadjuvant_treatment, 
                                     "'Yes, Radiation Prior to Resection' = 'Yes'") %>% 
             factor, 
           surgical_procedure = car::recode(surgical_procedure, 
                                            "'Other Method (please specify)' = 'Other'") %>% 
             factor(c('Whipple', 
                      'Distal Pancreatectomy', 
                      'Other')), 
           tumor_grade = factor(tumor_grade), 
           tumor_size = as.numeric(tumor_resected_max_dimension), 
           residual_tumor = factor(residual_tumor, 
                                   c('R0', 'R1', 'RX')), 
           pt_stage = factor(ajcc_tumor_pathologic_pt, 
                             c('T1', 'T2', 'T3', 'T4', 'TX')), 
           pn_stage = factor(ajcc_nodes_pathologic_pn, 
                             c('N0', 'N1', 'N1b', 'NX')), 
           pm_stage = factor(ajcc_metastasis_pathologic_pm), 
           p_stage = stri_replace(ajcc_pathologic_tumor_stage, 
                                  fixed = 'Stage ', 
                                  replacement = '') %>% 
             factor(c('I', 'IA', 'IB', 'II', 'IIA', 'IIB', 'III', 'IIIA', 'IIIB', 'IV')), 
           death = car::recode(vital_status, 
                               "'Alive' = 0; 
                               'Dead' = 1") %>% 
             as.numeric, 
           last_contact_days_to = as.numeric(last_contact_days_to), 
           death_days_to = as.numeric(death_days_to), 
           os_days = ifelse(vital_status == 'Alive', 
                            last_contact_days_to, 
                            death_days_to), 
           death_panc = ifelse(vital_status == 'Alive' | cause_of_death != 'Pancreatic Cancer', 
                               0, 
                               1), 
           history_chronic_pancreatitis = car::recode(history_chronic_pancreatitis, 
                                                      "'NO' = 'No'; 
                                                      'YES' = 'Yes'") %>% 
             factor(c('No', 'Yes')), 
           cancer_family_history = car::recode(family_history_cancer_indicator, 
                                               "'NO' = 'No'; 
                                               'YES' = 'Yes'") %>% 
             factor(c('No', 'Yes')), 
           treatment_outcome = car::recode(treatment_outcome_first_course, 
                                           "'Progressive Disease' = 'Progression'; 
                                           'Stable Disease' = 'Stable'; 
                                           'Partial Remission/Response' = 'Partial remission'; 
                                           'Complete Remission/Response' = 'Complete remission'") %>% 
             factor(c('Progression', 
                      'Stable', 
                      'Partial remission', 
                      'Complete remission')), 
           treatment_outcome_score = car::recode(treatment_outcome_first_course, 
                                                 "'Progressive Disease' = 0; 
                                           'Stable Disease' = 1; 
                                           'Partial Remission/Response' = 2; 
                                           'Complete Remission/Response' = 3") %>% 
             as.numeric, 
           new_tumor = car::recode(new_tumor_event_dx_indicator, 
                                   "'NO' = 'No'; 
                                   'YES' = 'Yes'") %>% 
             factor(c('No', 'Yes')), 
           relapse = car::recode(new_tumor_event_dx_indicator, 
                                 "'NO' = 0; 
                                 'YES' = 1") %>% 
             as.numeric, 
           age = as.numeric(age_at_initial_pathologic_diagnosis), 
           anatomic_subdivision = car::recode(anatomic_neoplasm_subdivision, 
                                              "'Head of Pancreas' = 'Head'; 
                                              'Body of Pancreas' = 'Body'; 
                                              'Tail of Pancreas' = 'Tail'; 
                                              'Other (please specify)' = 'Other'") %>% 
             factor(c('Head', 'Body', 'Tail', 'Other'))) %>% 
    select(patient_id, 
           histology, 
           histologic_diagnosis_other, 
           collection_type, 
           age, 
           sex, 
           race, 
           other_malignancy, 
           neoadjuvant, 
           surgical_procedure, 
           tumor_grade, 
           tumor_size, 
           residual_tumor, 
           pt_stage, 
           pn_stage, 
           pm_stage, 
           p_stage, 
           death, 
           os_days, 
           death_panc, 
           history_chronic_pancreatitis, 
           cancer_family_history, 
           treatment_outcome, 
           treatment_outcome_score, 
           new_tumor, 
           relapse, 
           anatomic_subdivision)
  
  
# Clearing the drug data ----
  
  insert_msg('Clearing the drugs data')
  
  tcga_clinics$drugs <- tcga_clinics$drugs[c(-1, -2), ] %>% 
    mutate(patient_id = bcr_patient_barcode, 
           drug_type = car::recode(pharmaceutical_therapy_drug_name, 
                                   "'[Unknown]' = NA;
                                   '[Not Available]' = NA;
                                   'Chemo, NOS' = NA;
                                   'Megace' = NA; ## no tumor specific-treatment
                                   '5 FU' = '5-FU'; 
                                   'Fluorouracil' = '5-FU'; 
                                   '5-fu' = '5-FU'; 
                                   'FU7' = '5-FU'; 
                                   '5-Fluorouracil' = '5-FU'; 
                                   'fluorouracil' = '5-FU'; 
                                   '5-fluorouracil' = '5-FU'; 
                                   '5-Fluorouracil?' = '5-FU'; 
                                   '5FU' = '5-FU'; 
                                   '5 FU' = '5-FU'; 
                                   'gemcitabine' = 'Gemcitabine'; 
                                   'Gemcitabine Injection' = 'Gemcitabine'; 
                                   'GEMCITABINE' = 'Gemcitabine'; 
                                   'Gemcitabine HCL' = 'Gemcitabine'; 
                                   'gemcitabine HCL' = 'Gemcitabine'; 
                                   'Gemcitibine' = 'Gemcitabine'; 
                                   'Gemzar' = 'Gemcitabine';
                                   'gemzar' = 'Gemcitabine'; 
                                   'CAPECITABINE' = 'Capecitabine'; 
                                   'Xeloda' = 'Capecitabine'; 
                                   'Camptosar' = 'Irinotecan'; 
                                   'Irinotecan Hydrochloride' = 'Irinotecan'; 
                                   'irinotecan' = 'Irinotecan'; 
                                   'Leucovorin Calcium' = 'FA'; 
                                   'Leucovorin calcium' = 'FA';
                                   'Leucovorin' = 'FA'; 
                                   'leucovorin' = 'FA'; 
                                   'folinic acid' = 'FA'; 
                                   'Folinic Acid' = 'FA'; 
                                   'Eloxatin' = 'Pt'; 
                                   'Cisplatin' = 'Pt'; 
                                   'oxaliplatin' = 'Pt'; 
                                   'cisplatin' = 'Pt'; 
                                   'Oxaliplatin' = 'Pt'; 
                                   'Carboplatin' = 'Pt'; 
                                   'Abraxane' = 'Tax'; 
                                   'Docetaxel' = 'Tax'; 
                                   'ABRAXANE' = 'Tax'; 
                                   'Tarceva' = 'other'; 
                                   'doxorubicin' = 'other';
                                   'cyclophosphamide' = 'other'; 
                                   'Dexamethasone' = 'other'") %>% 
             factor, 
           drug_response = car::recode(treatment_best_response, 
                                       "'Clinical Progressive Disease' = 'Progression'; 
                                       'Stable Disease' = 'Stable'; 
                                       'Partial Response' = 'Partial response'; 
                                       'Complete Response' = 'Complete response'") %>% 
             factor(c('Progression', 
                      'Stable', 
                      'Partial response', 
                      'Complete response')), 
           drug_response_score = car::recode(drug_response, 
                                             "'Clinical Progressive Disease' = 0; 
                                             'Stable Disease' = 1; 
                                             'Partial Response' = 2; 
                                             'Complete Response' = 3") %>% 
             as.numeric) %>% 
    select(patient_id, 
           drug_type, 
           drug_response, 
           drug_response_score)

# Clearing the radiation data -----
  
  insert_msg('Clearing the radiation data')
  
  tcga_clinics$radiation <- tcga_clinics$radiation[c(-1, -2), ] %>% 
    mutate(patient_id = bcr_patient_barcode, 
           radiation_dose = ifelse(radiation_adjuvant_units == 'cGy', 
                                   as.numeric(radiation_total_dose)/100, 
                                   as.numeric(radiation_total_dose)), 
           radiation_response = car::recode(treatment_best_response, 
                                       "'Radiographic Progressive Disease' = 'Progression'; 
                                       'Stable Disease' = 'Stable'; 
                                       'Partial Response' = 'Partial response'; 
                                       'Complete Response' = 'Complete response'") %>% 
             factor(c('Progression', 
                      'Stable', 
                      'Partial response', 
                      'Complete response')), 
           radiation_response_score = car::recode(treatment_best_response , 
                                             "'Radiographic Progressive Disease' = 0; 
                                             'Stable Disease' = 1; 
                                             'Partial Response' = 2; 
                                             'Complete Response' = 3") %>% 
             as.numeric) %>% 
    select(patient_id, 
           radiation_dose, 
           radiation_response, 
           radiation_response_score)
  
# Clearing the follow-up data -----
  
  insert_msg('Cleaning the follow-up data')
  
  tcga_clinics$fup1 <- tcga_clinics$fup1[c(-1, -2), ] %>% 
    mutate(patient_id = bcr_patient_barcode, 
           relapse_days = as.numeric(new_tumor_event_dx_days_to), 
           relapse_fup = car::recode(new_tumor_event_type, 
                                 "'[Not Available]' = NA; 
                                 'Distant Metastasis' = 1; 
                                 'Locoregional Recurrence' = 1; 
                                 'New Primary Tumor' = 0; 
                                 'Locoregional Recurrence|Distant Metastasis' = 1") %>% 
             as.numeric) %>% 
    select(patient_id, 
           relapse_days, 
           relapse_fup)
  
  tcga_clinics$fup2 <- tcga_clinics$fup2[c(-1, -2), ] %>% 
    mutate(patient_id = bcr_patient_barcode, 
           relapse_days = as.numeric(new_tumor_event_dx_days_to), 
           relapse_fup = car::recode(new_tumor_event_type, 
                                     "'[Not Available]' = NA; 
                                 'Distant Metastasis' = 1; 
                                 'Locoregional Recurrence' = 1; 
                                 'New Primary Tumor' = 0; 
                                 'Locoregional Recurrence|Distant Metastasis' = 1") %>% 
             as.numeric)  %>% 
    select(patient_id, 
           relapse_days, 
           relapse_fup)
  
  tcga_clinics$fup <- rbind(tcga_clinics$fup1, 
                            tcga_clinics$fup2)
  
  ## duplicate handling: the clinical record with the shortest time to relapse is kept
  
  tcga_clinics$fup <- tcga_clinics$fup %>% 
    filter(!is.na(relapse_days), 
           !is.na(relapse_fup), 
           relapse_fup == 1) %>% 
    ddply(.(patient_id), 
          filter, 
          relapse_days == min(relapse_days, na.rm = T)) %>% 
    ddply(.(patient_id), 
          function(x) x[1, ]) %>%
    as_tibble

# Merging the general clinical and fup data and clearing the relapse issue -----
  
  insert_msg('Merging and clearing the relapse issue')

  tcga_clinics$clinical <- tcga_clinics[c('clinical', 
                                          'fup')] %>% 
    reduce(left_join, 
           by = 'patient_id')
  
  tcga_clinics$clinical <- tcga_clinics$clinical %>% 
    mutate(relapse = ifelse(!is.na(relapse_fup), 
                            relapse_fup, 
                            relapse), 
           rfs_days = ifelse(relapse == 1, 
                             relapse_days, 
                             os_days)) %>% 
    select(- relapse_fup, 
           - relapse_days)

# END -----
  
  tcga_clinics <- tcga_clinics[c('clinical', 
                                 'drugs', 
                                 'radiation')]
  
  insert_tail()