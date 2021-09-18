# Packages ----

library('tidyverse')

# Source ----

load('Data/data_dfs.RData')

# For the dataset. Not using P-values, but removing contaminant proteins. 

cluster_data_df <- data_df %>%
  select(Accession, CTotal, OxTreatment, LDLTreatment, OxTotal, LDLTotal) %>%
  transmute(Accession = Accession, 
            Ox_treatment_ratio = OxTreatment / CTotal,
            LDL_treatment_ratio = LDLTreatment / CTotal,
            CTotal = CTotal,
            OxTotal = OxTotal,
            LDLTotal = LDLTotal) %>%
  filter(Ox_treatment_ratio <= 1.0 & LDL_treatment_ratio <= 1.0)

save(cluster_data_df, file = 'Data//cluster_data.Rdata')
