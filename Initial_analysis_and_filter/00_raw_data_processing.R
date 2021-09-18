# Packages ----

library('tidyverse')
library('readxl')

# Ensure that this file matches the required directory, or this will fail. 
Data_unedited <- read_excel("./Data/Data_unedited.xlsx")

# Note: To preserve the integrity of data that could be used in later
# publications, this data has been redacted by the author. 

# Global data dataframe ----

# calculate means by producing sums / n (total) trials
data_df <- data.frame(Description = Data_unedited$Description,
                      Accession = Data_unedited$Accession,
                      CA = rowSums(select(Data_unedited, c(3:6)),
                                   na.rm = T)/4,
                      CB = rowSums(select(Data_unedited, c(7:10)),
                                   na.rm = T)/4,
                      OXA = rowSums(select(Data_unedited, c(11:13)),
                                    na.rm = T)/3,
                      OXB = rowSums(select(Data_unedited, c(14:16)),
                                    na.rm = T)/3,
                      OxTreatment = rowSums(select(Data_unedited, c(17:19)),
                                            na.rm = T)/3,
                      LDLA = rowSums(select(Data_unedited, c(20:22)),
                                     na.rm = T)/3,
                      LDLB = rowSums(select(Data_unedited, c(23:25)),
                                     na.rm = T)/3,
                      LDLTreatment = rowSums(select(Data_unedited, c(26:28)),
                                             na.rm = T)/3
                      ) 
attach(data_df)

# Sample totals are useful for future calculations.  
data_df$CTotal <- CA + CB; data_df$OxTotal <- OXA + OXB; data_df$LDLTotal <- 
  LDLA + LDLB
attach(data_df)

# Quick facts about the data ----
data_df %>%
  count(OxTotal > 0) 
data_df %>%
  count(LDLTotal > 0)

# Ox-LDL-Specific Dataframe ----
# new dataframe containing just averages from OxLDL samples
oxLDLdata_df <- data.frame(Description = Description, 
                           Accession = Accession, 
                           CA = CA, 
                           CB = CB, 
                           CTotal = CTotal,
                           OXA = OXA,
                           OXB = OXB, 
                           OxTotal = OxTotal,
                           OxTreatment = OxTreatment)

# LDL-Specific Dataframe ----

# repeat for LDL samples
LDLdata_df <- data.frame(Description = Description, 
                         Accession = Accession, 
                         CA = CA, 
                         CB = CB,
                         CTotal = CTotal,
                         LDLA = LDLA, 
                         LDLB = LDLB, 
                         LDLTotal = LDLTotal,
                         LDLTreatment = LDLTreatment)

detach(data_df)