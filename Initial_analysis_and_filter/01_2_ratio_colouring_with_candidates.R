# Source ----

source('./Initial_analysis_and_filter/01_1_ratio_conditional_colouring.R')

# Concept ----

# Three candidate proteins were selected from the initial data set for further
# study; these were: CCN1 (O00622), FN1 (A0A024R462), and MMP1 (P03956).

# Expression of these proteins under Ox-LDL and LDL was tested using ELISA,
# which showed elevated expression from Ox-LDL. These proteins can help to serve
# as candidates to train our data on, to potentially reveal a characteristic 
# pattern from which other genuine expression changes could be detected. 

# Ox-LDL ----
# For closer inspection: 
# we can isolate the proteins of interest (poi) from the Ox- and LDL dataframes
ox_poi_df <- oxLDLdata_df %>% # new dataframe of only these accession #s 
  filter(Accession == 'O00622' # CCN1 
         | Accession == 'A0A024R462' # FN1 
         | Accession == 'P03956') # MMP1

# These can also be mapped directly on top of the other plot. 
oxLDL_treatment_ratio_colouring_ma_candidates <- 
  oxLDL_treatment_ratio_colouring_ma + 
  geom_point(data = subset(x = oxLDLdata_df, 
                           subset = Accession == 'O00622' # CCN1 
                           | Accession == 'A0A024R462' # FN1 
                           | Accession == 'P03956'), # MMP1
             size = 4,
             mapping = aes(shape = Accession)) + 
  labs(shape = "Candidate \nproteins")

# LDL ----

ldl_poi_df <- LDLdata_df %>% # new dataframe of only these accession 
  filter(Accession == 'O00622' # CCN1 
         | Accession == 'A0A024R462' # FN1 
         | Accession == 'P03956') # MMP1

LDL_treatment_ratio_colouring_ma_candidates <- 
  LDL_treatment_ratio_colouring_ma + 
  geom_point(data = subset(x = LDLdata_df, 
                           subset = Accession == 'O00622' # CCN1 
                           | Accession == 'A0A024R462' # FN1 
                           | Accession == 'P03956'), # MMP1
             size = 4,
             mapping = aes(shape = Accession)) + 
  labs(shape = "Candidate \nproteins")