# Packages ----

library('tidyverse')
library('svglite')
library('gridExtra')
library('xlsx')
library('UniprotR')

# Source ----

load('Data/data_dfs.RData')
load('Data/n_detections_dfs.RData')
load('Data/treatment_ratio_dfs.RData')
load('Data/P_value_dfs.RData')

# Concept ----

# The first role of the filter is to catch likely contaminants: 
# A suggested method is based on the ratio of abundance from the vehicle-only
# sample to the control. 
# In addition, if the Total Treated sample is very similar in abundance to the 
# vehicle-only sample, then a large amount of what is detected is likely just 
# from the treatment.

# As a basic filter - we can set thresholds on the characteristics discussed up
# to now. Since this is user-defined, it may be arbitrary, but it could still 
# produce a dataset that holds up for later clustering. 

# Sensitive thresholds: 
# Positive fold change (i.e. M > 0)
# Ratio < =1.00
# P-value < 0.05

# Note: Filtering by N.detections has been removed, as this is not independent 
# from P-value. 

# Ox-LDL sample filtering ----

# Adding all the filterable characteristics to the main Ox-LDL dataframe. 
oxLDLdata_df$Treatment_ratio <- oxLDL_treatment_ratio_ma_df$Treatment_ratio
oxLDLdata_df$ox_p_value <- oxLDL_p_ma_df$ox_p_value; 
oxLDLdata_df$Ox_detections_Mean <- oxLDL_N_detections_ma_df$Ox_detections_Mean

# Subset just the non-contaminants
oxLDL_non_contaminants_df <- subset(x = oxLDLdata_df, 
                                    subset = Treatment_ratio <= 1.0 &
                                      ox_p_value < 0.05 &
                                      M > 0.0)
# Plot on the MA
oxLDL_contaminant_filter_plot <- 
  ggplot(mapping = aes(x = A, y = M)) +
  # Using anti_join() removes proteins not considered contaminants from the 
  # main dataframe. 
  geom_point(data = anti_join(x = oxLDLdata_df, y = oxLDL_non_contaminants_df),
             aes(colour = "Filtered out"),
             size = 0.5) +  
  # non-contaminants are plotted as a separate level. 
  geom_point(data = oxLDL_non_contaminants_df, # Non-contaminants
             aes(colour = "Retained"),
             size = 0.5) +
  # Reference proteins are subsetted by Accession to change size and shape. 
  geom_point(data = subset(x = oxLDLdata_df, 
                          subset = Accession == 'O00622' # CCN1 
                          | Accession == 'A0A024R462' # FN1 
                          | Accession == 'P03956'), # MMP1
            size = 4,
            mapping = aes(shape = Accession,
                          colour = "Retained")) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12, 14), limits = c(0,14)) + 
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4, 6, 8, 10), limits = c(-4,10)) + 
  theme_bw() +
  theme(panel.grid = element_line(colour = 'grey', size = rel(0.3))) +
  scale_color_viridis_d(direction = -1) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  xlab('Log2 average intensity') + 
  ylab('Log2 fold change from control') +
  labs(colour = NULL) +
  ggtitle(label = 'Protein changes filtered from the dataset',
          subtitle = 'Treatment with Ox-LDL')

# LDL sample filtering ----

# Steps above repeated for the LDL datasets. Adding all the characteristics. 
LDLdata_df$Treatment_ratio <- LDL_treatment_ratio_ma_df$Treatment_ratio
LDLdata_df$ldl_p_value <- LDL_p_ma_df$ldl_p_value; 
LDLdata_df$LDL_detections_Mean <- LDL_N_detections_ma_df$LDL_detections_Mean

# Subsetting the non-contaminants. 
LDL_non_contaminants_df <- subset(x = LDLdata_df, 
                                  subset = Treatment_ratio <= 1.0 &
                                    ldl_p_value < 0.05 &
                                    M > 0.0)

LDL_contaminant_filter_plot <- 
  ggplot(mapping = aes(x = A, y = M)) +
  geom_point(data = anti_join(x = LDLdata_df, y = LDL_non_contaminants_df),
             aes(colour = "Filtered out"),
             size = 0.5) +  
  geom_point(data = subset(x = LDLdata_df, 
                           subset = Accession == 'A0A024R462' # FN1 
                           | Accession == 'P03956'), # CCN1
             size = 4,
             mapping = aes(shape = Accession,
                           colour = "Filtered out")) +
  geom_point(data = LDL_non_contaminants_df, # Non-contaminants
             aes(colour = "Retained"),
             size = 0.5) +
  geom_point(data = subset(x = LDLdata_df, 
                           subset = Accession == 'O00622'), # CCN1
             size = 4,
             mapping = aes(shape = Accession,
                           colour = "Retained")) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12, 14), limits = c(0,14)) + 
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4, 6, 8, 10), limits = c(-4,10)) + 
  scale_color_viridis_d(direction = -1) + 
  theme_bw() +
  theme(panel.grid = element_line(colour = 'grey', size = rel(0.3))) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  xlab('Log2 average intensity') + 
  ylab('Log2 fold change from control') +
  labs(colour = NULL) +
  ggtitle(label = 'Protein changes filtered from the dataset',
          subtitle = 'Treatment with LDL')

# Panels for figure ----

filter_visual_ma <- 
  arrangeGrob(nrow = 1, ncol = 2, 
             top = "Protein changes removed from each dataset",
             oxLDL_contaminant_filter_plot + 
               ggtitle(label = "Ox-LDL", subtitle = NULL) +
               theme(panel.grid = element_line(colour = 'grey', size = rel(1)),
                     legend.position = "none"), 
             LDL_contaminant_filter_plot + 
               ggtitle(label = "LDL", subtitle = NULL) +
               theme(panel.grid = element_line(colour = 'grey', size = rel(1)),
                     legend.position = "none"))

# plot(filter_visual_ma)

# Save this to file as part of Figure 6. 
ggsave(plot = filter_visual_ma,
       filename = './Images/filter_visual.svg', 
       width = 14, height = 9,
       units = "cm",
       device = "svg",
       dpi = "print")

# Prepare tables of filtered proteins ----

# I want a Table to use in a figure, showing just the retained proteins. 
# The table just needs the Accession, Description, Ox-Fold Change and 
# LDL Fold Change. Where the same protein has been filtered through in both 
# treated samples, these rows should be merged. 

table <- oxLDL_non_contaminants_df %>%
  full_join(LDL_non_contaminants_df, 
            by = c("Accession", "Description")) %>%
  select(c(Description, Accession)) %>%
  merge(oxLDLdata_df, by = c("Accession", "Description")) %>%
  merge(LDLdata_df, by = c("Accession", "Description")) %>%
  select(c("Accession", "Description", "CTotal.x", "OxTotal", "LDLTotal", 
           "ox_p_value", "ldl_p_value")) %>%
  rename("CTotal" = "CTotal.x") %>%
  transmute(Description = Description,
            Accession = Accession,
            'Ox-LDL Fold Change' = OxTotal / CTotal,
            'LDL Fold Change' = LDLTotal / CTotal,
            ox_p_value = ox_p_value,
            ldl_p_value = ldl_p_value) %>%
  # Renames the columns appropriately
  rename('Ox-LDL P-value' = 'ox_p_value',
         'LDL P-value' = 'ldl_p_value') %>%
  arrange(Description)

# Prevents manually having to remove this information. 
table$Description <- gsub(table$Description, 
                          replacement = "", 
                          pattern = "OS.*")

# Get subcellular location from UniProt to filter down to secreted proteins. 
subcell <- GetSubcellular_location(table$Accession)
# Add this to the table
table$Subcellular_location <- subcell$Subcellular.location..CC.

# Produce an Excel document for inspection. 
write.xlsx(x = table, file = 'Data/figure_table.xlsx')