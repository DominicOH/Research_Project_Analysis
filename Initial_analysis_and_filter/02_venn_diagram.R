# Packages ----

library('venneuler')
library('rJava')
# Please note: The venneuler package requires rJava to function. 
# rJava, in turn, requires Java software to be present on the OS. 
# This can be done here: https://www.java.com/en/ 
# Skip this script if unable to download Java. 

# Source ----

source('Initial_analysis_and_filter/00_raw_data_processing.R')

# Concept ----

# A Venn diagram helps get a good grasp of overlaps or unique datapoints within
# a set of data. 

# This could be useful to show how many proteins are present in each sample. 

# Proportional Venn diagram ----

# To show the shape of the data, it's first necessary to prepare a Venn diagram 
# showing approximate proportions. 

# Area vectors ----
attach(data_df) 
proportional_areas <- 
 c("A&B&C" = data_df %>% 
      tally(CTotal > 0 & 
              OxTotal > 0 &
              LDLTotal > 0),
    "A&B" = data_df %>%
      tally(CTotal > 0 &
               OxTotal > 0 &
               LDLTotal == 0),
    "A&C" = data_df %>%
      tally(CTotal > 0 &
               OxTotal == 0 & 
               LDLTotal > 0),
    "B&C" = data_df %>%
      tally(CTotal == 0 &
               OxTotal > 0 &
               LDLTotal > 0),
   A = data_df %>%
     tally(CTotal > 0 & 
               OxTotal == 0 &
               LDLTotal == 0),
   B = data_df %>%
     tally(CTotal == 0 & 
               OxTotal > 0 &
               LDLTotal == 0),
   C = data_df %>% 
     tally(CTotal == 0 & 
               OxTotal == 0 &
               LDLTotal > 0)
   )

detach(data_df)

# Prepare vector of areas for the Venn diagrams
proportional_areas <- c(A = proportional_areas[["A.n"]], # Control
           B = proportional_areas[["B.n"]], # Ox-LDL Treatment 
           C = proportional_areas[["C.n"]], # LDL Treatment 
           "A&B" = proportional_areas[["A&B.n"]], # Control u Ox-LDL
           "A&C" = proportional_areas[["A&C.n"]], # Control u LDL
           "B&C" = proportional_areas[["B&C.n"]], # Ox-LDL u LDL 
           "A&B&C" = proportional_areas[["A&B&C.n"]] # Control u Ox-LDL u LDL
)
proportional_venn <- venneuler(proportional_areas)

# Remove default labels. There is too much overlap for this to be clear. 
proportional_venn$labels <- c("", "", "")
proportional_venn$colors <- c(0.2, 0.5, 0.8)

proportional_venn_plot <- plot(proportional_venn)

size_vector <- 1.7

# Add text for each group size
proportional_venn_plot + 
  text(0.07, 0.7, "Ox-LDL\nTreated\n19", cex = size_vector) + # Ox-LDL only
  text(0.1, 0.23, "Control Only\n0", cex = size_vector) + # Control-only
  text(1, 0.58, "LDL\nTreated\n2", cex = size_vector) + # LDL-only 
  text(0.5, 0.5, "674", cex = size_vector) + # u All
  text(0.74, 0.35, "154", cex = size_vector) # Control u LDL 

# As a number would not directly fit on this space, we can add a line
# and a text element referring to it. Control u Ox-LDL 


proportional_venn_plot + 
  segments(y0 = 0.5, x0 = 0.15, 
           y1 = 0.5, x1 = 0.2) + 
  text(0.10, 0.5, "45", cex = size_vector) + 
  segments(y0 = 0.81, x0 = 0.48, 
           y1 = 0.82, x1 = 0.65) + 
  text(0.7, 0.82, "32", cex = size_vector)
