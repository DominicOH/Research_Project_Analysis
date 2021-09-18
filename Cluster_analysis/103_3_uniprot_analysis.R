# Packages ----

library('UniprotR')

# Source ----

source(file = 'Cluster_analysis/103_2_swing_dend_analysis.R')

# Subcellular location of cluster 4 ----

cluster_4_accessions <- cluster_4_df$Accession

# Only run when needed 
# cluster_4_subcel <- GetSubcellular_location(cluster_4_accessions)

# Cluster_10 ----
# Only run when needed
# cluster_10_subcellular <- GetSubcellular_location(cluster_10_df$Accession)

cluster_10_subcellular$Subcellular.location..CC. <- 
  gsub(cluster_10_subcellular$Subcellular.location..CC., 
                                  replacement = "", 
                                  pattern = "SUBCELLULAR_LOCATION: ")
cluster_10_df$Subcellular_location <- 
  cluster_10_subcellular$Subcellular.location..CC.
cluster_10_df$Subcellular_location <- gsub(cluster_10_subcellular$Subcellular.location..CC., 
                                           replacement = "", 
                                           pattern = "SUBCELLULAR LOCATION: ")
