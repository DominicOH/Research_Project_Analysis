# Packages ----

library('gridExtra')
library('svglite')
library('xlsx')
library('UniprotR')


# Source ----

source('Cluster_analysis/100_2_abundance_heatmap.R')

# Cluster analysis basics ----

# With the leaf colours I can visually compare clusters from the dendrogram 
# using the colour of the branches.

# In this tibble, I've got protein Accession, colour from the tree, and cluster 
# number. 
abundance_clusters_and_colours <- abundance_dend %>% {
  tibble(
    'Clusters' = cutree(., h = 0.1, order_clusters_as_data = F),
    'Colours' = leaf_colors(d = .),
    'Accession' = get_leaves_attr(., "label")
  )
}

# I can use this tibble to look at particular clusters based on colour, then 
# find the exact proteins I want to look at in more detail. Having done this, 
# I still need to find a function to pull the index of particular proteins, 
# so that I can look up the position of known proteins of interest. 

# The cutree function cuts the dendrogram at either a specified number of 
# clusters or a specified height. 

abundance_dend %>%
  cutree(h = 0.1) %>%
  max() # I can use max() to find how many clusters are extracted at this height

# get_subdendrograms allows me to cut the dendrogram at a specified number
# of clusters into separate dendrograms. I'm using the number produced above.
abundance_subdends_dendlist <- abundance_dend %>%
  get_subdendrograms(k = 7)

# And the indexing can be explored to find particular proteins or subtrees. 
abundance_subdends_dendlist[[2]][[1]]

# find_dendrogram finds the dendrogram that contains these exact items
abundance_subdends_dendlist %>%
  find_dendrogram(selected_labels = c("A0A087WUB9", "P20930", 
                                      "I4AY87", "C6GLV3"))

# I can also use cut_lower_fun to produce list objects with all the values
# below the height specified. 
abundance_subdends <- 
  cut_lower_fun(abundance_dend, 0.1)

# Grep can be used to search for the individual proteins, finding which cluster 
# they are in. 
as.integer(grep("P03956", abundance_subdends))

# I could incorporate this into a plot call 
plot(abundance_subdends_dendlist[[as.integer(grep("P03956", abundance_subdends))
]])

# Analysis of specific clusters ----

# Clusters 2, 3, and 4 are visually interesting, showing Ox-LDL-treatment-specific 
# changes in abundance. We can extract these to analyse them further. 

load('Data/P_value_dfs.RData')

# Cluster 2 ----
cluster_2 <- abundance_subdends_dendlist[[2]]

# Extract just the information from cluster 2.
cluster_2_df <- abundance_subdends[[2]] %>%
  tibble() %>%
  rename("Accession" = ".") %>%
  merge(oxLDLdata_df, by = "Accession") %>%
  merge(y = LDLdata_df, by = c("Accession", "Description", 
                               "CA", "CB", "CTotal"), 
        all.x = T) %>% 
  # Merging with the p-value dataframes to incorporate P-values
  merge(y = oxLDL_p_ma_df, by = c("Accession", "Description"),
        all.x = T) %>%
  select(Accession, Description, CTotal.x, OxTotal.x, ox_p_value) %>%
  merge(y = LDL_p_ma_df, by = c("Accession", "Description"), 
                                all.x = T) %>%
  select(Accession, Description, CTotal.x, OxTotal.x, LDLTotal, 
         ox_p_value, ldl_p_value) 

# Save an image of the cluster 
{svg(filename = 'Images/abundance_cluster_2.svg', width = 3.14, height = 3.14)
cluster_2 %>%
  set("labels_cex", 0.6) %>%
  plot()
dev.off()}

# Clusters 3 and 4 ---- 

# We can extract both clusters by making the cutree height higher 
larger_subdendlist <- abundance_dend %>%
  get_subdendrograms(k = 3)

larger_subdends <- abundance_dend %>%
  cut_lower_fun(0.6)

clusters_3_and_4_df <- larger_subdends[[2]] %>%
  tibble() %>%
  rename("Accession" = ".") %>%
  merge(y = LDLdata_df, by = "Accession", 
        all.x = T) %>% 
  merge(oxLDLdata_df, by = c("Accession", "Description", 
                             "CA", "CB", "CTotal"), 
        all.x = T) %>%
  
  # Merging with the p-value dataframes to incorporate P-values
  merge(y = oxLDL_p_ma_df, by = c("Accession", "Description"),
        all.x = T) %>%
  select(Accession, Description, CTotal.x, OxTotal.x, ox_p_value) %>%
  merge(y = LDL_p_ma_df, by = c("Accession", "Description"), 
        all.x = T) %>%
  select(Accession, Description, CTotal.x, OxTotal.x, LDLTotal, 
         ox_p_value, ldl_p_value) 

  
# Extract the dendlist that contains both subdendrograms
clusters_3_and_4_subdends <- larger_subdendlist[[2]]

write.xlsx(clusters_3_and_4_df, file = 'Data/clusters_3_4_df.xlsx')
write.xlsx(cluster_2_df, file = 'Data/clusters_2_df.xlsx')

# Save an image of the cluster

{svg(filename = 'Images/abundance_cluster_3_and_4.svg', width = 3.14, 
     height = 3.14)
clusters_3_and_4_subdends %>%
  set("labels_cex", 0.6) %>%
  plot()
dev.off()}

# Narrow down secreted proteins ---- 

c2_subcellular <- GetSubcellular_location(cluster_2_df$Accession)
cluster_2_df$Subcellular_location <- c2_subcellular$Subcellular.location..CC.
c2_secreted <- filter(cluster_2_df, grepl(Subcellular_location, 
                                          pattern = "Secreted"))
