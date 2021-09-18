# Packages -----

library('dplyr')
library('xlsx')
library('dendextend')

# Source ----

load('Data/swing_dend_data.RData')
load('Data/P_value_dfs.RData')
load('Data/data_dfs.RData')

# Analysis ----


# The fviz_cluster plot previously recommended three optimal clusters, which 
# I'll stick with for now. 

# Find the height of clustering. 
log2_subdends <- get_subdendrograms(swing_dend, k = 3) # 0.6478146
swing_fc_df$Accession <- row.names(swing_fc_df)

# Clusters ----
cluster_1_df <- cut_lower_fun(swing_dend, h = 0.6478146)[[1]] %>%
   tibble() %>%
   rename("Accession" = ".") %>%
  merge(oxLDLdata_df, by = "Accession", all.x = T) %>%
     merge(y = LDLdata_df, by = c("Accession", "Description", 
                               "CA", "CB", "CTotal"), 
        all.x = T) %>% 
  # Merging with the p-value dataframes to incorporate P-values
     merge(y = oxLDL_p_ma_df, by = c("Accession", "Description"),
        all.x = T) %>%
     select(Accession, Description, CTotal.x, OxTotal.x, ox_p_value, CA, CB,
            OXA, OXB, 
         LDLA, LDLB) %>%
     merge(y = LDL_p_ma_df, by = c("Accession", "Description"), 
        all.x = T) %>%
     select(Accession, Description, CTotal.x, OxTotal.x, LDLTotal, 
         ox_p_value, ldl_p_value, CA, CB, OXA, OXB, LDLA, LDLB) %>%
     rename(OxTotal = OxTotal.x, CTotal = CTotal.x, 
           `Ox-LDL P- Value` = ox_p_value, `LDL P-Value` = ldl_p_value) %>%
  merge(y = swing_fc_df, by = "Accession", all.x = T) %>%
  mutate(`Ox-LDL Fold Change` = OxTotal / CTotal,
         `LDL Fold Change` = LDLTotal / CTotal)

# Plot cluster 1

png(filename = "Images/log2_heatmap_1.png", width = 8, 
    height = 8, units = "cm", res = 200)

swing_dend[[1]] %>%
  set("labels_cex", 0.6) %>%
  plot()

dev.off()

cluster_2_df <- cut_lower_fun(swing_dend, h = 0.6478146)[[2]] %>%
  tibble() %>%
  rename("Accession" = ".") %>%
  merge(data_df, by = "Accession", all.x = T) %>%
  # Merging with the p-value dataframes to incorporate P-values
  merge(y = oxLDL_p_ma_df, by = c("Accession", "Description"),
        all.x = T) %>%
  select(Accession, Description, CTotal.x, OxTotal.x, ox_p_value, OXA, OXB, 
         LDLA, LDLB) %>%
  merge(y = LDL_p_ma_df, by = c("Accession", "Description"), 
        all.x = T) %>%
  select(Accession, Description, CTotal.x, OxTotal.x, LDLTotal, 
         ox_p_value, ldl_p_value, OXA, OXB, LDLA, LDLB) %>%
  rename(OxTotal = OxTotal.x, CTotal = CTotal.x, 
         `Ox-LDL P- Value` = ox_p_value, `LDL P-Value` = ldl_p_value) %>%
  mutate(`Ox-LDL Fold Change` = OxTotal / CTotal,
         `LDL Fold Change` = LDLTotal / CTotal)

cluster_3_df <- cut_lower_fun(swing_dend, h = 0.6478146)[[3]] %>%
  tibble() %>%
  rename("Accession" = ".") %>%
  merge(oxLDLdata_df, by = "Accession", all.x = T) %>%
  merge(y = LDLdata_df, by = c("Accession", "Description", 
                               "CA", "CB", "CTotal"), 
        all.x = T) %>% 
  # Merging with the p-value dataframes to incorporate P-values
  merge(y = oxLDL_p_ma_df, by = c("Accession", "Description"),
        all.x = T) %>%
  select(Accession, Description, CTotal.x, OxTotal.x, ox_p_value, OXA, OXB, 
         LDLA, LDLB) %>%
  merge(y = LDL_p_ma_df, by = c("Accession", "Description"), 
        all.x = T) %>%
  select(Accession, Description, CTotal.x, OxTotal.x, LDLTotal, 
         ox_p_value, ldl_p_value, OXA, OXB, LDLA, LDLB) %>%
  rename(OxTotal = OxTotal.x, CTotal = CTotal.x, 
         `Ox-LDL P- Value` = ox_p_value, `LDL P-Value` = ldl_p_value) %>%
  mutate(`Ox-LDL Fold Change` = OxTotal / CTotal,
         `LDL Fold Change` = LDLTotal / CTotal)

# Export this as an .xlsx

write.xlsx(cluster_1_df, file = 'Data/log2_heatmap_clusters.xlsx', 
           sheetName = "Cluster_1")
write.xlsx(cluster_2_df, file = 'Data/log2_heatmap_clusters.xlsx', 
           append = T, sheetName = "Cluster_2")
write.xlsx(cluster_3_df, file = 'Data/log2_heatmap_clusters.xlsx', 
           append = T, sheetName = "Cluster_3")


# I also want to extract the moderately large cluster in cluster 3, which seems
# to show higher Ox-LDL fold change and negative polarity swing (more apical)

# Cluster "10" ----

cluster_10_df <- cut_lower_fun(swing_dend, h = 8.232162e-02)[[10]] %>%
  tibble() %>%
  rename("Accession" = ".") %>%
  merge(data_df, by = "Accession", all.x = T) %>% 
  # Merging with the p-value dataframes to incorporate P-values
  merge(y = oxLDL_p_ma_df, by = c("Accession", "Description"),
        all.x = T) %>%
  select(Accession, Description, CTotal.x, OxTotal.x, ox_p_value, CA, CB, OXA,
         OXB, LDLA, LDLB) %>%
  merge(y = LDL_p_ma_df, by = c("Accession", "Description"), 
        all.x = T) %>%
  select(Accession, Description, CTotal.x, OxTotal.x, ox_p_value, CA, CB, OXA, 
         OXB, LDLA, LDLB) %>%
  merge(y = LDL_p_ma_df, by = c("Accession", "Description"), 
        all.x = T) %>%
  select(Accession, Description, CTotal.x, OxTotal.x, LDLTotal, 
         ox_p_value, ldl_p_value, CA, CB, OXA, OXB, LDLA, LDLB) %>%
  rename(OxTotal = OxTotal.x, CTotal = CTotal.x, 
         `Ox-LDL P- Value` = ox_p_value, `LDL P-Value` = ldl_p_value) %>%
  merge(y = swing_fc_df, by = "Accession", all.x = T) %>%
  mutate(`Ox-LDL Fold Change` = OxTotal / CTotal,
         `LDL Fold Change` = LDLTotal / CTotal) 

get_branches_heights(swing_dend)
cutree(swing_dend, h = 8.232162e-02)

# Plotting cluster "10" 
plot(get_subdendrograms(swing_dend, k = 14)[[10]])

png(filename = "Images/log2_heatmap_10.png", width = 12, 
    height = 8, units = "cm", res = 200)

get_subdendrograms(swing_dend, k = 14)[[10]] %>% 
  set("labels_cex", 0.2) %>%
  plot()

dev.off()

# Cluster "4" ----

# Plot cluster 4
png(filename = "Images/log2_heatmap_4.png", width = 10, 
    height = 8, units = "cm", res = 200)

get_subdendrograms(swing_dend, 5)[[4]] %>%
  set("labels_cex", 0.6) %>%
  plot()

dev.off()

cluster_4_df <- cut_lower_fun(swing_dend, h = 3.468447e-01)[[4]] %>%
  tibble() %>%
  rename("Accession" = ".") %>%
  merge(data_df, by = "Accession", all.x = T) %>% 
  # Merging with the p-value dataframes to incorporate P-values
  merge(y = oxLDL_p_ma_df, by = c("Accession", "Description"),
        all.x = T) %>%
  select(Accession, Description, CTotal.x, OxTotal.x, ox_p_value, CA, CB, OXA, 
         OXB, LDLA, LDLB) %>%
  merge(y = LDL_p_ma_df, by = c("Accession", "Description"), 
        all.x = T) %>%
  select(Accession, Description, CTotal.x, OxTotal.x, LDLTotal, 
         ox_p_value, ldl_p_value, CA, CB, OXA, OXB, LDLA, LDLB) %>%
  rename(OxTotal = OxTotal.x, CTotal = CTotal.x, 
         `Ox-LDL P- Value` = ox_p_value, `LDL P-Value` = ldl_p_value) %>%
  merge(y = swing_fc_df, by = "Accession", all.x = T) %>%
  mutate(`Ox-LDL Fold Change` = OxTotal / CTotal,
         `LDL Fold Change` = LDLTotal / CTotal) 

# Export this as an .xlsx

# Prevents manually having to remove this information. 
cluster_10_df$Description <- gsub(cluster_10_df$Description, 
                          replacement = "", 
                          pattern = "OS.*")

write.xlsx(cluster_1_df, file = 'Data/log2_heatmap_clusters.xlsx', 
           sheetName = "Cluster_1")
write.xlsx(cluster_4_df, file = 'Data/log2_heatmap_cluster.xlsx', 
           append = T, sheetName = "Cluster_4")
write.xlsx(cluster_10_df, file = 'Data/log2_heatmap_clusters.xlsx',
           append = T, sheetName = "Cluster_10")


