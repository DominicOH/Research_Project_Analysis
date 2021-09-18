# Packages ----

library('dendextend')
library('gplots')
library('dplyr')
library('viridis')
library('RColorBrewer')
library('factoextra')

# Source ----

load('Cluster_analysis/cluster_data.Rdata')
load('Data/data_dfs.RData')

# Dataframe ----

# Need to set zero values to a very small value (i.e. 0.01)
data_df[data_df == 0] <- 0.01

# Making a dataframe using the data_df
swing_fc_df <- data_df %>%
  select(Accession, OXA, OXB, LDLA, LDLB, CA, CB) %>%
  transmute(Accession = Accession,
            CA = CA,
            CB = CB,
            OXA = OXA,
            OXB = OXB,
            LDLA = LDLA,
            LDLB = LDLB) %>%
  mutate(CRatio = CA / (CA + CB),
         OxRatio = OXA / (OXA + OXB),
         LDLRatio = LDLA / (LDLA + LDLB)) %>%
  transmute(Accession = Accession, 
            Log2_Ox_LDL_FC = log2((OXA + OXB) / (CA + CB)),
            Ox_LDL_swing = OxRatio - CRatio,
            Log2_LDL_FC = log2((LDLA + LDLB) / (CA + CB)),
            LDL_swing = LDLRatio - CRatio) %>%
  semi_join(y = cluster_data_df, by = "Accession")

# Changing row names to Accession numbers
row.names(swing_fc_df) <- swing_fc_df$Accession
# Removing the Accession column to make it a completely numeric matrix
swing_fc_df$Accession <- NULL

# Clustering ----

swing_dend <- swing_fc_df %>%
  # Distances based on pearson correlation
  get_dist(method = "pearson") %>%
  # hierarchical clustering on average linkage method
  hclust(method = "average") %>%
  as.dendrogram()

swing_dend %>% plot()
# How many clusters can we see?

fviz_nbclust(swing_fc_df, FUNcluster = kmeans)
fviz_cluster(object = list(data = swing_fc_df, cluster = cutree(swing_dend,
                                                                k = 3)),
             geom = "point", 
             choose.vars = ) + 
  ggtitle(label = NULL) + 
  labs(colour = "Cluster", fill = "Cluster", shape = "Cluster")

# Based on the number of clusters, and appearance, the average linkage 
# algorithm appears best suited for this dataset. 

# Heatmap with simple scaling  ----

# By dividing the heatmapping log2 fold change values by 10, the dataset should appear 
# more visible in the heatmap. 

modified_swing_df <- swing_fc_df %>%
  transmute(Log2_Ox_LDL_FC = (Log2_Ox_LDL_FC / 10),
            Ox_LDL_swing = Ox_LDL_swing,
            Log2_LDL_FC = (Log2_LDL_FC / 10),
            LDL_swing = LDL_swing) %>%
  as.matrix()

# Save an image of the heatmap
png("Images/log2FC.png", width = 16, 
    height = 16, units = "cm", res = 200)

# Plot the heatmap 
heatmap.2(x = modified_swing_df,
          col = redgreen(20),
          Colv = F,
          labRow = "",
          labCol = "",
          Rowv = swing_dend, 
          dendrogram = "row",
          colsep = 2,
          cexCol = 0.5,
          cexRow = 0.5,
          scale = "none",
          breaks = seq(-1, 1, 0.1),
          key.xtickfun = function() {
            breaks = c(-1, 0, 1)
            list(at = parent.frame()$scale01(breaks),
                 labels = c("< -1", "0", "> 1"))
          })
dev.off()

# Save the dendrogram and modified data for deeper analysis
save(swing_dend, modified_swing_df, swing_fc_df, 
     file = 'Data/swing_dend_data.RData')
