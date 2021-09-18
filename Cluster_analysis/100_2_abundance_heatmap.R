# Packages ----

library('dendextend')
library('gplots')
library('viridis')
library('svglite')
library('RColorBrewer')
library('factoextra')
library('tidyverse')

# Source ----

load('Cluster_analysis/cluster_data.Rdata')
load('Data/data_dfs.RData')

# Concept ----

# Rather than using Euclidean distances of the total of each point, it's 
# informative to see how proteins are clustered with regard to the rate of 
# change or correlation. The same concept can be applied to apical / basal. 

# Going to use Euclidean averages as the default means of clustering and 
# measuring distance between points. 

# Dataframe ----

attach(cluster_data_df)
abundance_df <- data.frame(row.names = Accession,
                           CTotal = CTotal,
                           OxTotal = OxTotal, 
                           LDLTotal = LDLTotal)
detach(cluster_data_df)

# The purpose of this is to find the Pearson correlation of protein abundances 
# for each other protein in the dataframe.
# Code for this section is derived in part from Watson (2015).
# Performs a hierarchical cluster model from the matrix above. 

# dend_expend finds the optimal clustering and distance algorithm. 
# As I decided to use a Pearson correlation based approach, the most optimal
# algorithm turns out to be Euclidean Average.
# By contrast, in the base abundance_df, Euclidean Complete is marginally better

# Useful to save the hierarchical cluster data separately
abundance_dend <- abundance_df %>%
  get_dist(method = "pearson") %>%
  # Revised dendrogram using "average" method. 
  hclust(method = "average") %>%
  as.dendrogram()

test_dist <- abundance_df %>%
  get_dist(method = "pearson") %>%
  as.matrix()
abundance_dend %>% plot()

test <- abundance_df %>%
  t() %>%
  cor(method = "pearson") %>%
  as.matrix()
# Heatmap ----

# This confirms what we know from the statistical analysis, that
# the LDL-treated sample contains more abundant proteins than either other 
# sample. What can be compared are proteins showing consistent changes in 
# both treated samples. 
# pdf(file = 'Images/abundance_heatmap.pdf')
# Only run one image generator at a time
{png(filename = "Images/abundance_heatmap.png", width = 12, 
   height = 16, units = "cm", res = 200)

heatmap.2(x = as.matrix(abundance_df), 
          Rowv = abundance_dend,
          ylab = "",
          dendrogram = "row",
          col = plasma(n = 10, direction = -1),
          labCol = "",
          labRow = "",
          scale = "row")
dev.off()}

# Heatmap cluster 2 ----
png(filename = "Images/abundance_heatmap_2.png", width = 8, 
    height = 8, units = "cm", res = 200)

abundance_dend[[2]][[1]] %>%
  set("labels_cex", 0.6) %>%
  plot()

dev.off()

# Heatmap cluster 1 ----
png(filename = "Images/abundance_heatmap_1.png", width = 8, 
    height = 8, units = "cm", res = 200)

abundance_dend[[1]][[2]] %>%
  set("labels_cex", 0.6) %>%
  plot()

dev.off()
# Citation: 
# Watson, M. (2015). You probably don’t understand heatmaps [online]. Available 
# from: 
# http://www.opiniomics.org/you-probably-dont-understand-heatmaps/ [Accessed 
# 20.07.2021]. 

