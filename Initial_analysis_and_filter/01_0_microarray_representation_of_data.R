# Packages ----

library('ggplot2')

# Source ----

load('Data/data_dfs.RData')

# Concept ----

# An 'MA' plot, originally referring to a microarray, can be used to visualise 
# fold-difference between states, relative to the average abundance of either 
# state. This can later be overlaid with conditional colouring to demonstrate
# how fold-difference matches other factors. 

# M: Intensity ratio = (log2(Treated Sample/Control)) 
# A: Average intensity = (log2(Treated Sample) + log2(Control))/2

# Derived from:  
# Gatto, L., Breckels, L.M., Naake, L., and Gibb, S. (2015). 
# Visualization of proteomics data using R and Bioconductor. Proteomics, 15, 
# pp. 1375-1389. 

# Ox-LDL MA ----

# Defines a function to more quickly produce MA tables from a dataframe
ox_define_MA <- function(df){
  # Takes a dataframe object with columns 'OxTotal' and 'CTotal' and produces
  # corresponding MA columns for use in a MA plot. 
  
  x <- df$OxTotal
  y <- df$CTotal
  
  M <- log2(x/y) # binary intensity ratio 
  A <- (log2(x) + log2(y))/2 # average intensity
  
  ma_dataframe <- data.frame(M = M, A = A)
  
  return(ma_dataframe)
} 

# Define MA values from the Ox-LDL dataframe and add them to it. 
oxLDLdata_df$M <- ox_define_MA(oxLDLdata_df)$M
oxLDLdata_df$A <- ox_define_MA(oxLDLdata_df)$A

# For the purpose of graphing this data later on, it's useful to remove Inf and -Inf 
# values produced by the log transformation. 
oxLDLdata_df <- oxLDLdata_df %>% 
                  filter(is.finite(M) & is.finite(A))
# 719 data points remain

# LDL function ----

ldl_define_MA <- function(df){
  # Takes a dataframe object with columns 'OxTotal' and 'CTotal' and produces
  # corresponding MA columns for use in a MA plot. 
  
  x <- df$LDLTotal
  y <- df$CTotal
  
  M <- log2(x/y) # binary intensity ratio 
  A <- (log2(x) + log2(y))/2 # average intensity
  
  ma_dataframe <- data.frame(M = M, A = A)
  
  return(ma_dataframe)
} 

# Define MA values from the LDL dataframe and add them to it.
LDLdata_df$M <- ldl_define_MA(LDLdata_df)$M
LDLdata_df$A <- ldl_define_MA(LDLdata_df)$A

LDLdata_df <- LDLdata_df %>% 
  filter(is.finite(M) & is.finite(A))
# 828 data points remain.

# Blank microarrays ----

# Using the define function within a plot call to test it works correctly. 
png(filename = 'Images/blank_ox_ma.png', units = "cm", res = 200, 
    width = 8, height = 8)

ggplot(data = ox_define_MA(oxLDLdata_df), mapping = aes(x = A, y = M)) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(size = 0.6) + 
  xlab('Log2 average intensity') + 
  ylab('Log2 fold change from control') +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 14, 2), limits = c(0,14)) + 
  scale_y_continuous(breaks = seq(-4, 10, 2), limits = c(-4,10)) + 
  theme(panel.grid = element_line(colour = 'grey', size = rel(0.3))) +
  ggtitle(label = "Ox-LDL Treated")

dev.off()

png(filename = 'Images/blank_ldl_ma.png', units = "cm", res = 200, 
    width = 8, height = 8)

ggplot(data = ldl_define_MA(LDLdata_df), mapping = aes(x = A, y = M)) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(size = 0.6) + 
  xlab('Log2 average intensity') + 
  ylab('Log2 fold change from control') +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 14, 2), limits = c(0,14)) + 
  scale_y_continuous(breaks = seq(-4, 10, 2), limits = c(-4,10)) + 
  theme(panel.grid = element_line(colour = 'grey', size = rel(0.3))) +
  ggtitle(label = "LDL Treated")

dev.off()

# Save dataframes ----

save(oxLDLdata_df, LDLdata_df, data_df, file = 'Data/data_dfs.RData')
