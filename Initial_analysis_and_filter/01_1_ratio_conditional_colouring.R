# Packages ----

library('scales')
library('gridExtra')
library('svglite')

# Source ----

source('./Initial_analysis_and_filter/01_0_microarray_representation_of_data.R')

# Conditional colouring by ratio of Vehicle to Control abundance ---- 

# The Control sample aims to represent the endothelial secretome in its native
# and unaltered state. 

# The Vehicle-only sample represents proteins that are thought to be added by
# the Vehicle, causing contamination. 

# The ratio of Vehicle to Control abundance gives a good idea of proteins that 
# may be dangerously influenced by the Vehicle in the Treated Sample. 

# Ox-LDL ----
# Calculate the ratio of Treatment to Control and add to the dataframe. 
attach(oxLDLdata_df)
oxLDL_treatment_ratio_ma_df <- data.frame(Accession = Accession,
                                       Treatment_ratio = OxTreatment / CTotal,
                                       OxTreatment = OxTreatment,
                                       M = M,
                                       A = A)
detach(oxLDLdata_df)
# By separating subsets of the data, including one with Treatment value equal to 
# zero, it's easier to highlight proteins considered unlikely to be contaminants 
oxLDL_treatment_ratio_colouring_ma <- 
  ggplot(data = oxLDL_treatment_ratio_ma_df, 
         mapping = aes(x = A, y = M, colour = Treatment_ratio)) + 
  geom_point(data = subset(x = oxLDL_treatment_ratio_ma_df, 
                           subset = OxTreatment > 0),
                           aes(colour = Treatment_ratio),
             size = 0.6) +
  geom_point(data = subset(x = oxLDL_treatment_ratio_ma_df, 
                           subset = OxTreatment == 0),
             shape = 2, size = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  theme_bw() + 
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12, 14), limits = c(0,14)) + 
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4, 6, 8, 10), limits = c(-4,10)) + 
  theme(panel.grid = element_line(colour = 'grey', size = rel(0.3))) +
  scale_color_viridis_c(option = "viridis", limits = c(0, 1), 
                        breaks = c(0, 0.5, 1), oob = squish,
                        labels = c("0.0", 0.5, "1.0+")) + 
  xlab('Log2 average intensity') + 
  ylab('Log2 fold change from control') +
  labs(colour = "Medium / \nControl") +
  ggtitle(label = 'Abundance comparison from Control to Vehicle',
          subtitle = 'Treatment with Ox-LDL')
# LDL ----

# Repeating the steps for LDL 

attach(LDLdata_df)
LDL_treatment_ratio_ma_df <- 
  data.frame(Accession = Accession,
             Treatment_ratio = LDLTreatment / CTotal,
             LDLTreatment = LDLTreatment,
             M = M,
             A = A)
detach(LDLdata_df)

# Plot
LDL_treatment_ratio_colouring_ma <- 
  ggplot(data = LDL_treatment_ratio_ma_df, mapping = aes(x = A, y = M,
                                            colour = Treatment_ratio)) + 
  geom_point(data = subset(x = LDL_treatment_ratio_ma_df, 
                           subset = LDLTreatment > 0),
             aes(colour = Treatment_ratio),
             size = 0.6) +
  geom_point(data = subset(x = LDL_treatment_ratio_ma_df, 
                           subset = LDLTreatment == 0),
             shape = 2, size = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  theme_bw() + 
  # Use the same panel limits
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12, 14), limits = c(0,14)) + 
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4, 6, 8, 10), limits = c(-4,10)) + 
  theme(panel.grid = element_line(colour = 'grey', size = rel(0.3))) +
  scale_color_viridis_c(option = "viridis", limits = c(0, 1), 
                        breaks = c(0, 0.5, 1), oob = squish,
                        labels = c("0.0", 0.5, "1.0+")) + 
  xlab('Log2 average intensity') + 
  ylab('Log2 fold change from control') +
  labs(colour = "Vehicle / \ncontrol") +
  ggtitle(label = 'Abundance comparison from Control to Vehicle',
          subtitle = 'Treatment with LDL')


# Figure 3 ----

treatment_ratio_figure <- 
  arrangeGrob(nrow = 1, 
              ncol = 2, 
              oxLDL_treatment_ratio_colouring_ma + 
                theme(legend.position = "none",
                      panel.grid = element_line(colour = 'grey', 
                                                size = rel(1))) +  
                ggtitle(label = "Ox-LDL",
                        subtitle = NULL),
              LDL_treatment_ratio_colouring_ma + 
                theme(legend.position = "none", 
                      panel.grid = element_line(colour = 'grey', 
                                                size = rel(1))) +
                ggtitle(label = "LDL",
                        subtitle = NULL))

# Save this to file as part of Figure 3. 
ggsave(plot = treatment_ratio_figure,
       filename = './Images/treatment_ratio_figure.svg', 
       width = 14, height = 9,
       units = "cm",
       device = "svg",
       dpi = "print")

# Save dataframes ----

save(oxLDL_treatment_ratio_ma_df, LDL_treatment_ratio_ma_df, 
     file = 'Data/treatment_ratio_dfs.RData')

