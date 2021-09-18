# Source ----

source('./Initial_analysis_and_filter/01_0_microarray_representation_of_data.R')

# Concept ----

# Each protein has been detected a certain number of times in the data. 
# Some have been missed during mass spec. 
# The more detections, the more reliable the data point can be considered to be
# as fewer detections could fall vulnerable to bias by random variance. 

# Detection dataframe ----

# In the unedited dataset, undetected proteins are marked as NA.  
# Counting the number of cells that are not NA gives the number of detections. 
# This is repeated for each sample type. 

attach(data_df)
n_detections_df <- 
  data.frame(Accession = Accession,
             Cdetections = rowSums(!is.na(Data_unedited[3:10])),
             CTotal = CTotal,
             Ox_detections = rowSums(!is.na(Data_unedited[11:16])),
             Ox_Treat_detections = rowSums(!is.na(Data_unedited[17:19])),
             OxTotal = OxTotal,
             LDL_detections = rowSums(!is.na(Data_unedited[20:25])),
             LDL_Treat_detections = rowSums(!is.na(Data_unedited[26:28])),
             LDLTotal = LDLTotal)
detach(data_df)

# The average of these detection counts gives a scale, with higher meaning more
# overall detections, and lower meaning fewer overall detections. 

n_detections_df$Ox_detections_Mean <- rowMeans(n_detections_df[c(2, 4:5)])
n_detections_df$LDL_detections_Mean <- rowMeans(n_detections_df[c(2, 7:8)])

# Ox-LDL Plot ----

# Since the n_detection calculations were performed to the whole dataset, 
# MA calculations need to be repeated for Ox-LDL and LDL using the functions. 
oxLDL_N_detections_ma_df <- ox_define_MA(n_detections_df)
oxLDL_N_detections_ma_df$Ox_detections_Mean <- n_detections_df$Ox_detections_Mean
oxLDL_N_detections_ma_df$Accession <- data_df$Accession

# -Inf or Inf results are mapped inappropriately at corners of the plot. 
# These should be filtered out to avoid confusion. 
# Points of less than 1 detection don't appear on the plot either but they do 
# affect the scale, so these are filtered out as well. 

oxLDL_N_detections_ma_df <- oxLDL_N_detections_ma_df %>%
  filter(oxLDL_N_detections_ma_df$Ox_detections_Mean > 0.5 
         & is.finite(M) 
         & is.finite(A)) 

# Plotting as an MA
oxLDL_N_detections_ma <- 
  ggplot(data = oxLDL_N_detections_ma_df, 
       mapping = aes(x = A, y = M, 
                     colour = Ox_detections_Mean)) + 
  geom_point(size = 0.6) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  theme_bw() + 
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12, 14), limits = c(0,14)) + 
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4, 6, 8, 10), limits = c(-4,10)) + 
  theme(panel.grid = element_line(colour = 'grey', size = rel(0.3))) + 
  scale_color_viridis_c(option = "plasma", direction = -1) + 
  xlab('Log2 average intensity') + 
  ylab('Log2 fold change from control') +
  labs(colour = "Detections \n(mean)") +
  ggtitle(label = 'Number of detections', 
          subtitle = 'Treatment with Ox-LDL')

# Plotting as an MA with the candidate proteins
oxLDL_N_detections_ma_with_candidates <- 
  oxLDL_N_detections_ma + 
  geom_point(data = subset(x = oxLDL_N_detections_ma_df, 
                           subset = Accession == 'O00622' # CCN1 
                           | Accession == 'A0A024R462' # FN1 
                           | Accession == 'P03956'), # MMP1
             size = 4,
             mapping = aes(shape = Accession)) + 
  labs(shape = "Candidate \nproteins")

# LDL Plot ----

# Steps are repeated to LDL as above. 

LDL_N_detections_ma_df <- ldl_define_MA(n_detections_df)
LDL_N_detections_ma_df$LDL_detections_Mean <- n_detections_df$LDL_detections_Mean
LDL_N_detections_ma_df$Accession <- data_df$Accession

LDL_N_detections_ma_df <- LDL_N_detections_ma_df %>%
  filter(LDL_N_detections_ma_df$LDL_detections_Mean > 0.5 
         & is.finite(M) 
         & is.finite(A)) 

# Plotting as an MA

LDL_N_detections_ma <- 
  ggplot(data = LDL_N_detections_ma_df, 
         mapping = aes(x = A, y = M, 
                       colour = LDL_detections_Mean)) + 
  geom_point(size = 0.6) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  theme_bw() + 
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12, 14), limits = c(0,14)) + 
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4, 6, 8, 10), limits = c(-4,10)) + 
  theme(panel.grid = element_line(colour = 'grey', size = rel(0.3))) + 
  scale_color_viridis_c(option = "plasma", direction = -1) + 
  xlab('Log2 average intensity') + 
  ylab('Log2 fold change from control') +
  labs(colour = "Detections \n(mean)") +
  ggtitle(label = 'Number of detections', 
          subtitle = 'Treatment with LDL')

# Plotting as an MA with the candidate proteins
LDL_N_detections_ma_with_candidates <- 
  LDL_N_detections_ma + 
  geom_point(data = subset(x = LDL_N_detections_ma_df, 
                           subset = Accession == 'O00622' # CCN1 
                           | Accession == 'A0A024R462' # FN1 
                           | Accession == 'P03956'), # MMP1
             size = 4,
             mapping = aes(shape = Accession)) + 
  labs(shape = "Candidate \nproteins")

# Figure 4 ----

N_detections_figure <- 
  arrangeGrob(nrow = 1, 
              ncol = 2, 
              oxLDL_N_detections_ma + 
                theme(legend.position = "none",
                      panel.grid = element_line(colour = 'grey', 
                                                size = rel(1))) +  
                ggtitle(label = "Ox-LDL",
                        subtitle = NULL),
              LDL_N_detections_ma + 
                theme(legend.position = "none", 
                      panel.grid = element_line(colour = 'grey', 
                                                size = rel(1))) +
                ggtitle(label = "LDL",
                        subtitle = NULL))

# Save this to file as part of Figure 4. 
ggsave(plot = N_detections_figure,
       filename = './Images/N_detections_figure.svg', 
       width = 14, height = 9,
       units = "cm",
       device = "svg",
       dpi = "print")

# Save dataframe ----

save(oxLDL_N_detections_ma_df, LDL_N_detections_ma_df, 
     file = 'Data/n_detections_dfs.RData')
