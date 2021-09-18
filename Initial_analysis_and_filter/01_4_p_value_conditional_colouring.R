# Source ----

source('./Initial_analysis_and_filter/01_0_microarray_representation_of_data.R')

# Concept ----

# Fisher's Exact Test is a statistical tool in hypothesis testing. 
# The null hypothesis assumes that, between two states, there has been no 
# significant change beyond variance or noise. 

# A P-value gives the probability that a particular data point can be explained
# by the null hypothesis. A P-value closer to 1.0 means that a statistic could
# readily be the result of the null hypothesis, and any difference could just be
# variance from the expected value. 

# Two sets of P-values can be produced for the these data. The original
# state is represented by the Control, and Treated Samples reflect either test
# set. 

# P-value dataframe ---- 

# T-testing uses each or trial or repeat, taking two vectors of Control and 
# Test. This means the unedited dataset is needed. 
# For the sake of T-testing, it's much more convenient to reduce NA values to 0.
# NOTE: A value of 0 doesn't mean that the abundance was actually 0, but that 
# detection was uncertain to the extent that an abundance couldn't be produced. 
Data_unedited[is.na(Data_unedited)] <- 0

attach(Data_unedited)
# data frame combining apical and basolateral samples by trial 
p_data_df <- data.frame(Description = Description,
                        Accession = Accession,
                        C_1_Total = CA1 + CB1,
                        C_2_Total = CA2 + CB2,
                        C_3_Total = CA3 + CB3,
                        C_4_Total = CA4 + CB4,
                        OxTot_1_Total = OXA1 + OXB1,
                        OxTot_2_Total = OXA2 + OXB2,
                        OxTot_3_Total = OXA3 + OXB3,
                        OxTreatment = data_df$OxTreatment,
                        LDLTot_2_Total = LDLA2 + LDLB2,
                        LDLTot_3_Total = LDLA3 + LDLB3, 
                        LDLTot_4_Total = LDLA4 + LDLB4,
                        LDLTreatment = data_df$LDLTreatment)
detach(Data_unedited)

# in this case, the null hypothesis is that there will be no change from C to Ox
# i.e. H0: C = Ox
# the alternative hypothesis is that Ox will be different to C 
# i.e. Ha: C != Ox 

# T-Test Function ----

# By default, the t.test function takes COLUMNS of a dataframe as Control and
# Test statistics, and produces a LIST object containing the p-value. 
# For iterative use against a dataframe, this function performs the t.test and
# extracts only the P-value. 

t_test_p_value <- function(x, y) {
  # t.test requires a minimum of 2 'x' data values
  test <- t.test(x = x, 
                 y = y, 
                 paired = FALSE, 
                 alternative = "two.sided")
  
  return(test$p.value) # extracts just the p-value statistic
  
}

# Ox-LDL ----

# To avoid overwriting the original df, a new one is made for the Ox-LDL MA.  
# This is because a new P-value is needed for both Ox-LDL and LDL. 
oxLDL_p_ma_df <- p_data_df

# Since this is a dataframe. The function must be applied as a loop iterating 
# by row. 
attach(oxLDL_p_ma_df)
for (i in 1:nrow(oxLDL_p_ma_df)) {
  oxLDL_p_ma_df$ox_p_value[i] <- t_test_p_value(x = c(C_1_Total[i], C_2_Total[i], 
                                                     C_3_Total[i], C_4_Total[i]),
                                               y = c(OxTot_1_Total[i], 
                                                     OxTot_2_Total[i], 
                                                     OxTot_3_Total[i]))
}

detach(oxLDL_p_ma_df)

# LDL ----

# Above repeated for LDL. 

LDL_p_ma_df <- p_data_df
LDL_p_ma_df[is.na(LDL_p_ma_df) == TRUE] <- 0

attach(LDL_p_ma_df)
for (i in 1:nrow(LDL_p_ma_df)) {
  LDL_p_ma_df$ldl_p_value[i] <- t_test_p_value(x = c(C_1_Total[i], C_2_Total[i], 
                                                       C_3_Total[i], C_4_Total[i]),
                                                 y = c(LDLTot_2_Total[i], 
                                                       LDLTot_3_Total[i], 
                                                       LDLTot_4_Total[i]))
}
detach(LDL_p_ma_df)

# Ox-LDL MA ----

# This can now be used to overlay an MA plot as conditional colouring.
oxLDL_p_ma_df$OxTotal <- (rowSums(oxLDL_p_ma_df[7:9])/3)
oxLDL_p_ma_df$CTotal <- (rowSums(oxLDL_p_ma_df[3:6])/4)

# Using the define MA function to retrieve M and A values for the p-value set. 
ox_p_data_ma <- ox_define_MA(oxLDL_p_ma_df)
oxLDL_p_ma_df$M <- ox_p_data_ma$M; oxLDL_p_ma_df$A <- ox_p_data_ma$A

# As before, INF and -INF values caused by log transformation should be removed. 
oxLDL_p_ma_df <- oxLDL_p_ma_df %>%
  filter(is.finite(M) & is.finite(A)) 

oxLDL_P_value_ma <- 
  ggplot(data = oxLDL_p_ma_df,
                            aes(x = A, y = M,
                                colour = ox_p_value)) + 
  geom_point(size = 0.6) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  theme_bw() + 
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12, 14), limits = c(0,14)) + 
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4, 6, 8, 10), limits = c(-4,10)) + 
  theme(panel.grid = element_line(colour = 'grey', size = rel(0.3))) +
  scale_color_viridis_c(option = "plasma", direction = 1) + 
  xlab('Log2 average intensity') + 
  ylab('Log2 fold change from control') +
  labs(colour = "P-value") +
  ggtitle(label = 'Distribution of P-values',
          subtitle = 'Treatment with Ox-LDL')

# With the candidate proteins. 
oxLDL_P_value_ma_with_candidates <- 
  oxLDL_P_value_ma +
  geom_point(data = subset(x = oxLDL_p_ma_df, 
                           subset = Accession == 'O00622' # CCN1 
                           | Accession == 'A0A024R462' # FN1 
                           | Accession == 'P03956'), # MMP1
             size = 4,
             mapping = aes(shape = Accession)) + 
  labs(shape = "Candidate \nproteins")


# LDL MA ----

# Repeat of the above for LDL sample. 

LDL_p_ma_df$LDLTotal <- (rowSums(LDL_p_ma_df[11:13])/3)
LDL_p_ma_df$CTotal <- (rowSums(LDL_p_ma_df[3:6])/4)

ldl_p_data_ma <- ldl_define_MA(LDL_p_ma_df)
LDL_p_ma_df$M <- ldl_p_data_ma$M; LDL_p_ma_df$A <- ldl_p_data_ma$A

LDL_p_ma_df <- LDL_p_ma_df %>%
  filter(is.finite(M) 
         & is.finite(A)) 

# Plot
LDL_P_value_ma <- 
  ggplot(data = LDL_p_ma_df,
         aes(x = A, y = M,
             colour = ldl_p_value)) + 
  geom_point(size = 0.6) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  theme_bw() + 
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12, 14), limits = c(0,14)) + 
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4, 6, 8, 10), limits = c(-4,10)) + 
  theme(panel.grid = element_line(colour = 'grey', size = rel(0.3))) +
  scale_color_viridis_c(option = "plasma", direction = 1) + 
  xlab('Log2 average intensity') + 
  ylab('Log2 fold change from control') +
  labs(colour = "P-value") +
  ggtitle(label = 'Distribution of P-values',
          subtitle = 'Treatment with LDL')

# Plot with candidates
LDL_P_value_ma_with_candidates <- 
  LDL_P_value_ma +
  geom_point(data = subset(x = LDL_p_ma_df, 
                           subset = Accession == 'O00622' # CCN1 
                           | Accession == 'A0A024R462' # FN1 
                           | Accession == 'P03956'), # MMP1
             size = 4,
             mapping = aes(shape = Accession)) + 
  labs(shape = "Candidate \nproteins")

# Figure 5 ----

P_value_figure <- 
  arrangeGrob(nrow = 1, 
              ncol = 2, 
              oxLDL_P_value_ma + 
                theme(legend.position = "none",
                      panel.grid = element_line(colour = 'grey', 
                                                size = rel(1))) +  
                ggtitle(label = "Ox-LDL",
                        subtitle = NULL),
              LDL_P_value_ma + 
                theme(legend.position = "none", 
                      panel.grid = element_line(colour = 'grey', 
                                                size = rel(1))) +
                ggtitle(label = "LDL",
                        subtitle = NULL))

# Save this to file as part of Figure 4. 
ggsave(plot = P_value_figure,
       filename = './Images/P_value_figure.svg', 
       width = 14, height = 9,
       units = "cm",
       device = "svg",
       dpi = "print")

# Save dataframes ----

save(oxLDL_p_ma_df, LDL_p_ma_df, file = 'Data/P_value_dfs.RData')
