# Packages ----

library('ggplot2')

# Source ---- 

source('./Initial_analysis_and_filter/00_raw_data_processing.R')

# Concept ----

# A histogram helps represent basic facts about the shape of the data. 
# A log2 histogram is used to demonstrate that, in general, abundance values are
# higher in the Treated Samples than in the Control. 

# Single histogram ----

# First preparing a dataframe using only relevant columns of information. 
# Data are divided on sample, to help demonstrate the size difference. 

attach(data_df)
hist_df <- rbind(data.frame(Total = OxTotal, 
                            Log2 = log2(OxTotal),
                            Sample = "OxLDL"),
                 data.frame(Total = LDLTotal,
                            Log2 = log2(LDLTotal),
                            Sample = "LDL"),
                 data.frame(Total = CTotal, 
                            Log2 = log2(CTotal),
                            Sample = "Control"))
detach(data_df)

# A histogram of the data transformed by log2
# A log2 transformation helps concatenate the data into more visible categories
# This is necessary because it is so tailed. 

log2_histogram <- 
  ggplot(data = hist_df, aes(x = Log2, fill = Sample)) +
  geom_histogram(show.legend = F) + 
  scale_fill_brewer(palette = "Dark2") +
  theme_bw() +
  facet_wrap(~ Sample,
             scales = 'fixed') +
  theme_light() +
  theme(strip.text = element_text(size = 13),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13)) + 
  xlab('Total protein abundance (Log2)') +
  ylab('Frequency') 
# For reference, here is the histogram using only the Total abundance 
# Untransformed.

total_histogram <- 
  ggplot(data = hist_df, aes(x = Total, fill = Sample)) +
  geom_histogram(show.legend = F) + 
  scale_fill_brewer(palette = "Dark2") +
  theme_bw() +
  facet_wrap(~ Sample,
             scales = 'fixed') +
  theme_light() +
  ylab('Frequency') + 
  xlab('Total protein intensity') + 
  ggtitle('Total protein expression between sample types')

# Facetted histogram v1. ----

# Note: this section facets the two plots into a single space. 
# This can also be performed using grid.arrange from the gridExtra package. 

attach(data_df)
hist_df_2 <- rbind(data.frame(Total = OxTotal, 
                              Sample = "Ox-LDL",
                              Data = "Total"),
                   data.frame(Total = LDLTotal,
                              Sample = "LDL",
                              Data = "Total"),
                   data.frame(Total = CTotal, 
                              Sample = "Control",
                              Data = "Total"),
                   data.frame(Total = log2(OxTotal), 
                              Sample = "Ox-LDL",
                              Data = "Log2"),
                   data.frame(Total = log2(LDLTotal), 
                              Sample = "LDL",
                              Data = "Log2"),
                   data.frame(Total = log2(CTotal), 
                              Sample = "Control",
                              Data = "Log2")
                   )
detach(data_df)
# Use factor to change the order in which Total and Log2 are plotted.
hist_df_2$Data <- factor(hist_df_2$Data, levels = c("Total", "Log2"))

# Plot the facetted histogram. 
combined_histogram <- ggplot(data = hist_df_2, aes(x = Total)) +
  geom_histogram() + 
  theme_bw() +
  facet_wrap(facets = ~Data + Sample,
             nrow = 2,
             ncol = 3,
             scales = 'free') +
  theme_grey() + 
  ylab('Frequency') + 
  ggtitle('Total protein expression between sample types')