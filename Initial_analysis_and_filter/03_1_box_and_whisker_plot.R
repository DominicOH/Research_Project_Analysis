# Packages ----

library('ggplot2')
library('gridExtra')

# Source ----

source('./Initial_analysis_and_filter/03_0_protein_intensity_histogram.R')

# Concept ----

# Using the same dataframes needed to produce a histogram, it's also possible
# to present the data in the form of a boxplot. An advantage of this is making
# visual comparison between the median, mean, and percentile statistics of each 
# sample. Furthermore, this helps show the tailed-ness of the data better than 
# a histogram. 

# Plot ----
log2_boxplot <- 
  ggplot(data = hist_df, mapping = aes(x = Log2, 
                                       y = Sample,
                                       colour = Sample)) +
  geom_boxplot(outlier.shape = 1,
               outlier.size = 2,
               show.legend = F,
               lwd = 1,
               width = 0.4) + 
  scale_color_brewer(type = "qual", 
                     palette = "Dark2",
                     aesthetics = "colour") + 
  theme_light() + 
  theme(axis.title = element_text(size = 15),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text = element_text(size = 15)) +
  xlab(label = 'Total protein abundance (Log2)')

# Arranged with histogram ----

grid.arrange(bottom = "Log2 total abundance",
             log2_boxplot +
               xlab(label = NULL) +
               ylab(label = NULL), 
             log2_histogram +
               theme(panel.grid = element_line(size = 2)) +
               ggtitle(label = NULL) +
               xlab(label = NULL) +
               theme(strip.text.x = element_blank()),
             nrow = 2)
