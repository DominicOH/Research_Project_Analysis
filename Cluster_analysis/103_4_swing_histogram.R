# Packages ----

library('tidyverse')
library('factoextra')

# Source ----

load('Data/swing_dend_data.RData')

# Concept ----

# Making a histogram out of the swing data. Point is to show how much of the 
# samples are apical swinging. 

# Dataframe ----

# turning the swing dataframe with fold changes into just a dataframe with swing
# The name of the sample is also needed to split the two up. 
swing_hist_df <- swing_fc_df %>%
  select(Ox_LDL_swing, LDL_swing) %>%
  pivot_longer(cols = c(Ox_LDL_swing, LDL_swing),
               names_to = 'Treated sample',
               values_to = 'Swing')

# Histogram ----

swing_hist_df$`Treated sample`[swing_hist_df$`Treated sample` == "Ox_LDL_swing"] <-
  "Ox-LDL"

swing_hist_df$`Treated sample`[swing_hist_df$`Treated sample` == "LDL_swing"] <-
  "LDL"

png(filename = 'Images/swing_histogram.png', units = "cm", 
    width = 16, height = 8, res = 200)

ggplot(data = swing_hist_df, aes(fill = `Treated sample`)) + 
  geom_histogram(mapping = aes(x = Swing)) +
  scale_x_continuous(limits = c(-1, 1)) +
  facet_grid(~ `Treated sample`) +
  scale_fill_brewer(type = "div", palette = "Dark2") +
  theme_bw() +
  theme(legend.position = "none", ) + 
  labs(x = "Apical polarisation change", y = "Count") + 
  geom_vline(xintercept = 0, lty = "dashed")

  dev.off()

