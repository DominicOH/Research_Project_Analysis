# Packages ----
library('gridExtra')
library('svglite')

# Source ----

source('./Initial_analysis_and_filter/01_3_n_trials_conditional_colouring.R', 
       echo = F)
source('./Initial_analysis_and_filter/01_4_p_value_conditional_colouring.R', 
       echo = F)
# OxLDLdata_df is overwritten if the above two source files are run after the 
# first. 
source('./Initial_analysis_and_filter/01_2_ratio_colouring_with_candidates.R', 
       echo = F)

# Concept ----

# Combining all current conditional models to produce a single figure 
# helps to visualise the points likelier to be genuine changes, as well as
# those that are considered likely contaminants. 

# Final figure panels 


conditionals <- 
  arrangeGrob(nrow = 3, ncol = 2, 
              oxLDL_treatment_ratio_colouring_ma +
                ggtitle(label = NULL, subtitle = NULL) + 
                theme(legend.position = "none",
                      panel.grid = element_line(colour = 'grey', size = rel(1))), 
              LDL_treatment_ratio_colouring_ma +
                ggtitle(label = NULL, subtitle = NULL) +
                theme(legend.position = "none",
                      panel.grid = element_line(colour = 'grey', size = rel(1))),
              oxLDL_N_detections_ma + 
                ggtitle(label = NULL, subtitle = NULL) +
                theme(legend.position = "none", 
                      panel.grid = element_line(colour = 'grey', size = rel(1))),
              LDL_N_detections_ma +
                ggtitle(label = NULL, subtitle = NULL) +
                theme(legend.position = "none",
                      panel.grid = element_line(colour = 'grey', size = rel(1))),
              oxLDL_P_value_ma + 
                ggtitle(label = NULL, subtitle = NULL) +
                theme(legend.position = "none",
                      panel.grid = element_line(colour = 'grey', size = rel(1))),
              LDL_P_value_ma + 
                ggtitle(label = NULL, subtitle = NULL) +
                theme(legend.position = "none",
                      panel.grid = element_line(colour = 'grey', size = rel(1)))
  )
# Final figure panels with candidates ----

theme_update(axis.text = element_text(size = 13))

conditionals_with_candidates <- 
  arrangeGrob(nrow = 3, ncol = 2, 
             oxLDL_treatment_ratio_colouring_ma_candidates +
               ggtitle(label = NULL, subtitle = NULL) + 
               theme(legend.position = "none",
                     panel.grid = element_line(colour = 'grey', size = rel(1))), 
             LDL_treatment_ratio_colouring_ma_candidates +
               ggtitle(label = NULL, subtitle = NULL) +
               theme(legend.position = "none",
                     panel.grid = element_line(colour = 'grey', size = rel(1))),
             oxLDL_N_detections_ma_with_candidates + 
               ggtitle(label = NULL, subtitle = NULL) +
               theme(legend.position = "none", 
                     panel.grid = element_line(colour = 'grey', size = rel(1))),
             LDL_N_detections_ma_with_candidates +
               ggtitle(label = NULL, subtitle = NULL) +
               theme(legend.position = "none",
                     panel.grid = element_line(colour = 'grey', size = rel(1))),
             oxLDL_P_value_ma_with_candidates + 
               ggtitle(label = NULL, subtitle = NULL) +
               theme(legend.position = "none",
                     panel.grid = element_line(colour = 'grey', size = rel(1))),
             LDL_P_value_ma_with_candidates + 
               ggtitle(label = NULL, subtitle = NULL) +
               theme(legend.position = "none",
                     panel.grid = element_line(colour = 'grey', size = rel(1)))
)
ggsave(plot = conditionals,
                filename = './Images/conditionals.png', 
                width = 14, height = 22,
                units = "cm",
                device = "png",
                dpi = "print")

ggsave(plot = conditionals,
                filename = './Images/conditionals.svg', 
                width = 14, height = 22,
                units = "cm",
                device = "svg",
                dpi = "print")