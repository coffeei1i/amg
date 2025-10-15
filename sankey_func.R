
library(dplyr)
library(ggplot2)
library(ggalluvial)


create_sankey_diagram <- function(data, highlight_cell_types, sankey_colors, plot_title) {

  base_data <- data %>%
    mutate(group = interaction(cell_type, brain_region)) %>% 
    arrange(cell_type, brain_region)  
  
 
  all_data <- base_data %>%
    filter(!(cell_type %in% highlight_cell_types)) %>%
    mutate(alpha = 0.01, width = 0.5)
  

  highlight_data <- base_data %>%
    filter(cell_type %in% highlight_cell_types) %>%
    mutate(alpha = 1, width = 0.5)
  
  
  combined_data <- bind_rows(all_data, highlight_data)
  

  p <- ggplot(data = combined_data, aes(axis1 = cell_type, axis2 = brain_region, y = freq, group = group)) +
    geom_alluvium(aes(fill = cell_type, width = width, alpha = alpha)) +
    geom_stratum(width = 0.1, alpha = 1) + 
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, check_overlap = FALSE) +
    scale_fill_manual(values = sankey_colors) +  
    theme_minimal() +  
    theme(
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),  
      panel.background = element_blank(),  
      legend.position = "none"  
    ) +
    ggtitle(plot_title)  
  

  return(p)
}

