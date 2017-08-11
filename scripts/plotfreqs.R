library(tidyverse)

csv_to_frequency_plot <- function(filepath) {
  freqs <- read.csv(filepath)
  polymorphisms <- freqs %>% 
    unite(polymorphism, position, entry_type, mutation_detail, sep = ' ') %>%
    select(genome_id, polymorphism, generation, polymorphism_frequency) %>% 
    complete(nesting(polymorphism, genome_id), generation, fill = list(polymorphism_frequency = 0))
  write.csv(polymorphisms, file = paste(filepath, 'EDIT.csv', sep = '_'), row.names = FALSE, quote = FALSE)
  ggplot(data = polymorphisms, mapping = aes(x = generation, y = polymorphism_frequency)) + 
    geom_line(aes(group = polymorphism)) + 
    facet_wrap(~ genome_id, ncol = 1)
  ggsave(paste(filepath, 'plot.png', sep = '_'), last_plot())
}