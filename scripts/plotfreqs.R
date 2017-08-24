library(tidyverse)

get_polymorphisms <- function(filepath) {
  polymorphisms <- read.csv(filepath) %>% 
    unite(polymorphism, position, entry_type, mutation_detail, sep = ' ') %>%
    select(genome_id, polymorphism, generation, polymorphism_frequency) %>% 
    complete(nesting(polymorphism, genome_id), generation, fill = list(polymorphism_frequency = 0))
  return(polymorphisms)
}

plot_polymorphism_frequencies <- function(polymorphism_tibble, output_dir) {
  ggplot(data = polymorphism_tibble, mapping = aes(x = generation, y = polymorphism_frequency)) + 
    geom_line(aes(group = polymorphism)) + 
    facet_wrap(~ genome_id, ncol = 1)
  ggsave(paste(output_dir, 'polymorphism_plot.png', sep = ''), last_plot())
}