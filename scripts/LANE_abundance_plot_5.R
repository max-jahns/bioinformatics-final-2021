#Create Plot of Abundance of Bins in Metagenome

library(tidyverse)
library(magrittr)


#import bins and counts
data_1 = read.table('~/Downloads/total_bin_counts.tsv',
                  col.names = c("Bin", "Counts"))

#import bin taxonomy

genome_taxon_MAGS = read.table('~/Downloads/bin_taxonomy.tsv',
                               col.names = c("Bin", "Family")
)

#match bins/counts with family
data_2_a <- genome_taxon_MAGS %>%
  left_join(data_1, by = "Bin")


#sort data_2 on Family
data_2 <- data_2_a %>% arrange(Family)


#Make Donut Plot 

# Compute percentages
data_2$fraction = data_2$Counts / sum(data_2$Counts)

# Compute the cumulative percentages (top of each rectangle)
data_2$ymax = cumsum(data_2$fraction)

# Compute the bottom of each rectangle
data_2$ymin = c(0, head(data_2$ymax, n=-1))

data_2 <- data_2 %>% 
  mutate(final_family = paste0(Family, ' (' ,round(fraction*100, digits = 2), '%)'))

data_3 <- data_2 %>% 
  group_by(Family) %>% 
  mutate(total_fraction = sum(fraction)) %>% 
  ungroup() %>% 
  mutate(final_family = paste0(Family, ' (' ,round(total_fraction*100, digits = 2), '%)'))

#Plot it
ggplot(data_3, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=as_factor(final_family))) +
  geom_rect() +
  coord_polar(theta="y") +
  xlim(c(2, 4)) + 
  labs(fill = "Family", title = "Relative Abundance of Families in Metagenome") + 
  theme_void() +
  scale_fill_manual(values = c("steelblue2", "gray60", "orchid3", "orange2", "tomato2", "yellowgreen"))

