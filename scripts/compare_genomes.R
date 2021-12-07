stats_MAGS = read.table('~/Downloads/checkm_stats_MAGS.tsv',
                    col.names = c("Bin", "Marker_lineage", "Marker_lineage_code", "num_genomes", "num_markers", "num_marker_sets", "0", "1", "2", "3", "4", "5_plus", "Completeness", "Contamination", "Strain_heterogeneity"))


stats_NCBI = read.table('~/Downloads/checkm_stats_NCBI.tsv',
                    col.names = c("Bin", "Marker_lineage", "Marker_lineage_code", "num_genomes", "num_markers", "num_marker_sets", "0", "1", "2", "3", "4", "5_plus", "Completeness", "Contamination", "Strain_heterogeneity"))


genome_length_NCBI = read.table('~/Downloads/total_genome_length_NCBI.tsv.clean',
                    col.names = c("Bin", "Length"))

genome_length_MAGS = read.table('~/Downloads/total_genome_length_MAGS.tsv.clean',
                    col.names = c("Bin", "Length"))

genome_taxon_MAGS = read.table('~/Downloads/bin_taxonomy.tsv',
                                col.names = c("Bin", "Family"))

#copy the Marker_lineage in stats_NCBI into a Family column

#concatenate the two data frames
all_data <-stats_MAGS  %>% 
  rbind(stats_NCBI)

#concatenate genome length

#join genome length to all_data 
left_join(table1, by = 'Bin_ID')

#subset 
subset(all_data, select ="Bin", "Family", "Length", "Completeness", "Contamination")



#sort on Family
data_2 <- data_2_a %>% arrange(Family)

#output a table
color based on family?
scale_fill_manual(values = c("steelblue2", "gray60", "orchid3", "orange2", "tomato2", "yellowgreen"))
