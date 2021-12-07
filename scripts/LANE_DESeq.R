install.packages("BiocManager")
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("ggtern")
BiocManager::install("DESeq2")

library(DESeq2)
library(ggplot2)
library(gridExtra)
library(ggtern)

#data = as.matrix(read.csv("~/Downloads/deseq_input.tsv", sep=",", row.names = "Contig"))
data = read.csv("~/Downloads/deseq_input.tsv", row.names = "contigName" )

cdata = read.csv("~/Downloads/coldat.csv", row.names = "sample")

d <- as.matrix(data)
c <- as.matrix(cdata)

#Chose to run a deseq analysis across all three conditions: reseeded_24hrs, reseeded_05hrs, and control
dds = DESeqDataSetFromMatrix(countData = round(d), colData = c, design = ~ condition)
dds <- DESeq(dds)
res = results(dds, alpha = 0.05)
#so yeah, the factor column, then the “numerator” followed by “denominator”, denominator should be the control numerator the treatment
res_1 = results(dds, contrast = c('condition', 'reseeded_05', 'reseeded_24'), alpha = 0.05)
res_2 = results(dds, contrast = c('condition', 'control', 'reseeded_05'), alpha = 0.05)
res_3 = results(dds, contrast = c('condition', 'control', 'reseeded_24'), alpha = 0.05)


write.csv(res_1, file = "~/Downloads/DE_results_r05_r24.csv")
write.csv(res_2, file = "~/Downloads/DE_results_pairwise_c_r05.csv")
write.csv(res_3, file = "~/Downloads/DE_results_pairwise_c_r24.csv")



#deseq_res = results(dds, tidy=TRUE)

#Make the Triangles

#bring_back_deseq
stats = read.csv('~/Downloads/DE_results.csv')
r05_r24 = read.csv("~/Downloads/DE_results_r05_r24.csv")
c_r05 = read.csv("~/Downloads/DE_results_pairwise_c_r05.csv")
c_r24 = read.csv("~/Downloads/DE_results_pairwise_c_r24.csv")

#subset into dataframes for each Family
genome_taxon_MAGS = read.table('~/Downloads/bin_taxonomy.tsv',
                               col.names = c("Bin_ID", "Family"))
data_1 <- stats %>%
  separate(X, into = c('Bin_ID', 'scaffold'), sep = '_k')  #pipe it in instead of creating intermediate object, up to you to create intermediate opject

data_rr <- r05_r24 %>%
  separate(X, into = c('Bin_ID', 'scaffold'), sep = '_k')  #pipe it in instead of creating intermediate object, up to you to create intermediate opject

#data_rr <- data_rr %>% 
#  mutate(linearChange = (2^log2FoldChange) / 100)

data_cr5 <- c_r05 %>%
  separate(X, into = c('Bin_ID', 'scaffold'), sep = '_k')  #pipe it in instead of creating intermediate object, up to you to create intermediate opject

#data_cr5 <- data_cr5 %>% 
#  mutate(linearChange = (2^log2FoldChange) / 100)

data_cr24 <- c_r24 %>%
  separate("X", into = c("Bin_ID", "scaffold"), sep = '_k')  #pipe it in instead of creating intermediate object, up to you to create intermediate opject

#data_cr24 <- data_cr24 %>% 
#  mutate(linearChange = (2^log2FoldChange) / 100)



#add the family name
data_2 <- data_1 %>%
  left_join(genome_taxon_MAGS, by = 'Bin_ID')

data_rr <- data_rr %>%
  left_join(genome_taxon_MAGS, by = 'Bin_ID')

data_cr5 <- data_cr5 %>%
  left_join(genome_taxon_MAGS, by = 'Bin_ID')

data_cr24 <- data_cr24 %>%
  left_join(genome_taxon_MAGS, by = 'Bin_ID')


#colnames(data_1)[1] == colnames(genome_taxon_MAGS)[1]
#class(data_1)
#class(genome_taxon_MAGS)
#save(data_1, genome_taxon_MAGS, file = '~/Downloads/temp.rdata')

data_2 %>% 
  #filter(padj <= 0.05) %>%  # use  this instead of filtering pvalue to filter adjusted pvalues instead
  #filter(!is.na(Family)) %>%  #uncomment if you want to remove 'NA' Family rows
  ggplot(aes(log2FoldChange, padj)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~ Family) 

#data_2 %>% 
  #mutate(linearChange = (2^log2FoldChange) / 100) %>% 
  #filter(pvalue <= 0.05) %>%
  #filter(padj <= 0.05) %>%  # use  this instead of filtering pvalue to filter adjusted pvalues instead
  #filter(!is.na(Family)) %>%  #uncomment if you want to remove 'NA' Family rows
  #ggplot(aes(linearChange, pvalue) +
           #geom_point() +
           #theme_bw() +
           #facet_wrap(~ Family)

#turn the log2fold change into a range from 0 to 1, remove all NA families
data_3 <- data_2 %>% 
  mutate(linearChange = (2^log2FoldChange) / 100)
  #filter(!is.na(Family))  #uncomment if you want to remove 'NA' Family rows

#subset based on family
#data_z = data_3 %>%
#  as_tibble() %>% 
#  group_by(Family) %>% 
#  group_split()



#data_z = data_z %>% 
# set_names(nm = c(map_chr(data_z, ~ unique(.x$Family)) %>% 
#                     na.omit(), 'None'))


data_rr_A = subset(data_rr, Family == "Alteromonadaceae")
data_rr_A <- data_rr_A %>% 
  mutate(linearChange = (2^log2FoldChange) / 100)

data_cr5_A = subset(data_cr5, Family == "Alteromonadaceae")
data_cr5_A <- data_cr5_A %>% 
  mutate(linearChange = (2^log2FoldChange) / 100)

data_cr24_A = subset(data_cr24, Family == "Alteromonadaceae")
data_cr24_A <- data_cr24_A %>% 
  mutate(linearChange = (2^log2FoldChange) / 100)


data_rr_F = subset(data_rr, Family == "Flavobacteriaceae")
data_rr_F <- data_rr_F %>% 
  mutate(linearChange = (2^log2FoldChange) / 100)

data_cr5_F = subset(data_cr5, Family == "Flavobacteriaceae")
data_cr5_F <- data_cr5_F  %>% 
  mutate(linearChange = (2^log2FoldChange) / 100)

data_cr24_F = subset(data_cr24, Family == "Flavobacteriaceae")
data_cr24_F <- data_cr24_F %>% 
  mutate(linearChange = (2^log2FoldChange) / 100)


data_rr_L = subset(data_rr, Family == "Litoricolaceae")
data_rr_L <- data_rr_L %>% 
  mutate(linearChange = (2^log2FoldChange) / 100)

data_cr5_L = subset(data_cr5, Family == "Litoricolaceae")
data_cr5_L <- data_cr5_L %>% 
  mutate(linearChange = (2^log2FoldChange) / 100)
data_cr24_L = subset(data_cr24, Family == "Litoricolaceae")
data_cr24_L <- data_cr24_L %>% 
  mutate(linearChange = (2^log2FoldChange) / 100)

data_rr_N = subset(data_rr, Family == "Nitrincolaceae")
data_rr_N <- data_rr_N %>% 
  mutate(linearChange = (2^log2FoldChange) / 100)
data_cr5_N = subset(data_cr5, Family == "Nitrincolaceae")
data_cr5_N <- data_cr5_N %>% 
  mutate(linearChange = (2^log2FoldChange) / 100)
data_cr24_N = subset(data_cr24, Family == "Nitrincolaceae")
data_cr24_N <- data_cr24_N %>% 
  mutate(linearChange = (2^log2FoldChange) / 100)

data_rr_R = subset(data_rr, Family == "Rhodobacteraceae")
data_rr_R <- data_rr_R %>% 
  mutate(linearChange = (2^log2FoldChange) / 100)
data_cr5_R = subset(data_cr5, Family == "Rhodobacteraceae")
data_cr5_R <- data_cr5_R %>% 
  mutate(linearChange = (2^log2FoldChange) / 100)
data_cr24_R = subset(data_cr24, Family == "Rhodobacteraceae")
data_cr24_R <- data_cr24_R %>% 
  mutate(linearChange = (2^log2FoldChange) / 100)

data_rr_S = subset(data_rr, Family == "Sphingomonadaceae")
data_rr_S <- data_rr_S %>% 
  mutate(linearChange = (2^log2FoldChange) / 100)
data_cr5_S = subset(data_cr5, Family == "Sphingomonadaceae")
data_cr5_S <- data_cr5_S %>% 
  mutate(linearChange = (2^log2FoldChange) / 100)
data_cr24_S = subset(data_cr24, Family == "Sphingomonadaceae")
data_cr24_S <- data_cr24_S %>% 
  mutate(linearChange = (2^log2FoldChange) / 100)



scale_fill_manual(values = c("steelblue2", "gray60", "orchid3", "orange2", "tomato2", "yellowgreen"))

A = ggtern(data_rr_A, aes(linearChange,data_cr5_A$linearChange,data_cr24_A$linearChange)) +
  geom_point() +
  theme_hidetitles() +
  theme(panel.background = element_rect(fill = "steelblue2")) +
  labs(title = "Alteromonadaceae") +
  theme(plot.title=element_text(hjust=0.5))

  
  #labs(x = "Control",y="Reseeding 05hrs",z="Reseeding 24hrs")
A
  
F = ggtern(data_rr_F, aes(linearChange,data_cr5_F$linearChange,data_cr24_F$linearChange)) +
  geom_point() +
  theme(panel.background = element_rect(fill = "gray60")) +
  theme_hidetitles() +
  labs(title = "Flavobacteriaceae") +
  theme(plot.title=element_text(hjust=0.5))

  
F
L = ggtern(data_rr_L, aes(linearChange,data_cr5_L$linearChange,data_cr24_L$linearChange)) +
  geom_point() +
  theme(panel.background = element_rect(fill = "orchid3")) +
  theme_hidetitles() +
  labs(title = "Litoricolaceae") +
  theme(plot.title=element_text(hjust=0.5))
L

N = ggtern(data_rr_N, aes(linearChange,data_cr5_N$linearChange,data_cr24_N$linearChange))+
  geom_point() +
  theme(panel.background = element_rect(fill = "orange2")) +
  theme_hidetitles() +
  labs(title = "Nitrincolaceae") +
  theme(plot.title=element_text(hjust=0.5))
N

R = ggtern(data_rr_R, aes(linearChange,data_cr5_R$linearChange,data_cr24_R$linearChange)) +
  geom_point() +
  theme(panel.background = element_rect(fill = "tomato2")) +
  theme_hidetitles() +
  labs(title = "Rhodobacteraceae") +
  theme(plot.title=element_text(hjust=0.5))
R

S = ggtern(data_rr_S, aes(linearChange,data_cr5_S$linearChange,data_cr24_S$linearChange))+
  geom_point() +
  theme(panel.background = element_rect(fill = "yellowgreen")) +
  theme_hidetitles() +
  labs(title = "Sphingomonadaceae") +
  theme(plot.title=element_text(hjust=0.5))

S

grid.arrange(A,F,L,N,R,S,ncol=2)


  
  
                    
