# Environmental bioinformatics final project - Diatom metabolism

Team: Max Jahns, Kate Lane, Katie Halloran

Paper: Shibl, A.A., et al. Diatom modulation of select bacteria through use of two unique secondary metabolites. Proc. Nat. Acad. Sci. 117, 27445-27455 (2020). https://doi.org/10.1073/pnas.2012088117

## Introduction 

For this project, we addressed a paper by Shibl and colleagues which examined how phytoplankton modulate and control their bacterial community using excreted metabolites. In this paper, they used metagenomics, (meta)transcriptomics, and metabolomics to investigate interactions between the marine diatom _Asterionellopsis glacialis_ and its bacterial consortium. 

To reproduce their work, we proposed to asemble bacterial genomes from metagenomes, analyze differential gene expression during coculture using the metatranscriptomies, and complement the sequence-based -omics results with reanalysis of metabolite mass spectrometry data. This included replicating the data analysis and visualization underlying Figures 1 and 2. Figure 3 illustrates biochemical processes; figure 4 analyzes microscopy data, which is not bioinformatically relevant; figures 5a-c use experimental data, and we chose to skip reanalysis of figure 5d in order to keep our project within a reasonable scope. 

All genomic data is available on NCBI, and totaled about 250Gb of data. Mass spectral datasets are available via the MassIVE database. However, rather than analyzing chromatographic data directly, we chose to use the SI-provided lists of measured compounds (i.e., molecules defined by a mass-to-charge ratio and chromatographic retention time, or mzRTs), which felt more within the bounds of bioinformatics; it seemed that reanalysis of chromatograms would instead tend more towards instrumental chemistry. The relevant lists of mzRTs were contained in two tables, totaling 801 KB.

## Repository description 

This repository contains the following directories: 
 - directory: description
 - directory: description
 - directory: description

## How to run our analyses

### Type of analysis 1

### Type of analysis 2

### Metabolomics analysis

The metabolomic data used in this project were relatively small and manageable, and required analysis in R. Analysis was performed in R version 4.1.2.

Data underwent some manual manipulation prior to being used in R. "Clean" data reflects this manual manipulation: 
 - Metadata information was collapsed into a single identifier to simplify analysis
 - Table 4 was only provided in .pdf format by Shibl et al., so relevant data were manually transferred to a usable file format. Further manual manipulation enabled some early-stage data exploration: 
   - In Table 4, m/z were manually calculated based on the listed mass (M (Da)) and adduct.
   - Retention time (RT) in seconds was manually calculated in Excel from RT in minutes, which was given. 
 - Metadata (not collapsed into a single identifier) in Dataset_S01 was preserved in a separate file

R packages required for this analysis include tidyverse, BiocManager, rstatix, devtools, ComplexHeatmap, and pcaMethods. If not already installed, these packages can be installed as follows: 

```
install.packages("tidyverse")
install.packages("BiocManager")
install.packages("rstatix")
install.packages("devtools")
install_github("jokergoo/ComplexHeatmap")
BiocManager::install("pcaMethods")
```

Alternately, if it's easier to run a re-analysis locally than on Poseidon: data for re-running metabolomics analyses can be found in data/clean_data/metabolome; .Rmd scripts can be found in scripts/metabolomics; and output figures can be found in output/metabolomics. In addition, analysis/figure_2* contains the .Rmd files used to generate the relevant figures described by the directory name. Scripts will require editing the data path to reflect the local environment, if rerunning analysis is more done locally. 

In this project, we tried to replicate two kinds of metabolomic analysis reported by Shibl et al. - PCA 
