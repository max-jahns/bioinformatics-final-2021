# Environmental bioinformatics final project - Diatom metabolism

Team: Max Jahns, Kate Lane, Katie Halloran

Paper: Shibl, A.A., et al. Diatom modulation of select bacteria through use of two unique secondary metabolites. Proc. Nat. Acad. Sci. 117, 27445-27455 (2020). https://doi.org/10.1073/pnas.2012088117

## Introduction 

For this project, we addressed a paper by Shibl and colleagues which examined how phytoplankton modulate their bacterial community using excreted metabolites. In this paper, they used metagenomics, (meta)transcriptomics, and metabolomics to investigate interactions between the marine diatom _Asterionellopsis glacialis_ and its bacterial consortium. 

To reproduce their work, we proposed to asemble bacterial genomes from metagenomes, analyze differential gene expression during coculture using the metatranscriptomies, and complement the sequence-based -omics results with reanalysis of metabolite mass spectrometry data. This included replicating the data analysis and visualization underlying Figures 1 and 2. Figure 3 illustrates biochemical processes; figure 4 analyzes microscopy data, which is not bioinformatically relevant; figures 5a-c use experimental data, and we chose to skip reanalysis of figure 5d in order to keep our project within a reasonable scope. 

All genomic data is available on NCBI, and totaled about 250Gb of data. Mass spectral datasets are available via the MassIVE database. However, rather than analyzing chromatographic data directly, we chose to use the SI-provided lists of measured compounds (i.e., molecules defined by a mass-to-charge ratio and chromatographic retention time, or mzRTs), which felt more within the bounds of bioinformatics; it seemed that reanalysis of chromatograms would instead tend more towards instrumental chemistry. The relevant lists of mzRTs were contained in two tables, totaling 801 KB.

## Repository description 

This repository contains the following directories: 
 - data: raw data  organized into the analyses: diatom_genome, diatom_transcriptome, MAGS, metagenome, metatranscriptome, and metabolome
 - envs: yaml files that can be used to generate conda environments used in data analysis
 - jupyter-notebooks: Final jupyter notebooks for course and other notebooks for making figures with copies of files used as input for these notebooks. Note that our "final comparison" is split into three parts: 
      - "final-comparison-LANE.ipnyb," which reflects the final comparison for figure 1A
      - "final-comparison.ipnyb," which reflects the final comparison for figure 1B
      - "Metabolomic_analysis.ipynb," which reflects the final comparison for figure 2
 - logs: log files associated with organized in subdirectories based on analyses. Logs are numerous as jobs were often broken into parts in order to parallelize on poseidon HPC. Each directory contains a failed_logs directory of logs that did not run.
 - output: Directories for analyeses and for clean and processed data used in figure creation, organized into the analyses: diatom_genome, diatom_transcriptome,  Fig1A, Fig1B, Fig2, MAGS, metagenome, metatranscriptome, and metabolome
 - scripts: All slurm and R scripts used to analyze data and produce figures
 - tools: Accession number lists used for downloading data and directories for databases for tools used in this analysis

## How to run our analyses

### Metagenome, Metatranscriptome, and MAGS

All analyses are wrapped in slurm scripts, and numbers indicate the order in which they must be run. Downstream analyses' scripts used to generate /Fig1A are also present (and labeled in the order which they must be run), and supplemental files used in these analyses have also been copied there in directory /Fig1A_supps. These Fig1A_supps files are also located in the output/Fig1A directory. All wrappers include the activate on the conda environment to be used, and, if run in order, wrapper scripts will create the conda environments  needed from the /envs/ files. Possible specific modifications and instructions are included in each of the wrapper scripts.


### Diatom Transcriptome Analysis

Analysis of the diatom _Asterionellopsis glacialis_ transcriptome in reseeded and diatom_only cultures can be performed using the following pipeline:
 - _First_, all data aquistion, manipulation, and analysis of the transcriptome can be performed by running the 'diatom-transcriptome.qsub' script located in the 'scripts' directory
   - Please note, that running this script requires the diatom_transcriptome.yml conda environment to be active. You find the yml build in 'envs'.
   - The script is broken up into 'jobs' which are denoted by the :: '.job' :: header. Each job can be run seperately but MUST be run sequentially.
     - Each job lists its required inputs, including the path files where the job will look for the input.
     - Each job also lists the expected outputs. If running jobs seperately from the overall pipeline after the job has completed that the output can be found as listed, and that it is in the correct directory for the input of the next job.
   - The 'diatom-transcriptome.qsub' pipeline will produce all intermediate files, as if the jobs were run seperately.
   - The comments and outputted statement at the end of the 'diatom-transcriptome.qsub' will guide you through making sure the outputted data is ready for use in the diffrential expression step.
 - _Second_, after the pipeline has successfully completed, differential expression analysis can be performed using the DESeq2 package. The script for this process is provided in a jupyter notebook and can be found under 'jupyer-notebooks/differential-expression-analysis.ipnyb'.
     - To open this notebook, you can use the conda environment found in the 'r_jupyter.yml' file in the 'envs' directory.
     - This notebook requires four R packages to be installed which are listed at the top of the notebook. If these are not installed, uncomment the first chunk of the notebook to run instillation.
  - _Third_, to sumamrize all of this data and create figure1b, run the 'figure1b-creator-and-discussion.ipynb' jupyter notebook found in the 'jupyter-notebooks' directory.
     - To open this notebook, you can use the conda environment found in the 'r_jupyter.yml' file in the 'envs' directory.
     - This notebook requires three R packages to be installed which are listed at the top of the notebook. If these are not installed, uncomment the first chunk of the notebook t run instillation.

To summarize the pipeline is as follows... run 'scripts/diatom-transcriptome.qsub' -> open and run 'jupyer-notebooks/differential-expression-analysis.ipnyb' -> open and run 'figure1b-creator-and-discussion.ipynb'

Detailed analysis of the differences between our findings and the findings of Shibl et al. can be found in the 'figure1b-creator-and-discussion.ipynb' and 'final-comparision.ipnyb' jupyter notebooks.



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

Once packages are installed, metabolomic analysis can be replicated in the "Metabolomic analysis.ipynb" Jupyter notebook. 

Alternately, if it's easier to run a re-analysis locally than on Poseidon: data for re-running metabolomics analyses can be found in data/clean_data/metabolome; .Rmd scripts can be found in scripts/metabolomics; and output figures can be found in output/metabolomics. In addition, analysis/figure_2* contains the .Rmd files used to generate the relevant figures described by the directory name. Scripts will require editing the data path to reflect the local environment, if rerunning analysis is more done locally. 
