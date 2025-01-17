# Disentangling the mutational effects on protein stability and interaction of human MLH1

## Introduction
This respository contains all data (except from the raw FASTQ files, which are available at the NCBI Gene Expression Omnibus (GEO) repository (accession number: GSE273652)) and code to repeat the processing and analysis of the abundance and interaction data of MLH1 in Larsen-Ledet et al.: "Disentangling the mutational effects on protein stability and interaction of human MLH1".

## Overview of files
*Output files*
* **MLH1_data.csv** - Abundance and interaction scores as well as standard errors and epsilon scores for MLH1 variants.
* **Tile_1-7_not_rescaled.csv** - Abundance (not rescaled) and interaction scores for MLH1 variants in each tile.
* **Replicate_correlation.csv** - Correlation between three replicates of abundance and interaction scores for MLH1 variants.
* * **Human_MutLalpha_full_AF3.cif** - The best-scoring model of the full human MutLα complex predicted by AlphaFold 3. 
  
*Input files*
* **Rosetta_Gemme_data.csv** - Rosetta (ddG) and GEMME (ddE) calculations as well as experimental abundance and interacation scores for MLH1 variants.
* **ClinVar_gnomAD_data.csv** - ClinVar classifications and gnomAD allele frequencies as well as experimental abundance and interaction scores for MLH1 variants.
* **rSASA_Secondary.csv** - Relative solvent accessible surface area (rSASA) and secondary structure of each residue as well as median experimental abundance and interaction scores for each position.
* **Benchmark_Hinrichsen.csv** - Benchmarks our experimental abundance scores against MLH1 expression levels as reported by Hinrichsen et al. (2013).
* **Benchmark_Abildgaard.csv** - Benchmarks our experimental abundance scores against MLH1 expression levels as reported by Abildgaard et al. (2019).
* **Benchmark_Kosinski.csv** - Benchmarks our experimental interaction scores against PMS2 dimerization as reported by Kosinski et al. (2010).

*Excel files*
* **Primers_anneal.temp..xlsx** - Primers and annealing temperatures for the first PCR in amplicon preparation.
* **SupplementalFile1.xlsx** - All data files combined in a single Excel file.

## Analysis of raw sequencing data
The raw sequencing files were analyzed using Enrich2 to calculate abundance and interaction scores for MLH1 variants. The scores were calculated based on three replicates and normalized to synonymous wild-type variants. 

## Data analysis and plotting
The MLH1_data_analysis.R file is used to produce all plots in the main figures, and the MLH1_data_analysis_supplementary.R file is used to produce all plots in the supplementary figures.

## Preprint
https://doi.org/10.1101/2024.07.28.605491 
