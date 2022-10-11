Analysis of postmortem brain samples coming from healthy controls (CTL) and patients with schizophrenia (SCZ), bipolar disorder (BPD), and major depression disorder (MDD).

Following analyses will be performed: 

1. **Genome-Wide Association Study (GWAS)**
  - Pre-processing of genotype data
  - GWAS 
  - Generate Manhattan plots
2. **Polygenic Risk Score (PRS) calculation**
  - Pre-processing of genotype data
  - Pre-processing of cohort data 
  - Pre-processing of GWAS summary statistics files
  - Generate input files (for PRSice-2 and LDpred)
  - Run PRSice-2
  - Run Lassosum
  - Run LDpred
  - Evaluate and plot results
3. **GRScore (TRPS)**
  - Extract SNPs
  - Create score
  - Calculate score
  - Associate score
4. **Differential expression (DE) analysis**
  - Pre-processing of Affymetrix microarray data
  - QC and covariate evaluation/selection
  - Differential expression analysis (DEA)
  - Gene set enrichment analysis using MGSA/MSigDB and figure generation
  - Generate input files for FUMA pathway enrichment
5. **eQTL analysis**
  - Pre-processing of genotype data
  - Generate input files (genotype, expression and covariate files)
  - Run MatrixEQTL
  - Generate result tables
  - Clump eSNPs using PLINK (LD clumping) and generate their result tables
  - Generate GRanges objects (and save as .Rds files) 
  - Generate FUMA input files
  - SNP-based annotation to ChromHMM, VEP, GTEx, CMC and meQTLs using resampling (10,000 times)
  - Evaluate and plot results
6. **Overlap to rare variants of index-patients**
  - Preprocessing and GRanges object generation
  - Generate FUMA input files
  - Generate overlap to eeGenes and eQTS genes and plot results
