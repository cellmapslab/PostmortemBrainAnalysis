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
  - Gene set enrichment analysis usind MGSA and MSigDB and figure generation
5. **eQTL analysis**
  - Pre-processing of genotype data
  - Generate input files (genotype, expression and covariate files)
  - Run MatrixEQTL
  - Clump eSNPs using PLINK (LD clumping)
  - Generate GRanges objects (and save as .Rds files)
  - eSNP based annotation and resampling 
  - Evaluate and plot results
6. **Smoothing**
  - Generate trans-eQTLs
  - Filter, reweight and filter & reweight PRS using eQTLs or public data
7. **Co-expression analysis**
