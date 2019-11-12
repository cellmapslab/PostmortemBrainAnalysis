Analysis of postmortem brain samples coming from healthy controls (CTL) and patients with schizophrenia (SCZ), bipolar disorder (BPD), and major depression disorder (MDD).

Following analyses will be performed: 

1. **Genome-Wide Association Study (GWAS)**
  - Pre-processing of genotype data
  - GWAS 
  - Generate Manhattan plots
2. **Polygenic Risk Score (PRS) calculation**
  - Pre-processing of genotype data, cohort data and GWAS summary statistic files
  - Run LDpred
  - Run PRSice-2
  - Evaluate and plot results
3. **GRScore (TRPS)**
4. **Differential expression (DE) analysis**
  - Pre-processing of Affymetrix microarray data
  - QC and figures
  - Covariate evaluation and selection
  - DE analysis, evaluation and plotting of results
5. **eQTL analysis**
  - Pre-processing of genotype, gene expression and covariate data
  - Run MatrixEQTL
  - Evaluate and plot results
  - LD clumping
  - eQTL annotation and resampling
6. **Smoothing**
  - Generate trans-eQTLs
  - Filter, reweight and filter & reweight PRS using eQTLs or public data
7. **Co-expression analysis**
