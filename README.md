Analysis of postmortem brain samples coming from healthy controls (CTL) and patients with schizophrenia (SCZ), bipolar disorder (BPD), and major depression disorder (MDD).

Following analyses will be performed: 

1. **Genome-Wide Association Study (GWAS)**
2. **Polygenic Risk Score (PRS) calculation**
  - Pre-processing of cohort, genotype data and GWAS summary statistics
  - Run LDpred
  - Run PRSice-2
  - Evaluate and plot results
3. **GRScore**
4. **Differential expression (DE) analysis**
  - Pre-processing of Affymetrix microarray data
  - QC and figures
  - Covariate evaluation and selection
  - DE analysis
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
