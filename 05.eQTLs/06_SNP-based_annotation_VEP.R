library(MutationalPatterns)  # to load .vcf files as GRanges object

# Our data set consists of all human chromosomes including the variant consequences
# We create a vector containing the sample names:
sample.names <- c ( "homo_sapiens_incl_consequences-chr1", "homo_sapiens_incl_consequences-chr2", "homo_sapiens_incl_consequences-chr3", "homo_sapiens_incl_consequences-chr4", "homo_sapiens_incl_consequences-chr5", "homo_sapiens_incl_consequences-chr6", "homo_sapiens_incl_consequences-chr7", "homo_sapiens_incl_consequences-chr8", "homo_sapiens_incl_consequences-chr9", "homo_sapiens_incl_consequences-chr10", "homo_sapiens_incl_consequences-chr11", "homo_sapiens_incl_consequences-chr12", "homo_sapiens_incl_consequences-chr13", "homo_sapiens_incl_consequences-chr14", "homo_sapiens_incl_consequences-chr15", "homo_sapiens_incl_consequences-chr16", "homo_sapiens_incl_consequences-chr17", "homo_sapiens_incl_consequences-chr18", "homo_sapiens_incl_consequences-chr19", "homo_sapiens_incl_consequences-chr20", "homo_sapiens_incl_consequences-chr21", "homo_sapiens_incl_consequences-chr22", "homo_sapiens_incl_consequences-chrMT", "homo_sapiens_incl_consequences-chrX", "homo_sapiens_incl_consequences-chrY")

# We assemble a list of files we want to load, which match the sample names defined above
vcf.files <- list.files(path = "/storage/groups/ccm01/datasets/2019_Ensembl_GRCh37_variation",
                        pattern = ".vcf", full.names = TRUE)

# Get a reference genome BSgenome object.
library("BSgenome")
ref.genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref.genome, character.only = TRUE)

# This function loads the files as GRanges objects
vcfs <- read_vcfs_as_granges(vcf.files, sample.names, ref.genome)
