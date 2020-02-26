########################################################
### Change genotype file to matrix EQTL input format ###
########################################################
# Invoke the script in the commandline as follows: R CMD BATCH --slave '--args geno=<path/genotype-file> sample=<path/sample-file> ids=<path/ids-file> output=<path/output-file>" 01c_change_probs.R
# Change genotype probabilites file to matrixEQTL input file format (only one instead of three probabilities)

# add installed packages from home path to the library search path
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

# load packages
library(data.table)

# First read the arguments listed in the commandline, trailingOnly uses only the parameters after --args
args = (commandArgs(trailingOnly = TRUE))	

# go through the list of arguments
if(length(args) == 0){
    	cat("\nScript to change the genotype probabilities (.gen.gz) file into the matrixEQTL input file format with only one probability value instead of three.\n\nYou can invoke the script as follows:\nR CMD BATCH --slave '--args geno=\"/home/.../genotype-file.gen.gz\" sample=\"/home/.../file.sample\" ids=\"/home/.../ids.txt\" output=\"/home/.../out-file.gen\"' 01c_change_probs.R\nCommandline parameters:\n\ngeno\t\tpath and filename of the merged genotype probabilites (ending in gen.gz)\nsample\t\tpath and filename of the sample file corresponding to the genotype file\nids\t\tpath and filename, which includes the ids of overlapping samples from the gene expression and the genotype dataset (here: 169)\noutput\t\tpath and filename of the output (input file for matrixEQTL)\n\n")
} else {
    for(i in 1:length(args)){
    	eval(parse(text = args[[i]]))
   	}

   	# 1. Load sample file and save sample IDs
   	sam <- read.delim(sample, head = T, sep = "")
   	all_ids <- as.character(sam$ID_1)[-1]  # number of all genotype samples
   	ninds <- length(all_ids)  # number of individuals in genotype data (169)

   	# 2. Load sample IDs (overlapping between gene expression and genotype data) and save them
   	cases <- read.delim(ids, head = F, sep = "")
   	cases_ids <- as.character(cases$V1)
   	ncases <- length(cases_ids)  # number of individuals overlapping between gene expression data and genotype data (169)

   	# 3. Load genotype data and save in data table
	mat <- fread(geno, header = F)
	nc <- ncol(mat)  # number of columns in genotype data (650)
	nr <- nrow(mat)  # number of snps in genotype file (xxx)

	# 4. SNP annotation/location
	snps_anno <- as.data.frame(mat[,c(2,1,3)])  # snp chr pos 
	colnames(snps_anno) <- c("snp", "chr", "pos")
	#snps_anno$snp <- paste(snps_anno$snp, mat$V4, mat$V5, sep = ":")

	# 5. Calculate one probability/dosage from the three probabilites per SNP/patient combination
	m <- t(2 * mat[, seq(6, nc, by = 3), with = F] + mat[, seq(6, nc, by = 3) + 1, with = F])
	row.names(m) <- all_ids
	colnames(m) <- as.character(snps_anno$snp)

	# 6. Get only overlapping samples/ids and transpose matrix
	m_case <- m[row.names(m) %in% cases_ids, ]
	m_trans <- t(m_case)
	# print(head(m_trans))

	# 7. Save data to file (in chunks)
	out <- paste(output, "/pmb_info06_all_chr_maf0.05_merged_noDup_matrixEQTL_transformed.gen", sep = "")
	file.create(out)
	chunkssize <- round(nr / 100)
	for (j in 0:99){
		if(j == 0){
			print(paste(j * chunkssize + 1, ":", (j + 1) * chunkssize))
	    	start <- j * chunkssize + 1
	    	end <- (j + 1) * chunkssize
	    	write.table(m_trans[c(start:end), ], out, append = T, row.names = T, col.names = T, sep = "\t", quote = F)
	    	# ...write with column names = T
		} else if(j == 99){
		    print(paste(j * chunkssize + 1, ":", (j + 1) * chunkssize))
		    start <- j * chunkssize + 1
		    end <- nr  # number of rows in genotype file
		    write.table(m_trans[c(start:end), ], out, append = T, row.names = T, col.names = F, sep = "\t", quote = F)
	  	} else {
			print(paste(j * chunkssize + 1, ":", (j + 1) * chunkssize))
			start <- j * chunkssize + 1
			end <- (j + 1) * chunkssize
			write.table(m_trans[c(start:end), ], out, append = T, row.names = T, col.names = F, sep = "\t", quote = F)
	  	}
	}
}

