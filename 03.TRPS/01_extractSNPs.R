#use it on the 31.01.2019
### Änderung am 11.2018

##-incl-positions

#get GR-eQTL score SNPs
 score= read.delim("../../../TRPS/TRPS_n320_maxAbsFC_cor_withBetas.txt", head=T)
# id_mapping=read.delim("../../TRPS/rsid_update.txt", head=T)
# score$tag.SNP.eQTL.bin=as.character(score$tag.SNP.eQTL.bin)
# for(i in 1:length(id_mapping)){
#   j=grep(id_mapping$SNP[i], score$tag.SNP.eQTL.bin)
#   score$tag.SNP.eQTL.bin[j]= as.character(id_mapping$new_SNP_Id[i])
# }
# 
# #snipa: LD SNPs
# write.table(unique(score$tag.SNP.eQTL.bin), file="SNP_Ids.txt", row.names = F, col.names = F,quote = F,sep="\t")
# #http://snipa.helmholtz-muenchen.de/snipa3/index.php?task=proxy_search
# snipa=read.delim("SNIPA_proxySearch.results.EUR_v5.csv", head=T)
# snipa2=read.delim("SNIPA_proxySearch.results.EUR_v3.csv", head=T)

#generate list of SNPs to be extracted from genotype files
#1 based on pos.
#load hg19 postion and combine with LD
# load("../../../TRPS/eSNPs_hg19_296.rda")
# eqtlGR_hg19= as.data.frame(eqtlGR_hg19)
# pos=paste(eqtlGR_hg19$seqnames, eqtlGR_hg19$start, sep= ":")
# pos=as.character(sapply(pos,function(x)unlist(strsplit(as.character(x),"chr",""))[2]))
# 
# pos.ld=paste(snipa$CHR,snipa$POS2, sep= ":")
# pos.ld2=paste(snipa2$CHR,snipa2$POS2, sep= ":")
# length(unique(c(pos, pos.ld)))#6044
# length(unique(c(pos, pos.ld,pos.ld2)))#6153
 
##data=read.delim("01_eQTL_Score_Pos_withLD_EUR.txt",head=F)
##data_new=paste0("NA:",sapply(data$V1,function(x)unlist(strsplit(as.character(x),"\\:"))[2]))
##write.table(data_new, file="01_eQTL_Score_Pos_withLD_EUR_NA.txt",row.names = F,col.names = F,quote = F)


###working dir /storage/groups/ccm01/code/MPIP-Dean_postmortem_brain/

tools_path="/storage/groups/ccm01/tools/"
genotypes="/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/2018-12-20_post_mortem_brain_imputation/"
out_genotypes="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/genotype_data/4_TRPS/"

for((i=1; i<=22; i++)) do  ${tools_path}/qctool_v2.0.1/qctool -g ${genotypes}/pmb_info06_chr${i}.gen.gz -s ${genotypes}/pmb_info06_chr${i}.sample -snp-stats -osnp ${out_genotypes}/chr${i}.txt; done

#for((i=1; i<=22; i++)) do  ${tools_path}/qctool_v2.0.1/qctool -g ${genotypes}/pmb_info06_chr${i}.gen.gz -s ${genotypes}/pmb_info06_chr${i}.sample -incl-positions ${out_genotypes}/01_eQTL_Score_Pos_withLD_EUR_NA.txt -ofiletype vcf -og ${out_genotypes}/pmb_info06_chr${i}.vcf; done
