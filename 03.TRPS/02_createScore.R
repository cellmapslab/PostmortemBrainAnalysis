# #samples
# length(samples(header(vcf)) )
# #[1] 1556
# 
# #genomic postions
# head(rowRanges(vcf), 4)
# 
# #geno
# geno(vcf)
# 
# DS <-geno(vcf)$DS
# dim(DS)
# #[1]  136 1556
# #DS[1:3,1:4]
# 

##todo read summay stats filter maf und hwe! before!

# 
#read vcfs and get SNPs
#
#library("VariantAnnotation")
library(data.table)
eqtl=read.delim("01_eQTL_Score_Pos_withLD_EUR.txt", head=F)
#extraction postion based and check maf and hwe
snps_pos=c()
match=c()
for ( i in 1:22){
  stats=fread(paste0("geno_stats/chr",i, ".txt"), data.table = F, skip=9)
  stats$chrpos=paste0(stats$alternate_ids, ":", stats$position)
  stats.sub=stats[ stats$chrpos %in% eqtl$V1,]
  stats.sub.maf=stats.sub[ stats.sub$minor_allele_frequency>0.01,]
  stats.sub.maf.hwe = stats.sub.maf[ stats.sub.maf$HW_exact_p_value>=1e-6,]
  snps_pos= c(snps_pos,as.character(stats.sub.maf.hwe$chrpos))
  match= rbind(match,cbind(as.character(stats.sub.maf.hwe$chrpos), stats.sub.maf.hwe$rsid))
}
length(snps_pos)#6064
write.table(snps_pos, file="02_pmb_SNPs_Pos_withLD_EUR.txt", quote = F, col.names = F, row.names = F,sep="\t")
write.table(match, file="02_pmb_SNPs_Pos_withLD_EUR_match.txt", quote = F, col.names = F, row.names = F,sep="\t")

#GR data
load("../../../TRPS/eSNPs_hg19_296.rda")
snipa=read.delim("../../NCH/TRPS/SNIPA_proxySearch.results.EUR_v5.csv", head=T) #### use 2 snipa files
snipa2=read.delim("../../NCH/TRPS/SNIPA_proxySearch.results.EUR_v3.csv", head=T) #### use 2 snipa files

eqtlGR_hg19= as.data.frame(eqtlGR_hg19)
eqtlGR_hg19$SNP= row.names(eqtlGR_hg19)
tmp=paste(eqtlGR_hg19$seqnames, eqtlGR_hg19$start, sep= ":")
eqtlGR_hg19$pos_id= as.character(sapply(tmp,function(x)unlist(strsplit(as.character(x),"chr",""))[2]))
snipa$pos_id=paste(snipa$CHR,snipa$POS2, sep= ":")
snipa$v="v5"
snipa2$pos_id=paste(snipa2$CHR,snipa2$POS2, sep= ":")
snipa2$v="v3"
snipa=rbind(snipa, snipa2)

#delete proxy mit sich selbst
snipa= snipa[ as.character(snipa$QRSID) != as.character(snipa$RSID),]
snipa= snipa[ !duplicated(snipa$pos_id) ,] #5856

eqtlGR_hg19$pmb="no"
eqtlGR_hg19$pmb_SNP=""
eqtlGR_hg19$pmb_R2=""
eqtlGR_hg19$pmb_Proxy="no"
eqtlGR_hg19$pmb_Proxy_MAF=""
eqtlGR_hg19$pmb_Proxy_pos_id=""
eqtlGR_hg19$pmb_snipa=""
snps_pos=read.delim("02_pmb_SNPs_Pos_withLD_EUR.txt", head=F)

for(i in 1:dim(eqtlGR_hg19)[1]){
  #search for eSNP
  eSNP=grep(eqtlGR_hg19$pos_id[i],snps_pos$V1)
  if (length(eSNP)>0){
   # print(i)
    eqtlGR_hg19$pmb[i]="yes"
    eqtlGR_hg19$pmb_Proxy[i]="no"
  }
  
  #searh proxy
  if (length(eSNP)==0){
    print(i)
    proxys= snipa[ snipa$QRSID %in% eqtlGR_hg19$SNP[i],]
    proxys_pmb= proxys[ proxys$pos_id %in% snps_pos$V1,]
    proxy_snp= proxys_pmb[  which.max(proxys_pmb$R2),]
    if( dim(proxy_snp)[1]>0){
      eqtlGR_hg19$pmb[i]="yes"
      eqtlGR_hg19$pmb_SNP[i]=as.character(proxy_snp$RSID)
      eqtlGR_hg19$pmb_R2[i]=proxy_snp$R2
      eqtlGR_hg19$pmb_Proxy[i]="yes"
      eqtlGR_hg19$pmb_Proxy_MAF[i]=proxy_snp$MAF
      eqtlGR_hg19$pmb_Proxy_pos_id[i]=proxy_snp$pos_id 
      eqtlGR_hg19$pmb_snipa[i]=proxy_snp$v
      
    }
  }
  
}
write.table(eqtlGR_hg19, file="02_TRPS_pmb_overlap.txt", row.names = F,quote = F,sep="\t")

out= c(as.character(eqtlGR_hg19[ eqtlGR_hg19$pmb=="yes" & eqtlGR_hg19$pmb_Proxy=="yes",]$pmb_Proxy_pos_id),
as.character(eqtlGR_hg19[ eqtlGR_hg19$pmb=="yes" & eqtlGR_hg19$pmb_Proxy=="no",]$pos_id))
write.table(unique(out), file="02_TRPS_pmb_overlap_Pos.txt", row.names = F,quote = F,sep="\t", col.names = F)

length(unique(out)) #289
sum(eqtlGR_hg19$pmb=="no") #7 keine proxy 
sum(eqtlGR_hg19$pmb_Proxy=="yes") #2


##extract SNP and join/merge single data files from NHC to one combined file
#subset
#/storage/groups/ccm01/code/MPIP-Dean_postmortem_brain/
tools_path="/storage/groups/ccm01/tools/"
genotypes="/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/2018-12-20_post_mortem_brain_imputation/"
out_genotypes="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/genotype_data/4_TRPS/"

inclu=read.delim("02_TRPS_pmb_overlap.txt", head=T)
inclu = inclu[ inclu$pmb!="no",]
for(i in 1:22){

  inclu.chr=inclu[inclu$seqnames==paste0("chr",i),]
  inclu.chr.proxy= inclu.chr[inclu.chr$pmb_Proxy=="yes",]
  inclu.chr.noproxy= inclu.chr[inclu.chr$pmb_Proxy=="no",]
  
  if(dim(inclu.chr.proxy)[1]>0){
    write.table(c(paste0("NA:",inclu.chr.noproxy$start),paste0("NA:",unlist(strsplit(as.character(inclu.chr.proxy$pmb_Proxy_pos_id),":"))[2])) ,file="inlc_snps_chr.txt", row.names = F, col.names = F,quote = F)
  }else{
  write.table(paste0("NA:",inclu.chr$start),file="inlc_snps_chr.txt", row.names = F, col.names = F,quote = F)
  }
  system(paste0(tools_path,"qctool_v2.0.1/qctool -g ",genotypes,"pmb_info06_chr",i, ".gen.gz -s ", genotypes,"pmb_info06_chr",i, ".sample -incl-positions inlc_snps_chr.txt -og ",out_genotypes,"TRPS_chr",i,".gen -os ", out_genotypes,"TRPS_chr",i, ".sample"))

}

#merge 
command=c(paste0("cat "))
for (i in 1:22){
  tmp.chr = paste0(out_genotypes,"TRPS_chr",i,".gen ")
  command=paste0(command, as.character(tmp.chr))
}
system(paste0(command, "> ",out_genotypes,"TRPS_merged.gen"))
system(paste0("cp ",out_genotypes,"TRPS_chr1.sample ", out_genotypes,"TRPS_merged.sample"))


