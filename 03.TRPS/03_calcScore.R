#04 Scoring 
score= read.delim("../../../TRPS/TRPS_n320_maxAbsFC_cor_withBetas.txt", head=T)
score_new= read.delim("02_TRPS_pmb_overlap.txt", head=T)

score_final=merge(score_new, score[,c("P.ID","tag.SNP.eQTL.bin","BETA", "risk_allele", "Q.value")], by.x="SNP", by.y="tag.SNP.eQTL.bin")

#duplicates --> no handle snps with 2 expression probes
#check of cor in Gex
dupl=score_final[ score_final$SNP %in% score_final[ duplicated(score_final$pos_id),]$SNP,]
load("../../../CisAnalysis/ESet_afterBatch.rda")
Eset.dupl= ESetNorm[ featureNames(ESetNorm)%in% dupl$P.ID,]
dupl$excl=0 #0== exclude, 1== take and 2== take and calc separate

#if expression probes are correlated with r2>06. then just keep one
#else keep both
snps=unique(dupl$SNP)
for(i in 1: length(snps)){
  tmp=dupl[ dupl$SNP %in% snps[i],]
  Eset.tmp=ESetNorm[ featureNames(ESetNorm)%in% tmp$P.ID,]
  if( min(cor(t(exprs(Eset.tmp)))) >0.6){
    min.q=tmp[which.min(abs(tmp$Q.value)),]
    dupl[ dupl$pos_id %in% as.character(min.q$pos_id) & dupl$P.ID %in% as.character(min.q$P.ID),]$excl=1
    
  }else{ dupl[ dupl$pos_id %in% tmp$pos_id,]$excl=2}
}

score_final_1 = score_final[ !(score_final$SNP %in% dupl[ dupl$excl %in% c(0,2),]$SNP & score_final$P.ID %in% dupl[ dupl$excl %in% c(0,2),]$P.ID) ,]
score_final_1_out.proxy= score_final_1[ score_final_1$pmb=="yes" & score_final_1$pmb_Proxy=="yes",c("pmb_Proxy_pos_id", "risk_allele","BETA")]
colnames(score_final_1_out.proxy)[1] = "pos_id"
score_final_1_out.noproxy=score_final_1[ score_final_1$pmb=="yes" & score_final_1$pmb_Proxy=="no",c("pos_id", "risk_allele","BETA")]
score_final_1_out = rbind(score_final_1_out.noproxy, score_final_1_out.proxy)
#dim(score_final_1_out)#283   3

score_final_2 = score_final[ score_final$SNP %in% dupl[ dupl$excl %in% 2,]$SNP & score_final$P.ID %in% dupl[ dupl$excl %in% 2,]$P.ID ,]#10
#no proxys! so we can just take it
score_final_1_out=rbind(score_final_1_out, score_final_2[c(1,3,5,7,9),c("pos_id", "risk_allele","BETA")])
#dim(score_final_1_out)#288   3

score_final_2_out = score_final_2[c(2,4,6,8,10),c("pos_id", "risk_allele","BETA")]
#dim(score_final_2_out)#5   3

write.table(score_final_1_out, "03_TRPS_score_1.txt", col.names = F,row.names = F,quote = F,sep="\t")
write.table(score_final_2_out, "03_TRPS_score_2.txt", col.names = F,row.names = F,quote = F,sep="\t")



#calc score for the 2 sets as on eSNPs can have 2 betas == eQTLs
#remove NA in the first column
#awk '{$1="";print}' geno_stats/TRPS_merged.gen > geno_stats/TRPS_merged_noNA.gen
#cut -d' ' -f3 geno_stats/TRPS_merged_noNA.gen > geno_stats/TRPS_merged_noNA_ids
#chnage ids
#ids=read.delim("geno_stats/TRPS_merged_noNA_ids", head=F)
match=read.delim("02_pmb_SNPs_Pos_withLD_EUR_match.txt", head=F)
#length(intersect(ids$V1, match$V2))#289

score1= read.delim("03_TRPS_score_1.txt", head=F)
score1_match=merge(score1, match, by="V1")
write.table(score1_match[,c(4,2,3)], file="03_TRPS_score_1_match.txt", row.names = F,col.names = F,quote = F,sep="\t")
score2= read.delim("03_TRPS_score_2.txt", head=F)
score2_match=merge(score2, match, by="V1")
write.table(score2_match[,c(4,2,3)], file="03_TRPS_score_2_match.txt", row.names = F,col.names = F,quote = F,sep="\t")

system(paste0("/Applications/plink_mac/plink2 --gen geno_stats/TRPS_merged_noNA.gen --sample geno_stats/TRPS_merged.sample --out 03_testscore_1 --score 03_TRPS_score_1_match.txt"))
system(paste0("/Applications/plink_mac/plink2 --gen geno_stats/TRPS_merged_noNA.gen --sample geno_stats/TRPS_merged.sample --out 03_testscore_2 --score 03_TRPS_score_2_match.txt"))

##allele mis-match: n=7
nopred= read.delim("03_testscore_1.nopred", head=F, sep="")
score_final_1_out_match= merge(score_final_1_out, match, by.x="pos_id", by.y="V1")
for(i in 1:dim(nopred)[1]){
  risk_allele = score_final_1_out_match[ score_final_1_out_match$V2 %in% nopred$V2[i],]$risk_allele
  if(risk_allele=="T"){risk_allele="A"
  }else if(risk_allele=="A"){risk_allele="T"
  }else if(risk_allele=="C"){risk_allele="G"
  }else{risk_allele="C"}
  score_final_1_out_match[ score_final_1_out_match$V2 %in% nopred$V2[i],]$risk_allele =risk_allele
}
write.table(score_final_1_out_match[,c(4,2,3)], "03_TRPS_score_1_updateA.txt", col.names = F,row.names = F,quote = F,sep="\t")
system(paste0("/Applications/plink_mac/plink2 --gen geno_stats/TRPS_merged_noNA.gen --sample geno_stats/TRPS_merged.sample --out 03_testscore_1_update --score 03_TRPS_score_1_updateA.txt"))

### mean over all betas in score
score1=read.delim("03_testscore_1_update.profile", head=T, sep="")
score2=read.delim("03_testscore_2.profile", head=T, sep="")

score=merge(score1, score2, by=c("FID", "IID", "PHENO"))
colnames(score)[4:6]=c( "CNT_1" ,  "CNT2_1",  "SCORE_1")
colnames(score)[7:9]=c( "CNT_2" ,  "CNT2_2" , "SCORE_2")


score$gr_score= (score$SCORE_1 + score$SCORE_2) / 2
write.table(score[,-c(4:9)], "03_TRPS.txt",row.names = F,quote = F,sep="\t")


