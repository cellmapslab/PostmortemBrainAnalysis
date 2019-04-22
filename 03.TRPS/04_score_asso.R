
# IDs, race (1= black) and sex. 
pheno <- read.csv("../2018-10-23_sample_sheet_KW.csv", sep=";")
dim(pheno) #[1] 209

score <- read.delim( "03_TRPS.txt", head = T) 
#combine geno and pheno 
score.pheno<- merge( score, pheno, by.x="IID", by.y="IID")
dim(score.pheno)#[1]175

#preliminary asso
fit.sex = lm( score.pheno$gr_score ~ score.pheno$Sex)
summary(fit.sex)
# Call:
#   lm(formula = score.pheno$gr_score ~ score.pheno$Sex)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0041759 -0.0012062  0.0003999  0.0008242  0.0024145 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      -2.006e-03  1.847e-04 -10.865   <2e-16 ***
#   score.pheno$SexM  9.565e-05  2.185e-04   0.438    0.662    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.001306 on 173 degrees of freedom
# Multiple R-squared:  0.001106,	Adjusted R-squared:  -0.004668 
# F-statistic: 0.1916 on 1 and 173 DF,  p-value: 0.6621

fit.race = lm( score.pheno$gr_score ~ score.pheno$Ancestry)
summary(fit.race)
# Call:
#   lm(formula = score.pheno$gr_score ~ score.pheno$Ancestry)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0042287 -0.0011958  0.0005005  0.0008991  0.0024088 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 -0.0021096  0.0002014 -10.473   <2e-16 ***
#   score.pheno$AncestryCEU      0.0002516  0.0002326   1.082    0.281    
# score.pheno$AncestryCEU/CHB  0.0005140  0.0013209   0.389    0.698    
# score.pheno$AncestryCHB     -0.0003668  0.0005697  -0.644    0.521    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.001305 on 171 degrees of freedom
# Multiple R-squared:  0.01321,	Adjusted R-squared:  -0.004106 
# F-statistic: 0.7628 on 3 and 171 DF,  p-value: 0.5164

fit.dx = lm( score.pheno$gr_score ~ score.pheno$Dx)
summary(fit.dx)
# 
# Call:
#   lm(formula = score.pheno$gr_score ~ score.pheno$Dx)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0041848 -0.0012126  0.0003910  0.0008396  0.0025223 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       -0.0022232  0.0003179  -6.993 5.79e-11 ***
#   score.pheno$DxCTL  0.0003214  0.0003576   0.899    0.370    
# score.pheno$DxMDD  0.0003319  0.0004155   0.799    0.426    
# score.pheno$DxSCZ  0.0003052  0.0003544   0.861    0.390    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.001311 on 171 degrees of freedom
# Multiple R-squared:  0.005237,	Adjusted R-squared:  -0.01221 
# F-statistic: 0.3001 on 3 and 171 DF,  p-value: 0.8253

score.pheno$Dx_new=0
score.pheno [score.pheno$Dx %in% c("MDD","SCZ","BPD"),]$Dx_new=1
fit.dx = lm( score.pheno$gr_score ~ score.pheno$Dx_new)
# summary(fit.dx)
# 
# Call:
#   lm(formula = score.pheno$gr_score ~ score.pheno$Dx_new)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0041848 -0.0012025  0.0003910  0.0008495  0.0024056 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        -1.902e-03  1.633e-04 -11.647   <2e-16 ***
#   score.pheno$Dx_new -5.715e-05  2.050e-04  -0.279    0.781    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.001306 on 173 degrees of freedom
# Multiple R-squared:  0.0004489,	Adjusted R-squared:  -0.005329 
# F-statistic: 0.0777 on 1 and 173 DF,  p-value: 0.7808
sub=score.pheno[ score.pheno$Dx %in% c("MDD","CTL"),]
fit.dx = lm( sub$gr_score ~ as.factor(sub$Dx))
summary(fit.dx)

sub=score.pheno[ score.pheno$Dx %in% c("SCZ","CTL"),]
fit.dx = lm( sub$gr_score ~ as.factor(sub$Dx))
summary(fit.dx)

sub=score.pheno[ score.pheno$Dx %in% c("BPD","CTL"),]
fit.dx = lm( sub$gr_score ~ as.factor(sub$Dx))
summary(fit.dx)

sub=score.pheno[ score.pheno$Dx %in% c("MDD","BPD"),]
fit.dx = lm( sub$gr_score ~ as.factor(sub$Dx))
summary(fit.dx)
