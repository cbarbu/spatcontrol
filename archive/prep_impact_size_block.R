### get block size and check if possible cofactor
source("RanalysisFunctions.R")

#############
## Make data
#############
# base data
db.base <- read.csv("DB_simple_Pau_cyclo1_19Jul2011.csv",header=TRUE);
db.base$oanimal<-with(db.base,as.numeric(CO>0 | OV>0 | AV>0 |GA>0 | otros.animales != -1))

# number of house per block
housesPerBlock <- aggregate(cbind(rep(1,length(db.base$block_num)),db.base[,c("insp","pos","CU","PE","oanimal","I.NO","P.NO")]),by=list(db.base$block_num),sum)
names(housesPerBlock)[1:2]<-c("block_num","count")

# add number of house per block for each house
db.base2<-merge(db.base,housesPerBlock[,c("count","block_num")],by.x="block_num",by.y="block_num")

# plot of the numbers
plot_reel(db.base2$longitude,db.base2$latitude,db.base2$count,base=1,top=max(db.base2$count))
#=> looks well defined by block and consistent with size of blocks


#############
## Analysis 
#############

## univariate: positive rate versus size of block
# direct link positive rate by block and size of the block
housesPerBlock$ratePos<-housesPerBlock$pos/housesPerBlock$insp
plot(log(housesPerBlock$count),housesPerBlock$ratePos)
with(housesPerBlock,plot(count,ratePos))
#=> looks more going downward: the bigger the more difficult to have very high infestation rate

# glmanalysis:
## The BEST: zero-inflated poisson regression
library(pscl)
zi<-zeroinfl(pos~log(count+1)+offset(log(insp+1)),data=housesPerBlock)
# Call:
# zeroinfl(formula = housesPerBlock$pos ~ log(housesPerBlock$count + 1) + 
#     offset(log(housesPerBlock$insp + 1)))
# 
# Pearson residuals:
#     Min      1Q  Median      3Q     Max 
# -1.4854 -0.8126 -0.4682  0.4207  4.4899 
# 
# Count model coefficients (poisson with log link):
#                               Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                   -1.14306    0.16645  -6.867 6.55e-12 ***
# log(housesPerBlock$count + 1) -0.07672    0.05254  -1.460    0.144    
# 
# Zero-inflation model coefficients (binomial with logit link):
#                               Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                     1.4027     0.4264    3.29    0.001 ** 
# log(housesPerBlock$count + 1)  -1.5532     0.1472  -10.55   <2e-16 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
AIC(zi)
# AIC: 3216.168
#=> the size of the block very significantly decrease the probability of having no house
#   infested in the block, but once block infested doesn't change much the probability
#   for each house
housesPerBlock$predictResp<-predict(zi,type="response")
housesPerBlock$predictCount<-predict(zi,type="count")
housesPerBlock$predictZero<-predict(zi,type="zero")
par(mfrow=c(1,3))
with(housesPerBlock,plot(pos,predictResp))
abline(a=0,b=1)
with(housesPerBlock,plot(pos,predictZero))
abline(a=0,b=1)
with(housesPerBlock,plot(pos,predictCount))
abline(a=0,b=1)
dev.new()
with(housesPerBlock,plot(insp,pos))

## normal poisson regression
summary(glm(housesPerBlock$pos~log(housesPerBlock$count+1)+offset(log(housesPerBlock$insp+1))),family=poisson)
# Call:
# glm(formula = housesPerBlock$pos ~ log(housesPerBlock$count + 
#     1) + offset(log(housesPerBlock$insp + 1)))
# 
# Deviance Residuals: 
#     Min       1Q   Median       3Q      Max  
# -4.0652  -2.0471  -0.8336   1.0912  21.4688  
# 
# Coefficients:
#                               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                    -1.7683     0.3775  -4.685 3.26e-06 ***
# log(housesPerBlock$count + 1)   0.6309     0.1390   4.539 6.47e-06 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
# 
# (Dispersion parameter for gaussian family taken to be 10.09123)
# 
#     Null deviance: 8866.1  on 859  degrees of freedom
# Residual deviance: 8658.3  on 858  degrees of freedom
# AIC: 4432.6
#=> the size of the block increase the infestation

### more complete
housesPerBlock$rateCU<-housesPerBlock$CU/housesPerBlock$insp
housesPerBlock$ratePE<-housesPerBlock$PE/housesPerBlock$insp
housesPerBlock$rateI.NO<-housesPerBlock$I.NO/housesPerBlock$insp
housesPerBlock$rateP.NO<-housesPerBlock$P.NO/housesPerBlock$insp
## glm on rates
summary(glm(housesPerBlock$pos~log(housesPerBlock$count+1)+housesPerBlock$rateCU+housesPerBlock$ratePE+housesPerBlock$oanimal+housesPerBlock$rateI.NO+housesPerBlock$rateP.NO+offset(log(housesPerBlock$insp+1))),family=poisson)
# glm(formula = housesPerBlock$pos ~ log(housesPerBlock$count + 
#     1) + housesPerBlock$rateCU + housesPerBlock$ratePE + housesPerBlock$oanimal + 
#     housesPerBlock$rateI.NO + housesPerBlock$rateP.NO + offset(log(housesPerBlock$insp + 
#     1)))
# 
# Deviance Residuals: 
#     Min       1Q   Median       3Q      Max  
# -6.4801  -2.1340  -0.8863   1.3216  12.1050  
# 
# Coefficients:
#                                Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                   -1.895505   1.478094  -1.282   0.2009    
# log(housesPerBlock$count + 1)  0.008789   0.494023   0.018   0.9858    
# housesPerBlock$rateCU         -1.683466   1.515885  -1.111   0.2679    
# housesPerBlock$ratePE         -0.797161   1.243355  -0.641   0.5220    
# housesPerBlock$oanimal         0.386143   0.087085   4.434 1.41e-05 ***
# housesPerBlock$rateI.NO        0.652897   1.207614   0.541   0.5892    
# housesPerBlock$rateP.NO        2.096645   1.256520   1.669   0.0965 .  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
# 
# (Dispersion parameter for gaussian family taken to be 10.81201)
# 
#     Null deviance: 3062.2  on 247  degrees of freedom
# Residual deviance: 2605.7  on 241  degrees of freedom
#   (15 observations deleted due to missingness)
# AIC: 1303.1

# Number of Fisher Scoring iterations: 2

## glm on log(counts)
 summary(glm(housesPerBlock$pos~log(housesPerBlock$CU+1)+log(housesPerBlock$PE+1)+log(housesPerBlock$oanimal+1)+log(housesPerBlock$count+1)+log(housesPerBlock$I.NO+1)+log(housesPerBlock$P.NO+1)+offset(log(housesPerBlock$insp+1))),family=poisson)
# glm(formula = housesPerBlock$pos ~ log(housesPerBlock$CU + 1) + 
#     log(housesPerBlock$PE + 1) + log(housesPerBlock$oanimal + 
#     1) + log(housesPerBlock$count + 1) + log(housesPerBlock$I.NO + 
#     1) + log(housesPerBlock$P.NO + 1) + offset(log(housesPerBlock$insp + 
#     1)))
# 
# Deviance Residuals: 
#     Min       1Q   Median       3Q      Max  
# -5.1169  -2.2085  -0.8867   1.4167  12.0704  
# 
# Coefficients:
#                                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                     -0.25267    0.92181  -0.274 0.784224    
# log(housesPerBlock$CU + 1)      -0.01724    0.43660  -0.039 0.968529    
# log(housesPerBlock$PE + 1)      -0.04101    0.94750  -0.043 0.965512    
# log(housesPerBlock$oanimal + 1)  1.11094    0.54020   2.057 0.040743 *  
# log(housesPerBlock$count + 1)   -1.67927    0.89280  -1.881 0.061118 .  
# log(housesPerBlock$I.NO + 1)     0.73079    0.93662   0.780 0.435967    
# log(housesPerBlock$P.NO + 1)     1.40652    0.42005   3.349 0.000935 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
# 
# (Dispersion parameter for gaussian family taken to be 10.50181)
# 
#     Null deviance: 3062.7  on 262  degrees of freedom
# Residual deviance: 2688.5  on 256  degrees of freedom
# AIC: 1373.7

# Number of Fisher Scoring iterations: 2

### glm using bigblock
summary(glm(housesPerBlock$pos~log(housesPerBlock$CU+1)+log(housesPerBlock$PE+1)+log(housesPerBlock$oanimal+1)+(housesPerBlock$count>16)+log(housesPerBlock$I.NO+1)+log(housesPerBlock$P.NO+1)+offset(log(housesPerBlock$insp+1))),family=poisson)

 
##################
## Household level
##################

# basic glm on house level with the presence absence
model<-glm(pos~CU+PE+I.NO+P.NO+log(count),data=db.base2[which(db.base2$insp==1),],family=binomial)
# Call:
# glm(formula = pos ~ CU + PE + I.NO + P.NO + log(count), family = binomial, 
#     data = db.base2[which(db.base2$insp == 1), ])
# 
# Deviance Residuals: 
#     Min       1Q   Median       3Q      Max  
# -0.9768  -0.6796  -0.5868  -0.4524   2.2553  
# 
# Coefficients:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -2.22548    0.15906 -13.991  < 2e-16 ***
# CU           0.55532    0.05978   9.290  < 2e-16 ***
# PE           0.62300    0.05938  10.492  < 2e-16 ***
# I.NO        -0.21397    0.08004  -2.673  0.00751 ** 
# P.NO        -0.29344    0.06104  -4.807 1.53e-06 ***
# log(count)   0.15152    0.04953   3.059  0.00222 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
#     Null deviance: 9412.2  on 9848  degrees of freedom
# Residual deviance: 9154.0  on 9843  degrees of freedom
# AIC: 9166
# 
# Number of Fisher Scoring iterations: 4

# basic glm with the bigBlock
db.base2$bigBlock<-as.numeric(db.base2$count > median(housesPerBlock$count))
plot(db.base2$easting,db.base2$northing,col=(db.base2$count>21)+1,asp=1)
model<-glm(pos~CU+PE+I.NO+P.NO+oanimal+bigBlock,data=db.base2[which(db.base2$insp==1),],family=binomial)
# Call:
# glm(formula = pos ~ CU + PE + I.NO + P.NO + bigBlock, family = binomial, 
#     data = db.base2[which(db.base2$insp == 1), ])
# 
# Deviance Residuals: 
#     Min       1Q   Median       3Q      Max  
# -0.9725  -0.6970  -0.5930  -0.4555   2.2700  
# 
# Coefficients:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -1.96688    0.08979 -21.906  < 2e-16 ***
# CU           0.55277    0.05980   9.243  < 2e-16 ***
# PE           0.62699    0.05940  10.555  < 2e-16 ***
# I.NO        -0.23505    0.08000  -2.938   0.0033 ** 
# P.NO        -0.29553    0.06098  -4.846 1.26e-06 ***
# bigBlock     0.28389    0.06372   4.455 8.39e-06 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
#     Null deviance: 9412.2  on 9848  degrees of freedom
# Residual deviance: 9143.0  on 9843  degrees of freedom
# AIC: 9155
# 
# Number of Fisher Scoring iterations: 4


##########################
######## TRASH below
##########################
# probability of the block to be infected with the size
getGroup<-function(value,orderedLimits){
	classNum <- which(value<orderedLimits)[1]-1
	return(classNum)
}


housesPerBlock$binInsp<-as.numeric(housesPerBlock$insp>1)
housesPerBlock$binPos<-as.numeric(housesPerBlock$pos>1)
sizeGroupsBreaks<-c(1,2,4,8,16,32,64,128)
housesPerBlock$sizeGroup<-sapply(housesPerBlock$count,getGroup,sizeGroupsBreaks)
sel<-which(housesPerBlock$binInsp==1)

binPosSizeBlock<-t(table(housesPerBlock$binPos[sel],housesPerBlock$sizeGroup[sel]))
binPosSizeBlockRatePos<-binPosSizeBlock[,2]/(binPosSizeBlock[,1]+binPosSizeBlock[,2])
cbind(binPosSizeBlock,binPosSizeBlockRatePos)
#     0   1 binPosSizeBlockRatePos
# 2  13   2              0.1333333
# 3  65  16              0.1975309
# 4 173  81              0.3188976
# 5 181 167              0.4798851
# 6  26  39              0.6000000
# 7   0   3              1.0000000

chisq.test(housesPerBlock$binPos[sel],housesPerBlock$sizeGroup[sel])
#=> It's a bit stupid but ok the more houses in the block the more likely you have at least one 
#   infested house in the block

write.csv(db.base2,"DB_simple_Pau_cyclo1_19Jul2011_blockSize.csv");

