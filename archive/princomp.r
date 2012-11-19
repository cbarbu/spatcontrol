# principal component analysis in R on a dataset called data


#### usefull functions
## return the number of the column "namecol" in the dataset "dataset, namecol has to be entered with ""
numcol<-function(namecol,dataset){
	return(which(names(dataset)==namecol));
}

#### selection of the parameters in the dataset
# Please note: P.AD is not taken into account as it is a constant
params<-data[,c("FR_M","FR_A","No.RESID","I.num.amb","I.S","I.NO","I.LSE","I.AD","I.Gri","P.num.amb","P.S","P.NO","P.LSE","P.Gri","CU","CO","OV","PE","AV","GA","animal.techo","animal.patio","Residual.Rec","P.BLOQUETA","P.PIEDRA","P.panels")]

#### the principal components analysis:
params_pca<-prcomp(params,scale=TRUE)

summary(params_pca)

## let us be a little more selective
params<-data[,c("No.RESID","I.num.amb","I.S","I.NO","I.LSE","I.AD","P.num.amb","P.S","P.NO","P.LSE","CU","CO","OV","PE","AV","GA","P.BLOQUETA","P.PIEDRA")]
params_pca<-prcomp(params,scale=TRUE)
summary(params_pca)
# still at least 6 parameters

## on the data I had selected
params<-data[,c("I.Gri","P.num.amb","CU","PE")]
params_pca<-prcomp(params,scale=TRUE)
summary(params_pca)
# still 3 needed

## all the data selected by a lasso on everything
library(lasso2)
gl1_model<-gl1ce(pos~as.factor(FR_M)+FR_A+No.RESID+I.num.amb+I.S+I.NO+I.LSE+I.AD+I.Gri+I.panels+I.pir+P.num.amb+P.S+P.NO+P.LSE+P.panels+P.BLOQUETA+P.PIEDRA+P.Gri+CU+CO+OV+PE+AV+GA+animal.techo+animal.patio+Residual.Rec,data=data,family=binomial())
sort(coef(gl1_model)[abs(coef(gl1_model))>0.01])

# return all factors with a coef at least 10% of the higher coef. 
names(sort(coef(gl1_model)[abs(coef(gl1_model))>0.1*max(abs(coef(gl1_model)))]))

# that is (without the intercept):
 "Residual.Rec" "FR_A" "animal.patio" "P.num.amb" "P.Gri" "animal.techo" "as.factor(FR_M)9" "CU" "PE" "I.Gri" "as.factor(FR_M)10"

# or the same with "classical" risk factors
gl1_model<-gl1ce(pos~highseason+lowseason+No.RESID+I.num.amb+I.S+I.NO+I.LSE+I.AD+I.panels+I.pir+P.num.amb+P.S+P.NO+P.LSE+P.panels+P.BLOQUETA+P.PIEDRA+CU+CO+OV+PE+AV+GA+animal,data=data,family=binomial())
sort(coef(gl1_model))
names(sort(coef(gl1_model)[abs(coef(gl1_model))>0.0001*max(abs(coef(gl1_model)))]))
# only 13 parameters to start with
# then the pca
params<-data[,c("lowseason","highseason","I.NO","P.NO","P.BLOQUETA","No.RESID","OV","AV","P.num.amb","P.S","CO","CU","PE")]
params_pca<-prcomp(params,scale=TRUE)
# still  5 to 8 parameters to be kept (50-70%)

# or keep only high values (loose all materials)
coefs<-coef(gl1_model)[-1]
big_coefs<-names(sort(coefs[abs(coefs)>0.1*max(abs(coefs))]))
# theses are 11 params:"lowseason" "I.NO" "P.NO" "AV" "P.BLOQUETA" "P.num.amb" "P.S" "CO" "CU" "PE" "highseason"
params<-data[,big_coefs]
params_pca<-prcomp(params,scale=TRUE)
summary(params_pca)
# still 5 to 7 components to be kept, not that a big difference with keeping everything in a "readable format"

# oposant noble to everything else and CU,PE no animal to everything else, don't consider num amb, as supposedly already in No.RESID or animal presence
gl1_model<-gl1ce(pos~highseason+lowseason+No.RESID+I.NO+P.NO+CU+PE+noanimal,data=data,family=binomial())
# then everything is kept but Ro.RESID clearly weaker than the others so removed
# we stay with: highseason+lowseason+I.NO+P.NO+CU+PE+noanimal
coefs<-coef(gl1_model)[-1]
big_coefs<-names(sort(coefs[abs(coefs)>0.1*max(abs(coefs))]))
params<-data[,big_coefs]
params_pca<-prcomp(params,scale=TRUE)
summary(params_pca)
# then 3-4 is ok I can try to fit the 7 as is and look if I can find something to deal with this binary data

