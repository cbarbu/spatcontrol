# this is the gibbs/MH hasting sampler, on both generated or true data:
# with cofactors or not
graphics.off()
INTERMEDIARY<-FALSE

source("parameters_sampler.r")
source("functions_intercept.r")
source("pseudo_data_generation.r")
source("DeltaSampling.r")
source("prep_sampling.r");
source("init_sampler.r");

source("kernel_sampler.r")

source("visualize_sampler.r")
save.image(file="EndVisu.img")

source("getDIC.r")
cat(file="in_brief.txt","\nDIC (partial):",DIC,"ptheta:",ptheta,"meanD:",meanD,"Dmean:",Dmean,"\n",append=TRUE);
source("prior_post_sup.r")
source("visu_imp_risks.r")
source("proba_to_pos.r")

save.image(file="BeforePredictCheck.img")
source("predict_quality.r")
save.image(file="EndPredictCheck.img")
