### identify arms (based on site column)

numeric_id = as.numeric(gsub(" D.*","",genotypedata_latefailures$Sample.ID))
site_names = unique(genotypedata_latefailures$Site)

state_classification_all = c()
state_parameters_all = c()
ids_all = c()

for (site in site_names) {
	jobname = site
	genotypedata_RR = genotypedata_latefailures[genotypedata_latefailures$Site ==site,-c(2)]
	genotypedata_RR$Sample.ID = paste("_",genotypedata_RR$Sample.ID,sep="")
	additional_neutral =additional_genotypedata[additional_genotypedata$Site == site,-c(2)]
	if (dim(additional_neutral)[1] > 0) { additional_neutral$Sample.ID = paste("Additional_",1:dim(additional_neutral)[1],sep="")}
	source("mcmc_msp1.r")

	state_classification_all = rbind(state_classification_all,state_classification)
	state_parameters_all = rbind(state_parameters_all,state_parameters)
	ids_all = c(ids_all,ids)
}

rowMeans2 = function(x){
  if (length(dim(x)) == 0) {
    ret = mean(x)
  } else {
    ret = rowMeans(x)
  }
  ret
}

cbind2 = function(y,x){
  if (length(dim(x)) == 0) {
    ret = c(y,x)
  } else {
    ret = cbind(y,x)
  }
  ret
}

posterior_distribution_of_recrudescence = rbind(cbind2(ids_all,(state_classification_all)))
colnames(posterior_distribution_of_recrudescence)[1] = "ID"

probability_of_recrudescence = cbind(ids_all,rowMeans2(state_classification_all))

hist(rowMeans2(state_classification_all), breaks=10,main="Distribution of posterior probability of recrudescence", xlab="Posterior probability of recrudescence")

write.csv(probability_of_recrudescence, paste0(output_folder, "/", "probability_of_recrudescence_ALL.csv"))
