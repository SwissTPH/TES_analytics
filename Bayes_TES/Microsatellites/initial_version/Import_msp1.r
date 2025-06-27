


##### read in data
genotypedata_latefailures = read_excel(inputdata,sheet=1)



genotypedata_latefailures[genotypedata_latefailures == 0] = NA # missing data has to be coded as NA
genotypedata_latefailures[genotypedata_latefailures == "0"] = NA # missing data has to be coded as NA
genotypedata_latefailures[genotypedata_latefailures == "N/A"] = NA # missing data has to be coded as NA
genotypedata_latefailures[genotypedata_latefailures == "-"] = NA # missing data has to be coded as NA
genotypedata_latefailures[genotypedata_latefailures == "NA"] = NA # missing data has to be coded as NA

genotypedata_latefailures$Day[is.na(genotypedata_latefailures$Day)] = 0

sampleid = paste(genotypedata_latefailures$PatientID,"D",genotypedata_latefailures$Day,sep="")
genotypedata_latefailures = cbind(Sample.ID = sampleid,genotypedata_latefailures[,-c(1,2)])
### recode sample names so that each pair has a " Day 0" and a " Day Failure"
genotypedata_latefailures$Sample.ID = sub("D0$"," Day 0",genotypedata_latefailures$Sample.ID)
genotypedata_latefailures$Sample.ID = sub("D[0-9]+$"," Day Failure",genotypedata_latefailures$Sample.ID)


# each sample in genotypedata_RR has to have day 0 and day of Failure
ids = unique(unlist(strsplit(genotypedata_latefailures$Sample.ID[grepl("Day 0",genotypedata_latefailures$Sample.ID)]," Day 0")))
if (sum(!paste(ids, "Day Failure") %in% genotypedata_latefailures$Sample.ID) > 0) {
	print("Error - each sample must have day 0 and day of failure data")
}
ids = unique(unlist(strsplit(genotypedata_latefailures$Sample.ID[grepl("Day Failure",genotypedata_latefailures$Sample.ID)]," Day Failure")))
if (sum(!paste(ids, "Day 0") %in% genotypedata_latefailures$Sample.ID) > 0) {
	print("Error - each sample must have day 0 and day of failure data")
}

additional_genotypedata = read_excel(inputdata,sheet=2)
sampleid = paste(additional_genotypedata$PatientID,"D",additional_genotypedata$Day,sep="")
additional_genotypedata = cbind(Sample.ID = sampleid,additional_genotypedata[,-c(1,2)])
