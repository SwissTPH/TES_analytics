findposteriorfrequencies = function(x,tempdata) {
	data = tempdata[,c(1:maxMOI)+(x-1)*maxMOI];
	nalleles = frequencies_RR[[1]][x]
	freq_prior_alpha = rep(1,nalleles);
	freq_posterior_alpha = freq_prior_alpha + table(factor(c(data),levels=c(1:nalleles)));
	frequencies_RR[[2]][x,1:nalleles] <<- rdirichlet(1, freq_posterior_alpha);
}
