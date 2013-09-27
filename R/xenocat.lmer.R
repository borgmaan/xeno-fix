####
#
#	LMER
#	Functions for linear mixed-effects models
#
####

#
#	EM-algorithm for categorizing mixed-effects models (lmer)
#
xeno.EM = function(
	formula, 
	data, 
	max.iter=100, 
	loglh.difference=0.01, 
	responsename="Response", 
	treatmentname="Treatment", 
	tpname="Timepoint",  
	idname="Tumor_id", 
	discriminate=TRUE, 
	randomstart=FALSE, 
	return.iterations=FALSE,
	id=NULL,
	verbose=TRUE,
	...
){
	require(lme4)
	check_headers(data, name=c(responsename, treatmentname, tpname, idname), mer=FALSE, dataframe=TRUE)

	if(is.null(id) & verbose){
		type="print"
	}else{
		type="write"
	}

	growthname="Growth"
	if(discriminate & randomstart & verbose){
		print("Warning! Discriminate and Random start should not be used simultaneously, as this can lead to non-feasible local maxima of categorization.")
		print("Suggested combinations are Discriminate with No random start or No discrimination with Random start.")
	}

	# Predicting tumor volumes for given growth category using fixed effects of the formula
	predict.tumor = function(fitted_formula, growth){

		modelmatrix = model.matrix(fitted_formula)
		
		# Testing which fixed effects hold Growth-labels
		if(length(grep("Growth", names(fixef(fitted_formula))))>0 &(loc = grep("Growth", names(fixef(fitted_formula))))[1]!=0){
			modelmatrix[,loc] = modelmatrix[,loc]*growth
		}
		
		result = as.vector(modelmatrix %*% fixef(fitted_formula))
		result

	}

	# Expectation step estimation of theta	
	expected.theta = function(fit, data){
		pred.1 = predict.tumor(fitted_formula=fit, growth=1)
		pred.0 = predict.tumor(fitted_formula=fit, growth=0)

		## equal priors on both classes cancel in Bayes Theorem

		# If there are missing values, cannot compare to original response although we can do the prediction

		exp_theta = vector(length=length(unique(data[[idname]])))
		
		for(i in 1:length(unique(data[[idname]]))){
			included_response = (!is.na(data[[responsename]])) & unique(data[[idname]])[i]==data[[idname]]
			included_pred = (unique(data[[idname]])[i]==data[[idname]][!is.na(data[[responsename]])])
			exp_theta[i] = prod(dnorm(data[[responsename]][included_response], pred.1[included_pred], sd=sd(data[[responsename]][!is.na(data[[responsename]])]))) / (prod(dnorm(data[[responsename]][included_response], pred.0[included_pred], sd=sd(data[[responsename]][!is.na(data[[responsename]])]))) + prod(dnorm(data[[responsename]][included_response], pred.1[included_pred], sd=sd(data[[responsename]][!is.na(data[[responsename]])]))))
		}


		# Replicating theta to be same for each of the measurement values within a specific tumor
		exp_theta = rep(exp_theta, times=tps_per_tumor)
		exp_theta

	}

	tps_per_tumor = vector(length=0)
	uniques = unique(data[[idname]])
	
	for(i in 1:length(unique(data[[idname]]))){
		sub = subset(data, data[[idname]]==uniques[i])
		tps_per_tumor = append(tps_per_tumor, length(unique(sub[[tpname]])))
	}

	# Whether we're using random start or prior that all tumors are growing
	if(randomstart){
		if(verbose){print("Randomizing theta")}
		theta = rbinom(length(unique(data[[idname]])), 1, 0.5)
		# Replicating randomly starting theta to be same within a specific tumor
		theta = rep(theta, times=tps_per_tumor)
		# Cannot allow all theta to be == 0 as this would result in a downdated formula matrix
		while(sum(theta)==0 | sum(theta)==length(data[[tpname]])){
			theta = rbinom(length(unique(data[[idname]])), 1, 0.5)
			#theta = rep(theta, each=length(unique(data[[tpname]])))
			theta = rep(theta, times=tps_per_tumor)
		}
		updated_data = data.frame(data, Growth=theta)	
	}else{
		if(verbose){print("Set first step theta to 1")}
		updated_data = data.frame(data, Growth=1)	
	}
	
	fit = lmer(formula, data=updated_data, ...)
	
	likelihood = logLik(fit)
	
	
	iterations = list()
	datas = list()
	iter=1
	old.likelihood = likelihood - 1
	if(verbose){print("Starting iterations...")}
	while ( abs(likelihood-old.likelihood)> loglh.difference & iter < max.iter){
		#if(iter%%50==0){ print(paste("Iter:",iter, "of", max.iter, "/ Current logLik:",likelihood))} 
		if(!is.null(id) | (verbose & iter%%25==0)){
			save_progress(
				paste("EM iteration ",iter," of maximum ",max.iter, " (current log-lik: ",likelihood,")", sep=""),
				id=id,
				type=type
				)
		}

		iterations[[iter]] = fit
		datas[[iter]] = updated_data
		
		## Expectation step
		exp_theta = expected.theta(fit=fit, data=updated_data)
		# Discriminating expected theta to be between 0 and 1
		if(discriminate){exp_theta = (exp_theta>0.5)*1}

		updated_data["Growth"] = exp_theta

		## Maximization step, where we assume that exp_theta is theta
		fit = lmer(formula, data=updated_data, ...)

		old.likelihood = likelihood
		likelihood = logLik(fit)
		
		iter = iter + 1
	}
	
	# Testing the estimated random effects
	for(i in 1:length(ranef(fit))){
		rans = ranef(fit)[[i]]
		dims = dim(ranef(fit)[[i]])
		zeros = rep(0, times=dims[1])
		for(j in 1:dims[2]){
			if(identical(rans[,j], zeros)){
				stop(paste("Model random effects estimation for component '", names(rans)[j],
					"' resulted in false convergence, please consider adjusting the model formulation."))
			}
		}
	}
	
	# Saving/printing progress
	if(!is.null(id) | verbose){
		save_progress(
			paste("EM-algorithm finished. Total amount of iterations:",iter),
			id=id,
			type=type
			)
	}
	
	# Returning solution(s)
	if(return.iterations){
		list(fit, updated_data, iterations)
	}else{
		updated_data
	}
}

#
#	Diagnostics function for testing global optimum and coefficients of the model
#
xeno.diag = function(
	x,
	model = Response ~ 1 + Treatment + Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0 + Timepoint|Tumor_id),
	starts=1000,
	discriminate=TRUE,
	idname="Tumor_id",
	tpname="Timepoint",
	treatmentname="Treatment",
	responsename="Response",
	verbose = FALSE
){
	growthname="Growth"
	
	if(class(x)=="data.frame"){
		dat = x	
	}else if(class(x)=="mer"){
		dat = x@frame
	}else{
		stop("Parameter x should be either a data.frame or a lme4 mer-object.")
	}
	
	
	if(growthname %in% names(dat)){
		dat = dat[,-which(growthname==names(dat))]
	}
	

	datas = list()
	fits = list()
	logs = vector()
	
	datas[[1]] = xeno.EM(model, dat, discriminate=discriminate, randomstart=FALSE, responsename=responsename, tpname=tpname, idname=idname, treatmentname=treatmentname, verbose=verbose)	
	fits[[1]] = lmer(model, datas[[1]])
	logs[1] = logLik(fits[[1]])	
	
	for(i in 1:starts){
		datas[[i+1]] = try(xeno.EM(model, dat, discriminate=discriminate, randomstart=TRUE, responsename=responsename, tpname=tpname, idname=idname, treatmentname=treatmentname, verbose=verbose), silent=TRUE)
		while(class(datas[[i+1]])=="try-error"){
			datas[[i+1]] = try(xeno.EM(model, dat, discriminate=discriminate, randomstart=TRUE, responsename=responsename, tpname=tpname, idname=idname, treatmentname=treatmentname, verbose=verbose), silent=TRUE)		
		}
		fits[[i+1]] = lmer(model, datas[[i+1]])
		logs = append(logs, logLik(fits[[i+1]]))
	}
	
	logliks = as.vector(unique(logs))
	count = 0
	for(i in 1:(length(logliks)-1)){
		if(logliks[i+1]>logliks[1]){
			count = count + 1
		}
	}
	# Global optimum
	if(count==0){
		print("Global optimum for categories was identified using a prior start.");
	}else{
		print("Warning!")
		print(paste("Global optimum for categories identified with random start, amount of random start solutions better than a prior solution:", count));
	}
	
	# Growth inhibition |b4/b3|
	if(paste(treatmentname,":",tpname,":",growthname, sep="") %in% names(fixef(fits[[1]]))
		&& paste(tpname,":",growthname, sep="") %in% names(fixef(fits[[1]]))){
		est_b3 = fixef(fits[[1]])[which(names(fixef(fits[[1]]))==paste(tpname,":",growthname, sep=""))]
		est_b4 = fixef(fits[[1]])[which(names(fixef(fits[[1]]))==paste(treatmentname,":",tpname,":",growthname, sep=""))]
		print(paste("Growth inhibition ratio (|b4/b3|):", round(abs(est_b4/est_b3),3)*100, "%"))
		if(abs(est_b4/est_b3)>1){
			print("Warning! Ratio |b4/b3| above 100%, possible violation of model assumptions.")
		}
	}else if(paste(treatmentname,":",tpname, sep="") %in% names(fixef(fits[[1]]))
		&& tpname %in% names(fixef(fits[[1]]))){
		est_b3 = fixef(fits[[1]])[which(names(fixef(fits[[1]]))==tpname)]
		est_b4 = fixef(fits[[1]])[which(names(fixef(fits[[1]]))==paste(treatmentname,":",tpname, sep=""))]
		print(paste("Growth inhibition ratio (|b4/b3|):", round(abs(est_b4/est_b3),3)*100, "%"))
		if(abs(est_b4/est_b3)>1){
			print("Warning! Ratio |b4/b3| above 100%, possible violation of model assumptions.")
		}
	}else{
		print("Could not find b3 and b4 terms from model formulation.")
	}
	
	# Statistical significance of model terms
	if(class(x)=="mer"){
		fixs = fixef(x)
		if(length(fixs)==4){
			names(fixs) = paste("b", 1:length(fixs), sep="")
		}else if(length(fixs)==2){
			print("Two fixed effects detected, assuming a model formulation without a target size")
			names(fixs) = paste("b", 3:4, sep=="")
		}else{
			print("Custom amount of fixed effects detected!")
		}
		sign = xeno.fixef.pvals(x)
		print("Statistically significant model fixed effects estimates (according to default MCMC):")
		print(fixs[sign[,6]<0.05])
		if("(Intercept)" %in% names(fixef(x)[sign[,6]<0.05])){
			print("b1 detected as statistically significant - confirm coherence with the tumor inclusion criteria.")
		}
		print("Statistically insignificant:")
		print(fixs[sign[,6]>=0.05])
	}
}


#
#	Function for collecting fixed effect statistical significance through MCMC etc.
#
xeno.fixef.pvals = function(
	fit, 
	LRT=FALSE, 
	MCMC=TRUE, 
	MCMCsim=10000, 
	digits=3,
	id=NULL,
	verbose=FALSE)
{
	require(lme4)
	if(!(class(fit) == "mer" || class(fit) == "lmerMod")){
		stop("Input object should be a mer-object")
	}
	
	# Namespace issues with lme4 in R 2.14
	#summarycoef = try(summary(fit)@coefs, silent=TRUE)
	summarycoef = try(extract_str(fit=fit, type="coefs"), silent=TRUE)
	# Problem in newer versions of R (apparently since 2.13) with the @coefs
	if(class(summarycoef)=="try-error"){
		print("Error extracting @coefs from the mer-object, this may be an issue with the newer versions of R.")
		print("To extract the standard deviations and t values of the object please use e.g. R version 2.11")
		summarycoef = matrix(nrow=length(fixef(fit)), ncol=3)
		summarycoef[,1] = fixef(fit)
		colnames(summarycoef) = c("Estimate", "Std. Error", "t value")
		rownames(summarycoef) = names(fixef(fit))
	}

	if(MCMC){
		if(is.null(id)){
			MCMCpvals = xeno.fixef.MCMC(fit=fit, MCMCsim=MCMCsim, digits=digits, draw=FALSE, verbose=verbose)
		}else{
			MCMCpvals = xeno.fixef.MCMC(fit=fit, MCMCsim=MCMCsim, digits=digits, draw=FALSE, id=id, verbose=verbose)
		}
		summarycoef = cbind(summarycoef, MCMCpvals)
	}
	if(LRT){
		LRTpvals = round(xeno.fixef.LRT(fit), digits)
		names(LRTpvals) = names(fixef(fit))
		summarycoef = cbind(summarycoef, LRTpvals)
	}
	
	summarycoef
}

#
#	Fixed effects statistical significance by MCMC-simulation
#
xeno.fixef.MCMC = function(
	fit,
	MCMCsim=10000,
	digits = 6,
	draw=TRUE,
	burn.draw = FALSE,
	burn.cut = 0,
	id=NULL,
	verbose=FALSE,
	time=FALSE)
{

	draw.burn.func = function(
		samps){
		
		fixsper2 = round((dim(samps@fixef)[1])/2,0)		

		par(mfrow=c(fixsper2, 2))
		for(i in 1:dim(samps@fixef)[1]){
			plot(samps@fixef[i,], type="l", xlab="Sample index", ylab=rownames(samps@fixef)[i])		
		}
		
	}

	require(lme4)
	if(!(class(fit) == "mer" || class(fit) == "lmerMod")){
		stop("Input object should be a mer-object")
	}
	if(is.null(id) & verbose){
		type="print"
	}else{
		type="write"
	}
	
	if(!is.null(id) | verbose){
		test = system.time(mcmcsamp(fit, 1000))[3]
		real = (MCMCsim/1000)*test
		mins = round(real/60,0)
		secs = round(real%%60,0)
		time_est = paste(mins,"m ",secs,"s",sep="")
		if(!time){
			save_progress(
				paste("Generating MCMC-samples (Estimated time: ",time_est,")",sep=""),
				id=id,
				type=type
			)
		}else{
			return(paste("Estimated time for generating the MCMC samples: ",time_est, sep=""))		
		}
	}
	if(!time){
		samples = mcmcsamp(fit, MCMCsim)
		# Visualizing MCMC samples according to their index for detecting any possible shift or convergence
		if(burn.draw){
			draw.burn.func(samples)	
		}
		# Possibly dropping out initial samples
		if(!burn.cut==0 && burn.cut>0){
			samples@fixef = samples@fixef[,burn.cut:length(samples@fixef[1,])]	
		}
		if(!is.null(id) | verbose){
			save_progress(
				"Computing HPD-intervals and p-values",
				id=id,
				type=type
			)		
		}
		intervals = HPDinterval(samples)
		ncols = ncol(samples@fixef)
		abovezero = rowSums(1*(samples@fixef > 0))/ncols
		MCMCpvals = round(2*pmax(0.5/ncols, pmin(abovezero, 1-abovezero)), digits)
		names(MCMCpvals) = names(fixef(fit))
		bounds = intervals$fixef[,1:2]
		colnames(bounds) = c("HPD95%lower", "HPD95%higher")
		if(draw){
			xeno.draw.HPDfixef(fit=fit, MCMCsim=MCMCsim, samples=samples)
		}
		if(!is.null(id) | verbose){
			save_progress(
				"MCMC-simulations finished",
				id=id,
				type=type
			)
		}
		cbind(bounds, MCMCpvals)
	}
}

#
#	Visualizing for the Highest Posterior Density, MCMC-samples, mean of MCMC and fixed effect estimates
#
xeno.draw.HPDfixef = function(
	fit, 
	MCMCsim=10000,
	samples=NULL
){
	require(lattice)
	require(lme4)
	if(!(class(fit) == "mer" || class(fit) == "lmerMod")){
		stop("Input object should be a mer-object")
	}
	cols=c("red","green","orange")
	ltys=c(2,3,4)
	
	if(is.null(samples)){
		samples = mcmcsamp(fit, MCMCsim)
	}else if(!class(samples)=="merMCMC"){
		stop("samples-parameter has to be a merMCMC-object generated with the mcmcsamp-function")
	}
	intervals = HPDinterval(samples)$fixef
	fixefs = length(fixef(fit))
	sampleplots = data.frame(sim=as.vector(t(samples@fixef)), fix=rep(names(fixef(fit)), each=MCMCsim))
	
	print(densityplot( ~ sim | fix,
		data = sampleplots,
		xlab="Value",
		layout=c(round(fixefs/2,0),2),
		pch=".",
		subscripts=TRUE,
		key = list(
			text = list(c("HPD 95%","MCMC mean","Estimate")),
			lines = list(
				lty = c(ltys[1],ltys[2],ltys[3]),
				col = c(cols[1],cols[2],cols[3]),
				lwd = c(2,2,2)
			)
		),
		panel = function(x,subscripts){
			index = round(subscripts[1]/MCMCsim,0)+1
			panel.densityplot(x, pch=".")
			panel.abline(v=intervals[index,1], type="l", lty=ltys[1], lwd=2, col=cols[1])
			panel.abline(v=intervals[index,2], type="l", lty=ltys[1], lwd=2, col=cols[1])
			panel.abline(v=mean(x), type="l", lty=ltys[2], lwd=2, col=cols[2])
			panel.abline(v=fixef(fit)[index], type="l", lty=ltys[3], lwd=2, col=cols[3])
			index=index+1
		}
	))
}

#
#	Function for extracting fixed effect significance by fitting the nested model with ML with omitted term(s) and computing significance according to the LRT
#
xeno.fixef.LRT = function(
	fit,
	digits=6) 
{
	if(!(class(fit) == "mer" || class(fit) == "lmerMod")){
		stop("Input object should be a mer-object")
	}
	model = as.formula(fit@call)
	data.lmer = data.frame(model.matrix(fit))
	
	names(data.lmer) = names(fixef(fit))
	names(data.lmer) = gsub(pattern=":", x=names(data.lmer), replacement=".", fixed=TRUE)
	names(data.lmer) = ifelse(names(data.lmer)=="(Intercept)", "Intercept", names(data.lmer))
	string.call = strsplit(x=as.character(model), split=" + (", fixed=TRUE)
	var.dep = unlist(strsplit(x=unlist(string.call)[2], " ~ ", fixed=TRUE))[1]
	vars.fixef = names(data.lmer)
	
	formula.ranef = paste("+ (", string.call[[3]][-1], sep="")
	formula.ranef = paste(formula.ranef, collapse=" ")
	formula.full = as.formula(paste(var.dep, "~ -1 +", paste(vars.fixef, collapse=" + "), 
	formula.ranef))
	
	data.ranef = data.frame(fit@frame[,which(names(fit@frame) %in% names(ranef(fit)))])
	names(data.ranef) = names(ranef(fit))
	data.lmer = data.frame(fit@frame[, 1], data.lmer, data.ranef)
	names(data.lmer)[1] = var.dep
	for(i in 1:length(names(fit@frame))){
		if(!(names(fit@frame)[i] %in% data.lmer)){
			data.lmer[names(fit@frame)[i]] = fit@frame[[names(fit@frame)[i]]]
		}
	}
	out.full = lmer(formula.full, data=data.lmer, REML=FALSE)
	p.value.LRT = vector(length=length(vars.fixef))
	for(i in 1:length(vars.fixef)) {
		formula.reduced = as.formula(paste(var.dep, "~ -1 +", paste(vars.fixef[-i], 
		collapse=" + "), formula.ranef))
		out.reduced = lmer(formula.reduced, data=data.lmer, REML=FALSE)
		out.LRT = data.frame(anova(out.full, out.reduced))
		p.value.LRT[i] = round(out.LRT[2, 7], digits)
	}

	p.value.LRT
}

#
#	Function for permutating the tumor labels for assessing statistical significance of treatment fixed effects
#
xeno.fixef.perm = function(
	fit,
	perms = 10000,
	idname="Tumor_id",
	tpname="Timepoint",
	treatmentname="Treatment",
	responsename="Response",
	id=NULL,
	verbose=FALSE,
	time=FALSE
){
	require(lme4)
	if(!(class(fit) == "mer" || class(fit) == "lmerMod")){
		stop("Input object should be a mer-object")
	}
	check_headers(fit, name=c(responsename, treatmentname, tpname, idname), mer=TRUE, dataframe=FALSE)
	
	if(is.null(id) & verbose){
		type="print"
	}else{
		type="write"
	}
	
	growthname="Growth"
	form = as.formula(fit@call)
	
	# Original data frame
	orig_frame = fit@frame
	uniques = as.vector(unique(orig_frame[[idname]]))
	# Original t staticics
	#orig_t = try(summary(asS4(fit))@coefs[,3], silent=TRUE)
	orig_t = try(extract_str(fit=fit, type="coefs")[,3], silent=TRUE)
	if(class(orig_t)=="try-error"){
		stop("Error extracting coefficient estimates, std. dev. and t values from the mer-object summary. Try reverting to R version 2.12 or older for this functionality. (@coefs for mer-object)")
	}
	
	if(!is.null(id) | verbose){
		test = system.time(xeno.fixef.perm(fit=fit, perms=10, idname=idname, tpname=tpname,treatmentname=treatmentname,responsename=responsename))[3]
		real = (perms/10)*test
		mins = round(real/60,0)
		secs = round(real%%60,0)
		time_est = paste(mins,"m ",secs,"s", sep="")
		if(!time){
			save_progress(
				paste("Permutating 0% (Estimated time: ",time_est,")",sep=""),
				id=id,
				type=type
			)
		}else{
			return(paste("Estimated time for permutations: ",time_est, sep=""))
			#save_progress(
			#	paste("Estimated time for permutations: ",time_est, sep=""),
			#	id=id,
			#	type=type
			#)
		}
	}
	if(!time){
		perm_t = matrix(nrow=perms, ncol=length(fixef(fit)))
		# Permutating treatment-labels while conserving growth-category if using categorizing approach
		for(i in 1:perms){
			if((!is.null(id) | verbose )& i%%25==0){
				save_progress(
					paste("Permutating ",round(i/perms*100,0),"% (Estimated time: ",time_est,")",sep=""),
					id=id,
					type=type				
				)	
			}
			perm_frame = data.frame()
			samp = sample(length(uniques), replace=FALSE)

			for(j in 1:length(uniques)){
				subdata = subset(orig_frame, orig_frame[idname]==uniques[samp[j]])
				ids = paste("Perm_",j)
				treat = subdata[[treatmentname]][1]
				if(growthname %in% names(orig_frame)){
					grow = subdata[[growthname]][1]
				}
				realdata = subset(orig_frame, orig_frame[idname]==uniques[j])
				resp = realdata[[responsename]]
				tps = realdata[[tpname]]
				if(growthname %in% names(orig_frame)){			
					perm = data.frame(resp, treat, grow, tps, ids)
					names(perm) = c(responsename, treatmentname, growthname, tpname, idname)
				}else{
					perm = data.frame(resp, treat, tps, ids)
					names(perm) = c(responsename, treatmentname, tpname, idname)
				}

				perm_frame = rbind(perm_frame, perm)
			}
			fit_perm = lmer(form, data=perm_frame)
			#perm_t[i,] = try(summary(asS4(fit_perm))@coefs[,3], silent=TRUE)
			perm_t[i,] = try(extract_str(fit=fit_perm, type="coefs")[,3], silent=TRUE)		
			if(class(perm_t)=="try-error"){
				stop("Error extracting coefficients estimates, std. dev. and t values from the mer-object summary. Try reverting to R version 2.12 or older for this functionality. (@coefs for mer-object)")
			}
		}
		list(orig_t, perm_t)
		pvals = vector(length=0)
		for(i in 1:length(orig_t)){
			pvals = append(pvals, sum(abs(orig_t[i])<abs(perm_t[,i]))/perms)
		}
		names(pvals) = names(fixef(fit))
		if(!is.null(id) | verbose){
			save_progress(
				"Permutations finished",
				id=id,
				type=type
			)			
			#prog = try(write.table("Permutations finished",file=paste("progress_",id,".txt",sep=""), row.names=FALSE, col.names=FALSE), silent=TRUE)
			#if(class(prog)=="try-error"){
			#	print("Error writing current progress to a text file...")
			#}
		}
		pvals
	}
}

#
#	Computation for within-tumor autocorrelation of residuals
#
xeno.test.autocor = function(
	fit, 
	idname="Tumor_id", 
	tpname="Timepoint", 
	draw=FALSE, 
	maintitle="Autocorrelation"
){
	require(lme4)
	check_headers(fit, name=c(tpname, idname), mer=TRUE, dataframe=FALSE)

	# Checking residual autocorrelation so that only residuals within the same tumor will be considered for the autocorrelation
	resids = list()
	minlag = 0
	tumors = unique(fit@frame[[idname]])
	maxlag = 0

	for(i in 1:length(tumors)){
		if(length(fit@resid[which(fit@frame[[idname]]==tumors[i])])-1>maxlag){
			maxlag = length(fit@resid[which(fit@frame[[idname]]==tumors[i])])-1
		}
		
		
		if(length(which(fit@frame[[idname]]==tumors[i])>0)){	
			resids[[i]] = fit@resid[which(fit@frame[[idname]]==tumors[i])]
		}else{
			resids[[i]] = NULL
		}
	}
	
	lags_a = list()
	lags_b = list()
	for(i in minlag:maxlag){
		temp_a = vector(length=0)
		temp_b = vector(length=0)
		for(j in 1:length(tumors)){
			resids_a = resids[[j]]
			resids_b = resids[[j]]
			if(i<length(resids_a)){
				if(i>0){
					resids_a = resids_a[-c((length(resids_a)-i+(1-minlag)):length(resids_a))]
					resids_b = resids_b[-c(1:i)]
				}
				temp_a = append(temp_a, resids_a)
				temp_b = append(temp_b, resids_b)
			}
		}
		lags_a[[i+1]] = temp_a
		lags_b[[i+1]] = temp_b
	}

	
	autocors = vector(length=(maxlag+1))
	for(i in 1:length(autocors)){
		autocors[i] = cor(lags_a[[i]], lags_b[[i]])
	}
	
	if(draw){
		plot(0,0,col="white",xlim=c(minlag,maxlag), ylim=c(min(autocors, na.rm=TRUE),max(autocors, na.rm=TRUE)), xlab="Lag", ylab="Autocorrelation of residuals (within ID label)", main=maintitle)
		for(i in 1:length(autocors)){
			lines(c(i-1, i-1), c(0,autocors[i]))
		}
	}
	labs = paste("Lag_",0:maxlag, sep="")
	names(autocors) = labs
	autocors
}

#
# 	Power for chosen N and T according to the simulated datasets
#
xeno.sim.power = function(
	fit, 
	file="POW_ANALYSIS.txt", 
	nsim=10000, 
	Ns = c(10,15,20), 
	pvalMCMC=TRUE, 
	MCMCsim=10000, 
	alpha=0.05,
	responsename="Response",
	treatmentname="Treatment",
	tpname="Timepoint",
	idname="Tumor_id",
	T=max(fit@frame[[tpname]])
){
	require(lme4)
	check_headers(fit, name=c(responsename, treatmentname, tpname, idname), mer=TRUE, dataframe=FALSE)
	
	growthname="Growth"	
	formula = as.formula(fit@call)
	pow = matrix(nrow=length(Ns),ncol=length(names(fixef(fit))))
	if(growthname %in% names(fit@frame)){
		categories = TRUE
		print("Categories found from fit frame, using categorized power analysis...")
		if(!length(unique(fit@frame[[growthname]]))==2){
			stop("The current power analysis procedure only supports binary categorization")
		}
	}else{
		categories = FALSE
		print("Categories not found from fit frame, using non-categorized power analysis...")
	}
	
	colnames(pow)=names(fixef(fit))
	rownames(pow)=Ns
	for(nindex in 1:length(Ns)){
		print(paste("Simulating for N:",Ns[nindex]))
		powers = matrix(nrow=nsim,ncol=length(names(fixef(fit))))
		for(i in 1:nsim){
			simulated_data = xeno.sim.tumor(fit=fit, N=Ns[nindex], T=T, responsename=responsename, treatmentname=treatmentname, tpname=tpname, idname=idname)

			if(i%%10==0){print(paste("Time:",Sys.time(), " - i",i,"of",nsim))}
			sim_fit = try(lmer(data=simulated_data, formula=formula), silent=TRUE)
			# Simulation resulted in unfeasible dataset:
			while(class(sim_fit)=="try-error"){
				simulated_data = xeno.sim.tumor(fit=fit, N=Ns[nindex], T=T, responsename=responsename, treatmentname=treatmentname, tpname=tpname, idname=idname)
				sim_fit = try(lmer(data=simulated_data, formula=formula), silent=TRUE)
			}			

			if(pvalMCMC){
				# Computing significance according to the MCMC simulations		
				pval = xeno.fixef.pvals(
						fit=sim_fit,
						LRT=FALSE,
						MCMC=TRUE,
						MCMCsim=MCMCsim,
						digits=4)[,6]
			}else{
				# Computing significance according to the LRT testing
				pval = xeno.fixef.LRT(sim_fit)
			}
			powers[i,] = 1*(pval<alpha)
			for(j in 1:length(names(fixef(fit)))){
				pow[nindex,j] = mean(powers[,j])
			}
		}
		print(pow)
		write.table(pow, file=file)
	}
	pow
}

#
#	Data simulation according to Gelman & Hill
#
xeno.sim.tumor = function(
	fit, 
	responsename="Response", 
	treatmentname="Treatment", 
	tpname="Timepoint",  
	idname="Tumor_id" ,
	N=10, 
	T=max(fit@frame[[tpname]])
){
	require(lme4)
	check_headers(fit, name=c(responsename, treatmentname, tpname, idname), mer=TRUE, dataframe=FALSE)	
	
	growthname="Growth"
	# Information extracted from original data
	mintp = min(fit@frame[[tpname]])
	
	# Categories
	if(growthname %in% names(fit@frame)){
		if(!length(unique(fit@frame[[growthname]]))==2){
			stop("The current power analysis procedure for categorized profiles supports only binary categorization")
		}
		orig_treat_0_tp0 = subset(fit@frame, fit@frame[[treatmentname]]==0 & fit@frame[[tpname]]==mintp)
		orig_treat_1_tp0 = subset(fit@frame, fit@frame[[treatmentname]]==1 & fit@frame[[tpname]]==mintp)
		orig_treat_0_growthprop = mean(orig_treat_0_tp0[[growthname]])
		orig_treat_1_growthprop = mean(orig_treat_1_tp0[[growthname]])
	# No categories, all tumor curves treated equal
	}else{
		orig_treat_0_growthprop = 1
		orig_treat_1_growthprop = 1
	}
	# Equal amount of tumors for each treatment group
	N=N*2

	treat_0_growths = rbinom(N/2, 1 ,orig_treat_0_growthprop)
	treat_1_growths = rbinom(N/2, 1 ,orig_treat_1_growthprop)
	
	timepoints = rep(seq(mintp,T),N)
	treatments = c(rep(0,N/2, each=(T+1-mintp)), rep(1,N/2,each=(T+1-mintp)))
	id = paste("ID_",rep(1:N, each=(T+1-mintp)),sep="")
	growths = rep(c(treat_0_growths, treat_1_growths), each=(T+1-mintp))
	sim_data = data.frame(treatments, timepoints,id,growths)
	names(sim_data)=c(treatmentname,tpname,idname,growthname)

	# variances for the fixed effects
	fix = fixef(fit)
	fixnames = names(fixef(fit))
	
	# Random effects matrix
	#reff_mat = try(as.matrix(summary(asS4(fit))@REmat), silent=TRUE)
	reff_mat = try(as.matrix(extract_str(fit=fit, type="REmat")), silent=TRUE)
	# In newest versions of R some extractor functions do not work with lme4
	if(class(reff_mat)=="try-error"){
		stop("Error extracting random effects matrix from the mer-object summary. Try reverting to R version 2.12 or older for this functionality. (@REmat for mer-object)")
	}
	
	rannames = names(ranef(fit))
	
	fixes = list()

	
	#Intercept
	intercept_fixpart = 0
	if("(Intercept)" %in% names(fix)){
		intercept_fixpart = fix["(Intercept)"]
		fix = fix[-which(names(fix)=="(Intercept)")]	
	}else{
		intercept_fixpart = 0
	}

	if("(Intercept)" %in% reff_mat[,"Name"]){
		fixes[[1]] = rep(rnorm(N, intercept_fixpart, sd=as.numeric(reff_mat[which(reff_mat[,"Name"]=="(Intercept)"),"Std.Dev."])), each=(T+1-mintp))	
	}else{
		fixes[[1]] = rep(intercept_fixpart, each = N, times = (T + 1 - mintp))
	}
	
	## ASSUMING ONLY TUMOR ID IS GIVEN AS A GROUPING FACTOR
	#Other fixed effects that are affected by random effects
	for(i in 2:(length(reff_mat[,"Name"])-1)){
		name = reff_mat[i,"Name"]
		namesplit = strsplit(name, ":")[[1]]
		product=0
		if(length(which(namesplit==tpname))>0){
			fixes[[i]] = rep(rnorm(N, product, sd=as.numeric(reff_mat[which(reff_mat[,"Name"]==name),"Std.Dev."])), each=(T+1-mintp), times=1)*sim_data[[tpname]]
		}else{
			fixes[[i]] = rep(rnorm(N, product, sd=as.numeric(reff_mat[which(reff_mat[,"Name"]==name),"Std.Dev."])),each=(T+1-mintp))
		}
	}
	
	#Fixed effects not affected by random effects
	for(i in 1:length(fix)){
		name = names(fix)[i]
		namesplit = strsplit(name, ":")[[1]]
		product = 1
		for(j in 1:length(namesplit)){
			if(namesplit[j] %in% names(sim_data)){
				product = product * sim_data[,which(namesplit[j]==names(sim_data))]
			}
		}
		# Adding the coefficient
		fixes[[length(fixes)+1]] = product * fix[name]
	}
	
	result = fixes[[1]]
	for(i in 2:length(fixes)){
		result = result + fixes[[i]]
	}
	
	response = rnorm(N*(T+1-mintp), result, sd=as.numeric(reff_mat[which(reff_mat[,"Groups"]=="Residual"),"Std.Dev."]))
	sim_data = data.frame(response, sim_data)
	names(sim_data)[1] = responsename
	sim_data
}

#
#	Function for testing a component extracted from a fit against a data.frame column
#
xeno.test.marker.fit = function(
	fit,
	orig_data, 
	componentname="Growth",
	other="", 
	idname="Tumor_id",
	tpname="Timepoint",
	valueat=1,
	verbose=TRUE,
	rm.na=TRUE
){
	require(lme4)
	check_headers(fit, name=c(tpname, idname), mer=TRUE, dataframe=FALSE)
	
	if(!other %in% names(orig_data)){
		stop(paste("Could not find column from the original data:",other))
	}else if(!componentname %in% names(fit@frame)){
		stop(paste("Could not find column from the fit frame:",componentname))
	}
	
	fit_data = fit@frame
	
	tumors = as.vector(unique(fit_data[[idname]]))
	components = unique(fit_data[[componentname]])
	components = components[order(components)]
	
	if(verbose) cat("\nComponents found:",paste(components, sep=" , "),"\n")
	
	markers = orig_data[[other]]
	non.na.markers = markers[!is.na(markers)]
	if(verbose) cat("\nData values found (amount):",length(non.na.markers),"\n")
	
	component_values = list()
	for(i in 1:length(components)){
		component_values[[i]] = vector(length=0)
	}
	
	
	# Looping over the different component values
	for(i in 1:length(tumors)){
		tumor_fit = subset(fit_data, fit_data[idname]==tumors[i])
		id = which(tumor_fit[[componentname]][1] == components)
		tumor_orig = subset(orig_data, orig_data[idname]==tumors[i])
		other_value = tumor_orig[[other]][valueat]
		component_values[[id]] = append(component_values[[id]], other_value)
		if(rm.na) component_values[[id]] = component_values[[id]][!is.na(component_values[[id]])]
	}
	
	if(verbose){
		for(i in 1:length(components)){
			cat("\nComponent (",componentname,") values for component: ", components[i],"\n")
			cat("\nValues (",other,"):",paste(component_values[[i]],sep=" , "),"\n")
		}
	}
	cat("\n\n")

	return(component_values)
}

#
#	Function for extracting a component from the data.frame and comparing it to a column from the same data.frame
#
xeno.test.marker.data = function(
	orig_data, 
	componentname="Growth",
	other="", 
	idname="Tumor_id",
	tpname="Timepoint",
	valueat=1,
	verbose=TRUE,
	rm.na=TRUE
){
	require(lme4)
	check_headers(orig_data, name=c(tpname, idname), mer=FALSE, dataframe=TRUE)
	if(!other %in% names(orig_data)){
		stop(paste("Could not find column from the data frame:",other))
	}else if(!componentname %in% names(orig_data)){
		stop(paste("Could not find column from the data frame:",componentname))
	}
	
	tumors = as.vector(unique(orig_data[[idname]]))
	components = unique(orig_data[[componentname]])
	components = components[order(components)]
	
	if(verbose) cat("\nComponents found:",paste(components, sep=" , "),"\n")
	
	markers = orig_data[[other]]
	non.na.markers = markers[!is.na(markers)]
	
	component_values = list()
	for(i in 1:length(components)){
		component_values[[i]] = vector(length=0)
	}
	
	
	# Looping over the different component values
	for(i in 1:length(tumors)){
		tumor_orig = subset(orig_data, orig_data[idname]==tumors[i])
		id = which(tumor_orig[[componentname]][1] == components)
		other_value = tumor_orig[[other]][valueat]
		component_values[[id]] = append(component_values[[id]], other_value)
		if(rm.na) component_values[[id]] = component_values[[id]][!is.na(component_values[[id]])]
	}
	
	if(verbose){
		for(i in 1:length(components)){
			cat("\nComponent (",componentname,") values for other component (", other, "): ", components[i],"\n")
			cat("\nValues (",other,"):",paste(component_values[[i]],sep=" , "),"\n")
		}
	}
	cat("\n\n")

	component_values
}

#
#	Function for comparing a particular random effect with a column of the original data
#
xeno.test.ran = function(
	fit, 
	orig_data, 
	rand_component="Timepoint", 
	other="", 
	treatmentname="Treatment", 
	idname="Tumor_id", 
	valueat=1, 
	cormethod="pearson"
){
	require(lme4)
	growthname="Growth"
	if(!other %in% names(orig_data)){
		stop(paste("Could not find column from the original data:",other))
	}
	ran_efs = ranef(fit)[[1]][rand_component]
	labels = rownames(ran_efs)
	treatments = unique(orig_data[[treatmentname]])
	treatments = treatments[order(treatments)]

	tumors_tr0 = unique(subset(orig_data, orig_data[treatmentname]==treatments[1])[[idname]])
	tumors_tr1 = unique(subset(orig_data, orig_data[treatmentname]==treatments[2])[[idname]])
	
	other_values = vector(length=0)
	for(i in 1:length(labels)){
		other_values = append(other_values, orig_data[which(orig_data[[idname]]==labels[i]),other][valueat])
	}
	print(paste("Correlation",cormethod,"for random effect",rand_component,"and data.frame column",other,"(whole data)"))
	print(cor(ran_efs[,1], other_values, use="pairwise.complete.obs", method=cormethod))
	a = cbind(ran_efs, other=other_values)
	rownames(a) = labels
	
	tr0 = a[which(rownames(a) %in% tumors_tr0),]
	tr1 = a[which(rownames(a) %in% tumors_tr1),]

	print(paste("Correlation",cormethod,"for random effect",rand_component,"and data.frame column",other,"(",treatmentname,":",treatments[1],")"))
	print(cor(tr0[,1], tr0[,2], use="pairwise.complete.obs", method=cormethod))
	print(paste("Correlation",cormethod,"for random effect",rand_component,"and data.frame column",other,"(",treatmentname,":",treatments[2],")"))
	print(cor(tr1[,1], tr1[,2], use="pairwise.complete.obs", method=cormethod))
	
	if(growthname %in% names(fit@frame) & length(unique(fit@frame[[growthname]]))==2){
		growths = unique(fit@frame[[growthname]])
		growths = growths[order(growths)]
		growths_0 = unique(subset(orig_data, fit@frame[growthname]==growths[1])[[idname]])
		growths_1 = unique(subset(orig_data, fit@frame[growthname]==growths[2])[[idname]])
	
		grwth0 = a[which(rownames(a) %in% growths_0),]
		grwth1 = a[which(rownames(a) %in% growths_1),]

		print(paste("Correlation",cormethod,"for random effect",rand_component,"and data.frame column",other,"(",growthname,":",growths[1],")"))
		print(cor(grwth0[,1], grwth0[,2], use="pairwise.complete.obs", method=cormethod))
		print(paste("Correlation",cormethod,"for random effect",rand_component,"and data.frame column",other,"(",growthname,":",growths[2],")"))
		print(cor(grwth1[,1], grwth1[,2], use="pairwise.complete.obs", method=cormethod))

		list(a, tr0, tr1, grwth0, grwth1)			
	}else{
		list(a, tr0, tr1)
	}
}	


#
#	Function for summarizing characteristics of the data
#
xeno.summary = function(
	orig_data,
	responsename="Response", 
	tpname="Timepoint", 
	treatmentname="Treatment", 
	idname="Tumor_id"	
){
	check_headers(orig_data, name=c(responsename, treatmentname, tpname, idname), mer=FALSE, dataframe=TRUE)	

	growthname="Growth"
	groups = unique(orig_data[[treatmentname]])
	cat("\norig_data fields available in orig_data frame:", paste(names(orig_data),sep=","))
	cat(paste("\nAmount of groups '",treatmentname, "' in orig_data frame:",length(groups)))
	cat("\nGroup indicators:", paste(groups))
	group_subsets = list()
	tumors = list()
	tumoramount = list()
	obsN = list()
	maxN = list()
	missingvalues = list()
	for(i in 1:length(groups)){
		group_subsets[[i]] = subset(orig_data, orig_data[treatmentname]==groups[i])
		tumors[[i]] = unique(group_subsets[[i]][[idname]])
		cat(paste("\n\nGroup",groups[i]))
		tumoramount[[i]] = length(tumors[[i]])
		cat(paste("\nUnique unit amount:", tumoramount[[i]]))
		cat("\nLabels:", paste(tumors[[i]]))
		cat(paste("\nFirst measured TP:", min(group_subsets[[i]][[tpname]], na.rm=TRUE)))
		cat(paste("\nLast measured TP:", max(group_subsets[[i]][[tpname]], na.rm=TRUE)))
		tpamount = length(min(group_subsets[[i]][[tpname]]):max(group_subsets[[i]][[tpname]]))
		cat(paste("\nTP amount:", tpamount))
		obsN[[i]] = 0
		for(j in 1:length(group_subsets[[i]][[responsename]])){
			if(!is.na(group_subsets[[i]][[responsename]][j])){
				obsN[[i]] = obsN[[i]] + 1
			}
		}
		cat(paste("\nObservations N:", obsN[[i]]))
		maxN[[i]] = tumoramount[[i]]*tpamount
		missingvalues[[i]] = maxN[[i]] - obsN[[i]]
		cat(paste("\nMissing values:", missingvalues[[i]]))
		missingprop = round((missingvalues[[i]]/(maxN[[i]]))*100, 2)
		cat(paste("\nMissing value percentage of all possible observations:", missingprop, "%"))
	}
	cat("\n\nAll groups combined:")
	cat("\nTotal amount of unique units:", sum(unlist(tumoramount)))
	cat("\nTotal amount of observations:", sum(unlist(obsN)))
	cat("\nTotal amount of missing observations:", sum(unlist(missingvalues)))	
	cat("\nMissing values percentage of all possible observations:", round((sum(unlist(missingvalues))/sum(unlist(maxN)))*100, 2), "%")
	cat("\n")
}

#
#	Function for extracting vector of the categorization for individual tumor units
#
xeno.cat = function(
	fit,
	tpname="Timepoint", 
	idname="Tumor_id"
){
	check_headers(fit, name=c(tpname, idname), mer=TRUE, dataframe=TRUE)	

	growthname="Growth"
	
	if(class(fit)=="mer"){
		data=  fit@frame
	}else if(class(fit)=="data.frame"){
		data = fit
	}else{
		stop("Invalid input type: should be either mer or data.frame")
	}
	
	firsttp = min(data[[tpname]])
	cut = subset(data, data[tpname]==firsttp)
	cats = as.vector(cut[[growthname]])
	if(!is.null(cats)){
		names(cats) = as.vector(cut[[idname]])
		cats = cats[unique(names(cats))]
	}
	cats
}

#
#	Function for extracting fuzzy categorization of tumors according to fitting of Gaussian mixture models with mixtools-package
#
xeno.cat.fuzzy = function(
	fit,
	tpname="Timepoint", 
	idname="Tumor_id",
	responsename="Response",
	mu = c(0.1,0.5,0.9),
	sigma = c(0.1,0.1,0.1),
	draw=TRUE,
	legendposition="top",
	rand=FALSE,
	randeff="Timepoint",
	...
){
	require(mixtools)
	require(lme4)
	check_headers(fit, name=c(responsename, tpname, idname), mer=TRUE, dataframe=FALSE)	

	
	growthname="Growth"
	if(!rand){
		if(identical(round(xeno.cat(fit, tpname=tpname, idname=idname),0),xeno.cat(fit, tpname=tpname, idname=idname))){
			stop("Fuzzy categorization ought to be used only for probabilistic Growth-covariate")
		}
	}
	# Using mixtools to fit the Gaussian mixture model for extracting fuzzy categories
	library(mixtools)
	
	colors = c("blue", "red", "green", "purple", rainbow(10))
	if(!rand){
		# Growth-covariate
		probcats = xeno.cat(fit=fit, tpname=tpname, idname=idname)
	}else{
		# Chosen random effect
		probcats = ranef(fit)[[1]][,randeff]
		names(probcats) = rownames(ranef(fit)[[1]])
	}
	mix = normalmixEM(x=probcats, mu = mu, sigma = sigma)
	ndist = length(mix$mu)
	mus = mix$mu
	sigmas = mix$sigma
	# Ordering the groups according to the mu in x-axis
	sigmas = sigmas[order(mus)]
	mus = mus[order(mus)]
	print("Mu:")
	print(mus)
	print("Sigma:")
	print(sigmas)
	
	groups = vector(length=length(probcats))
	for(i in 1:length(groups)){
		groups[i] = which(max(dnorm(probcats[i], mus, sigmas), na.rm=TRUE)==dnorm(probcats[i],mus,sigmas))
	}
	
	if(draw){
		if(!rand){
			maintitle="Probabilistic categorization"
		}else{
			maintitle="Random effect categorization"
		}
		ymin = 0
		ymax = max(dnorm(mus,mus,sigmas), rm.na=TRUE)
		if(min(probcats)>=0){
			xmin = 0
		}else{
			xmin = min(probcats, na.rm=TRUE)
		}
		if(max(probcats)<=1){
			xmax = 1
		}else{
			xmax = max(probcats, na.rm=TRUE)
		}
		plot.new()
		plot.window(xlim=c(xmin,xmax), ylim=c(ymin, ymax))
		if(!rand){
			title(maintitle, xlab="Probabilistic Growth-covariate", ylab="Density")
		}else{
			title(maintitle, xlab=paste("Random effect:",randeff), ylab="Density")
		}
		axis(1, at=seq(xmin,xmax,0.1))
		axis(2, las=1)	
		# Fitted Gaussian distributions
		for(i in 1:ndist){
			points(seq(xmin,xmax,0.01), dnorm(seq(xmin,xmax,0.01), mus[i], sigmas[i]), type="l", col=colors[i], lwd=2)
		}
		# Probabilistic categories
		points(probcats, abs(runif(length(probcats))*ymax/25), pch=20, col=colors[groups])
		if(!is.null(legendposition)){
			legend(legendposition, 
				c(paste("Subgroup",0:(ndist-1)),"Probabilistic\ncategorization"),
				lwd=c(rep(2,times=ndist),0),
				lty=c(rep(1,times=ndist),0),
				pch=c(rep(NA_integer_,times=ndist),20),
				col=c(colors[1:ndist],"black"),
				bty="n"
				)
		}
	}
	
	frame = fit@frame
	
	tumors = unique(frame[[idname]])
	for(i in 1:length(tumors)){
		frame[which(frame[[idname]]==names(probcats[i])),growthname] = groups[i]-1
	}
	frame
}

#
#	Function for drawing power curves from data generated with xeno.sim.power
#
xeno.draw.power = function(
	powermats=list(), 
	powerlabs=c(), 
	minpw=0, 
	maxpw=1, 
	minN=14, 
	maxN=50, 
	maintitle="Power curves", 
	drawcrit=0.8, 
	legendposition="right"
	)
{
	cols=c(1,2,3,4) 	
	plot.new()
	plot.window(xlim=c(minN,maxN), ylim=c(minpw, maxpw))
	title(maintitle, xlab="Tumor amount per group", ylab="Power")
	axis(1, at=as.numeric(rownames(powermats[[1]])))
	axis(2, las=1)
	
	colors = rainbow(length(cols), start=0.1, end=0.9)
	pchs = 1:length(powermats)
	
	for(i in 1:length(powermats)){
		for(j in 1:length(cols)){
			points(as.numeric(rownames(powermats[[i]])), powermats[[i]][,cols[j]], col=colors[j], pch=(20+pchs[i]), type="b")
		}
	}
	if(!is.na(drawcrit)){
		abline(col="red", h=drawcrit, lty=2)
	}
	if(!length(powerlabs)==length(powermats)){
		powerlabs = paste("Data #", 1:length(powermats))
	}
	if(!is.null(legendposition)){
		legend(legendposition, 
			c(powerlabs, colnames(powermats[[1]])[cols], "Power 0.8 criteria"), 
			lty=c(rep(0, times=length(powerlabs)), rep(1, times=length(cols)), 2), 
			pch=c(20+pchs, rep(NA_integer_, times=length(cols)), NA_integer_), 
			col=c(rep("black", times=length(powerlabs)), colors[1:length(cols)], "red"),
			cex=1,
			bty="n"
			)
	}
}

#
#	Function for drawing precision curves from data generated with functions like xeno.test.precision.fit
#
xeno.draw.precision = function(
	fits = list(), 
	testedK = list(), 
	Klabels = c(), 
	fitlabels = c(), 
	responsename = "Response", 
	tpname = "Timepoint", 
	treatmentname = "Treatment", 
	idname = "Tumor_id", 
	verbose = TRUE, 
	legendposition = "topright"
	)
{
	growthname="Growth"
	if(length(fits)==0 | !class(fits)=="list"){
		stop("Please supply a mer-object fitted with the lme4-package for the parameter fit")
	}
	if(!class(fits[[1]])=="mer"){
		stop("Please supply a mer-object fitted with the lme4-package for the parameter fit")
	}
	if(length(testedK)==0 | !class(testedK)=="list"){
		stop("Please supply a list of vectors for parameter testedK")
	}
	tpwise = list()
	Kwise = list()
	
	for(i in 1:length(fits)){
		# Fifth element is the precision
		tpwise[[i]] = xeno.test.precision.tps(fit=fits[[i]], responsename=responsename, tpname=tpname, treatmentname=treatmentname, idname=idname)
		Kwise[[i]] = list()
		for(j in 1:length(testedK)){
			Kwise[[i]][[j]] = xeno.test.precision.fit(fit=fits[[i]], K=testedK[[j]], responsename=responsename, tpname=tpname, treatmentname=treatmentname, idname=idname)
		}
	
	}
	# The two first rows are intercept and growth term
	
	tpmatrix = matrix(nrow=(length(tpwise[[1]])-2), ncol=length(tpwise))
	for(i in 1:(length(tpwise[[1]])-2)){
		for(j in 1:length(fits)){
			tpmatrix[i,j] = as.numeric(tpwise[[j]][[i]][5])
		}
	}
	if(!length(fitlabels)==length(fits)){
		fitlabels = paste("Fit #",1:length(fits),sep="")
	}
	Kmatrix = matrix(nrow=length(Kwise[[1]][[1]]), ncol=length(fits)*length(testedK))
	
	rows = vector(length=0)
	for(i in 1:length(Kwise[[1]][[1]])){
		rows = append(rows, paste("Omitted",Kwise[[1]][[1]][[i]][1],"TP"))
	}
	rownames(Kmatrix) = rows
	index=1
	for(i in 1:length(fits)){
		for(j in 1:length(testedK)){
			for(k in 1:length(Kwise[[1]][[1]])){
				Kmatrix[k,index] = Kwise[[i]][[j]][[k]][5]
			}
		index = index + 1
		}
	}
	
	par(mfrow=c(2,1))
	
	plot.new()
	plot.window(xlim=c(0,length(Kmatrix[,1])-1), ylim=c(as.numeric(min(Kmatrix, na.rm=TRUE)), as.numeric(max(Kmatrix, na.rm=TRUE))))
	colors = rainbow(length(testedK), start=0.1, end=0.9)
	pchs = 20 + 1:length(fits)
	
	title(main="Precisions", xlab="Time points omitted from end", ylab="Precision")
	axis(1, at=0:(length(Kmatrix[,1])-1), labels=(length(Kmatrix[,1])-1):0)
	axis(2, las=1)
	index = 1
	for(i in 1:length(fits)){
		for(j in 1:length(testedK)){
			points((length(Kmatrix[,index]):1)-1, Kmatrix[,index], type="b", col=colors[j], pch=pchs[i])
			index = index + 1
		}
		
	}
	
	if(!is.null(legendposition)){
		legend(legendposition, 
			append(fitlabels, Klabels), 
			lty=append(rep(0, times=length(fitlabels)), 
			rep(1, times=length(Klabels))), 
			pch=append(pchs, rep(NA_integer_, times=length(Klabels))), 
			col=append(rep("black", times=length(fitlabels)), colors[1:length(Klabels)]),
			bty="n"
		)
	}
	
	plot.new()
	plot.window(xlim=c(0,nrow(tpmatrix)-1), ylim=c(as.numeric(min(tpmatrix, na.rm=TRUE)), as.numeric(max(tpmatrix, na.rm=TRUE))), log="y")
	par("ylog")
	title(main="Precisions", xlab="Treatment-term specific to the time point", ylab="Precision (log-scale)")
	axis(1, at=0:(nrow(tpmatrix)-1))
	axis(2, las=1)	
	colors = c("purple")
	pchs = 20 + 1:length(fits)
	for(i in 1:length(fits)){	
		points(as.numeric(1:nrow(tpmatrix))-1, as.numeric(tpmatrix[,i]), type="b", col=colors[1], pch=pchs[i])
	}

	if(!is.null(legendposition)){
		legend(legendposition, 
			fitlabels, 
			pch=pchs, 
			col=rep(colors[1], times=length(fits)),
			bty="n"
		)
	}
	
	list(tpmatrix, Kmatrix)
	
	
}

#
#	Function for testing the precision of a fitted mixed effects model with a vector K for the tested hypothesis
#
xeno.test.precision.fit = function(
	fit, 
	K = c(0,0,0,1), 
	responsename="Response", 
	tpname="Timepoint", 
	treatmentname="Treatment", 
	idname="Tumor_id"
){
	require(lme4)
	check_headers(fit, name=c(responsename, treatmentname, tpname, idname), mer=TRUE, dataframe=FALSE)	

	growthname="Growth"
	mintp = min(fit@frame[[tpname]])
	maxtp = max(fit@frame[[tpname]])
	fit_data = fit@frame

	fit_formula = as.formula(fit@call)
	results = list()
	
	for(i in 0:(length(mintp:maxtp)-2)){
	
		new_data = subset(fit_data, fit_data[tpname]<=(maxtp-i))
		fitted_subset = lmer(fit_formula, data=new_data)

		results[[i+1]] = append(i, xeno.test.precision(fit=fitted_subset, K=K, responsename=responsename))
	}
	
	results
}


#
#	Function for extracting fitted model's power, lambda and precision according to Stroup
#
xeno.test.precision = function(
	fit, 
	K=c(0,1,0,0), 
	alpha=0.05, 
	responsename="Response"
){
	require(lme4)
	check_headers(fit, name=c(responsename), mer=TRUE, dataframe=FALSE)		
	
	# Internal function
	fitsigma =extract_str(fit, type="sigma")
	# Namespace issues in R 2.14
	#V = (crossprod(as.matrix(fit@A)) + diag(1, ncol(fit@A)))*(summary(fit)@sigma^2)
	V = (crossprod(as.matrix(fit@A)) + diag(1, ncol(fit@A)))*(fitsigma^2)
	Vinv = solve(V)
	b = as.matrix(fixef(fit))
	X = fit@X
	Xt = t(X)
	Kt = t(K)
	n = length(fit@frame[[responsename]][!is.na(fit@frame[[responsename]])])
	rankK = qr(K)$rank
	rankX = qr(X)$rank
	Fstat = as.numeric((t(Kt %*% b) %*% solve( Kt %*% solve( Xt %*% Vinv %*% X ) %*% K ) %*% (Kt %*% b))/rankK)
	# Non centrality parameter
	lambda = as.numeric(t(Kt %*% b) %*% solve( Kt %*% solve( Xt %*% Vinv %*% X ) %*% K ) %*% (Kt %*% b))
	Fcrit = qf(1-alpha, rankK, n-rankX)
	power = 1-pf(Fcrit, rankK, n-rankX, ncp=lambda)
	varKB = Kt %*% solve(Xt %*% Vinv %*% X) %*% K
	precision = as.numeric(1/varKB)
	result = vector(length=4)
	result = c(lambda,Fcrit,power,precision)
	names(result)=c("Non centrality parameter lambda", paste("F-crit, alpha=",alpha), "Power", "Precision")
	result
}


#
#	Function for extracting the timepoint-specific treatment term precision
#
xeno.test.precision.tps = function(
	fit, 
	responsename="Response", 
	tpname="Timepoint", 
	treatmentname="Treatment", 
	idname="Tumor_id"
){
	growthname="Growth"
	check_headers(fit, name=c(responsename, treatmentname, tpname, idname), mer=TRUE, dataframe=FALSE)	
	
	timepoint_indicators = function(data, tpname="Timepoint"){
		mintp = min(data[[tpname]], na.rm=TRUE)
		maxtp = max(data[[tpname]], na.rm=TRUE)
		tps = mintp:maxtp
		for(i in 1:length(tps)){
			indicator_vector = 1*(data[[tpname]]==tps[i])
			data[paste("TP",tps[i],sep="")] = indicator_vector
		}
		data
	}
	
	mintp = min(fit@frame[[tpname]])
	maxtp = max(fit@frame[[tpname]])
	fit_data = fit@frame
	
	fit_data = timepoint_indicators(data=fit_data, tpname=tpname)
	
	if(!growthname %in% names(fit@frame)){
		fit_fix_base = paste(responsename, " ~ 1 + ", tpname, sep="")
	}else{
		fit_fix_base = paste(responsename, " ~ 1 + ", tpname, ":", growthname, sep="")
	}
	fit_random_base = paste(" +  (1 | ",idname, ") + (0 + ", tpname, " | ", idname, " )", sep="")
	for(i in mintp:maxtp){
		fit_fix_base = paste(fit_fix_base, " + TP", i, ":", treatmentname, sep="")
	}
	fit_formula = as.formula(paste(fit_fix_base, fit_random_base))
	#print(paste("Using formula: ", fit_formula))
	
	
	results = list()
	fitted_tp = lmer(fit_formula, data=fit_data)
	for(i in 1:length(fixef(fitted_tp))){
		Kvector = 0*vector(length=length(fixef(fitted_tp)))
		Kvector[i] = 1
		results[[i]] = append(paste("Term:",names(fixef(fitted_tp))[i]), xeno.test.precision(fit=fitted_tp, K=Kvector, responsename=responsename))
	}
	results
}


#
#	Function for extracting a 2x2 contingency table for testing with fisher's test
#
xeno.test.cat = function(
	x,
	treatmentname="Treatment",
	tpname="Timepoint",
	idname="Tumor_id"
){
	check_headers(x, name=c(treatmentname, tpname, idname), mer=TRUE, dataframe=TRUE)	

	growthname="Growth"
	# Defining the matrix
	res = matrix(nrow=2,ncol=2)
	rownames(res)=c("Growing","Poorly growing")
	colnames(res)=c("Control","Treatment")
	# Extracting the table depending on if it's the data.frame or a lme4-x
	if(class(x)=="data.frame"){
		res[1,1] = length(unique(subset(x, x[treatmentname]==0 & x[growthname]==1)[[idname]]))
		res[2,1] = length(unique(subset(x, x[treatmentname]==0 & x[growthname]==0)[[idname]]))
		res[1,2] = length(unique(subset(x, x[treatmentname]==1 & x[growthname]==1)[[idname]]))
		res[2,2] = length(unique(subset(x, x[treatmentname]==1 & x[growthname]==0)[[idname]]))
	}else if(class(x)=="mer"){
		frame = x@frame
		IDs = length(unique(frame[[idname]]))
		tps = length(unique(frame[[tpname]]))
		# Check whether it's a non-linear fit, where each fixed effect causes the frame to multiple rows in dimension
		if(length(frame[[idname]])>IDs*tps){
			frame = frame[names(fixef(x))[1]]
		}
		res[1,1] = length(unique(subset(x@frame, x@frame[treatmentname]==0 & x@frame[growthname]==1)[[idname]]))
		res[2,1] = length(unique(subset(x@frame, x@frame[treatmentname]==0 & x@frame[growthname]==0)[[idname]]))
		res[1,2] = length(unique(subset(x@frame, x@frame[treatmentname]==1 & x@frame[growthname]==1)[[idname]]))
		res[2,2] = length(unique(subset(x@frame, x@frame[treatmentname]==1 & x@frame[growthname]==0)[[idname]]))
	# Error, neither of the required classes
	}else{
		stop("x must be a data.frame or a mer (lme4) object")
	}
	# Returning the 2x2 contingency table
	res
}


#
#	Drawing functions: data from a fitted model's frame
#
xeno.draw.data = function(
	x,
	responsename="Response", 
	treatmentname="Treatment", 
	tpname="Timepoint", 
	idname="Tumor_id",
	maintitle="Data",
	legendposition="topleft"
){
	check_headers(x, name=c(responsename, treatmentname, tpname, idname), mer=TRUE, dataframe=TRUE)	

	turn_y_labs=TRUE
	#mar=c(5, 4, 1, 0.5)
	growthname="Growth"

	colours = c("black", "green")

	if(class(x)=="data.frame"){
		draw_data = x
	}else if(class(x)=="mer"){
		draw_data = x@frame
	}else{
		stop("Invalid x: should be either a data.frame or fitted lme4 (mer) object")
	}
	tumors = unique(draw_data[[idname]])
	ymax = max(draw_data[[responsename]], na.rm=TRUE)
	ymin = min(draw_data[[responsename]], na.rm=TRUE)
	tpmin = min(draw_data[[tpname]], na.rm=TRUE)
	tpmax = max(draw_data[[tpname]], na.rm=TRUE)
	
	#par(mar=mar)
	plot.new()
	plot.window(xlim=c(tpmin,tpmax), ylim=c(ymin,ymax))
	title(maintitle, xlab=tpname, ylab=responsename)
	axis(1, at=tpmin:tpmax)
	if(turn_y_labs){
		axis(2, las=1)
	}else{
		axis(2)
	}
	for(i in 1:length(tumors)){
		current_tumor = subset(draw_data, draw_data[[idname]]==tumors[i])
		#points(unique(current_tumor[[tpname]]), current_tumor[[responsename]], type="l", lwd=2, col=colours[current_tumor[[treatmentname]][1]+1], pch=21)
		points(current_tumor[[tpname]], current_tumor[[responsename]], type="l", lwd=2, col=colours[current_tumor[[treatmentname]][1]+1], pch=21)
	}

	if(!is.null(legendposition)){
		legend(legendposition, c("Control", "Treatment"), col=c(colours[1], colours[2]), pch=c(NA_integer_, NA_integer_), lwd=2, box.col=NA)
	}
}

#
#	Drawing functions: fixed fit from a fitted model
#
xeno.draw.fixed = function(
	fit,
	orig_data,
	responsename="Response", 
	treatmentname="Treatment", 
	tpname="Timepoint", 
	idname="Tumor_id",
	maintitle="Fixed effect fit",
	legendposition="topleft",
	drawgroups=TRUE,
	drawcat=FALSE,
	colorgroups=FALSE
){
	require(lme4)
	check_headers(fit, name=c(responsename, treatmentname, idname), mer=TRUE, dataframe=FALSE)		

	turn_y_labs=TRUE
	#mar=c(5, 4, 1, 0.5)
	growthname="Growth"

	data = fit@frame
	data[["Fixed_fit"]] = as.vector((fit@X) %*% fixef(fit))
	data[["Whole_fit"]] = as.vector((fit@X) %*% fixef(fit) + t(as.matrix(fit@Zt)) %*% unlist(ranef(fit)))
	formula = as.formula(fit@call)

	colours = c("black", "green")

	if(!(tpname %in% names(data))){
		data[[tpname]] = orig_data[[tpname]][!is.na(orig_data[[responsename]])]
	}	

	tumors = unique(data[[idname]])
	ymax = max(data[[responsename]], na.rm=TRUE)
	ymin = min(data[[responsename]], na.rm=TRUE)
	tpmin = min(data[[tpname]], na.rm=TRUE)
	tpmax = max(data[[tpname]], na.rm=TRUE)
		
	#par(mar=mar)
	plot.new()
	plot.window(xlim=c(tpmin,tpmax), ylim=c(ymin,ymax))
	title(maintitle, xlab=tpname, ylab=responsename)
	axis(1, at=tpmin:tpmax)
	if(turn_y_labs){
		axis(2, las=1)
	}else{
		axis(2)
	}


	if(drawgroups && !drawcat){
		for(i in 1:length(tumors)){
			current_tumor = subset(data, data[[idname]]==tumors[i])
			if(!colorgroups){
				points(unique(current_tumor[[tpname]]), current_tumor[[responsename]], type="l", lwd=1, col="grey")
			}else{
				if(current_tumor[[treatmentname]][1]==0){
					points(unique(current_tumor[[tpname]]), current_tumor[[responsename]], type="l", lwd=1, col="black")
				}else{
					points(unique(current_tumor[[tpname]]), current_tumor[[responsename]], type="l", lwd=1, col="green")
				}
			}
		}
	}else if(drawcat && growthname %in% names(data)){
		# Fixed color palettes
		discrete.colours = c("blue", "red", "green", "purple")
		continuous.colours= rainbow(1000, start=0.7, end=0.95)
	
		# Drawing discriminated categories with blue/red/etc
		if(identical(round(xeno.cat(fit, tpname=tpname, idname=idname),0),xeno.cat(fit, tpname=tpname, idname=idname))){
			cols = discrete.colours
			for(i in 1:length(tumors)){
				current_tumor = subset(data, data[[idname]]==tumors[i])
				points(unique(current_tumor[[tpname]]), current_tumor[[responsename]][1:length(unique(current_tumor[[tpname]]))], type="l", lwd=2, col=cols[current_tumor[[growthname]][1]+1])
			}
		}else {
		# Drawing a continuous range of categories with a colour scale from blue to red
			cols = continuous.colours
			for(i in 1:length(tumors)){
				current_tumor = subset(data, data[[idname]]==tumors[i])
				points(unique(current_tumor[[tpname]]), current_tumor[[responsename]][1:length(unique(current_tumor[[tpname]]))], type="l", lwd=2, col=cols[round((999*current_tumor[[growthname]][1]+1),0)])
			}
		}
	}


	for(i in 1:length(tumors)){
		current_tumor = subset(data, data[[idname]]==tumors[i])
		if(!treatmentname %in% names(current_tumor)){
			treatmentgroup = 0
		}else{
			treatmentgroup = current_tumor[[treatmentname]][1]
		}
		if(growthname %in% names(current_tumor)){
			if(current_tumor[[growthname]][1]==0){
				points(unique(current_tumor[[tpname]]), current_tumor[["Fixed_fit"]], type="l", lwd=3, lty=3, col=colours[treatmentgroup+1], pch=21)
			}else{
				points(unique(current_tumor[[tpname]]), current_tumor[["Fixed_fit"]], type="l", lwd=3, lty="dashed", col=colours[treatmentgroup+1], pch=22)
			}
		}else{
			points(unique(current_tumor[[tpname]]), current_tumor[["Fixed_fit"]], type="l", lwd=3, col=colours[treatmentgroup+1], pch=21)
		}
	}
	
	if(!is.null(legendposition)){
		if(!drawgroups){
			legend(legendposition, c("Growing fixed effects","Poorly growing fixed effects", "Control", "Treatment", "Data"), col=c("black", "black", colours[1], colours[2], "grey"), pch=c(NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_), lty=c(2,3,1,1,1), lwd=c(3,3,3,3,1), box.col=NA)
		}else{
			legend(legendposition, c("Growing fixed effects","Poorly growing fixed effects", "Control", "Treatment"), col=c("black", "black", colours[1], colours[2]), pch=c(NA_integer_, NA_integer_, NA_integer_, NA_integer_), lty=c(2,3,1,1), lwd=c(3,3,3,3), box.col=NA)
		}
	}

}

#
#	Drawing functions: Fixed + random effects (full fit) from a fitted model's frame
#
xeno.draw.fit = function(
	fit,
	orig_data,
	responsename="Response", 
	treatmentname="Treatment", 
	tpname="Timepoint", 
	idname="Tumor_id",
	maintitle="Fixed + Random effect fit",
	legendposition="topleft",
	drawgroups=TRUE,
	drawcat=FALSE,
	colorgroups=FALSE
){
	require(lme4)
	check_headers(fit, name=c(responsename, treatmentname, idname), mer=TRUE, dataframe=FALSE)	

	turn_y_labs=TRUE
	#mar=c(5, 4, 1, 0.5)
	growthname="Growth"
	
	data = fit@frame
	data[["Fixed_fit"]] = as.vector((fit@X) %*% fixef(fit))
	data[["Whole_fit"]] = as.vector((fit@X) %*% fixef(fit) + t(as.matrix(fit@Zt)) %*% unlist(ranef(fit)))
	formula = as.formula(fit@call)

	colours = c("black", "green")

	if(!(tpname %in% names(data))){
		data[[tpname]] = orig_data[[tpname]][!is.na(orig_data[[responsename]])]
	}	

	tumors = unique(data[[idname]])
	ymax = max(data[[responsename]], na.rm=TRUE)
	ymin = min(data[[responsename]], na.rm=TRUE)
	tpmin = min(data[[tpname]], na.rm=TRUE)
	tpmax = max(data[[tpname]], na.rm=TRUE)
		
	#par(mar=mar)
	plot.new()
	plot.window(xlim=c(tpmin,tpmax), ylim=c(ymin,ymax))
	title(maintitle, xlab=tpname, ylab=responsename)
	axis(1, at=tpmin:tpmax)
	if(turn_y_labs){
		axis(2, las=1)
	}else{
		axis(2)
	}

	if(drawgroups && !drawcat){
		for(i in 1:length(tumors)){
			current_tumor = subset(data, data[[idname]]==tumors[i])
			if(!colorgroups){
				points(unique(current_tumor[[tpname]]), current_tumor[[responsename]], type="l", lwd=1, col="grey")
			}else{
				if(current_tumor[[treatmentname]][1]==0){
					points(unique(current_tumor[[tpname]]), current_tumor[[responsename]], type="l", lwd=1, col="black")
				}else{
					points(unique(current_tumor[[tpname]]), current_tumor[[responsename]], type="l", lwd=1, col="green")
				}
			}
		}
	}else if(drawcat && growthname %in% names(data)){
		# Fixed color palettes
		discrete.colours = c("blue", "red", "green", "purple")
		continuous.colours= rainbow(1000, start=0.7, end=0.95)
	
		# Drawing discriminated categories with blue/red/etc
		if(identical(round(xeno.cat(fit, tpname=tpname, idname=idname),0),xeno.cat(fit, tpname=tpname, idname=idname))){
			cols = discrete.colours
			for(i in 1:length(tumors)){
				current_tumor = subset(data, data[[idname]]==tumors[i])
				points(unique(current_tumor[[tpname]]), current_tumor[[responsename]][1:length(unique(current_tumor[[tpname]]))], type="l", lwd=2, col=cols[current_tumor[[growthname]][1]+1])
			}
		}else {
		# Drawing a continuous range of categories with a colour scale from blue to red
			cols = continuous.colours
			for(i in 1:length(tumors)){
				current_tumor = subset(data, data[[idname]]==tumors[i])
				points(unique(current_tumor[[tpname]]), current_tumor[[responsename]][1:length(unique(current_tumor[[tpname]]))], type="l", lwd=2, col=cols[round((999*current_tumor[[growthname]][1]+1),0)])
			}
		}
	}

	for(i in 1:length(tumors)){
		current_tumor = subset(data, data[[idname]]==tumors[i])
		if(!treatmentname %in% names(current_tumor)){
			treatmentgroup = 0
		}else{
			treatmentgroup = current_tumor[[treatmentname]][1]
		}
		if(growthname %in% names(current_tumor)){
			if(current_tumor[[growthname]][1]==0){
				points(unique(current_tumor[[tpname]]), current_tumor[["Whole_fit"]], type="l", lwd=2, lty=3, col=colours[treatmentgroup+1], pch=21)
			}else{
				points(unique(current_tumor[[tpname]]), current_tumor[["Whole_fit"]], type="l", lwd=2, lty="dashed", col=colours[treatmentgroup+1], pch=22)
			}
		}else{
			points(unique(current_tumor[[tpname]]), current_tumor[["Whole_fit"]], type="l", lwd=2, col=colours[treatmentgroup+1], pch=21)
		}
	}

	if(!is.null(legendposition)){
		legend(legendposition, c("Growing fit","Poorly growing fit", "Control", "Treatment", "Data"), col=c("black", "black", colours[1], colours[2], "grey"), pch=c(NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_), lty=c(2,3,1,1,1), lwd=c(2,2,2,2,1), box.col=NA)
	}
}

#
#	Drawing functions: identified categories by EM
#
xeno.draw.cat = function(
	fit,
	orig_data,
	responsename="Response", 
	treatmentname="Treatment", 
	tpname="Timepoint", 
	idname="Tumor_id",
	maintitle="Found categorization",
	legendposition="topleft",
	ymax=NULL,
	ymin=NULL,
	scaled.growth=FALSE,
	discrete.colours = c("blue", "red", "green", "purple"),
	#continuous.colours= rainbow(1000, start=0.65, end=1)
	continuous.colours= rainbow(1000, start=0.7, end=0.95)
){
	require(lme4)
	check_headers(fit, name=c(responsename, treatmentname, tpname, idname), mer=TRUE, dataframe=TRUE)	
	
	turn_y_labs=TRUE
	#mar=c(5, 4, 1, 0.5)
	growthname="Growth"

	if(class(fit)=="mer"){
		data = fit@frame
	}else if(class(fit)=="data.frame"){
		data = fit
	}else{
		stop("Invalid input type.")
	}
	#data[["Fixed_fit"]] = as.vector((fit@X) %*% fixef(fit))
	#data[["Whole_fit"]] = as.vector((fit@X) %*% fixef(fit) + t(as.matrix(fit@Zt)) %*% unlist(ranef(fit)))
	formula = as.formula(fit@call)

	tumors = unique(data[[idname]])
	if(is.null(ymax)){
		ymax = max(data[[responsename]], na.rm=TRUE)
	}
	if(is.null(ymin)){
		ymin = min(data[[responsename]], na.rm=TRUE)
	}

	if(!(tpname %in% names(data))){
		data[[tpname]] = orig_data[[tpname]][!is.na(orig_data[[responsename]])]
	}	

	tpmin = min(data[[tpname]], na.rm=TRUE)
	tpmax = max(data[[tpname]], na.rm=TRUE)
	
	#par(mar=mar)
	plot.new()
	plot.window(xlim=c(tpmin,tpmax), ylim=c(ymin,ymax))
	title(maintitle, xlab=tpname, ylab=responsename)
	axis(1, at=tpmin:tpmax)
	if(turn_y_labs){
		axis(2, las=1)
	}else{
		axis(2)
	}
	
	if(scaled.growth & growthname %in% names(data)){
		mingrowth = min(data[[growthname]])
		maxgrowth = max(data[[growthname]])
		data[growthname] = (data[growthname]-mingrowth)/(maxgrowth-mingrowth)
	}
	
	if(growthname %in% names(data)){
		# Drawing discriminated categories with blue/red/etc
		if(identical(round(xeno.cat(fit, tpname=tpname, idname=idname),0),xeno.cat(fit, tpname=tpname, idname=idname))){
			cols = discrete.colours
			for(i in 1:length(tumors)){
				current_tumor = subset(data, data[[idname]]==tumors[i])
				points(unique(current_tumor[[tpname]]), current_tumor[[responsename]][1:length(unique(current_tumor[[tpname]]))], type="l", lwd=2, col=cols[current_tumor[[growthname]][1]+1])
			}
		}else {
		# Drawing a continuous range of categories with a colour scale from blue to red
			cols = continuous.colours
			for(i in 1:length(tumors)){
				current_tumor = subset(data, data[[idname]]==tumors[i])
				points(unique(current_tumor[[tpname]]), current_tumor[[responsename]][1:length(unique(current_tumor[[tpname]]))], type="l", lwd=2, col=cols[round((999*current_tumor[[growthname]][1]+1),0)])
			}
		}
	}else{
		for(i in 1:length(tumors)){
			current_tumor = subset(data, data[[idname]]==tumors[i])
			points(unique(current_tumor[[tpname]]), current_tumor[[responsename]][1:length(unique(current_tumor[[tpname]]))], type="l", lwd=2, col="orange", pch=23)
		}
	}
	
	if(!is.null(legendposition) & length(unique(data[[growthname]]))==2){
		legend(legendposition, c("Growth=0", "Growth=1"), col=c(cols[1], cols[2]), lwd=c(2,2), box.col=NA)
	}else if(!growthname %in% names(data)){
		legend(legendposition, 
			c("No growth category detected"), 
			col=c("orange"), 
			lwd=c(2), 
			box.col=NA)
	}else if(identical(round(xeno.cat(data, tpname=tpname, idname=idname),0),xeno.cat(data, tpname=tpname, idname=idname))){
		legend(legendposition, 
			c(paste("Growth=",min(data[[growthname]], na.rm=TRUE):max(data[[growthname]], na.rm=TRUE))), 
			col=c(cols[1:length(unique(data[[growthname]]))]), 
			lwd=c(rep(2,times=length(unique(data[[growthname]])))), 
			box.col=NA)
	}else if(!is.null(legendposition)){
		yrange = ymax - ymin
		ylow = 0.85*yrange
		xrange = tpmax - tpmin
		xleft = 0.025*xrange
		xright = 0.025*xrange + 0.075*xrange
		for(i in 1:length(continuous.colours)){
			ychange = (i/length(continuous.colours))*0.1*yrange
			lines(c(xleft, xright), c(ylow + ychange, ylow + ychange), col=continuous.colours[i])
		}
		text(c(xright+0.0*xrange, xright+0.0*xrange), 
			c(ylow, ylow+0.1*yrange), 
			labels=c("Growth=0", "Growth=1"), 
			pos=4,
			col=c(continuous.colours[1],continuous.colours[length(continuous.colours)]))
	}
}

#
#	Drawing functions: residuals
#
xeno.draw.res = function(
	fit,
	responsename="Response", 
	maintitle="Residuals",
	drawline=FALSE,
	drawquad=FALSE
){
	check_headers(fit, name=c(responsename), mer=TRUE, dataframe=FALSE)	

	#mar=c(5, 4, 1, 0.5)
	turn_y_labs=TRUE
	
	data = fit@frame[1:length(resid(fit)),]
	
	colours = c("black", "green")

	resxmin = min(data[[responsename]], na.rm=TRUE)
	resxmax = max(data[[responsename]], na.rm=TRUE)
	resymin = min(resid(fit), na.rm=TRUE)
	resymax = max(resid(fit), na.rm=TRUE)

	resfit = lm(resid(fit) ~ data[[responsename]])

	#par(mar=mar)
	plot.new()
	plot.window(xlim=c(resxmin,resxmax), ylim=c(resymin,resymax))
	box()
	title(maintitle, xlab="Observed values", ylab="Residuals")
	axis(1)
	if(turn_y_labs){
		axis(2, las=1)
	}else{
		axis(2)
	}
	
	points(data[[responsename]], resid(fit), pch=20, cex=1)
	if(drawline){
		abline(coef(resfit), col="orange", lwd=2)	
	}
	if(drawquad){
		xresp = data[[responsename]]
		xresp2 = xresp^2 
		quadfit = lm(resid(fit) ~ xresp + xresp2)
		a = coef(quadfit)[1]
		b = coef(quadfit)[2]
		c = coef(quadfit)[3]
		curve(a + b*x + c*x*x, col="purple", lwd=2, add=TRUE, 
			from=(resxmin-(0.05*resxmax)), to=(resxmax+(0.05*resxmax)))
	}
}

#
#	Drawing functions: proportional residuals
#
xeno.draw.propres = function(
	fit,
	orig_data,
	responsename="Response", 
	treatmentname="Treatment", 
	tpname="Timepoint", 
	idname="Tumor_id",
	maintitle="Proportional Residuals",
	drawline=FALSE,
	highlight=TRUE,
	times_sd=2,
	drawhighlightfit=FALSE
){
	check_headers(fit, name=c(responsename, treatmentname, tpname, idname), mer=TRUE, dataframe=FALSE)	

	growthname="Growth"
	#mar=c(5, 4, 1, 0.5)
	turn_y_labs=FALSE
	
	data = fit@frame
	data[["Fixed_fit"]] = as.vector((fit@X) %*% fixef(fit))
	data[["Whole_fit"]] = as.vector((fit@X) %*% fixef(fit) + t(as.matrix(fit@Zt)) %*% unlist(ranef(fit)))
	formula = as.formula(fit@call)
	
	colours = c("black", "green")

	if(!(tpname %in% names(data))){
		data[[tpname]] = orig_data[[tpname]][!is.na(orig_data[[responsename]])]
	}	

	tumors = unique(data[[idname]])
	ymax = max(data[[responsename]], na.rm=TRUE)
	ymin = min(data[[responsename]], na.rm=TRUE)
	tpmin = min(data[[tpname]], na.rm=TRUE)
	tpmax = max(data[[tpname]], na.rm=TRUE)

	prop_res = resid(fit)
	# Removing zeros
	prop_res = prop_res[!data[[responsename]]==0]
	obs = data[[responsename]][!data[[responsename]]==0]

	prop_res = prop_res/obs

	resxmin = min(obs, na.rm=TRUE)
	resxmax = max(obs, na.rm=TRUE)
	resymin = min(prop_res, na.rm=TRUE)
	resymax = max(prop_res, na.rm=TRUE)

	resfit = lm(prop_res ~ obs)
	
	if(!drawhighlightfit){
		#par(mar=mar)
		plot.new()
		plot.window(xlim=c(resxmin,resxmax), ylim=c(resymin,resymax))
		title(maintitle, xlab="Observed values", ylab="Proportional Residuals")
		box()
		axis(1)
		if(turn_y_labs){
			axis(2, las=1)
		}else{
			axis(2)
		}
		if(!highlight){
			points(obs, prop_res, pch=20, cex=1)
		}else{
			abovetwosd = abs(mean(prop_res)-prop_res)>abs(mean(prop_res)-times_sd*sd(prop_res))
			points(obs[abovetwosd], prop_res[abovetwosd], pch=20, cex=1, col="red")
			points(obs[!abovetwosd], prop_res[!abovetwosd], pch=20, cex=1, col="black")
		}
		if(drawline){
			abline(coef(resfit), col="black", lwd=2)	
		}
	}else{
		xeno.draw.fit(
			fit=fit, 
			orig_data=orig_data, 
			responsename=responsename, 
			treatmentname=treatmentname,
			idname=idname, 
			legendposition=NULL)
		abovetwosd = abs(mean(prop_res)-prop_res)>abs(mean(prop_res)-times_sd*sd(prop_res))
		points(fit@frame[[tpname]][!data[[responsename]]==0][abovetwosd], obs[abovetwosd], cex=2, col="red", pch=20)	
	}
}

#
#	Drawing functions: within-tumor autocorrelation
#
xeno.draw.autocor = function(
	fit,
	orig_data,
	responsename="Response", 
	treatmentname="Treatment", 
	tpname="Timepoint", 
	idname="Tumor_id",
	maintitle="Autocorrelation"
){
	check_headers(fit, name=c(responsename, treatmentname, tpname, idname), mer=TRUE, dataframe=FALSE)	

	growthname="Growth"

	data = fit@frame
	#data[["Fixed_fit"]] = as.vector((fit@X) %*% fixef(fit))
	#data[["Whole_fit"]] = as.vector((fit@X) %*% fixef(fit) + t(as.matrix(fit@Zt)) %*% unlist(ranef(fit)))
	formula = as.formula(fit@call)

	tumors = unique(data[[idname]])
	ymax = max(data[[responsename]], na.rm=TRUE)
	ymin = min(data[[responsename]], na.rm=TRUE)
	tpmin = min(data[[tpname]], na.rm=TRUE)
	tpmax = max(data[[tpname]], na.rm=TRUE)
	
	colours = c("black", "green")

	if(!(tpname %in% names(data))){
		data[[tpname]] = orig_data[[tpname]][!is.na(orig_data[[responsename]])]
	}	
	
	autocors = xeno.test.autocor(fit=fit, draw=TRUE, idname=idname, tpname=tpname, maintitle=maintitle)
}

#
#	Drawing functions: visualizing random effects of a fit model object
#
xeno.draw.ranef = function(
	fit,
	nclasses = 30
){
	if(!(class(fit) == "mer" || class(fit) == "lmerMod")){
		stop("Input object should be a mer-object")
	}

	# Number of random effects components in the model
	rands = 0
	for(i in 1:length(ranef(fit))){
		for(j in 1:(dim(ranef(fit)[[i]])[2])){
			rands = rands + 1
		}
	}

	# Drawing histograms
	par(mfrow=c(round(rands/2, 0), 2))	
	for(i in 1:length(ranef(fit))){
		for(j in 1:(dim(ranef(fit)[[i]])[2])){
			xlims = c(-abs(round(max(ranef(fit)[[i]][,j]),0)),abs(round(max(ranef(fit)[[i]][,j]),0)))
			hist(ranef(fit)[[i]][,j], xlim=xlims, nclass=nclasses, main=paste("Random effect",colnames(ranef(fit)[[i]])[j],"grouped by",names(ranef(fit))[i]), xlab=colnames(ranef(fit)[[i]])[j])
		}
	}
}

#
#	Drawing functions: steps of the EM-algorithm
#
xeno.draw.EM = function(
	formula,
	data,
	mfrow=NULL,
	responsename="Response",
	tpname="Timepoint",
	idname="Tumor_id",
	treatmentname="Treatment",
	legendposition="topleft"
){
	require(lme4)
	check_headers(data, name=c(responsename, treatmentname, tpname, idname), mer=FALSE, dataframe=TRUE)	

	growthname="Growth"
	mar=c(5, 4, 4, 2) + 0.1
	turn_y_labs=TRUE
	
	fits = xeno.EM(data=data, formula=formula, max.iter=100, return.iterations=TRUE, responsename=responsename, tpname=tpname, idname=idname, treatmentname=treatmentname)[[3]]

	timepoints = min(data[[tpname]], na.rm=TRUE):max(data[[tpname]], na.rm=TRUE)
	y = data[[responsename]]
	tumors=unique(data[[idname]])
	x = rep(timepoints, length(tumors))
	colours = c("black", "green")
	
	# Removing datapoints with missing response value
	x = x[!is.na(y)]
	labels = data[[idname]]
	labels = labels[!is.na(y)]
	y = y[!is.na(y)]

	if(!is.null(mfrow)){
		par(mfrow=mfrow, mar=mar)
	}else{
		#par(mar=mar)
	}

	# Plotting iterative steps
	for(iter in 1:length(fits)){
		if(iter == 1 | iter == 2 | iter == 3 | iter == 4 | iter == 5 | iter == length(iter)){
			EM_data = fits[[iter]]@frame
			EM_data[["Fixed_fit"]] = as.vector((fits[[iter]]@X) %*% fixef(fits[[iter]]))
			EM_data[["Whole_fit"]] = as.vector((fits[[iter]]@X) %*% fixef(fits[[iter]]) + t(as.matrix(fits[[iter]]@Zt)) %*% unlist(ranef(fits[[iter]])))
			ymax = max(EM_data[[responsename]], na.rm=TRUE)
			ymin = min(EM_data[[responsename]], na.rm=TRUE)
			tpmin = min(EM_data[[tpname]], na.rm=TRUE)
			tpmax = max(EM_data[[tpname]], na.rm=TRUE)
			title = paste("EM iteration",iter)
			plot.new()
			plot.window(xlim=c(tpmin,tpmax), ylim=c(ymin,ymax))
			title(title, xlab=tpname, ylab=responsename)
			axis(1, at=tpmin:tpmax)
			if(turn_y_labs){
				axis(2, las=1)
			}else{
				axis(2)
			}
			tumors = unique(EM_data[[idname]])
			# INDIVIDUAL TUMORS AND THEIR CURRENT CATEGORIZATION
			for(i in 1:length(tumors)){
				current_tumor = subset(EM_data, EM_data[[idname]]==tumors[i])
				if(growthname %in% names(current_tumor)){
					#points(unique(current_tumor[[tpname]]), current_tumor[[responsename]], type="l", col="grey")
					if(current_tumor[[growthname]][1]==0){
						points(unique(current_tumor[[tpname]]), current_tumor[[responsename]], col="blue", type="l", lwd=1, pch=20, cex=1)
					}else{
						points(unique(current_tumor[[tpname]]), current_tumor[[responsename]], col="red", type="l", lwd=1, pch=20, cex=1)
					}
				}else{
					points(unique(current_tumor[[tpname]]), current_tumor[["Fixed_fit"]], type="l", lwd=3, col=colours[current_tumor[[treatmentname]][1]+1], pch=21)
				}
			}
			# FIXED EFFECTS
			for(i in 1:length(tumors)){
				current_tumor = subset(EM_data, EM_data[[idname]]==tumors[i])
				if(growthname %in% names(current_tumor)){
					if(current_tumor[[growthname]][1]==0){
						points(unique(current_tumor[[tpname]]), current_tumor[["Fixed_fit"]], type="l", lwd=3, lty=3, col=colours[current_tumor[[treatmentname]][1]+1], pch=21)
					}else{
						points(unique(current_tumor[[tpname]]), current_tumor[["Fixed_fit"]], type="l", lwd=3, lty="dashed", col=colours[current_tumor[[treatmentname]][1]+1], pch=22)
					}
				}else{
					points(unique(current_tumor[[tpname]]), current_tumor[["Fixed_fit"]], type="l", lwd=3, col=colours[current_tumor[[treatmentname]][1]+1], pch=21)
				}
			}
		}
		if(iter == 1 & !is.null(legendposition)){
			legend(legendposition, c("Growing", "Poorly growing"), col=c("red","blue"), lty=c(1,1), lwd=c(2,2), box.col=NA)
		}
		if(iter == 2 & !is.null(legendposition)){
			legend(legendposition, c(expression(paste("Growing fixed effects")),expression(paste("Poorly growing fixed effects"))), col=c("black", "black"), pch=c(NA_integer_, NA_integer_), lty=c(2,3), lwd=c(3,3), box.col=NA)
		}
	}
	if(!is.null(legendposition)){
		legend(legendposition, c("Control", "Treatment", "Data"), col=c(colours[1], colours[2], "grey"), pch=c(NA_integer_, NA_integer_, NA_integer_), lty=c(1,1,1), lwd=c(2,2,1), box.col=NA)
	}
}
	
#
#	Function for normalizing the response values by dividing/subtracting the starting values for tumors and then possibly taking a logarithm
#
xeno.norm.start = function(
	data,
	newfield="LogResponse",
	responsename="Response",
	tpname="Timepoint",
	idname="Tumor_id",
	log=TRUE,
	start=TRUE,
	type="divide"
){
	check_headers(data, name=c(responsename, tpname, idname), mer=FALSE, dataframe=TRUE)	

	responses = data[[responsename]]
	uniques = as.vector(unique(data[[idname]]))
	if(!log & !start){
		stop("No normalization method specified (log nor start TRUE)")
	}
	data[,newfield] = 0
	for(i in 1:length(uniques)){
		data_subset = subset(data, data[idname]==uniques[i])
		if(start==TRUE){
				if(type=="divide"){
					response_vector = data_subset[[responsename]] / data_subset[[responsename]][1]
				}else if(type=="subtract"){
					response_vector = data_subset[[responsename]] - data_subset[[responsename]][1]
				}else{
					stop("Incorrect type-parameter: should be either divide or subtract")
				}
		}else{
			response_vector = data_subset[[responsename]]
		}
		if(log){
			response_vector = log(response_vector)
		}
		data[which(data[[idname]]==uniques[i]),newfield] = response_vector
	}
	data
}

#
#	Function for extracting first and last measurements for each of the tumors and testing for their correlation
#
xeno.test.endcor = function(
	x,
	responsename="Response",
	tpname="Timepoint",
	idname="Tumor_id",
	method="pearson",
	subgroup="Treatment",
	draw=TRUE, 
	legendposition="topright",
	maintitle="End-point measurements"
){
	check_headers(x, name=c(responsename, tpname, idname), mer=TRUE, dataframe=TRUE)	

	colors = c("black", "green", "blue", "red")
	first_vec = vector(length=0)
	last_vec = vector(length=0)
	groups_vec = vector(length=0)
	
	if(class(x)=="data.frame"){
		uniques = as.vector(unique(x[[idname]]))
		for(i in 1:length(uniques)){
			id = subset(x, x[idname]==uniques[i])
			resp = id[[responsename]]
			resp = resp[!is.na(resp)]
			first = resp[1]
			last = resp[length(resp)]
			first_vec = append(first_vec, first)
			last_vec = append(last_vec, last)
			if(!is.null(subgroup)){
				groups = as.vector(unique(x[[subgroup]]))
				groups = groups[order(groups)]
			}
			groups_vec = append(groups_vec, which(id[[subgroup]][1]==groups))
		}
		
	}else if(class(x)=="mer"){
		uniques = as.vector(unique(x@frame[[idname]]))
		for(i in 1:length(uniques)){
			id = subset(x@frame, x@frame[idname]==uniques[i])
			resp = id[[responsename]]
			resp = resp[!is.na(resp)]
			first = resp[1]
			last = resp[length(resp)]
			first_vec = append(first_vec, first)
			last_vec = append(last_vec, last)
			if(!is.null(subgroup)){
				groups = as.vector(unique(x@frame[[subgroup]]))
				groups = groups[order(groups)]
			}
			groups_vec = append(groups_vec, which(id[[subgroup]][1]==groups))
		}
	}else{
		stop("Invalid x: should be either a data.frame or fitted lme4 (mer) object")
	}
	mat = matrix(ncol=2, nrow=length(uniques))
	rownames(mat) = uniques
	colnames(mat) = c("First","Last")
	mat[,1] = first_vec
	mat[,2] = last_vec
	#print(cor(mat[,1], mat[,2], method=method))
	print(cor.test(mat[,1], mat[,2], method=method))
	if(draw){
		plot.new()
		plot.window(xlim=c(min(mat[,1], na.rm=TRUE),max(mat[,1], na.rm=TRUE)), ylim=c(min(mat[,2], na.rm=TRUE),max(mat[,2], na.rm=TRUE)))
		title(maintitle, xlab="First measured value", ylab="Last measured value")
		box()
		axis(1)
		axis(2)
		points(mat[,1], mat[,2], pch=20, col=colors[groups_vec])
		if(!is.null(legendposition)){
			legend(legendposition, c(paste(subgroup, ":", groups)), col=colors[1:length(groups)], pch=20, bty="n")
		}
	}
	mat
}

#
# Helper function for constructing the model fit formula
#
xeno.formula = function(
	# Was there a target size in the experiment design - if true includes intercept (b1) and offset (b2) in model
	target = TRUE,
	# Should categories be included in the model - if true the covariate Growth is included in the b3 and b4 terms
	cat = TRUE,
	# Response column name
	responsename = "Response",
	# Treatment column name
	treatmentname = "Treatment",
	# Time point column name
	tpname = "Timepoint",
	# Unique tumor ID column names
	idname = "Tumor_id",
	# Boolean values whether to include fixed effects (b1,2,3,4) and/or random effects (u1,2)
	terms = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
){
	if(!terms[5] & !terms[6]){
		stop("There has to be random effects term in the model formulation to fit a mixed-effects model.")
	}
	
	# Response
	frml = paste(responsename, " ~ ", sep="")

	# Intercept
	if(!target | !terms[1]){
		frml = paste(frml, " 0 + ", sep="")
	}
	if(target & terms[1]){
		frml = paste(frml, " 1 + ", sep="")
	}

	# Offset	
	if(target & terms[2]){
		frml = paste(frml, treatmentname, " + ", sep="")
	}
	
	
	# Overall growth
	if(cat & terms[3]){
		frml = paste(frml, tpname, ":Growth + ", sep="")
	}else if(!cat & terms[3]){
		frml = paste(frml, tpname, " + ", sep="")
	}
	
	# Slope effect
	if(cat & terms[4]){
		frml = paste(frml, treatmentname, ":", tpname, ":Growth + ", sep="")
	}else if(!cat & terms[4]){
		frml = paste(frml, treatmentname, ":", tpname, " + ", sep="")
	}
	
	# Random intercept & slope
	if(terms[5] & terms[6]){
		frml = paste(frml, "(1|", idname, ") + (0 + ", tpname, "|", idname, ")", sep="")
	# Just random intercept
	}else if(terms[5] & !terms[6]){
		frml = paste(frml, "(1|", idname, ")", sep="")
	# Just random slope
	}else{
		frml = paste(frml, "(0 + ", tpname, "|", idname, ")", sep="")
	}
	
	as.formula(frml)
}
