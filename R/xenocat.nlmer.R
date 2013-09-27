####
#
#	NLMER
#	Functions for non-linear mixed-effects models
#
####

#
#	Function for the non-linear version of the categorizing EM
#
xeno.nlmer.EM = function(
	formula,
	data, 
	Model,
	ModelGradient,
	start,
	max.iter=100, 
	loglh.difference=0.01, 
	responsename="Response",
	treatmentname="Treatment",
	tpname="Timepoint", 
	idname="Tumor_id", 
	verbose=FALSE, 
	discriminate=TRUE, 
	randomstart=FALSE,
	return.iterations=FALSE,
	...
){
	require(lme4)
	check_headers(data, name=c(responsename, treatmentname, tpname, idname), mer=FALSE, dataframe=TRUE)	
	
	
	# Predicting tumor volumes for given growth category using fixed effects of the Model
	get.fixed = function(
		fit,
		data,
		Model,
		growth
	){
		if(verbose){print(paste("Predicting",growth,"..."))}
		growth_vector = as.vector(rep(growth, times=length(data[[tpname]])))
		est = as.data.frame(list(Growth=growth_vector))

		est_l = length(est[,1])

		fixest = rbind(fixef(fit))[rep(1, times=est_l),]
		est = as.data.frame(cbind(est, fixest))
		
		for(i in 1:length(colnames(fit@frame))){
			if(!colnames(fit@frame)[i] %in% names(est)){
				est = cbind(est, fit@frame[,i][1:length(data[[tpname]])])
				names(est)[length(names(est))] = colnames(fit@frame)[i]
			}
		}		

		attach(est, warn.conflicts=FALSE)
		fixes = as.vector(as.formula(body(Model)[[2]]))
		detach(est)
		
		fixes = fixes[1:length(data[[tpname]])]
		fixes
	}


	# Expectation step estimation of theta	
	expected.theta = function(fit, data){
		# LINEAR-PREDICTION
		#pred.1 = predict.tumor(fitted_Model=fit, growth=1)
		#pred.0 = predict.tumor(fitted_Model=fit, growth=0)
	
		# NON-LINEAR LME4
		pred_1 = get.fixed(fit=fit, data=data, Model=Model, growth=1)
		pred_0 = get.fixed(fit=fit, data=data, Model=Model, growth=0)

                # If there are missing values, cannot compare to original response although we can do the prediction
                exp_theta = vector(length=length(unique(data[[idname]])))

                for(i in 1:length(unique(data[[idname]]))){
                        
			included_response = (!is.na(data[[responsename]])) & 
			unique(data[[idname]])[i]==data[[idname]]

			included_pred = 
			(unique(data[[idname]])[i]==data[[idname]][!is.na(data[[responsename]])])

			dnorms1 = dnorm(data[[responsename]][included_response], pred_1[included_pred], sd=sd(data[[responsename]][!is.na(data[[responsename]])]))
			dnorms0 = dnorm(data[[responsename]][included_response], pred_0[included_pred], sd=sd(data[[responsename]][!is.na(data[[responsename]])]))

			logdnorms1 = log(dnorms1)
			logdnorms0 = log(dnorms0)

			logtrans1 = logdnorms1 - max(logdnorms1)   
			logtrans0 = logdnorms0 - max(logdnorms0)

			translogprod1 = sum(logtrans1)
			translogprod0 = sum(logtrans0)

			exp_theta[i] = exp(translogprod1) / (exp(translogprod0) + exp(translogprod1))

                }

                # Replicating theta to be same for each of the measurement values within a specific tumor
                exp_theta = rep(exp_theta, each=length(unique(data[[tpname]])), times=1)
                exp_theta

        }

	# Whether we're using random start or prior that all tumors are growing
	if(randomstart){
		theta=rbinom(length(unique(data[[idname]])), 1, 0.5)
		updated_data = data.frame(data, Growth=theta)	
	}else{
		updated_data = data.frame(data, Growth=1)	
	}
	
	fit = nlmer(formula, data=updated_data, start=start, ...)
	
	if(verbose){ print("First fit complete:")
		print(fit)
	}

	likelihood = logLik(fit)
	
	iterations = list()
	datas = list()
	iter=1
	old.likelihood = likelihood - 1
	while ( abs(likelihood-old.likelihood)> loglh.difference & iter < max.iter){
		if(iter%%1==0){ print(paste("Iter:",iter, "of", max.iter, "/ Current logLik:",likelihood))} 

		iterations[[iter]] = fit
		datas[[iter]] = updated_data
		
		## Expectation step
		exp_theta = expected.theta(fit=fit, data=updated_data)
		# Discriminating expected theta to be between 0 and 1
		if(discriminate){exp_theta = (exp_theta>0.5)*1}

		updated_data["Growth"] = exp_theta
		
		## Maximization step, where we pretend that exp_theta is theta
		fit = nlmer(formula, data=updated_data, start=start, ...)

		if(verbose){
			print(paste("Iteration",iter))
			print(fit)
			#print(xeno.test.cat(fit, responsename=responsename, treatmentname=treatmentname,tpname=tpname, idname=idname))
		}
		old.likelihood = likelihood
		likelihood = logLik(fit)
		
		iter = iter + 1
	}
	
	if(verbose) print(paste("EM finished - iterations:",iter,", logLik in end:",logLik(fit)))
	
	if(return.iterations){
		list(fit, updated_data, iterations)
	}else{
		updated_data
	}
}



#
#	Drawing functions (nlmer): drawing fixed effect fit
#
xeno.nlmer.draw.fixed = function(
	fit,
	orig_data,
	Model,
	responsename="Response",
	treatmentname="Treatment",
	idname="Tumor_id",
	tpname="Timepoint",
	maintitle="Fixed effect fit",
	ymax=NULL,
	draw_orig=TRUE
	){
	require(lme4)
	check_headers(fit, name=c(responsename, treatmentname, tpname, idname), mer=TRUE, dataframe=FALSE)	
	growthname="Growth"

	get.fixed = function(
		fit,
		orig_data,
		Model,
		growth
	){
		if(length(growth)==1){
			growth_vector = as.vector(rep(growth, times=length(orig_data[[tpname]])))
		}else{
			growth_vector = growth
		}

		est = as.data.frame(list(Growth=growth_vector))

		est_l = length(est[,1])

		fixest = rbind(fixef(fit))[rep(1, times=est_l),]
		est = as.data.frame(cbind(est, fixest))

		for(i in 1:length(colnames(fit@frame))){
			if(!colnames(fit@frame)[i] %in% names(est)){
				est = cbind(est, fit@frame[,i])
				names(est)[length(names(est))] = colnames(fit@frame)[i]
			}
		}		

		est = est[1:length(orig_data[[tpname]][!is.na(orig_data[[responsename]])]),]

		attach(est, warn.conflicts=FALSE)
		fixes = as.vector(as.formula(body(Model)[[2]]))
		detach(est)

		fixes = fixes[1:length(orig_data[[tpname]])]
		fixes = fixes[!is.na(fixes)]
		fixes
	}
	growthname="Growth"

	Time = orig_data[[tpname]][!is.na(orig_data[[responsename]])]
	IDs = orig_data[[idname]][!is.na(orig_data[[responsename]])]
	Treatment = orig_data[[treatmentname]][!is.na(orig_data[[responsename]])]

	if(growthname %in% names(fit@frame)){
		Growth = fit@frame[[growthname]][1:length(orig_data[[tpname]][!is.na(orig_data[[responsename]])])]
	}else{
		Growth = rep(0, times=length(orig_data[[tpname]][!is.na(orig_data[[responsename]])]))
	}
	draw_data = data.frame(Treatment, Time, IDs)

	individuals = unique(IDs)

	colours = c("black", "green")

	Fixed_fit = get.fixed(fit=fit, orig_data=orig_data, Model=Model, growth=Growth)
	draw_data = data.frame(draw_data, Fixed=Fixed_fit)


	plot.new()
	if(is.null(ymax)){
		plot.window(xlim=c(min(Time,na.rm=TRUE), max(Time,na.rm=TRUE)), ylim=c(min(c(Fixed_fit,orig_data[[responsename]]),na.rm=TRUE), max(c(Fixed_fit,orig_data[[responsename]]),na.rm=TRUE)))
	}else{
		plot.window(xlim=c(min(Time,na.rm=TRUE), max(Time,na.rm=TRUE)), ylim=c(min(c(Fixed_fit,orig_data[[responsename]]),na.rm=TRUE), ymax))
	}
	title(maintitle, xlab=tpname, ylab=responsename)
	axis(1)
	axis(2)
	if(draw_orig){
		for(i in 1:length(individuals)){
			draw = subset(orig_data, orig_data[[idname]]==individuals[i])
			tps = unique(draw[[tpname]])
			points(tps, draw[[responsename]], type="l", col="grey", lwd=1)
		}
	}
	for(i in 1:length(individuals)){
		draw = subset(draw_data, draw_data[["IDs"]]==individuals[i])

		if(!treatmentname %in% names(fit@frame)){
			treatmentgroup = 0
		}else{
			treatmentgroup = draw[[treatmentname]][1]
		}
		if(growthname %in% names(draw)){
			if(draw[[growthname]][1]==0){
				points(draw[["Time"]][1:(length(unique(draw[["Time"]])))], draw[["Fixed"]][1:(length(unique(draw[["Time"]])))], type="l", lwd=2, lty=3, col=colours[treatmentgroup+1], pch=21)
			}else{
				points(draw[["Time"]][1:(length(unique(draw[["Time"]])))], draw[["Fixed"]][1:(length(unique(draw[["Time"]])))], type="l", lwd=2, lty="dashed", col=colours[treatmentgroup+1], pch=22)
			}
		}else{
			points(draw[["Time"]][1:(length(unique(draw[["Time"]])))], draw[["Fixed"]][1:(length(unique(draw[["Time"]])))], type="l", lwd=2, col=colours[treatmentgroup+1], pch=21)
		}
	}
	
}

#
#	Drawing functions (nlmer): drawing fixed + random effect (full) fit
#
xeno.nlmer.draw.fit = function(
	fit,
	orig_data,
	responsename="Response",
	treatmentname="Treatment",
	idname="Tumor_id",
	tpname="Timepoint",
	draw_orig=TRUE,
	maintitle="Full fit",
	ymax=NULL,
	per_individual=FALSE
	){	
	require(lme4)
	check_headers(fit, name=c(responsename, treatmentname, tpname, idname), mer=TRUE, dataframe=FALSE)	
	growthname="Growth"

	Whole_fit = fit@eta
	Time = orig_data[[tpname]][!is.na(orig_data[[responsename]])]
	IDs = orig_data[[idname]][!is.na(orig_data[[responsename]])]
	Treatment = orig_data[[treatmentname]][!is.na(orig_data[[responsename]])]
	draw_data = data.frame(Whole_fit, Treatment, Time, IDs)
	individuals = unique(IDs)

	if(growthname %in% names(fit@frame)){
		draw_data = data.frame(draw_data, Growth=fit@frame[[growthname]][1:length(orig_data[[tpname]])][!is.na(orig_data[[responsename]])])
	}else{
		draw_data = data.frame(draw_data, Growth=0)
	}

	colours = c("black", "green")

	if(!per_individual){
		plot.new()
		if(is.null(ymax)){
			plot.window(xlim=c(min(Time,na.rm=TRUE), max(Time,na.rm=TRUE)), ylim=c(min(c(Whole_fit,orig_data[[responsename]]),na.rm=TRUE), max(c(Whole_fit,orig_data[[responsename]]),na.rm=TRUE)))
		}else{
			plot.window(xlim=c(min(Time,na.rm=TRUE), max(Time,na.rm=TRUE)), ylim=c(min(c(Whole_fit,orig_data[[responsename]]),na.rm=TRUE), ymax))
		}
		title(maintitle, xlab=tpname, ylab=responsename)
		axis(1)
		axis(2)
		if(draw_orig){
			for(i in 1:length(individuals)){
				draw = subset(orig_data, orig_data[[idname]]==individuals[i])
				tps = unique(draw[[tpname]])
				points(tps, draw[[responsename]], type="l", col="grey", lwd=1)
			}
		}
		for(i in 1:length(individuals)){
			draw = subset(draw_data, draw_data[["IDs"]]==individuals[i])

			if(!treatmentname %in% names(fit@frame)){
				treatmentgroup = 0
			}else{
				treatmentgroup = draw[[treatmentname]][1]
			}
			if(growthname %in% names(draw)){
				if(draw[[growthname]][1]==0){
					points(draw[["Time"]][1:(length(unique(draw[["Time"]])))], draw[["Whole_fit"]][1:(length(unique(draw[["Time"]])))], type="l", lwd=2, lty=3, col=colours[treatmentgroup+1], pch=21)
				}else{
					points(draw[["Time"]][1:(length(unique(draw[["Time"]])))], draw[["Whole_fit"]][1:(length(unique(draw[["Time"]])))], type="l", lwd=2, lty="dashed", col=colours[treatmentgroup+1], pch=22)
				}
			}else{
				points(draw[["Time"]][1:(length(unique(draw[["Time"]])))], draw[["Whole_fit"]][1:(length(unique(draw[["Time"]])))], type="l", lwd=2, col=colours[treatmentgroup+1], pch=21)
			}

		}
	}else{

		for(i in 1:length(individuals)){
			plot.new()
			if(is.null(ymax)){
				plot.window(xlim=c(min(Time,na.rm=TRUE), max(Time,na.rm=TRUE)), ylim=c(min(c(Whole_fit,orig_data[[responsename]]),na.rm=TRUE), max(c(Whole_fit,orig_data[[responsename]]),na.rm=TRUE)))
			}else{
				plot.window(xlim=c(min(Time,na.rm=TRUE), max(Time,na.rm=TRUE)), ylim=c(min(c(Whole_fit,orig_data[[responsename]]),na.rm=TRUE), ymax))
			}
			box()
			axis(1)
			axis(2)

			draw = subset(orig_data, orig_data[[idname]]==individuals[i])
			tps = unique(draw[[tpname]])
			points(tps, draw[[responsename]], type="b", col="grey", lwd=1)
			draw = subset(draw_data, draw_data[["IDs"]]==individuals[i])
			points(draw[["Time"]], draw[["Whole_fit"]], type="b", lwd=1, col="black")				

			title(paste("Individual",draw[["IDs"]][1]), xlab="Time", ylab="Response")
			legend("topleft", lwd=c(1,1), lty=c(1,1), col=c("grey", "black"), legend=c("orig_data", "Full fit"))
		}
	}
}

