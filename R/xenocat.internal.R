####
#
#	XenoCat
#	Internal functions not intended for ordinary users
#
####

#
# Saved as a .Rdata to avoid namespace issues in newer R versions (2.14 or so)
#

#extractstr = function(
#	fit,
#	type="coefs"
#){
#	require(lme4, quietly = TRUE)
#	fit = asS4(fit)
#	if(!class(fit)=="mer"){
#		stop("Invalid fit-parameter for extract_str-function, should be lme4 mer-object")
#	}
#	
#	if(type=="coefs"){
#		ret = summary(fit)@coefs
#	}else if(type=="REmat"){
#		ret = summary(fit)@REmat
#	}else if(type=="sigma"){
#		ret = summary(fit)@sigma
#	}else{
#		stop("Invalid type-parameter for extract_str-function")
#	}
#	return(ret)
#}

extract_str = function(
	fit, 
	type="coefs"
){
	data(extractstr)
	return(extractstr(fit=fit, type=type))
}


extract_progress = function(
	id
){
	progress = try(read.table(paste("progress_",id,".txt",sep=""), header=FALSE), silent=TRUE)
	if(!class(progress)=="try-error"){
		progress = as.character(progress[1,1])
	}else{
		progress = "Error retrieving current progress or the process has not yet started."
	}
	progress
}


save_progress = function(
	progress,
	id,
	type="write"
){
	if(type=="write"){
		prog = try(write.table(progress,file=paste("progress_",id,".txt",sep=""), row.names=FALSE, col.names=FALSE), silent=TRUE)
		if(class(prog)=="try-error"){
			print(paste("Error writing current progress to a text file (id: ",id,") (progress: ", progress, ")", sep=""))
		}
	}else if(type=="print"){
		print(progress)
	}else{
		stop("Invalid type-parameter for save_progress-function")
	}	
}

check_headers = function(
	x,
	name,
	mer=TRUE,
	dataframe=TRUE
){
	for(i in 1:length(name)){
		if(class(x)=="mer" & mer){
			for(i in 1:length(name)){
				if(!name[i] %in% names(x@frame)){
					stop(paste("Could not find column header (", name[i], ") from the input mer-object frame. Please check that column headers are specified correctly in the function call (parameters ending in -name)"))
				}
			}
		}else if(class(x)=="data.frame" & dataframe){
			for(i in 1:length(name)){
				if(!name[i] %in% names(x)){
					stop(paste("Could not find column header (", name[i], ") from the input data.frame. Please check that column headers are specified correctly in the function call (parameters ending in -name)"))
				}
			}
		}else{
			if(mer & !dataframe){
				stop(paste("Invalid object class for the input object (",class(x), "), should be data.frame"))
			}else if(!mer & dataframe){
				stop(paste("Invalid object class for the input object (",class(x), "), should be mer"))
			}else if(mer & dataframe){
				stop(paste("Invalid object class for the input object (",class(x), "), should be mer or data.frame"))
			}else{
				stop("Invalid object class or check_headers-function call")
			}
		}
	}
}

