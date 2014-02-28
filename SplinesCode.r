require(gss)
require(pracma)


##permTestNew function computes the 1000 different permutations for class

permTestNew<-function(dat,t=0,B=1000){
  counts = dat$resp;
  time   = dat$posi;
  statusP= dat$cl;
  mat = array(dat$cl,dim=c(length(data$cl),B))
  tmp = strsplit(rownames(dat),split="")
  #tmp2 = sapply(tmp,function(i){
    #x = i[[1]][1]
      #for(j in 2:4){ x=paste(x,i[j],sep="")}
    #return(x)})
  tmp2=dat$inde
  
  k = match(unique(tmp2),tmp2)
  cl  = statusP[k+3];
  
  for(l in 1:B){
    clp = sample(cl,length(cl));

    for(j in 1:(length(k)-1)){
      statusP[k[j]:k[j+1]]=clp[j]
    }
    statusP[k[j+1]:length(statusP)]=clp[j+1]
    
    mat[,l] = statusP;
  }
  dat$statusp = mat;
  
  return(dat)
}


####SSRegionFinder uses SSANOVA to find regions of difference. 
###Function takes in data with response, class, Individual ID, position, and perm # for to perform
newSSobject<-function(class,position,response){
	obj = data.frame(class = class, position = position, response = response)
	return(obj)
}

extractRegions<- function(interval_index,sign="",control){
	verbose = control$verbose

	# stopping conditions
	if(length(interval_index) <= 0 ){ 
		if(verbose) show(sprintf("There are no %s regions.",sign))
		return()
	}

	# initialize values
	diff_list = list()
	
	n_diff_intervals = sum( diff( interval_index ) > 1 ) < 1 

	# There are only one differing time points
	if(n_diff_intervals){
		diff_list[[1]] = interval_index
	} else{

	# There are multiple differing time points
	# This code should be sped up.
		ind = which(diff(interval_index) > 1)
		i=1
		while(i <= length(ind)){
			diff_list[[i]]=interval_index[1:ind[i]]
			interval_index=interval_index[-(1:ind[i])]
			i=i+1;
		}
		i=length(diff_list)+1
		diff_list[[i]]=interval_index
	}
	return(diff_list)
}

ssControl <- function(verbose=TRUE,seed=1){
	list(verbose=verbose,seed=seed)
}

SSRegionFinder2<- function(obj,formula,terms,offset = 0,permMat,B = 1000,control=ssControl()){
	# initialize
	position =obj$position
	verbose  =control$verbose
	set.seed(control$seed)
	# fit
	fit   = ssanova(formula, data = obj)
	x_seq = seq(min(position), max(position), by=1)
	df    = data.frame(position = x_seq, class = factor(1))
	pred  = predict(fit, df, se = TRUE, include = terms)

	upper_interval = 2*pred$fit + (2 * 1.645 * pred$se)
	lower_interval = 2*pred$fit - (2 * 1.645 * pred$se)

	positive_regions = which( upper_interval < 0 ) # Which intervals are less than zero. 
	negative_regions = which( lower_interval > 0 ) # Which lower intervals are greater than zero.

	posR = extractRegions(positive_regions,"positive",control)
	negR = extractRegions(negative_regions,"negative",control)

	regions = c(posR, negR)
	nregs = length(regions)

	if(nregs > 0){
		# initializing the permutations
		permutationMatrix = permMat #### THIS WILL need to change to be a function
		B = dim(permutationMatrix)[2]
		permuteData = obj
		
		starts = sapply(regions,min)
		ends   = sapply(regions,max)

		areaPSp=matrix(NA, B, nregs)
		# run permutations
		for (i in 1:B){
			permuteData$class = factor(permutationMatrix[,i])
			fit_i             = ssanova(formula, data=permuteData)
			pred_i            = predict(fit_i, df, se = T, include = terms)

			# This next part can be simplified
			areas_i = sapply(1:nregs,function(j){
				mins = starts[j]  
				maxs = ends[j]    
				trapz(x = x_seq[starts[j]:ends[j]], y = abs(2*pred_i$fit[starts[j]:ends[j]]))
			})
			areaPSp[i,] = areas_i
			if( (verbose & i%%50==0) ) show(i);
		}

		# Write the summary table now
		areas = sapply(1:length(regions),function(i){
			mins = starts[i]
			maxs = ends[i]
			trapz(x=x_seq[mins:maxs], y=abs(2*pred$fit[mins:maxs]))
		})
		
		results = matrix(NA,length(regions),4)
		colnames(results) = c("Start","End","Area","P-Value")
		results[,1] = starts - offset
		results[,2] = ends - offset
		results[,3] = areas
		# this next step can be simplified to rowMeans/colMeans (check which)
		for (i in 1:nregs){
			results[i,4]=mean(results[i,3]<areaPSp[,i]) # Is this less than or less than and equal to
		}
		res = list(result = results, permutated_areas = areaPSp,call=match.call())
		return(res)
	} else{
		return("No interesting regions found.")
	}
}