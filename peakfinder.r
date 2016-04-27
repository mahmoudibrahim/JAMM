########################################################################
# JAMMv1.0.7rev2 is a peak finder for joint analysis of NGS replicates.
# Copyright (C) 2014-2015  Mahmoud Ibrahim
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact: mahmoud.ibrahim@mdc-berlin.de
########################################################################




# ======================= 
# User-defined variables
# =======================
samplingSeed = 1011414 #the seed makes the random sampling that JAMM does "deterministic". You can change this to your integer seed of choice to reproduce exact same results or randomize it this way: ceiling(rnorm(1,50,300)). 
reportNoClust = "n" #report windows for which clustering failed? Reported windows will be marked by "NoClust" in the peak name (4th column). 
presched = TRUE #Set preschedule in mclapply to TRUE or FALSE
cutoff = NA #To enforce an SNR cutoff ratio for bin enrichment calling, delete NA and enter the number you want.
strict = 1 #To make bin enrichment calling more / less strict, increase / decrease this number.
meanAdjust = "n" #Adjust the initialization mean vector for each window before clustering? If you want to do so, change this to "y".
options(warn = -1, scipen = 1000) #R will not report any warnings (warn = -1), R will not use scientific notation (scipen = 1000).
#=======================> DONE! 




# ================================ 
# Required Libraries check & load
# ================================ 
if ((is.element('mclust', installed.packages()[,1]) == FALSE) || (is.element('signal', installed.packages()[,1]) == FALSE) || (is.element('parallel', installed.packages()[,1]) == FALSE)) {
	stop("R package 'mclust', 'signal' and 'parallel' are required. Please install them!")
}
if ((is.element('Rsamtools', installed.packages()[,1]) == FALSE)) {
  stop("R package'Rsamtools' is required. Please install it!")
}
if ((is.element('GenomicRanges', installed.packages()[,1]) == FALSE)) {
  stop("R package'GenomicRanges' is required. Please install it!")
}
if ((is.element('GenomicAlignments', installed.packages()[,1]) == FALSE)) {
  stop("R package'GenomicAlignments' is required. Please install it!")
}
suppressPackageStartupMessages(library("mclust"))
suppressPackageStartupMessages(library("signal"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("Rsamtools"))
suppressPackageStartupMessages(library("GenomicAlignments"))
suppressPackageStartupMessages(library("GenomicRanges"))
#=======================> DONE! 





# ================= 
# Custom Functions
# =================
makeCountList <- function(actlevels, userlevels, userlevellength, binchromv)
{
  countlist <- list()
  userlevels <- tolower(userlevels)
  userlevels <- gsub("chr","",userlevels)
  userlevels <- gsub("seq","",userlevels)
  if(identical(class(binchromv),"character")) 
  {
    chromv <- tolower(binchromv)
    chromv <- gsub("chr","",chromv)
    chromv <- gsub("seq","",chromv)
    st <- match(userlevels,chromv)
    for(i in 1:length(st))
    {
      if(is.na(st[[i]])) {binchromv[[i]] <- 0} 
      else {binchromv[[i]] <- 1}
    }
  }
  tempactlevels <- tolower(actlevels)
  tempactlevels <- gsub("chr","",tempactlevels)
  tempactlevels <- gsub("seq","",tempactlevels)
  temp <- match(userlevels,tempactlevels)
  temp2 <- match(tempactlevels,userlevels)
  #false if 1)what's in the file is provided and what's provided is in the file (or more)
  if (any(is.na(temp))||any(is.na(temp2)))
  {
    actlevels <- actlevels[which(!is.na(temp2))]
    #what's in the file is provided, but there's more provided, so remove the unnecessary information
    userlevels <- userlevels[which(!is.na(temp))]
    binchromv <- binchromv[which(!is.na(temp))]
    userlevellength <- userlevellength[which(!is.na(temp))]
    temp <- Filter(Negate(function(x) is.na(unlist(x))), temp)
    if(any(is.na(temp))){
      warning("There are less chromosomes in the bam file than the sizes file indicates")
    }
  }
  newbinchromv <- binchromv[order(temp)]
  newuserlevellength <- userlevellength[order(temp)]
  
  actlevels <- actlevels[which(newbinchromv==1)]
  newuserlevellength <- newuserlevellength[which(newbinchromv==1)]
  length(countlist) <- length(actlevels)
  names(countlist) <- actlevels
  for (i in 1:length(actlevels))
  {
    countlist[[actlevels[[i]]]] <- newuserlevellength[[i]]
  }
  return(countlist)
  #mapply(defineLength,actlevels,newuserlevellength)
  #return()
}

defineLength <- function(chromname,chromlength)
{
  #countlist[[chromname]] <<- chromlength
  #return()
  countlist[[chromname]] <- chromlength
  return(countlist)
}


makeVectors <- function(start,end,value)
{
  return(c(rep(value,times=end-start+1)))
}

removeOutOfRangeReads <- function(start,end,seqname,strand,chromlength)
{
  if (start>0 && end<=chromlength) 
  { 
    output <- list(start,end,seqname,strand)
    return (output) 
  }
}


callCountReads <- function(cornum,type,bamfiles,index,reads,frag,chromsize,filelist,chrcount,kpduplicates,RCV)
{
  if (type == "single") {
    if (cornum > 1) {
      counts = mclapply(bamfiles,countreads,index=index,reads=reads,frag=frag,chromsize=chromsize,filelist=filelist,chrcount=chrcount,mc.cores=cornum,mc.preschedule=presched,kpduplicates=kpduplicates,RCV=RCV)
    } else {
      counts = lapply(bamfiles,countreads,index=index,reads=reads,frag=frag,chromsize=chromsize,filelist=filelist,chrcount=chrcount,kpduplicates=kpduplicates,RCV=RCV)
    }
  }  
  if (type == "paired") {
    if (cornum > 1) {
      counts = mclapply(bamfiles,countreadspe,index=index,reads=reads,frag=frag,chromsize=chromsize,filelist=filelist,chrcount=chrcount,mc.cores=cornum,mc.preschedule=presched,RCV=RCV)
    } else {
      counts = lapply(bamfiles,countreadspe,index=index,reads=reads,frag=frag,chromsize=chromsize,filelist=filelist,chrcount=chrcount,RCV=RCV)
    }
  }
  return (counts)
}




#Get per-row Geometric mean (takes list, returns vectors, not lists!)
geomeanL <- function(mat){
	n = length(mat)
	if (n > 1) {
		mult = (mat[[1]])*(mat[[2]])
		if (n > 2) {
			for (i in 3:n) {
				mult = mult*(mat[[i]])
			}
		}
		mult = mult^(1/n)
		mat = mult
		return(mat)
	} else {
		return(mat[[1]])
	}
}



#Get per-row Geometric mean (takes matrix, returns vectors, not lists!)
geomean <- function(mat){
	n = NCOL(mat)
	if (n > 1) {
		mult = (mat[,1])*(mat[,2])
		if (n > 2) {
			for (i in 3:n) {
				mult = mult*(mat[,i])
			}
		}
	mult = mult^(1/n)
	mat = mult
	}
return(mat)
}



#Read in bed(st2) file
parsein = function(bedfile) {
  l = read.table(bedfile, header = FALSE)[[1]]
	l = l + 1
	return(l)
}


#Read in bedpe(st2) file
parseinpe = function(bedfile) {
	l = read.table(bedfile, header = FALSE)
	l = cbind((l[[1]] + 1), l[[2]])
	return(l)
}


#Produces normalized extended read counts (takes output of parsein(), return a vector of floats)
countreads = function(bamfile, index, reads, frag, chromsize, filelist, chrcount, kpduplicates, RCV=ReadChromVector) {
	
  #####
  #print(paste("filelist",filelist))
  #print(paste("bamfile",bamfile))
  #print(paste("chrcount",chrcount))
  #print(paste("frag",frag))
  
  #o ist die Nummer des Chromosoms
  ###########COMMENT FROM MAHMOUD################
  ## In my implementation, o was not the chromosome, it was which replicate. Just pointing it out in case you didn't set it up here on purpose
  ###############################################
	o = which(filelist == bamfile)
	#print(paste("o",o,class(o)))
  
	#if (reads[[o]][length(reads[[o]])] > chromsize) {
	#	message(paste0(chromName, ", Warning: Read alignments do not match chromosome length, Skipped!"))
	#	
	#	if (chrcount == 1) {
	#		message(paste0(chromName, ", ERROR: The first chromosome in the analysis was skipped. I can not calculate normalization factors. You can either delete this chromosome from your chromosome size file or fix the previous warning!"))
	#		quit()
	#	} else {
	#		quit()
	#	}
	#}
  cat("\n", "peakfinder: Counting",bamfile)
  if (o > length(frag))
  {
    FragLength <- frag[length(frag)]
  } else {
    FragLength <- frag[o]
  }
	
  als <- NA
  #o = which(filelist == bamfile)
  bamfile <- BamFile(bamfile)
  param <- ScanBamParam()
  #countlist <<- list()
  countlist <- list()
  #creates a list for all chromosomes to call peaks at
  countlist <- makeCountList(seqnames(seqinfo(bamfile)),chromName,chromSize,RCV)
  chromnames <- names(countlist)
  starting=TRUE
  curnum <- 1
  for (element in chromnames)
  {
    if(curnum == 1) {
      cat("\n","Counting:",length(chromnames),"elements","\n")
      cat(paste0(element,", ",curnum,"; "))
    } else {
      cat(paste0(element,", ",curnum,"; "))
    }
    #print(paste("Counting",element))
    #chromlength <- length(countlist[[element]])
    chromlength <- as.integer(countlist[[element]])
    
    #####is it important to skip chromosome if length doesn't fit?
	###########COMMENT FROM MAHMOUD################
	## I think so, yes. Would be nice to get an error here and skip
	###############################################
    
    als <- readGAlignments(bamfile,index=indexfile,param=ScanBamParam(which=GRanges(element,IRanges(1,as.integer(chromlength)))))
   	###########COMMENT FROM MAHMOUD################
	## what is 536870912? You call the same function twice but with that difference
	###############################################
    alsover <- readGAlignments(bamfile,index=indexfile,param=ScanBamParam(which=GRanges(element,IRanges(as.integer(chromlength),536870912))))
    if(!!length(alsover)){
      #message(paste0(chromName, ", Warning: Chromosome ",element,"has no reads, Skipped!"))
      message(paste0(element, ", Warning: Chromosome ",element,"has no reads, Skipped!"))
      #if (chrcount == 1) {
      if (starting) {
        #message(paste0(chromName, ", ERROR: The first chromosome in the analysis was skipped. I can not calculate normalization factors. You can either delete this chromosome from your chromosome size file or fix the previous warning!"))
        message(paste0(element, ", ERROR: The first chromosome in the analysis was skipped. I can not calculate normalization factors. You can either delete this chromosome from your chromosome size file or fix the previous warning!"))
        quit(status=1)
      } else {
        quit()
      }
    }
    rm(alsover)
    gc()
    if (!length(als))
    { 
      #return(NULL)
      #countlist[[element]] <<- NULL
      countlist[[element]] <- NULL
      next()
    }
    als <- GRanges(als)
    #extends/truncates reads to user specified length and deletes the ones that don't match chromosome size anymore
    als <- resize(als, FragLength)
    #print(paste("FragLength is",FragLength))
    starts<-start(ranges(als));ends<-end(ranges(als));seqnames<-as.character(seqnames(als));strands<-as.character(strand(als))
    als <- mapply(removeOutOfRangeReads,starts,ends,seqnames,strands,chromlength)
    unl <- unlist(als)
    starts <- as.integer(unl[c(TRUE,FALSE,FALSE,FALSE)]);ends <- as.integer(unl[c(FALSE,TRUE,FALSE,FALSE)])
    seqnames <- as.character(unl[c(FALSE,FALSE,TRUE,FALSE)]);strands <- as.character(unl[c(FALSE,FALSE,FALSE,TRUE)])
    als <- GRanges(seqnames=Rle(seqnames),strand=Rle(strands),ranges=IRanges(start=starts,end=ends))
    if (!kpduplicates) {
      als <- unique(als)
    }
    #calculate the number of reads on each position
    #starts <- NA; ends <- NA; seqnames <- NA; strands <- NA
    rm(starts,ends,seqnames,strands)
    gc()
    cov <- coverage(als)
    curcov <- cov[element]
    curcov <- as(curcov, "GRanges")
    values <- score(curcov)
    st <- start(ranges(curcov))
    en <- end(ranges(curcov))
    rm(curcov,cov)
    gc()
    ###COMMENT FROM MAHMOUD########
    ## I couldn't follow this part right before here
    ## Maybe you can explain to me later
    #####################################
    #decompress counts to a chromosome size long vector
    reslist <- mapply(makeVectors,st,en,values,SIMPLIFY=FALSE)
    curvector <- unlist(reslist)
    #normalize the read counts
    rm(reslist,st,en)
    gc()
    #reslist <- NA
    #st <- NA; en <- NA
    #curvector <- curvector/mean(curvector)
    
    
    mCount = mean(curvector)
    
    #####continue here!
    #if (chrcount == 1) {
    if (starting) {
      curvector = curvector/mCount
      #####adapt for different chromosomes, as I do all at once
      write(paste(mCount), file = paste0(out, "/norma.", o, ".info"))
      #write(paste(mCount), file = paste0(out, "/norma.", element, ".info"))
    } else {
      meanCounts = mean(as.numeric(read.table(paste0(out, "/norma.", o, ".info"))[[1]]))
      if ((mCount >  (5*meanCounts)) || (mCount <  (0.2*meanCounts))) {
        mCount = meanCounts
      } else {
        write(paste(mCount), file = paste0(out, "/norma.", o, ".info"), append = TRUE)
      }
      curvector = curvector/mCount
    }
    
    ###COMMENT FROM MAHMOUD########
    ## like in "xcorr", we shouldn't construct this "countlist" variable. 
    ## If I got it right, this has all the counts for the whole genome. If so, it will for sure be huge, regardless of whether you can delete the  previous "reslist", "st"...etc.
    ## The alternative to be able to directly call 
    #####################################
    countlist[[element]] <- curvector
    rm(als,curvector)
    gc()
    gc()
    #als <- NA
    starting <- FALSE
    #curvector <- NA
    curnum=curnum+1
  }
  return (countlist)
  
  
  #####
  #print(paste("Anzahl reads M",length(reads[[o]])))
	#counts = vector(mode = "numeric", length = chromsize)
	#for (j in 1:length(reads[[o]])) {
	#	if ((reads[[o]][j]+frag[o]-1) <= chromsize) {
	#		counts[(reads[[o]][j]):(reads[[o]][j]+frag[o]-1)] =  counts[(reads[[o]][j]):(reads[[o]][j]+frag[o]-1)] + 1
	#	}
	#}
	
  #####
  #counts2 <<- counts
  
	#mCount = mean(counts)
	
	#if (chrcount == 1) {
	#	counts = counts/mCount
	#	write(paste(mCount), file = paste0(out, "/norma.", o, ".info"))
	#} else {
	#	meanCounts = mean(as.numeric(read.table(paste0(out, "/norma.", o, ".info"))[[1]]))
	#	if ((mCount >  (5*meanCounts)) || (mCount <  (0.2*meanCounts))) {
	#		mCount = meanCounts
	#	} else {
	#		write(paste(mCount), file = paste0(out, "/norma.", o, ".info"), append = TRUE)
	#	}
	#	counts = counts/mCount
	#}

	#return(counts)
}



#Produces normalized extended read counts (takes output of parsein(), return a vector of floats)
countreadspe = function(bedfile, reads, chromsize, filelist, chrcount) {
	
	o = which(filelist == bedfile)
	
	
	if (reads[[o]][length(reads[[o]]),2] > chromsize) {
		message(paste0(chromName, ", Warning: Read alignments do not match chromosome length, Skipped!"))
		
		if (chrcount == 1) {
			message(paste0(chromName, ", ERROR: The first chromosome in the analysis was skipped. I can not calculate normalization factors. You can either delete this chromosome from your chromosome size file or fix the previous warning!"))
			quit()
			system("exit 1")
		} else {
			quit()
		}
	}
	
	counts = vector(mode = "numeric", length = chromsize)
	for (j in 1:length(reads[[o]][,1])) {
		counts[(reads[[o]][j,1]):(reads[[o]][j,2])] = counts[(reads[[o]][j,1]):(reads[[o]][j,2])] + 1
	}
	
	mCount = mean(counts)
	
	if (chrcount == 1) {
		counts = counts/mCount
		write(paste(mCount), file = paste0(out, "/norma.", o, ".info"))
	} else {
		meanCounts = mean(as.numeric(read.table(paste0(out, "/norma.", o, ".info"))[[1]]))
		if ((mCount >  (5*meanCounts)) || (mCount <  (0.2*meanCounts))) {
			mCount = meanCounts
		} else {
			write(paste(mCount), file = paste0(out, "/norma.", o, ".info"), append = TRUE)
		}
		counts = counts/mCount
	}

	return(counts)
}



#find enriched bins
pickbins = function(winStart, counts, binSize, chromSize, numdup, C, cutoff, strict, mCs, dCs, bkgd) {

	if ((winStart + binSize) <= chromSize) {
		winEnd = winStart + binSize
	} else {
		winEnd = chromSize
	}
	binSizeTemp = winEnd - winStart
	tempend = winEnd - 1

	#extract subset of the background
	if (bkgd != "None") {
		Cs = counts[[numdup+1]][winStart:tempend]
		mCs = mean(Cs)
		dCs = sd(Cs)
	}

	go = rep(0, numdup)
	for (g in 1:numdup) {
		mS = (mean(counts[[g]][winStart:tempend])) 
		ratio = mS/dCs
		if ((mS > (mCs * strict)) && (ratio > cutoff)) {
			go[g] = 1
		}
	}
	veep = sum(go)
	return(veep)
}




#find enriched wins
pickwins = function(winStart, coffeeshopSud, counts, numdup, startlist, winSize) {
	
	plz = which(startlist == winStart)
	winEnd = coffeeshopSud[plz]
	rWinSize = winEnd - winStart + 1


	if(rWinSize >= winSize) {
		mS = rep(0, numdup)
		for (g in 1:numdup) {
			mS[g] = (mean(counts[[g]][winStart:winEnd])) 
		}
		veep = mean(mS)
	} else {
		veep = FALSE
	}
	
	return(veep)
}



#score windows for fast analysis
scorewindow = function(winStart, coffeeshopSud, numdup, C, bkgd, counts, startlist) {
	
	plz = which(startlist == winStart)
	winEnd = coffeeshopSud[plz]
	
	#will store peak information
	writethis = list()

	rWinSizeTemp = winEnd - winStart + 1
	
	#extract subset of the IP
	Rs = matrix(nrow = rWinSizeTemp, ncol = numdup)
	Rsr = Rs
	for (j in 1:numdup) {
		Rsr[,j] = counts[[j]][winStart:winEnd]
		Rs[,j] = filtfilt(rep(1,80)/80,1,Rsr[,j])
	}
	#extract subset of the background
	if (bkgd != "None") {
		Cs = counts[[numdup+1]]
		Cmin = min(Cs[Cs > 0])
		Cs = Cs[winStart:winEnd]
		Cs = filtfilt(rep(1,80)/80,1,Cs) + Cmin #gets rid of Inf in the fold change
	} else {
		set.seed(samplingSeed)
		Cs = sample(C, rWinSizeTemp, replace = TRUE)
		Cs = filtfilt(rep(1,80)/80,1,Cs)
	}
			
	#start scoring
	signal = (geomean(Rs))
	cairo = (mean(signal)) / (mean(Cs))
	return(cairo)
}


#Initialize MClust clustering parameters
smoothcounts = function(winStart, coffeeshopSud, numdup, counts, startlist) { #helper function1

	plz = which(startlist == winStart)
	winEnd = coffeeshopSud[plz]
	
	#extract subset of the IP
	Rs = matrix(0, nrow = (winEnd - winStart + 1), ncol = numdup)
	for (j in 1:numdup) {
		Rs[,j] = counts[[j]][winStart:winEnd]
	}
	#smooth extended read counts
	for (j in 1:numdup) {
		Rs[,j] = filtfilt(rep(1,80)/80,1,Rs[,j])
	}	
	return(Rs)
}
cluster = function(model, sig, init, clustnummer, noise) { #helper function2
	set.seed(samplingSeed)
	noisy = sample(noise, length(sig[,1]), replace = TRUE)
	clust = me(model, sig+noisy, init)
	bicc =  bic(model, clust$loglik, length(sig[,1]), length(sig[1,]), clustnummer)
	out = list(bicc = bicc, param = clust$parameters)
	return(out)
}
initparam = function(coffeeshopNord, coffeeshopSud, numdup, counts, cornum, clustnummer, modelnames, noise) { #main function
	
	n = length(coffeeshopNord)
	#smooth extended read counts
	if (cornum > 1) {
		sig = mclapply(coffeeshopNord, smoothcounts, coffeeshopSud, numdup, counts, startlist = coffeeshopNord, mc.cores = cornum, mc.preschedule = presched)
	} else {
		sig = lapply(coffeeshopNord, smoothcounts, coffeeshopSud, numdup, counts, startlist = coffeeshopNord)
	}
	sig = do.call(rbind, sig) 

	#kmeans initialization
	set.seed(samplingSeed)
	init = kmeans(sig, clustnummer, nstart = 20)
	init = unmap(init$cluster)
	

	if (cornum > 1) {
		param = mclapply(modelnames, cluster, sig, init, clustnummer, noise, mc.cores = cornum, mc.preschedule = presched)
	} else {
		param = lapply(modelnames, cluster, sig, init, clustnummer, noise)
	}
	
	bicc = vector(mode = "numeric", length = length(modelnames))
	for (i in 1:length(modelnames)) {
		bicc[i] = as.numeric(param[[i]]$bicc)
	}
	bicc = which.max(bicc)
	
	out = list(initparam = param[[bicc]]$param, modelname = modelnames[bicc])
	return(out) 
}




#find peaks
findpeak = function(winStart, coffeeshopSud, numdup, C, param, bkgd, resol, counts, noise, startlist, meanAdjust, clustnummer) {
	
		
	plz = which(startlist == winStart)
	winEnd = coffeeshopSud[plz]
	
	
	#will store peak information
	writethis = list()
	ccx = 1 #default is clustering didNOT work

	
	rWinSizeTemp = winEnd - winStart + 1
	
	#extract subset of the IP
	Rs = matrix(nrow = rWinSizeTemp, ncol = numdup)
	Rsr = Rs
	if (meanAdjust == "y") {
		for (j in 1:numdup) {
			Rsr[,j] = counts[[j]][winStart:winEnd]
			Rs[,j] = filtfilt(rep(1,80)/80,1,Rsr[,j])
			kabel = which.max(param$init$mean[j,])
			param$initparam$mean[j,kabel] = mean(Rs[,j])
		}
	} else {
		for (j in 1:numdup) {
			Rsr[,j] = counts[[j]][winStart:winEnd]
			Rs[,j] = filtfilt(rep(1,80)/80,1,Rsr[,j])
		}
	}
	
	if (resol != "window") {
		
		#clustering (take 1)
		take = 1
		set.seed(samplingSeed)
		noisy = sample(noise, rWinSizeTemp, replace = TRUE)
		clust = em(param$modelname, Rs+noisy, param$initparam)
		clust$classification = map(clust$z)
		if (!((any(diff(clust$classification)) != 0) && (!(any(is.na(clust$classification)))))) { #clustering didn't work, take1
			
			#repeat clustering from scratch, take 2!
			set.seed(samplingSeed)
			init = kmeans(Rs, clustnummer, nstart = 20)
			init = unmap(init$cluster)
			set.seed(samplingSeed)
			noisy = sample(noise, rWinSizeTemp, replace = TRUE)
			clust = me(param$modelname, Rs+noisy, init)
			clust$classification = map(clust$z)
			if ((any(diff(clust$classification)) != 0) && (!(any(is.na(clust$classification))))) {
				ccx = 0 #clustering worked, take2
				take = 2
			}
		} else {ccx = 0}  #clustering worked, take1
		 
		if (ccx != 1) { #clustering worked either in take1 or take2
					
			if (numdup > 1) { #check whether all components replicates agreed on clustering assignments
				cc = vector(mode = "numeric", length = numdup)
				for (g in 1:numdup) {
					cc[g] = which.max(clust$parameters$mean[g,]) #which cluster has the largest mean (this is the peak cluster, hopefully!)
				}
				ccx = sum(diff(cc))
				cluster = cc[1]
				rm(cc)
			
				if ((ccx != 0) && (take == 1)) { #not all replicates agreed? Repeat the clustering with from scratch if not already done!
					set.seed(samplingSeed)
					init = kmeans(Rs, clustnummer, nstart = 20)
					init = unmap(init$cluster)
					set.seed(samplingSeed)
					noisy = sample(noise, rWinSizeTemp, replace = TRUE)
					clust = me(param$modelname, Rs+noisy, init)
					clust$classification = map(clust$z)
					if ((any(diff(clust$classification)) != 0) && (!(any(is.na(clust$classification))))) { #clustering worked? check whether replicates agreed take 3
						cc = vector(mode = "numeric", length = numdup)
						for (g in 1:numdup) {
							cc[g] = which.max(clust$parameters$mean[g,]) #which cluster has the largest mean (this is the peak cluster, hopefully!)
						}
						ccx = sum(diff(cc))
						cluster = cc[1]
						rm(cc)
						take = 3
					}
				}
			} else { #no replicates!
				cluster = which.max(clust$parameters$mean) #which cluster has the largest mean (this is the peak cluster, hopefully!)
			}
		}
		
		if ((ccx != 0) && (reportNoClust=="y")) { resol = "window" } #clustering did not work and windows should be reported
			
		
		if (ccx == 0) { #clustering worked and all replicates agree on the cluster assignments
			
			#extract subset of the background
			if (bkgd != "None") {
				Cs = counts[[numdup+1]][winStart:winEnd]
				Cs = filtfilt(rep(1,80)/80,1,Cs)
			} else {
				set.seed(samplingSeed)
				Cs = sample(C, rWinSizeTemp, replace = TRUE)
				Cs = filtfilt(rep(1,80)/80,1,Cs)
			}
			
			#find region boundaries
			loc = 1:length(clust$classification)
			gmclass = cbind(loc, clust$classification)
			locPeak = gmclass[gmclass[,2] == cluster,,drop=FALSE]
			rStart = locPeak[1] #start position of the region
			rEnd = locPeak[length(locPeak[,1]),1] #end position of the region
		
			#peak resolution check
			if (resol == "region") {
				pSize = rEnd - rStart
				signal = (geomean(Rs[rStart:rEnd,])) 
				signal2 = (signal) - (Cs[rStart:rEnd])
				gm = mean(signal2)
				summit = which.max(geomean(Rsr[rStart:rEnd,])) - 1
				will2k = wilcox.test(signal, Cs[rStart:rEnd])
				
				#Is there signal in the region above background
				if (gm > 0) {
					writethis[[1]] = rStart + winStart - 1
					writethis[[2]] = rEnd + winStart
					#writethis[[3]] = paste0(chromName, ".", rStart+winStart -1)
					writethis[[3]] = paste0(element, ".", rStart+winStart -1)
					writethis[[4]] = "1000"
					writethis[[5]] = "."
					writethis[[6]] = gm
					writethis[[7]] = will2k$p.value
					writethis[[8]] = "-1"
					writethis[[9]] = summit
				}
			} else if (resol == "peak") {
				#find out where separate peaks are
				d = diff(locPeak[,1]) 
				d[length(d)+1] = 0		
				locPeak = cbind(locPeak, d)
				bound1 = which(locPeak[,3] > 1, arr.in=TRUE)
				bound2 = bound1 + 1
				bound = locPeak[sort(c(bound1,bound2))]
				bound = c(rStart, bound, rEnd)
				w = 1
				warum = 0
				while (w < length(bound)) {
					pStart = bound[w] + winStart - 1 
					pEnd = bound[w+1] + winStart
					pSize = pEnd - pStart
					signal = (geomean(Rs[(bound[w]):(bound[w+1]),])) 
					signal2 = (signal) - (Cs[bound[w]:bound[w+1]])
					gm = mean(signal2)
					summit = which.max(geomean(Rsr[(bound[w]):(bound[w+1]),])) - 1
					will2k = wilcox.test(signal, Cs[(bound[w]):(bound[w+1])])
										
					weil = warum * 9
					#Is there signal in the region above background
					if (gm > 0) {
						writethis[[1+weil]] = pStart
						writethis[[2+weil]] = pEnd
						#writethis[[3+weil]] = paste0(chromName, ".", pStart)
						writethis[[3+weil]] = paste0(element, ".", pStart)
						writethis[[4+weil]] = "1000"
						writethis[[5+weil]] = "."
						writethis[[6+weil]] = gm
						writethis[[7+weil]] = will2k$p.value
						writethis[[8+weil]] = "-1"
						writethis[[9+weil]] = summit
					}
					w = w + 2
					warum = warum + 1			
				}
			} #peak resolution check
		} #clustering worked and all replicates agree on clustering assignments?
	} #window resolution check
	
	if (resol == "window") {
	
		#extract subset of the background
		if (bkgd != "None") {
			Cs = counts[[numdup+1]][winStart:winEnd]
			Cs = filtfilt(rep(1,80)/80,1,Cs)
		} else {
			set.seed(samplingSeed)
			Cs = sample(C, rWinSizeTemp, replace = TRUE)
			Cs = filtfilt(rep(1,80)/80,1,Cs)
		}
		
		#calculate scores
		pSize = rWinSizeTemp
		signal = geomean(Rs) 
		signal2 = (signal) - (Cs)
		gm = mean(signal2)
		summit = which.max(geomean(Rsr)) - 1
		will2k = wilcox.test(signal, Cs)
		
		#Is there signal in the region above background
		if (gm > 0) {
			writethis[[1]] = winStart - 1
			writethis[[2]] = winEnd
			#writethis[[3]] = paste0(chromName, ".", winStart -1, ".NoClust")
			writethis[[3]] = paste0(element, ".", winStart -1, ".NoClust")
			writethis[[4]] = "1000"
			writethis[[5]] = "."
			writethis[[6]] = gm
			writethis[[7]] = will2k$p.value
			writethis[[8]] = "-1"
			writethis[[9]] = summit
		}
	} #window reporting

return(writethis)
}


#filter return value of findpeak()
processPeaks = function(peaks) {
	peaks = matrix(unlist(peaks), ncol=9, byrow=TRUE)
	peaks = peaks[peaks[,1] != FALSE,,drop=FALSE]
	peaks = data.frame(peaks)
	return(peaks)
}
#=======================> DONE!









# ========================== 
# Parse-in System Variables
# ==========================
args = commandArgs(trailingOnly = TRUE) # Read Arguments from command line

#Parsing arguments and storing values
for (each.arg in args) {
	#chormosome size file
	if (grepl('-sfile=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				size.file <- arg.split[2]
		} else {
			stop('No genome size file')
		} 
	}
	#bed files directory
	if (grepl('-frag=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				frag <- arg.split[2]
		} else {
			stop('No fragment length given')
		} 
	}	
	#background files directory
	if (grepl('-bkgd=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				bkgd <- arg.split[2]
		} else {
			message("No reads found in one or more replicates!")
			quit()
		} 
	}	
	#background files directory
	if (grepl('-out=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				out <- arg.split[2]
		} else {
			stop('No output directory given')
		} 
	}
	#Cluster number
	if (grepl('-clustnummer=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				clustnummer <- as.numeric(arg.split[2])
		} 
	}
	#resolution
	if (grepl('-resolution=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				resol <- arg.split[2]
		} 
	}
	#processor cores
	if (grepl('-p=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				cornum <- as.numeric(arg.split[2])
		} 
	}
	#minimum window size
	if (grepl('-window=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			winSize <- arg.split[2]
		} 
	}
	#window size
	if (grepl('-bin=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			binsize <- arg.split[2]
		} 
	}
	#type (paired / single)
	if (grepl('-type=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			type <- arg.split[2]
		} 
	}
	#chromosome number
	if (grepl('-chrcount=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			chrcount <- as.numeric(arg.split[2])
		} 
	}
	#bed file names
	if (grepl('-bednames=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				bednames <- arg.split[2]
		} else {
			message("No reads found in one or more replicates!")
			#if (chrcount == 1) {
			#	message(paste0("ERROR: The first chromosome in the analysis has no reads in one or more of your BED files (-s). You can either delete this chromosome from your chromosome size file!"))
			#	quit(status = 1)
			#} else {
			#	quit()
			#}
		} 
	}
	#index files
	if (grepl('-iindex=',each.arg)) {
	  arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
	  if (! is.na(arg.split[2]) ) {
	    iindex <- arg.split[2]
	  } else {
	    stop('No index file for one or multiple bam files given')
	  } 
	}	
	#window enrichment cutoff
	if (grepl('-windowe=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			windowe <- arg.split[2]
		} 
	}
	#initialize
	if (grepl('-initModel=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			initialize <- arg.split[2]
		} 
	}
	#keep duplicates or not
	if (grepl('-uniq=',each.arg)) {
	  arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
	  if (! is.na(arg.split[2]) ) {
	    uniq <- as.character(arg.split[2])
	    if(uniq=="y"){
	      uniq=TRUE
	    } else {
	      uniq=FALSE
	    }
	    print(paste("uniq",uniq))
	  } 
	}
	#number of replicates
	if (grepl('-nreps=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				nreps <- as.numeric(arg.split[2])
		} 
	}
  
	#number of background files
	if (grepl('-nbkgd=',each.arg)) {
	  arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
	  if (! is.na(arg.split[2]) ) {
	    nbkgd <- as.numeric(arg.split[2])
	  } 
	}
	
}



##Parse in variables
chromosomes = read.table(size.file, header=FALSE)
chromName = as.character(chromosomes$V1); #which chromosome
chromSize = as.numeric(chromosomes$V2); #chromosomes size
rm(chromosomes)
ReadChromVector <- chromName
iindex = strsplit(iindex, ",", fixed = TRUE)[[1]]
#print(paste("iindex",iindex))

#if (chrcount == 1) {
#	write(paste(samplingSeed), file = paste0(out, "/seed.info"))
#} else {
#	samplingSeed = as.numeric(read.table(paste0(out, "/seed.info"), header = FALSE)[[1]])
#}
readsFiles = as.list(strsplit(bednames, ",", fixed = TRUE)[[1]])
numdup = length(readsFiles) #number of replicates

#if (numdup != nreps) {
if (numdup != nreps+nbkgd) {
  print(paste("numdup, Anzahl der samplefiles:",numdup,"nreps:",nreps))
	message("No reads found in one or more replicates!")
	quit()
}
#if (bkgd != "None") {
#	readsFiles[[numdup+1]] = bkgd
#}

winSize = as.numeric(winSize)
binSize = as.numeric(binsize)
winSize = binSize * winSize

if (type == "single") {
	frags = as.numeric(strsplit(frag, ",", fixed = TRUE)[[1]])
}
rm(bednames)

if (numdup > 1) {
	modelnames = c("VVV","VEV")
} else {
	modelnames = "V"
}

if (windowe != "auto") {
	windowe = as.numeric(windowe)
}

nothing = FALSE #default is we found some peaks
options(stringsAsFactors = FALSE)
#=======================> DONE!





# ======================= 
# Some preliminary stuff
# =======================
#if (type == "single") {
#	if (cornum > 1) {
#		datain = mclapply(readsFiles, parsein, mc.cores = cornum, mc.preschedule = presched) #read in all bed files (samples and control)
#	} else {
#		datain = lapply(readsFiles, parsein) #read in all bed files (samples and control)
#	}
#}

#if (type == "paired") {
#	if (cornum > 1) {
#		datain = mclapply(readsFiles, parseinpe, mc.cores = cornum, mc.preschedule = presched) #read in all bed files (samples and control)
#	} else {
#		datain = lapply(readsFiles, parseinpe) #read in all bed files (samples and control)
#	}
#}

#minimum peak size (only a recommendation)
minpeak = floor(binSize / 4)

#make bins vector
#bins = seq(from = 1, to = (chromSize - 1), by = binSize)
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
bins = seq2(from = 1, to = (chromSize - 1), by = binSize)
#print(bins[[1]]);print(bins[[20]])
#=======================> DONE!


#####
#BamFile <- BamFile ("/home/arettig/Documents/Files/RG5_mm_Aligned.out.coordSrt.rmdup.uniqNH.chrY.bam")
#ChromSizes <- read.table("/home/arettig/Documents/chrSizes/mm10.txt")
#ReadChromVector <- c("chrY")
#ReadLength <- 100

# =============== 
# Counting Reads
# ===============



areinbkgd <- c()

if(nbkgd!=0)
{
  ibkgd <- readsFiles[(length(readsFiles)-nbkgd+1):length(readsFiles)]
  #print(paste("ibkgd",ibkgd))
  ibkgdindex <- iindex[(length(iindex)-nbkgd+1):length(iindex)]
  #print(paste("length iindex",length(iindex)))
  counts=callCountReads(cornum,type,bamfiles=ibkgd,index=ibkgdindex,reads=datain,frag=frags,chromsize=chromSize,filelist=readsFiles,chrcount=chrcount,kpduplicates=uniq,RCV=ReadChromVector)
  if(!length(counts[[1]])){
    quit()
  }
  for (i in 1:length(counts)) 
  {
    areinbkgd <- unique(c(areinbkgd,names(counts[[i]])))
  }
  #print(paste("areinbkgd",areinbkgd))
  #areinbkgd <- names(counts[[1]])
  ifrgd <- readsFiles[1:(length(readsFiles)-nbkgd)]
  #print(paste("length iindex",length(iindex)))
  ifrgdindex <- iindex[1:(length(iindex)-nbkgd)]
  countsf=callCountReads(cornum,type,bamfiles=ifrgd,index=ifrgdindex,reads=datain,frag=frags,chromsize=chromSize,filelist=readsFiles,chrcount=chrcount,kpduplicates=uniq,RCV=as.character(areinbkgd))
  #countsf = mclapply(ifrgd, countup, indexfile=ifrgdindex, filelist = ifrgd, cornum = cornum, presched = presched,RCV=as.character(areinbkgd))
  counts = c(countsf,counts)
} else {
  counts = callCountReads(cornum,type,bamfiles=readsFiles,index=iindex,reads=datain,frag=frags,chromsize=chromSize,filelist=readsFiles,chrcount=chrcount,kpduplicates=uniq,RCV=ReadChromVector)
}


#if (type == "single") {
#	if (cornum > 1) {
#		counts = mclapply(readsFiles, countreads, index=iindex, reads = datain, frag = frags, chromsize = chromSize, filelist = readsFiles, chrcount = chrcount, mc.cores = cornum, mc.preschedule = presched, kpduplicates=uniq)
#	} else {
#		counts = lapply(readsFiles, countreads, index=iindex, reads = datain, frag = frags, chromsize = chromSize, filelist = readsFiles, chrcount = chrcount, kpduplicates=uniq)
#	}
#}

#if (type == "paired") {
#	if (cornum > 1) {
#		counts = mclapply(readsFiles, countreadspe, index=iindex, reads = datain, chromsize = chromSize, filelist = readsFiles, chrcount = chrcount, mc.cores = cornum, mc.preschedule = presched)
#	} else {
#		counts = lapply(readsFiles, countreadspe, index=iindex, reads = datain, chromsize = chromSize, filelist = readsFiles, chrcount = chrcount)
#	}
#}

if(nbkgd > 1)
{
  chrsinbkgd <- names(counts[[nreps+1]])
  for (i in (nreps+2):numdup)#nreps+1 would be the first background file, the one the others will be added to
  {
    chrsini <- names(counts[[i]])
    #for (j in 1:length(counts[[nreps+1]]))
    for (j in 1:length(counts[[i]]))
    {
      if (chrsini[[j]] %in% chrsinbkgd){
        counts[[nreps+1]][[chrsini[[j]]]] <- counts[[nreps+1]][[chrsini[[j]]]] + counts[[i]][[chrsini[[j]]]]
      } else {
        chrsinbkgd <- c(chrsinbkgd,chrsini[[j]])
        counts[[nreps+1]][[chrsini[[j]]]] <- counts[[i]][[chrsini[[j]]]]
      }
      #counts[[nreps]][[j]][[1]] = datain[[nreps]][[j]][[1]] + datain[[i]][[j]][[1]]
      #counts[[nreps]][[j]][[2]] = datain[[nreps]][[j]][[2]] + datain[[i]][[j]][[2]]
    }
  }
}

#print(names(counts[[3]]))
#print(class(counts))
#print(counts)
#print(names(counts))
#print(length(counts))
#print(length(counts[[1]]))
#counts2 <- counts[[1]]


#counts <- Jamm.prep(BamFile,ChromSizes,ReadChromVector,ReadLength=ReadLength)[["chrY"]]
#counts3 <- counts
#print(length(counts3))
#print(all(counts2==counts3))
#print(length(which(counts2!=counts3)))
#write((which(counts2!=counts3)), file="/home/arettig/Documents/test3.txt", append=TRUE)
#print(length(counts2))
#print(length(counts3))
#print(counts2[299376])
#print(counts3[299376])

#####
#o = which(readsFiles == readsFiles)
#mCount = mean(counts)
#if (length(names(counts)) == 1) {
#  counts = counts/mCount
#  write(paste(mCount), file = paste0(out, "/norma.", o, ".info"))
#} else {
#  meanCounts = mean(as.numeric(read.table(paste0(out, "/norma.", o, ".info"))[[1]]))
#  if ((mCount >  (5*meanCounts)) || (mCount <  (0.2*meanCounts))) {
#    mCount = meanCounts
#  } else {
#    write(paste(mCount), file = paste0(out, "/norma.", o, ".info"), append = TRUE)
#  }
#  counts = counts/mCount
#}


numdup=nreps
rm(datain)

#=======================> DONE!


if (is.null(areinbkgd)){
  chrs <- c()
  for(i in 1:length(counts))
  {
    chrs <- unique(c(chrs,names(counts[[i]])))
  }
} else {
  chrs <- areinbkgd
}

chrcount=1
datapool=counts
counts <- NA
for (element in chrs)
{
  print(paste("Calling peaks of chromosome",element))
  if (chrcount == 1) {
    write(paste(samplingSeed), file = paste0(out, "/seed.info"))
  } else {
    samplingSeed = as.numeric(read.table(paste0(out, "/seed.info"), header = FALSE)[[1]])
  }
  
  counts=list()
  for (i in 1:length(datapool))
  {
    counts[[i]]=datapool[[i]][[element]]
  }
  
  
  
  # ============================ 
  # Estimating Background Model
  # ============================ 
  if (chrcount == 1){ #first chromosome, estimate bkgd (includes SNR cutoff)
    if (is.na(cutoff)) {
      if (bkgd != "None") {
        cutoff = vector(length = numdup)
        sdC = sd(counts[[numdup+1]])
        for (x in 1:numdup) {
          cutoff[x] = (mean(counts[[x]]))/(sdC)
        }
        cutoff = max(cutoff)
        C = NULL
        mCs = NULL
        write(paste(c(cutoff,NA,NA)), file = paste0(out, "/bkgd.info"), append = TRUE)
      } else {
        cutoff = vector(length = numdup)
        mmV = var(geomeanL(counts))
        mmM = mean(geomeanL(counts))
        sigma = log(1+((mmV) / ((mmM)^2)))
        mu = (log(mmM)) - (0.5 * (sigma))
        set.seed(samplingSeed)
        C = rlnorm(100000, mu, sqrt(sigma))
        for (x in 1:numdup) {
          cutoff[x] = (mean(counts[[x]]))/(sd(C))
        }
        cutoff = max(cutoff)
        set.seed(samplingSeed)
        snow = sample(C, binSize*5, replace = TRUE)
        mCs = mean(snow)
        dCs = sd(snow)
        write(paste(c(cutoff,sigma,mu)), file = paste0(out, "/bkgd.info"), append = TRUE)
      }
    }
  } else { #bkgd estiamted from before
    bkgdInfo = read.table(paste0(out, "/bkgd.info"), header = FALSE)
    if (is.na(cutoff)) {
      cutoff = as.numeric(bkgdInfo[[1]][1])
    }
    
    if (bkgd != "None") {
      C = NULL
      mCs = NULL
    } else {
      sigma = as.numeric(bkgdInfo[[1]][2])
      mu = as.numeric(bkgdInfo[[1]][3])
      set.seed(samplingSeed)
      C = rlnorm(100000, mu, sqrt(sigma))
      set.seed(samplingSeed)
      snow = sample(C, binSize*5, replace = TRUE)
      mCs = mean(snow)
      dCs = sd(snow)
    }
  }
  print("estimated background model")
  #=======================> DONE!
  
  
  ###COMMENT FROM MAHMOUD####
  ## Starting here, stuff would be groupped in one function, that just gets called directly inside the readcounting function"..that way peaks are called "on the fly" while getting reads from a file.
  ###########################
  
  # ======================== 
  # Picking Enriched Windows
  # ========================
  if (cornum > 1) {
    coffeeshop = mclapply(bins, pickbins, counts, binSize, chromSize, numdup, C, cutoff, strict, mCs, dCs, bkgd, mc.cores = cornum, mc.preschedule = presched)
  } else {
    coffeeshop = lapply(bins, pickbins, counts, binSize, chromSize, numdup, C, cutoff, strict, mCs, dCs, bkgd)
  }
  coffeeshop = as.numeric(unlist(coffeeshop))
  coffeeshop[coffeeshop != numdup] = 0
  
  if (sum(coffeeshop) != 0) { #Any enriched bins?
    coffeeshop = c(0, diff(coffeeshop))
    coffeeshop = cbind(coffeeshop, bins)
    coffeeshopNord = coffeeshop[coffeeshop[,1] == numdup,,drop=FALSE]
    coffeeshopSud = coffeeshop[coffeeshop[,1] == -numdup,,drop=FALSE]
    if (nrow(coffeeshopNord) != 0) { #Any enriched bins?
      coffeeshopNord = coffeeshopNord[,2]
      coffeeshopSud = coffeeshopSud[,2] - 1
      if (length(coffeeshopSud) < length(coffeeshopNord)) {
        coffeeshopSud = c(coffeeshopSud, chromSize) 
      } else if (length(coffeeshopSud) > length(coffeeshopNord)) {
        coffeeshopNord = c(1, coffeeshopNord)
      }
      if (coffeeshopSud[length(coffeeshopSud)] > chromSize) {
        coffeeshopSud[length(coffeeshopSud)] = chromSize
      }
      
      
      if (cornum > 1) {
        coffeeshop = mclapply(coffeeshopNord, pickwins, coffeeshopSud, counts, numdup, startlist = coffeeshopNord, winSize, mc.cores = cornum, mc.preschedule = presched)
      } else {
        coffeeshop = lapply(coffeeshopNord, pickwins, coffeeshopSud, counts, numdup, startlist = coffeeshopNord, winSize)
      }
      coffeeshop = as.numeric(unlist(coffeeshop))
      coffeeshop = cbind(coffeeshopNord, coffeeshopSud, coffeeshop)
      coffeeshop = coffeeshop[coffeeshop[,3] != FALSE,,drop=FALSE]
      if (nrow(coffeeshop) > 1) {
        coffeeshop = coffeeshop[order(coffeeshop[,3], decreasing = TRUE),]
      }
      rm(bins)
      print("Nach enriched windows")
      #=======================> DONE!
      
      
      
      
      # =================================== 
      # Initializing Clustering Parameters
      # ===================================
      if (length(coffeeshop[,1]) > 0) { #any enriched windows detected?
        if (initialize == "deterministic") {
          yummy = ceiling(length(coffeeshop[,1]) / 1000)
          if (yummy == 0) {
            yummy = 1
          }
        }
        if (initialize == "stochastic") {
          yummy = ceiling(length(coffeeshop[,1]) / 4)
          if (yummy > 20) {
            set.seed(samplingSeed)
            yummy = sample(1:yummy, 20)
          } else if (yummy > 0) {
            yummy = 1:yummy
          } else {
            yummy = 1
          }
        }
        coffeeshopNord = coffeeshop[yummy,1]
        coffeeshopSud = coffeeshop[yummy,2]
        set.seed(samplingSeed)
        noise = rnorm(100000, mean=0, sd=0.1)
        param = initparam(coffeeshopNord, coffeeshopSud, numdup, counts, cornum, clustnummer, modelnames, noise)
        print("Clustering parameters")
        #=======================> DONE!
        
        
        
        
        
        # ========================== 
        # Enriched Window Filtering
        # ==========================
        if (windowe != 1) { #do it only if window fold enrichment filtering is required
          if (cornum > 1) {
            scores = mclapply(coffeeshop[,1], scorewindow, coffeeshop[,2], numdup, C, bkgd, counts, startlist = coffeeshop[,1], mc.cores = cornum, mc.preschedule = presched)
          } else {
            scores = lapply(coffeeshop[,1], scorewindow, coffeeshop[,2], numdup, C, bkgd, counts, startlist = coffeeshop[,1])
          }
          scores = unlist(scores)
          
          if (windowe == "auto") {
            lscores = log(scores)
            if (length(scores) > 0) {
              if (chrcount == 1) {
                cutthisTEMP = ((mean(lscores)) + (sd(lscores)*1))
                write(paste(cutthisTEMP), file = paste0(out, "/bkgd.info"), append = TRUE)
              } else {
                cutthisTEMP = as.numeric(bkgdInfo[[1]][4])
              }
              finalwins = which(lscores > cutthisTEMP)
              cutthisW = min(scores[finalwins])
              coffeeshop = cbind(coffeeshop[finalwins,1], coffeeshop[finalwins,2], coffeeshop[finalwins,3])
            } else {
              if (chrcount == 1) {
                cutthisTEMP = 0
                cutthisW = "Not Applicable, All Windows Analyzed!"
                write(paste(cutthisTEMP), file = paste0(out, "/bkgd.info"), append = TRUE)
              }
            }
          } else {
            cutthisW = windowe
            if (length(scores) > 0) {
              finalwins = which(scores >= windowe)
              coffeeshop = cbind(coffeeshop[finalwins,1], coffeeshop[finalwins,2], coffeeshop[finalwins,3])
            }
          }
        } else { cutthisW = 1 }
        print("Enriched window filtering")
        #=======================> DONE!
        
        
        
        
        
        # ============== 
        # Finding Peaks
        # ==============
        if (nrow(coffeeshop) > 0) { #any enriched windows left after filtering?
          coffeeshop = cbind(coffeeshop[,1], coffeeshop[,2])
          if (cornum > 1) {
            peaks = mclapply(coffeeshop[,1], findpeak, coffeeshop[,2], numdup, C, param, bkgd, resol, counts, noise, startlist = coffeeshop[,1], meanAdjust, clustnummer, mc.cores = cornum, mc.preschedule = presched)
          } else {
            peaks = lapply(coffeeshop[,1], findpeak, coffeeshop[,2], numdup, C, param, bkgd, resol, counts, noise, startlist = coffeeshop[,1], meanAdjust, clustnummer)
          }
          if (!(is.null(peaks))) { #any peaks discovered?
            writethis = processPeaks(peaks)
            print("Finding peaks")
            #=======================> DONE!
            
            
            
            
            
            # =========================
            # Writing Peak Information
            # =========================
          } else { nothing = TRUE } #no peaks
        } else { nothing = TRUE } #no enriched windows left after filtering
      } else { nothing = TRUE; cutthisW = windowe } #no enriched widnows discovered
    } else { nothing = TRUE; cutthisW = windowe } #no enriched bins discovered
  } else { nothing = TRUE; cutthisW = windowe } #no enriched bins discovered
  
  if (isTRUE(nothing)) {
    #file.create(paste0(out, "/", chromName, ".peaks.bed"))
    file.create(paste0(out, "/", element, ".peaks.bed"))
    print(paste0("file is: ",out, "/", element, ".peaks.bed"))
    #write(paste(chromName, minpeak, sep = "	"), file = paste0(out, "/min.peaksize"), append=TRUE)
    write(paste(element, minpeak, sep = "  "), file = paste0(out, "/min.peaksize"), append=TRUE)
    
    if (chrcount == 1) {
      message(paste0("No peaks found! - Window Fold Enrichment: ", cutthisW, " - Seed: ", samplingSeed))
    } else {
      message("No peaks found!")
    }
  } else {
    #write(paste(chromName, writethis$X1, writethis$X2, writethis$X3, writethis$X4, writethis$X5, writethis$X6, writethis$X7, writethis$X8, writethis$X9, minpeak, sep = "	"), file = paste0(out, "/", chromName, ".peaks.bed"), ncolumns = 1)
    write(paste(element, writethis$X1, writethis$X2, writethis$X3, writethis$X4, writethis$X5, writethis$X6, writethis$X7, writethis$X8, writethis$X9, minpeak, sep = "  "), file = paste0(out, "/", element, ".peaks.bed"), ncolumns = 1)
    #write(paste(chromName, minpeak, sep = "	"), file = paste0(out, "/min.peaksize"), append=TRUE)
    write(paste(element, minpeak, sep = "  "), file = paste0(out, "/min.peaksize"), append=TRUE)
    
    
    if (chrcount == 1) {
      message(paste0("Done! - Window Fold Enrichment: ", cutthisW, " - Seed: ", samplingSeed))
    } else {
      message("Done!")
    }
  }
  print("end")
  #=======================> DONE!
  
  
  
  
    ###COMMENT FROM MAHMOUD########
    ## If we can call peaks from with the counting function, we wouldn't need this big for loop over all chromosomes.
    #####################################
  chrcount=chrcount+1
}
