########################################################################
# JAMMv1.0.7rev4 is a peak finder for joint analysis of NGS replicates.
# Copyright (C) 2014-2016  Mahmoud Ibrahim
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
presched = TRUE #Set preschedule in mclapply to TRUE or FALSE
options(warn = -1, scipen = 1000) #R will not report any warnings (warn = -1), R will not use scientific notation (scipen = 1000).
#=======================> DONE! 








# ========================= 
# Required Libraries check 
# ========================= 
if ((is.element('parallel', installed.packages()[,1]) == FALSE)) {
	stop("R package'parallel' is required. Please install it!")
}
suppressPackageStartupMessages(library("parallel"))
#=======================> DONE! 




# ============================================== 
# Parsing Arguments (source: phantom SPP script)
# ============================================== 
args = commandArgs(trailingOnly = TRUE) # Read Arguments from command line


#Set arguments to default values
ibed = NA # input bed file
sFile = NA # chromosome size
storeFile = NA # file to store result

for (each.arg in args) {
	if (grepl('^-s=',each.arg)) {			
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			sFile <- arg.split[2] 
		} else {
			stop('No chromosome size file')
		}
		###
		chromosomes = read.table(sFile, header=FALSE)
		chromName = chromosomes$V1 #which chromosome
		chromSize = as.numeric(chromosomes$V2) #chromosome size
		rm(chromosomes)
		###
	}

	#input bed file
	if (grepl('^-ibed=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			ibed <- arg.split[2]
		} else {
			message(paste0(chromName, ": No reads found in one or more replicates!"))
			quit()
		}
	}
	#bakcground files directory
	if (grepl('-bkgd=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				bkgd <- arg.split[2]
		} else {
			stop('No background file')
		} 
	}	
	if (grepl('^-rl=',each.arg)) {			
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			rl <- arg.split[2] 
		} else {
			stop('Read length missing')
		}
	}
	if (grepl('^-d=',each.arg)) {			
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			storeFile <- arg.split[2] 
		} else {
			stop('No file to store result')
		}
	}
	if (grepl('^-p=',each.arg)) {			
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			cornum <- as.numeric(arg.split[2]) 
		} else {
			stop('No number of cores given')
		}
	}
	#number of replicates
	if (grepl('-nreps=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			nreps <- as.numeric(arg.split[2])
		} 
	}
}

#Read in variables
ibed = strsplit(ibed, ",", fixed = TRUE)[[1]]
numdup = length(ibed) #number of replicates
if (bkgd != "None") {
	nreps = nreps + 1
}
if (numdup != nreps) {
	message(paste0(chromName, ": No reads found in one or more replicates!"))
	quit()
}
rl = as.numeric(strsplit(rl, ",", fixed = TRUE)[[1]])
#=======================> DONE! 


counting = function(readst, br){
	if (max(readst) > max(br)) {
		return(rep(0, length = length(br)))
	} else {
		return(hist(readst, breaks = br, plot = FALSE)$counts)
	}
}


countingreads = function(bedfile, readlen, filelist, gv) {
	reads = read.table(bedfile, header = FALSE)
	reads = cbind(as.character(reads[[1]]), as.character(reads[[2]]))
	o = which(filelist == bedfile)
	tri = list(as.numeric(reads[reads[,2] == "+",,drop = FALSE][,1]), (as.numeric(reads[reads[,2] == "-",,drop = FALSE][,1])) + readlen[o] - 1)
	readcounts = mclapply(tri, counting, br = gv, mc.cores = cornum, mc.preschedule = presched)
	return(readcounts)
}

xc = function(countlist) {
	if (sum(unlist(countlist)) != 0) {
		crossCorrelation = ccf(countlist[[2]], countlist[[1]], plot = FALSE); #xcorr
		crossCorrelation$lag = crossCorrelation$lag * 20; #correct lag for counts window
		maxCorr = which.max(crossCorrelation$acf);
		maxCorr = abs(crossCorrelation$lag[maxCorr]);
		return(maxCorr)
	} else {
		return(0)
	}
}


# ================== 
# Read Start Count
# ==================
#make the window vector (a vector with start positions of chromosome windows)
genomevec = seq(0, chromSize, by = 20);
if (max(genomevec) < chromSize) {
	genomevec = c(genomevec, chromSize);
}
datain = mclapply(ibed, countingreads, readlen = rl, filelist = ibed, gv = genomevec, mc.cores = cornum, mc.preschedule = presched)
if (sum(unlist(datain)) == 0) {
	message(paste0(chromName, ", Warning: Read alignments do not match chromosome length, Skipped!"))
	quit()
}
#=======================> DONE! 


# ================== 
# Cross Correlation
# ==================
xcorr = unlist(mclapply(datain, xc, mc.cores = cornum, mc.preschedule = presched))
#=======================> DONE! 


# ===================== 
# Write result to File
# =====================
for (i in 1:length(xcorr)) {
	filename = strsplit(ibed[i], "/", fixed = TRUE)[[1]]
	filename = filename[length(filename)]
	filename = strsplit(filename, ".", fixed = TRUE)[[1]]
	filename = filename[3]
	if ((xcorr[i] <= 500) && (xcorr[i] >= 50)) { #only write this to file if xcorr value was plausible
		message(paste0(chromName, ", ", filename, ": Ok!"))
		write(paste(chromName, xcorr[i], sep = "\t"), file = paste0(storeFile, "/xc.", filename, ".tab"), append = TRUE)
	} else {
		if (xcorr[i]) {
		} else {
			message(paste0(chromName, ", ", filename, ": Value Not Used!"))
		}
	}
}
#=======================> DONE!
