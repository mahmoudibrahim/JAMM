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


if(is.na(bednamesS)){
  stop("No bed file names")
}



# ============================================== 
# Parsing Arguments (source: phantom SPP script)
# ============================================== 
#args = commandArgs(trailingOnly = TRUE) # Read Arguments from command line
#
#
##Parsing arguments and storing values
#for (each.arg in args) {
#  #bed file names
#  if (grepl('^-filelist=',each.arg)) {
#    arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
#    if (! is.na(arg.split[2]) ) {
#      bednames <- arg.split[2]
#    } else {
#      stop('No bed file names')
#    } 
#  }
#  if (grepl('^-dest=',each.arg)) {
#    arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
#    if (! is.na(arg.split[2]) ) {
#      destfile <- arg.split[2]
#    } else {
#      stop()
#    } 
#  }
#  if (grepl('^-abovezero=',each.arg)) {
#    arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
#    if (! is.na(arg.split[2]) ) {
#      abovezerofile <- arg.split[2]
#    } else {
#      stop()
#    } 
#  }
#}
#options(stringsAsFactors = FALSE)
#=======================> DONE! 



# ================= 
# Custom Functions
# =================
helpfilter <- function(chr,start,end,name,score,strand,signal,pvalue,qvalue,summit,minpeak,geom,size){
  if(as.numeric(size)>=as.numeric(minpeak)){
    if(as.numeric(signal)>=as.numeric(geom)){
      summit = gsub("\r","",summit)
      summit = gsub("\n","",summit)
      return(list(chr=chr,start=start,end=end,name=name,score=score,strand=strand,signal=as.numeric(signal),pvalue=pvalue,qvalue=qvalue,summit=as.numeric(summit)))
    }
    else {
      return(list(chr=NA,start=NA,end=NA,name=NA,score=NA,strand=NA,signal=NA,pvalue=NA,qvalue=NA,summit=NA))
    }
  }
  else {
    return(list(chr=NA,start=NA,end=NA,name=NA,score=NA,strand=NA,signal=NA,pvalue=NA,qvalue=NA,summit=NA))
  }
}

writeout = function(bedfile) {
  
  writethis = read.table(bedfile, header=FALSE)
  #print(head(writethis))
#  abovezero<-as.matrix(writethis)
  abovezero<-writethis
#  abovezero[7,]<-as.numeric(abovezero[7,])
  abovezero[,7]<-as.numeric(abovezero[,7])
#  if(!!length(which(abovezero[7,]<=0))){
  if(!!length(which(abovezero[,7]<=0))){
    abovezero=abovezero[-which(abovezero[,7]<=0),]
  }
  #print(head(abovezero))
#  abovezero=abovezero[ ,order( -unlist(abovezero[7,])) ]
  abovezero=abovezero[ order( -unlist(abovezero[,7])), ]
  #print(head(abovezero))
#  write(paste(abovezero[1,], abovezero[2,], abovezero[3,], abovezero[4,], abovezero[5,], abovezero[6,], abovezero[7,], abovezero[8,], abovezero[9,], abovezero[10,], sep = "\t"), file = abovezerofile, ncolumns = 1)
  write(paste(abovezero[,1], abovezero[,2], abovezero[,3], abovezero[,4], abovezero[,5], abovezero[,6], abovezero[,7], abovezero[,8], abovezero[,9], abovezero[,10], sep = "\t"), file = abovezerofileS, ncolumns = 1)
  
  sizes = writethis$V3 - writethis$V2
  #print(head(sizes))
  writethis = mapply(helpfilter,writethis$V1,writethis$V2,writethis$V3,writethis$V4,writethis$V5,writethis$V6,writethis$V7,writethis$V8,writethis$V9,writethis$V10,writethis$V11,writethis$V12,sizes)
  dimnames(writethis)[[2]]=NULL
  delcol=c()
  for (i in 1:ncol(writethis)){
    if (all(is.na(writethis[,i]))){
      delcol<-c(delcol,i)
    }
  }
  if(!!length(delcol)){
    writethis<-writethis[,-delcol]
  }
  writethis=writethis[ ,order( -unlist(writethis[7,])) ]
  write(paste(writethis[1,], writethis[2,], writethis[3,], writethis[4,], writethis[5,], writethis[6,], writethis[7,], writethis[8,], writethis[9,], writethis[10,], sep = "\t"), file = destfileS, ncolumns = 1)
  #writethis<-writethis[,-which(writethis[7,]<=0)]
  #write(paste(writethis[1,], writethis[2,], writethis[3,], writethis[4,], writethis[5,], writethis[6,], writethis[7,], writethis[8,], writethis[9,], writethis[10,], sep = "\t"), file = abovezerofile, ncolumns = 1)
  return(bedfile)
}
#=======================> DONE! 




bedfiles = writeout(bednamesS) #filter all bed files (samples and control)
#=======================> DONE!




