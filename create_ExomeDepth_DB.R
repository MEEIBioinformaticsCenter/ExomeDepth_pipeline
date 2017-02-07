#!/usr/bin/env Rscript

# This script creates a ED.RData file with counts to be used later for CNV calls
# ./create_ExomeDepth_DB.R -b bamlist.txt -o /output/path/ -v

start_time = Sys.time()

require(optparse) # http://cran.r-project.org/web/packages/optparse/optparse.pdf
require(ExomeDepth)
data(exons.hg19) # from ExomeDepth package

options(stringsAsFactors=FALSE) # crucial for handling BAM filenames as strings

option_list = list(
  make_option(c("-b", "--bamlist"), action="store", default='', 
              type='character', help="Path to list of BAMs"),
  make_option(c("-o", "--outdir"), action="store", default='./',
              type='character', help="Output directory [default %default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print verbose output [default %default]"),
  make_option(c("-q", "--quiet"), action="store_false", dest="verbose",
              help="Do not print verbose output (this is the default)")
)
opt = parse_args(OptionParser(option_list=option_list))

print(opt,file=stderr())

# read list of BAMs
# to avoid writing a tryCatch I use file.exists, and set default to '', which is
# a file that never exists.
if (file.exists(opt$bamlist)) {
    # read bam list directly into a vector (note use of $V1)
    bams = read.table(opt$bamlist,header=FALSE)$V1
} else {
    cat("You need to specify a valid BAM list using -b.\n",file=stderr())
    cat(paste("The filename you specified was '",opt$bamlist,"'.",sep=''),file=stderr())
    stop()
}

# read output directory
# note right now if not specified, I stop execution. 
# an alternative is to create the dir.
# see http://stackoverflow.com/questions/4216753/check-existence-of-directory-and-create-if-doesnt-exist
if (file.exists(opt$outdir)) {
    setwd(opt$outdir)
} else {
    cat("You need to specify a valid output directory using -o.\n",file=stderr())
    cat(paste("The directory you specified was '",opt$outdir,"'.",sep=''),file=stderr())
    stop()
}

if (opt$verbose) {
    cat(paste("Read BAM list from ",opt$bamlist,"\n",sep=''),file=stdout())
}


counts=getBamCounts(bed.frame = exons.hg19, bam.files = bams)

if (opt$verbose) {
    cat(paste("Calculated counts\n",sep=''),file=stdout())
}

#####
# If desired, at this point you can save the counts, then have a second script
# which re-loads them. To do that, uncomment this part and split accordingly.

save(counts,file="ED.RData")

if (opt$verbose) {
	cat(paste("Wrote counts to ",getwd(),"\n",sep=''),file=stdout())
}


