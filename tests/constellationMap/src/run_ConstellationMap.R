# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# The MIT License (MIT)

# Copyright (c) 2015 The Broad Institute of Harvard and MIT
#   
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

args <- commandArgs(trailingOnly=TRUE)

vers <- "3.0"            # R version
libdir <- args[1]
server.dir <- args[2]
patch.dir <- args[3]
javaexec <- args[4]

print("server ")
print(server.dir)
print("patch ")
print(patch.dir)
print("java ")
print(javaexec)


source(file.path(libdir, "loadRLibrary.R"))
load.packages(libdir, patch.dir, server.dir, vers)

option_list <- list(
  make_option("--input.gct.file", dest="input.gct.file"),
  make_option("--input.cls.file", dest="input.cls.file"),
  make_option("--gene.sets.database", dest="gene.sets.database"),
  make_option("--gene.sets.file", dest="gene.sets.file", default=NULL),
  make_option("--top.n", dest="top.n", default=NULL),
  make_option("--direction", dest="direction", default=NULL),
  make_option("--image.format", dest="image.format", default=NULL),
  make_option("--target.class", dest="target.class", default=NULL),
  make_option("--jaccard.threshold", dest="jaccard.threshold", default=NULL)
)

opt <- parse_args(OptionParser(option_list=option_list), positional_arguments=TRUE, args=args)
print(opt)
opts <- opt$options

# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

# Returns NULL if the value is blank or a numeric if the value parses as numeric.
# Otherwise it stops execution with an error.
require.numeric.or.null <- function(value, paramName) {
  if (is.null(value)) { return(NULL) }
  trimmedValue = trim(value)
  if (trimmedValue == "") { return(NULL) }
  suppressWarnings(val <- as.numeric(trimmedValue))
  if (is.na(val)) {
    stop(paste("The parameter '", paramName, "' must be numeric or blank.  Received '", value, "'.", sep=""))  
  }
  return(val)
}

# Source the setup code
source(file.path(libdir, "ConstellationMap_functions.R"))

gct.filename <- opts$input.gct.file
cls.filename <- opts$input.cls.file
gs.database <- opts$gene.sets.database
gs.url <- opts$gene.sets.file
top.n <- require.numeric.or.null(opts$top.n, "top.n")   # how many genesets to include in the constellationmap
direction <- opts$direction   # positive or negative
image.format <- opts$image.format
target.class <- opts$target.class
jaccard.threshold <- require.numeric.or.null(opts$jaccard.threshold, "jaccard.threshold")

sessionInfo()

# Convert direction to UP (+) or DOWN (-)
if(direction == "positive") {
  direction <- "UP"
} else if(direction == "negative") {
  direction <- "DOWN"
}

# If gene.set.database.file given, parse file; else fetch gene.sets.database from MSigDB
if(!is.null(gs.url)) {
  gs.url.split <- unlist(strsplit(gs.url, "\\."))
  if (tail(gs.url.split, n=1) == "gmt") {
    # is gmt formatted
    GSDB <- Read.GeneSets.gmt(gs.url, thres.min = 2, thres.max = 2000, gene.names = NULL)
    gs.basename <- basename(gs.url)
    gs.prefix <- sub(".gmt", "", gs.basename)
  } else {
    # is gmx formatted
    GSDB <- Read.GeneSets.gmx(gs.url, thres.min = 2, thres.max = 2000)
    gs.basename <- basename(gs.url)
    gs.prefix <- sub(".gmx", "", gs.basename)
  }
} else {
  tmpGSDBFile <- tempfile()

  Get.GeneSets.db(javaexec, libdir, gs.database, tmpGSDBFile)  # only gmt-formatted files are downloaded from MSigDB
  GSDB <- Read.GeneSets.gmt(tmpGSDBFile, thres.min=2, thres.max=2000, gene.names=NULL)
  unlink(tmpGSDBFile)
  gs.basename <- basename(gs.database)
  gs.prefix <- sub(".gmt", "", gs.basename)
}

# Error if top.n<3
if(top.n<3) {stop("top.n must be greater than 2")}

# If target.class is not given, take first class listed in .cls file as target.class
CLS <- MSIG.ReadPhenFile.2(file=cls.filename)
if(is.null(target.class)) {
  target.class <- CLS$phen.list[1]
} else if(!(target.class %in% CLS$phen.list)) {
  err.target.class <- paste("No class", target.class, "found in", cls.filename, sep=" ")
  stop(err.target.class)
}

# Default jaccard.threshold to 0.1
if(is.null(jaccard.threshold)) {
  jaccard.threshold = 0.1
}

# Error if jaccard not [0,1]
if(jaccard.threshold > 1 || jaccard.threshold < 0) {
  err.jaccard.threshold <- "jaccard.threshold must be greater than or equal to 0 and less than or equal to 1."
  stop(err.jaccard.threshold)
}

# Read the enrichment data
dataset <- MSIG.Gct2Frame( gct.filename )

# Remove gene sets that are too large/small from the enrichment dataset
if (!is.null(GSDB$gs.removed)){
  rem.lgcl <- !(dataset$row.names %in% GSDB$gs.removed)
  dataset$ds <- dataset$ds[rem.lgcl,]
  dataset$row.names <- dataset$row.names[rem.lgcl]
  dataset$descs <- dataset$descs[rem.lgcl]
}

expr <- as.matrix( dataset$ds )

dataset.genesets <- dataset$row.names
GSDB.genesets <- GSDB$gs.names

# Error if gene sets in GCT file and GMT/GMX file do not match
if(!all(dataset.genesets %in% GSDB.genesets)) {
  err.geneset.match <- paste("Failed to match all gene set names from", basename(gct.filename), "with corresponding gene sets from", gs.basename, sep=' ')
  stop(err.geneset.match)
}

# Calculate mutual information
suppressWarnings(MI.rank(
  input.ds = gct.filename,
  dataset = dataset,
  parsed.cls = CLS,
  gs.prefix = gs.prefix,
  signatures = "ALL",
  phenotype = NULL,
  target.class = target.class,
  target.type = "discrete",
  weight = 0.75,
  statistic = "area.under.RES",
  top.n = top.n,
  direction = direction,
  image.format = image.format
))

report.filename <- output.txt <- file.path(paste(sub(".gct", "", basename(gct.filename)), ".", gs.prefix, ".", target.class, ".REPORT.txt", sep=""))

# Read the annotation file
# CLS <- CLS
phens <- CLS$phen.list
target.ind <- match(target.class, phens)
if(CLS$is.continuous) {
  class.v <- as.vector(as.numeric(CLS$class.v[target.ind,]))
} else {
  class.v <- as.numeric( CLS$class.v )
  class.v <- ifelse(class.v==target.ind, 0, -1)
}
# 
# for(i in 1:length(class.v)){
#   if(class.v[i] == target.ind){
#     class.v[i] = 1
#   } else {
#     class.v[i] = 2
#   }
# }
# print(class.v)

# Read the mutual information report and select the top.n genesets
gs.up <- as.character( extract.gene.set( report.filename, top.n, direction) )
gs <- gs.up
print(gs)

output.prefix <- sub(".gct", "", basename(gct.filename))

# Create output filenames (PDF/PNG)
if(image.format=="PDF") {
  output.filename <- paste(output.prefix, gs.prefix, target.class, "CONSTELLATION_MAP.pdf", sep=".")
  print(output.filename)
  pdf(file.path(output.filename), width=11, height=8.5)
} else if(image.format=="PNG") {
  output.filename <- paste(output.prefix, gs.prefix, target.class, "CONSTELLATION_MAP.png", sep=".")
  print(output.filename)
  png(filename=file.path(output.filename), width=10, height=7.5, units="in", res=120)
}

layout(matrix(c(1,2,3,4), 2, 2, byrow=T), heights=c(7,1.5), widths=c(7, 4))

# Plot constellation map
res <- suppressWarnings(constellation.map( expr, as.vector( class.v ), gs, GSDB, jaccard.threshold, direction, target.class))
# res <- suppressWarnings(constellation.map( expr, as.vector(1 - class.v ), gs, GSDB, jaccard.threshold, target.class))

plot.data2odf.file(x=res$x,
                   y=res$y,
                   gs=gs,
                   gs.prefix=gs.prefix,
                   GSDB=GSDB,
                   jaccard.threshold=jaccard.threshold,
                   target.class=target.class,
                   gct.filename=basename(gct.filename),
                   parsed.cls=CLS,
                   cls.filename=basename(cls.filename),
                   direction=direction,
                   edgepairs=res$edgepairs,
                   jcoef=res$connectivity)

plot.new()
legend("left", paste(1:length(gs), gs, sep=": "), cex=0.55)
dev.off()

# Copy html file so it appears in job results
html.filename <- "Visualizer.html"
copy.html.file(libdir, html.filename)
