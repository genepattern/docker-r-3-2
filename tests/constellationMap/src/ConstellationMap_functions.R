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

# Auxiliary functions and definitions

MSIG.Gct2Frame <- function(filename = "NULL") { 
  #
  # Reads a gene expression dataset in GCT format and converts it into an R data frame
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  ds <- read.delim(filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T, na.strings = "")
  descs <- ds[,1]
  ds <- ds[-1]
  row.names <- row.names(ds)
  names <- names(ds)
  return(list(ds = ds, row.names = row.names, descs = descs, names = names))
} # End MSIG.Gct2Frame

#######
# MSIG.ReadPhenFile.2 <- function(file = "NULL") { 
#   #
#   # Reads a matrix of class vectors from a CLS file and defines phenotype and class labels vectors
#   #  (numeric and character) for the samples in a gene expression file (RES or GCT format)
#   #
#   # The Broad Institute
#   # SOFTWARE COPYRIGHT NOTICE AGREEMENT
#   # This software and its documentation are copyright 2003 by the
#   # Broad Institute/Massachusetts Institute of Technology.
#   # All rights are reserved.
#   #
#   # This software is supplied without any warranty or guaranteed support
#   # whatsoever. Neither the Broad Institute nor MIT can be responsible for
#   # its use, misuse, or functionality.
#   
#   cls.cont <- readLines(file)
#   num.lines <- length(cls.cont)
#   
#   # check for categorical (i.e., not continuous)
#   if (cls.cont[[1]] == "#numeric") { stop(paste0("Continuous phenotype labels detected in file ", basename(file), ". ConstellationMap currently only accepts categorical phenotype labels.")) }
#   
#   error.msg <- paste0("Failed to parse input cls file ", basename(file), ".")
#   
#   # check for 3 lines
#   if (num.lines != 3) { stop(paste0(error.msg, " Expected 3 lines. Instead saw ", num.lines, " lines.")) }
#   
#   # Determine delimiter or error out
#   line1 <- unlist(strsplit(cls.cont[[1]], " "))
#   line1 <- line1[line1!=""]
#   if (length(line1) == 3) {
#     delim <- " "
#   } else {
#     line1 <- unlist(strsplit(cls.cont[[1]], "\t"))
#     line1 <- line1[line1!=""]
#     if (length(line1) == 3) {
#       delim <- "\t"
#     } else {
#       stop(paste0(error.msg, " Improper header line. CLS files should be space- or tab-delimited."))
#     }
#   }
#   
#   line1 <- as.numeric(line1) # convert to numeric
#   
#   # Extract phen.names, col.phen
#   if (length(line1) == 3) {
#     phen.names <- NULL
#     col.phen <- NULL
#   } else {
#     l.phen.names <- match("phen.names:", line1)
#     l.col.phen <- match("col.phen:", line1)
#     phen.names <- line1[(l.phen.names + 1):(l.col.phen - 1)]
#     col.phen <- line1[(l.col.phen + 1):length(line1)]
#   }
#   
#   # Extract phen.list (second line)
#   temp <- unlist(strsplit(cls.cont[[2]], delim))
#   temp <- temp[temp!=""]
#   
#   # Check if number of phenotypes agrees with line1
#   if (line1[2] != length(temp)-1) { stop(paste0(error.msg, " Expected ", line1[2], " phenotypes. Instead saw ", length(temp)-1, " phenotypes.")) }
#   phen.list <- temp[2:length(temp)]
#   
#   # Extract class.list
#   phen <- NULL
#   temp <- unlist(strsplit(cls.cont[[3]], delim))
#   len <- length(temp)
#   # Check if number of samples agrees with line1
#   if (line1[1] != len) { stop(paste0(error.msg, " Expected ", line1[1], " samples. Instead saw ", len, " samples.")) }
#   class.list <- temp
#   classes <- unique(temp)
#   # Check if number of unique labels matches number of phenotypes
#   if (length(classes) != line1[2]) { stop(paste0(error.msg, " Expected ", line1[2], " labels. Instead saw ", length(classes), " labels.")) }
#   class.v <- match(temp, classes)
#   phen <- c(phen, classes)
#   
# #   phen <- NULL
# #   for (k in 1:(num.lines - 2)) {
# #     temp <- unlist(strsplit(cls.cont[[k + 2]], delim))
# #     if (k == 1) {
# #       len <- length(temp)
# #       class.list <- matrix(0, nrow = num.lines - 2, ncol = len)
# #       class.v <- matrix(0, nrow = num.lines - 2, ncol = len)
# #       #           phen <- NULL
# #     }
# #     class.list[k, ] <- temp
# #     classes <- unique(temp)
# #     class.v[k, ] <- match(temp, classes)
# #     #        phen[[k]] <- classes
# #     phen <- c(phen, classes)
# #   }
# #   if (num.lines == 3) {
# #     class.list <- as.vector(class.list)
# #     class.v <- as.vector(class.v)
# #     #         phen <- unlist(phen)
# #   }
#   return(list(phen.list = phen.list, phen = phen, phen.names = phen.names, col.phen = col.phen,
#               class.v = class.v, class.list = class.list))
# } # End MSIG.ReadPhenFile.2
######

MSIG.ReadPhenFile.2 <- function(file = "NULL") { 
  #
  # Reads a matrix of class vectors from a CLS file and defines phenotype and class labels vectors
  #  (numeric and character) for the samples in a gene expression file (RES or GCT format)
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  cls.cont <- readLines(file)
  num.lines <- length(cls.cont)
  file.name <-  basename(file)
  
  error.msg <- paste0("Failed to parse input cls file ", file.name, ".")
  
  # check for categorical (i.e., not continuous)
  if (cls.cont[[1]] == "#numeric") {
    # Continuous (not discrete/categorical)
    is.continuous <- T
    phen <- NULL # phen is NULL if file is continuous
    
    # Error check
    if (num.lines == 1) {
      stop(paste0(error.msg, " No phenotype labels found."))
    } else if (num.lines %% 2 != 1) { # check for odd number of lines
      stop(paste0(error.msg, " Unmatched phenotype or continuous label set."))
    }
    
    phen.lines <- seq(2, num.lines, 2)
    phen.lines <- unlist(cls.cont[phen.lines])
    
    label.lines <- seq(3, num.lines, 2)
    label.lines <- unlist(cls.cont[label.lines])
    
    phen.list <- NULL
    class.v <- NULL
    class.list <- NULL
    
    # Determine delimiter or error out
    delim <- " "
    label.tmp <- unlist(strsplit(label.lines[1], delim))
    label.tmp.str <- label.tmp[label.tmp!=""]
    label.tmp.num <- as.numeric(label.tmp.str)
    if (any(is.na(label.tmp.num))) {
      delim <- "\t"
      label.tmp <- unlist(strsplit(label.lines[1], delim))
      label.tmp.str <- label.tmp[label.tmp!=""]
      label.tmp.num <- as.numeric(label.tmp.str)
      if (any(is.na(label.tmp.num))) {
        stop(paste0(error.msg, " Non-numeric values found in line 3. CLS files should be space- or tab-delimited."))
      } else {
        warning(paste0("Tab-delimiters detected in ", file.name, ". Space-delimited CLS files are preferred."))
      }
    }
    
    n.row <- length(label.lines)
    n.col <- length(label.tmp.num)
    
    for (i in 1:length(label.lines)) {
      # Parse phenotype lines
      phen.tmp <- phen.lines[i]
      if (substr(phen.tmp, 1, 1) != "#") {
        stop(paste0(error.msg, " In continuous type CLS files phenotypes must be preceeded by a pound/hash (#)."))
      } else {
        phen.tmp <- substr(phen.tmp, 2, nchar(phen.tmp))
        phen.tmp <- gsub(delim, "", phen.tmp) # remove trailing delimiters
        phen.list <- c(phen.list, phen.tmp)
      }
      
      if(i > 1) {
        label.tmp <- unlist(strsplit(label.lines[i], delim))
        label.tmp.str <- label.tmp[label.tmp!=""]
        label.tmp.num <- as.numeric(label.tmp.str)
        
        if (any(is.na(label.tmp.num))) { stop(paste0(error.msg, " Non-numeric values found in line ", i, ".")) }
        if (length(label.tmp.num) != n.col) { stop(paste0(error.msg, " Missing values in line ", i, ". Expected ", n.col, " values.")) }
      }
      
      class.v <- c(class.v, label.tmp.num)
      class.list <- c(class.list, label.tmp.str)
    }
    # Make class.v, class.list matrices
    class.v <- matrix(class.v, nrow=n.row, ncol=n.col, byrow=T)
    class.list <- matrix(class.list, nrow=n.row, ncol=n.col, byrow=T)
  } else {
    # Discrete/categorical (not continuous)
    is.continuous <- F
    
    # Check for 3 lines
    if (num.lines != 3) { stop(paste0(error.msg, " Expected 3 lines in categorical type CLS file. Instead saw ", num.lines, " lines.")) }
    
    # Determine delimiter or error out
    delim <- " "
    line1 <- unlist(strsplit(cls.cont[[1]], delim))
    line1 <- line1[line1!=""]
    if (length(line1) != 3) {
      delim <- "\t"
      line1 <- unlist(strsplit(cls.cont[[1]], delim))
      line1 <- line1[line1!=""]
      if (length(line1) == 3) {
        warning(paste0("Tab-delimiters detected in ", file.name, ". Space-delimited CLS files are preferred."))
      } else {
        stop(paste0(error.msg, " Improper header line. CLS files should be space- or tab-delimited."))
      }
    }
    
    line1 <- as.numeric(line1) # convert to numeric
    
    # Extract phen.names, col.phen
    if (length(line1) == 3) {
      phen.names <- NULL
      col.phen <- NULL
    } else {
      l.phen.names <- match("phen.names:", line1)
      l.col.phen <- match("col.phen:", line1)
      phen.names <- line1[(l.phen.names + 1):(l.col.phen - 1)]
      col.phen <- line1[(l.col.phen + 1):length(line1)]
    }
    
    # Extract phen.list (second line)
    temp <- unlist(strsplit(cls.cont[[2]], delim))
    temp <- temp[temp!=""]
    
    # Check if number of phenotypes agrees with line1
    if (line1[2] != length(temp)-1) { stop(paste0(error.msg, " Expected ", line1[2], " phenotypes. Instead saw ", length(temp)-1, " phenotypes.")) }
    phen.list <- temp[2:length(temp)]
    
    # Extract class.list
    phen <- NULL
    temp <- unlist(strsplit(cls.cont[[3]], delim))
    temp <- temp[temp!=""]
    len <- length(temp)
    # Check if number of samples agrees with line1
    if (line1[1] != len) { stop(paste0(error.msg, " Expected ", line1[1], " samples. Instead saw ", len, " samples.")) }
    class.list <- temp
    classes <- unique(temp)
    # Check if number of unique labels matches number of phenotypes
    if (length(classes) != line1[2]) { stop(paste0(error.msg, " Expected ", line1[2], " labels. Instead saw ", length(classes), " labels.")) }
    class.v <- match(temp, classes)
    phen <- c(phen, classes)
  }
  
  return(list(phen.list = phen.list, phen = phen, class.v = class.v, class.list = class.list, is.continuous = is.continuous))
} # End MSIG.ReadPhenFile.2

# FUNCTIONS FOR READING GENE SET FILES (GMT, GMX)
#----------------------------------------------------------------------

# Read.GeneSets.gmt <- function(
#   gs.db,
#   thres.min = 2,
#   thres.max = 2000,
#   gene.names = NULL) {
#   
#   temp <- readLines(gs.db)
#   max.Ng <- length(temp)
#   temp.size.G <- vector(length = max.Ng, mode = "numeric") 
#   for (i in 1:max.Ng) {
#     temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
#   }
#   max.size.G <- max(temp.size.G)      
#   gs <- matrix(rep("null", max.Ng*max.size.G), nrow=max.Ng, ncol= max.size.G)
#   temp.names <- vector(length = max.Ng, mode = "character")
#   temp.desc <- vector(length = max.Ng, mode = "character")
#   gs.count <- 1
#   for (i in 1:max.Ng) {
#     gene.set.size <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
#     gs.line <- noquote(unlist(strsplit(temp[[i]], "\t")))
#     gene.set.name <- gs.line[1] 
#     gene.set.desc <- gs.line[2] 
#     gene.set.tags <- vector(length = gene.set.size, mode = "character")
#     for (j in 1:gene.set.size) {
#       gene.set.tags[j] <- gs.line[j + 2]
#     }
#     if (is.null(gene.names)) {
#       existing.set <- rep(TRUE, length(gene.set.tags))
#     } else {
#       existing.set <- is.element(gene.set.tags, gene.names)
#     }
#     set.size <- length(existing.set[existing.set == T])
#     if ((set.size < thres.min) || (set.size > thres.max)) next
#     temp.size.G[gs.count] <- set.size
#     gs[gs.count,] <- c(gene.set.tags[existing.set], rep("null", max.size.G - temp.size.G[gs.count]))
#     temp.names[gs.count] <- gene.set.name
#     temp.desc[gs.count] <- gene.set.desc
#     gs.count <- gs.count + 1
#   }
#   Ng <- gs.count - 1
#   gs.names <- vector(length = Ng, mode = "character")
#   gs.desc <- vector(length = Ng, mode = "character")
#   size.G <- vector(length = Ng, mode = "numeric") 
#   
#   gs.names <- temp.names[1:Ng]
#   gs.desc <- temp.desc[1:Ng]
#   size.G <- temp.size.G[1:Ng]
#   
#   return(list(N.gs = Ng, gs = gs, gs.names = gs.names, gs.desc = gs.desc, size.G = size.G, max.N.gs = max.Ng))
# } # End Read.GeneSets.gmt

Read.GeneSets.gmt <- function(
  gs.db,
  thres.min = 2,
  thres.max = 2000,
  gene.names = NULL) {
  
  temp <- readLines(gs.db)
  max.Ng <- length(temp)
  temp.size.G <- vector(length = max.Ng, mode = "numeric") 
  for (i in 1:max.Ng) {
    temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
  }
  max.size.G <- max(temp.size.G)      
  gs <- matrix(rep("null", max.Ng*max.size.G), nrow=max.Ng, ncol= max.size.G)
  temp.names <- vector(length = max.Ng, mode = "character")
  temp.desc <- vector(length = max.Ng, mode = "character")
  gs.removed <- vector(length = max.Ng, mode = "character")
  gs.count <- 1
  gs.removed.count <- 1
  for (i in 1:max.Ng) {
    gene.set.size <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
    gs.line <- noquote(unlist(strsplit(temp[[i]], "\t")))
    gene.set.name <- gs.line[1] 
    gene.set.desc <- gs.line[2] 
    gene.set.tags <- vector(length = gene.set.size, mode = "character")
    for (j in 1:gene.set.size) {
      gene.set.tags[j] <- gs.line[j + 2]
    }
    if (is.null(gene.names)) {
      existing.set <- rep(TRUE, length(gene.set.tags))
    } else {
      existing.set <- is.element(gene.set.tags, gene.names)
    }
    set.size <- length(existing.set[existing.set == T])
    if ((set.size < thres.min) || (set.size > thres.max)) {
      gs.removed[gs.removed.count] <- gene.set.name
      gs.removed.count <- gs.removed.count + 1
      next
    }
    temp.size.G[gs.count] <- set.size
    gs[gs.count,] <- c(gene.set.tags[existing.set], rep("null", max.size.G - temp.size.G[gs.count]))
    temp.names[gs.count] <- gene.set.name
    temp.desc[gs.count] <- gene.set.desc
    gs.count <- gs.count + 1
  }
  Ng <- gs.count - 1
  gs.names <- vector(length = Ng, mode = "character")
  gs.desc <- vector(length = Ng, mode = "character")
  size.G <- vector(length = Ng, mode = "numeric") 
  
  gs.names <- temp.names[1:Ng]
  gs.desc <- temp.desc[1:Ng]
  size.G <- temp.size.G[1:Ng]
  if (gs.removed.count == 1) {
    gs.removed <- NULL
  } else{
    gs.removed <- gs.removed[1:(gs.removed.count-1)]
  }
  
  # Remove "null"s
  rem.lgcl <- apply(gs, 1, function(x) !all(x == "null"))
  gs <- gs[rem.lgcl,]
  
  return(list(N.gs = Ng, gs = gs, gs.names = gs.names, gs.desc = gs.desc, size.G = size.G, max.N.gs = max.Ng, gs.removed = gs.removed))
} # End Read.GeneSets.gmt

Read.GeneSets.gmx <- function(
  gs.gmx,
  thres.min = 2,
  thres.max = 2000) {
  
  # gmx file defines gene sets column-wise
  # thus, length of temp (number of lines/rows of gmx file) is the size of the largest gene set in collection
  
  df.temp <- t(read.table(gs.gmx, header=TRUE, sep="\t"))
  all.gs.names <- row.names(df.temp)
  all.gs.desc <- as.vector(df.temp[,1])
  all.gs <- as.matrix(df.temp[,-1])
  all.gs.sizes <- as.vector(rowSums(all.gs != ""))
  pass.thresholds <- (all.gs.sizes >= thres.min & all.gs.sizes <= thres.max)  
  gs.names <- all.gs.names[pass.thresholds]
  gs.removed <- all.gs.names[!pass.thresholds]
  if (length(gs.removed) == 0) { gs.removed <- NULL }
  gs.desc <- all.gs.desc[pass.thresholds]
  gs.sizes <- all.gs.sizes[pass.thresholds]
  gs <- all.gs[pass.thresholds,]
  max.Ng <- length(all.gs.names)
  Ng <- length(gs.names)
  
  # N.gs = number of gene sets defined in gmx file that satisfy the min and max thresholds
  # gs = matrix containing gene set collections, one per line, satisfying min/max thresholds
  # gs.names = vector of names of gene sets (of length N.gs)
  # gs.desc = vector of descriptions of gene sets (of length N.gs)
  # size.G = vector with sizes of each gene set (of length N.gs)
  # max.N.gs = total number of gene sets defined in gmx file; includes those that do not satisfy min/max thresholds
  return (list(N.gs = Ng, gs = gs, gs.names = gs.names, gs.desc = gs.desc, size.G = gs.sizes, max.N.gs = max.Ng, gs.removed = gs.removed))
  
} # End of Read.GeneSets.gmx

# OPAM Functions
#----------------------------------------------------------------------

OPAM.Evaluate.Results.2 <- function(
  input.ds,
  dataset,
  parsed.cls,
  target.class = NULL,
  sort.results = T,
  display.top.n = 20,
  output.txt,
  output.img,
  direction,
  image.format) {
  
  if(image.format=="PDF") {
    pdf(file=output.img, height=8.5, width=11)
  } else if(image.format=="PNG") {
    png(filename=output.img, height=7.5, width=11, units="in", res=120)
  }
  
#   dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
  m <- data.matrix(dataset$ds)
  N <- dim(m)[1]
  model.names <- dataset$row.names
  model.descs <- dataset$descs
  Ns <- length(as.matrix(m)[1,])
  for (i in 1:N) {
    if (sd(as.matrix(m)[i,]) == 0) {
      val <- as.matrix(m)[i, 1]
      m[i,] <- as.matrix(m)[i,] + runif(n=Ns, min= val - 0.001, max=val + 0.001)  # add small noise to flat profiles
    }
  }
  
  sample.names <- dataset$names
  
  # Get phenotypic target vector
  #   CLS <- MSIG.ReadPhenFile.2(file = input.cls) # Read phenotype file (CLS format)
  #   cls.labels <- CLS$class.v
  #   cls.phen <- CLS$phen
  #   cls.list <- CLS$class.list 
  if (parsed.cls$is.continuous) {
    phen.loc <- match(target.class, parsed.cls$phen.list)
    target <- parsed.cls$class.v[phen.loc,]
  } else {
    target.symbol <- parsed.cls$phen[match(target.class, parsed.cls$phen.list)]
    target <- ifelse(parsed.cls$class.list == target.symbol, 1, 0)
  }
  
  #   # Added some code to accept class names and convert to symbols
  #   target.class <- CLS$phen[match(target.class, CLS$phen.list)]
  #   
  #   if (is.null(phenotype)) {
  #     phen.loc <- 1
  #   } else {
  #     phen.loc <- match(phenotype, CLS$phen.list)
  #   }
  #   if (is.vector(CLS$class.list)) {
  #     target.vec <- CLS$class.list
  #   } else {
  #     target.vec <- CLS$class.list[phen.loc,]
  #   }
  #   if (target.type == "continuous") {
  #     target <- target.vec
  #   } else if (target.type == "discrete") {
  #     target <- ifelse(target.vec == target.class, 1, 0)    
  #   }
  
  ind <- order(target)
  target <- target[ind]
#   target.vec <- target.vec[ind]
  m <- as.matrix(m)[, ind]
  sample.names <- sample.names[ind]
  #   class.v <- CLS$class.v
  #   if (is.vector(class.v)) {
  #     class.v <- class.v[ind]
  #   } else {
  #     class.v <- class.v[, ind]
  #   }
  annot <- MI <- AUC <- AUC.pval <- t.stat <- t.pval <- vector(length=N, mode="numeric")
  
  NMI.ref <- mutual.inf.P(x = target, y = target, n.grid=100)$NMI
  
  for (i in 1:N) {
    if (N == 1) {
      feature <- m
    } else {
      feature <- m[i,]
    }
    MI[i] <- signif(mutual.inf.P(target, feature, n.grid=100)$NMI/NMI.ref, 4)
    if (parsed.cls$is.continuous) {
      #     if (target.type == "continuous") {
      AUC[i] <- AUC.pval[i] <- t.stat[i] <- t.pval[i] <- "-"
    } else {
      #     } else if (target.type == "discrete") {
      feature.norm <- (feature - min(feature))/(max(feature) - min(feature))
      perf.auc <- roc.area(target, feature.norm)
      AUC[i] <- ifelse(perf.auc$A < 0.5, -(1 - perf.auc$A), perf.auc$A)
      AUC[i] <- signif(AUC[i], digits=4)
      p.val <- perf.auc$p.value
      p.val <- ifelse(p.val > 0.5, 1 - p.val, p.val)
      AUC.pval[i] <- signif(p.val, digits=4)
      temp <- split(feature, target)
      x <- temp$'1'
      y <- temp$'0'
      t.stat[i] <- signif(t.test(x=x, y=y)$statistic, digits=4)
      p.val <- t.test(x=x, y=y)$p.value
      p.val <- ifelse(p.val > 0.5, 1 - p.val, p.val)
      t.pval[i] <- signif(p.val, digits=4)
    }
    annot[i] <- paste(MI[i], "     ", AUC[i], " (", AUC.pval[i], ")    ", t.stat[i], " (", t.pval[i], ")", sep="")
  }
  
  if ((N > 1) & (sort.results == T)) {
    direction.isup <- ifelse(direction=="UP", T, F)
    MI.order <- order(MI, decreasing=direction.isup)
    MI <- MI[MI.order]
    AUC <- AUC[MI.order]
    AUC.pval <- AUC.pval[MI.order]
    t.stat <- t.stat[MI.order]
    t.pval <= t.pval[MI.order]
    m <- as.matrix(m)[MI.order,]
    annot <- annot[MI.order]
  }
  
  mycol <- vector(length=512, mode = "numeric")
  for (k in 1:256) mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
  for (k in 257:512) mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
  mycol <- rev(mycol)
  cex.axis = 1
  ncolors <- length(mycol)
  
  nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), 1, c(4, 10), FALSE)
  par(mar = c(1, 15, 5, 15))
  par(oma = c(0, 5, 0, 0))
  max.v <- max(max(target), -min(target))
  V1 <- target
  if (parsed.cls$is.continuous) {
    #   if (target.type == "continuous") {
    target.min <- min(target)
    target.max <- max(target)
    
    col.matrix <- colorRamp(c("yellow","purple"))(1:100 / 100)
    col.matrix <- round(col.matrix)
    col.vector <- apply(col.matrix, 1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
    image(1:length(target), 1:1, as.matrix(V1), zlim = c(target.min, target.max), col=col.vector, axes=FALSE, main="", sub = "", xlab= "", ylab="")
    axis(3, at=c(1, round(length(target)/2), length(target)), labels=c(target.min, target.class, target.max), adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
  } else {
    #   } else if (target.type == "discrete") {
    image(1:length(target), 1:1, as.matrix(V1), zlim = c(0, 1), col=c("yellow", "purple"), axes=FALSE, main="", sub = "", xlab= "", ylab="")
    axis(3, at=length(target):length(target), labels=target.class, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
  }
  
  #Yan made some changes here
  #   axis(3, at=length(target):length(target), labels=target.class, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
  #   axis(3, at=1:1, labels=paste(phenotype, target.class), adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
  axis(4, at=1:1, labels="NMI     AUC (p-val)     t-test (p-val)", adj= 0.5, tick=FALSE, 
       las = 1, cex.axis=0.80, font.axis=1, line=-1) 
  par(mar = c(5, 15, 1, 15))
  
  if (direction.isup) {
    N.max <- sum(MI > 0)
  } else {
    N.max <- sum(MI <= 0)
  }
  
  if (display.top.n > N.max) {display.top.n <- N.max}
  
  if (N == 1) {
    V1 <- m
    V1 <- (V1 - mean(V1))/sd(V1)
  } else {
    V1 <- m[1:display.top.n, ]
    for (i in 1:display.top.n) V1[i,] <- (V1[i,] - mean(V1[i,]))/sd(V1[i,])
  }
  
  max.v <- max(max(V1), -min(V1))
  V1 <- ceiling(ncolors * (V1 - (- max.v))/(1.001*(max.v - (- max.v))))
  
  if (N > 1) {
    V1 <- apply(V1, MARGIN=2, FUN=rev)
    image(1:Ns, 1:display.top.n, t(V1), zlim = c(0, ncolors), col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
    axis(2, at=1:display.top.n, labels=row.names(V1), adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
    axis(4, at=1:display.top.n, labels=rev(annot[1:display.top.n]), adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
    axis(1, at=1:Ns, labels=sample.names, adj= 0.5, tick=FALSE, las = 3, cex.axis=0.60, font.axis=1, line=-1)
  } else {
    image(1:Ns, 1:1, as.matrix(V1), zlim = c(0, ncolors), col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
    axis(2, at=1:1, labels=model.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
    axis(4, at=1:1, labels=annot[1], adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
    axis(1, at=1:Ns, labels=sample.names, adj= 0.5, tick=FALSE, las = 3, cex.axis=0.60, font.axis=1, line=-1)
  }
  dev.off()
  NMI <- MI # Express normalized mutual information with proper acronym
  annot2 <- data.frame(cbind(NMI, AUC, AUC.pval, t.stat, t.pval))
  row.names(annot2) <- row.names(m)
  write(paste(c("gene set ", noquote(colnames(annot2))), collapse="\t"), file = output.txt, append = F, 
        ncolumns = length(colnames(annot2)))
  write.table(annot2, file=output.txt, append=T, quote=F, sep="\t", eol="\n", col.names=F, row.names=T)
  
} # End OPAM.Evaluate.Results.2

mutual.inf.P <- function(x, y, n.grid=100) {
  # for definitions of mutual information and the universal metric (NMI) see the 
  # definition of "Mutual Information" in wikipedia and Thomas and Cover's book
  
  #   kde2d.xy <- kde2d(x, y, n = n.grid, h = c(width.SJ(x, method="dpi"), width.SJ(y, method="dpi")))
  kde2d.xy <- kde2d(x, y, n = n.grid, h = c(bcv(x), bcv(y)))
  X <- kde2d.xy$x
  Y <- kde2d.xy$y  
  #   PXY <- kde2d.xy$z/sum(kde2d.xy$z)
  PXY <- kde2d.xy$z/sum(kde2d.xy$z) + .Machine$double.eps
  
  #   PX <- apply(PXY, MARGIN=1, sum)
  PX <- rowSums(PXY)
  PX <- PX/sum(PX)
  HX <- -sum(PX * log2(PX))
  PX <- matrix(PX, nrow=n.grid, ncol=n.grid)
  
  #   PY <- apply(PXY, MARGIN=2, sum)
  PY <- colSums(PXY)
  PY <- PY/sum(PY)
  HY <- -sum(PY * log2(PY))
  PY <- matrix(PY, byrow = TRUE, nrow=n.grid, ncol=n.grid)
  
  MI <- sum(PXY * log2(PXY/(PX*PY)))
  HXY <- -sum(PXY * log2(PXY))
  NMI <- sign(cor(x, y)) * ((HX + HY)/HXY - 1)  # use pearson correlation the get the sign (directionality)
  
  return(list(MI=MI, HXY=HXY, HX=HX, HY=HY, NMI=NMI))
} # End mutual.inf.P


# Other functions

MI.rank <- function( 
  input.ds,   # the input filename with directory information
  dataset, # parsed input file
  parsed.cls,   # the input class file with directory information
  gs.prefix,   # the basename of the gene set file
  #                     output.dir, # the directory to save output files
  #                     res.sub = "Test",   # the sub folder to hold the data, default as Test folder
  signatures = "ALL",
  phenotype = NULL,
  target.class =  "TRUE",
  target.type  =  "discrete",
  weight       =  0.75,
  statistic    =  "area.under.RES",
  top.n,
  direction,
  image.format
) {
  # Get the base name of the input dataset file
  input.ds.name <- basename(input.ds)
  output.prefix <- sub(".gct", "", input.ds.name)
  
  # Generate the output name with the directory information
  output.txt <- file.path(paste(output.prefix, ".", gs.prefix, ".", target.class, ".REPORT.txt", sep=""))
  if(image.format=="PDF") {
    output.img <- file.path(paste(output.prefix, ".", gs.prefix, ".", target.class, ".HEATMAP.pdf", sep=""))
  } else if(image.format=="PNG") {
    output.img <- file.path(paste(output.prefix, ".", gs.prefix, ".", target.class, ".HEATMAP.png", sep=""))
  }
  
  OPAM.Evaluate.Results.2(
    input.ds           = input.ds,
    dataset            = dataset,
    parsed.cls         = parsed.cls,
    target.class       = target.class,
    sort.results       = T,
    display.top.n      = top.n,
    output.txt         = output.txt,
    output.img         = output.img,
    direction          = direction,
    image.format       = image.format)
  
} # End MI.rank

# function to generate constellation map
constellation.map <- function( expr, class.v, gs, GSDB, jaccard.threshold, direction, target.class) {
  # expr is the enrichment matrix where samples are in columns
  # class.list is a list of phenotypes corresponding to each sample in the expr matrix
  # gs is the gene sets whose distances will be calculated
  
  ### First calculate the pairwise mutual information
  row.ind <- match( as.character( gs ), rownames( expr ) )
  df <- as.matrix( rbind( class.v, expr[ row.ind, ] ) )
  rownames( df )[ 1 ] <- "Phenotype"
  n.fea <- dim(df)[1]
  NMI.dist.list <- NMI.dist.mat(df)
  dist.matrix <- NMI.dist.list$dist.matrix
  phen.NMI <- NMI.dist.list$phen.NMI

  # Normalize the distance matrix to 0 - 1
  dist.matrix <- (dist.matrix - min(dist.matrix)) / (max(dist.matrix) - min(dist.matrix))
  diag(dist.matrix) <- 0
  
  # Calculate the variation within all the features, didn't use it in the current version
  max.2.phen <- max(dist.matrix[1, ])
  max.bw.gs <- max(dist.matrix[ -1, -1])
  pvar <- max.bw.gs / (2 * max.2.phen)   # percent variation
  #print(pvar)
  
  ### MDS for the distance matrix
  fea.dist.matrix <- dist.matrix[2:n.fea, 2:n.fea]
  # smacof.map <- smacofSphere.primal(fea.dist.matrix, ndim = 2, weightmat = NULL, init = NULL,
  #                         metric = TRUE, ties = "primary", verbose = FALSE, modulus = 1, itmax = 1000, eps = 1e-6)

  smacof.map <- smacofSphere(fea.dist.matrix, algorithm="primal", ndim = 2, weightmat = NULL, init = NULL, ties = "primary", verbose = FALSE, modulus = 1, itmax=1000, eps = 1e-6)
  # Calculate the coordinates
  x0 <- smacof.map$conf[,1]
  y0 <- smacof.map$conf[,2]
  r <- sqrt(x0*x0 + y0*y0)
  radius <- dist.matrix[1, 2:n.fea]
  if(direction == "DOWN") { radius <- -radius }

  # Normalize the radius distance, convert the radius into relative distance
  radius <- (radius - min(radius)) / (max(radius) - min(radius))	+ 0.5
  x <- x0 * radius / r
  y <- y0 * radius / r
  angles <- atan2(y, x)
  
  ### Plot the constellation map 
  
  # Here calculate the connectivity of genesets
  npoints <- length(gs)
  connectivity.mat <- matrix( 0, nrow = npoints, ncol = npoints )
  size.G <- gs.size(GSDB, gs)
  jcoef <- jaccard.coef(GSDB, gs)
  
  const.plot.list <- constellation.plot(x, y, radius, phen.NMI, angles, jcoef, jaccard.threshold, target.class)
  edgepairs <- const.plot.list$edgepairs
  NMI.annot.df <- const.plot.list$NMI.annot.df
  
  return(list(x=x, y=y, radius=radius, angles=angles, connectivity=jcoef, edgepairs=edgepairs, NMI.annot.df=NMI.annot.df, phen.NMI=phen.NMI))
} # End constellation.map

constellation.plot <- function(x, y, radius, phen.NMI, angles, connectivity.mat, jaccard.threshold, target.class){
  
  # Plot the outline of the constellation plot
  plot(x, y, pch=20, bty="n", xaxt='n', axes = FALSE, type="n", xlab="", ylab="",
       xlim=1.5*c(-max(radius), max(radius)), ylim=1.5*c(-max(radius), max(radius)))
  
  # Plot the outline of the radial plot
  radial.col <- "gray80"
  line.angle <- seq(0, 2*pi-0.001, 0.001)
  
  # Plot circles with annotations
  phen.NMI <- phen.NMI[which(!grepl("Phenotype", row.names(phen.NMI))), 'Phenotype'] #Remove 'Phenotype' from phen.NMI
  radius.unq <- unique(radius)
  tmp.ln <- length(radius.unq)
  NMI.annot <- rep(0, tmp.ln)
  max.rings <- 5 # Plot only up to max.rings
  if(tmp.ln <= max.rings) {
    indx.to.print <- 1:max.rings
  } else { # If more radii than max.rings, then choose rings to plot
    indx.to.print <- c(1, tmp.ln)
    tmp.int <- floor((tmp.ln-2)/(max.rings-2) * 1:(max.rings-2)) + 1
    indx.to.print <- c(indx.to.print, tmp.int)
  }
  prev.NMI <- -999
  for(i in 1:tmp.ln) {
    cur.NMI.range <- as.numeric(phen.NMI[radius %in% radius.unq[i]])
    NMI.annot[i] <- signif(mean(cur.NMI.range), 2)
    if((i %in% indx.to.print) && (NMI.annot[i] != prev.NMI)) {
      line.max.x <- radius.unq[i] * cos(line.angle)
      line.max.y <- radius.unq[i] * sin(line.angle)
      points(line.max.x, line.max.y, type="l", col=radial.col, lwd=1)
      text(0, radius.unq[i], NMI.annot[i], cex=0.5, col=radial.col, pos=3, offset=0.1)
      prev.NMI <- NMI.annot[i]
    }
  }
  NMI.annot.df <- data.frame('radius'=radius.unq, 'NMI.annot'=NMI.annot, 'in.plot'=(1:tmp.ln %in% indx.to.print))

  # Plot the features
  points(x, y, pch=1, col="darkblue", cex=2)
  
  # Plot the phenotype
  points(0, 0, pch = 1, col = "red", cex = 1.5)
  text(0, 0, labels = target.class, cex = 1, col = "red", pos = 1)
  
  n.fea <- length(x)
  textxy(x, y, paste(1:n.fea), cex = 0.65, offset=0)
  
  # Plot the connectivity lines
  edgepairs <- matrix(nrow=length(x)*(length(x)-1)/2, ncol=2, NA)
  counter <- 1
  npoints <- dim(connectivity.mat)[1]
  for (i in 1:(npoints - 1)){
    for(j in (i + 1):npoints){
      coef <- connectivity.mat[ i, j ] 
      if (coef  >= jaccard.threshold){ # Jaccard coef threshold (default 0.1)
        lines(c(x[i], x[j]), c(y[i], y[j]), col = "darkgreen", lwd = 15 * coef)
        edgepairs[counter,] <- c(i,j)
        counter <- counter + 1
      }
    }
  }
  
  # Remove NAs from edgepairs
  if(all(is.na(edgepairs))) {
    edgepairs <- NULL
  } else {
    edgepairs <- edgepairs[!is.na(edgepairs[,1]),]
    # coerce to matrix if removing NAs gives vector
    if(!is.matrix(edgepairs)) {edgepairs <- t(as.matrix(edgepairs))}
  }
  
  return(list('edgepairs'=edgepairs, 'NMI.annot.df'=NMI.annot.df))
} # End constellation.plot

# Calculate the pairwise mutual information of a data matrix, each row is one feature here
NMI.dist.mat <- function(df){
  
  selected.features <- row.names( df )
  n.fea <- length( selected.features )
  
  # Normalized mutual information matrix
  NMI.matrix <- matrix( 0, nrow = n.fea, ncol = n.fea )
  
  # mutual information transformed distance matrix
  dist.matrix <- matrix(1, nrow = n.fea, ncol = n.fea)
  
  row.names( NMI.matrix ) <- selected.features
  colnames( NMI.matrix ) <- selected.features
  rownames( dist.matrix ) <- selected.features
  colnames( dist.matrix ) <- selected.features
  
  # M: number of samples
  M <- ncol( df )
  
  for ( i in 1:n.fea ){
    for ( j in 1:n.fea ){
      x <- df[ selected.features[ i ], ]
      y <- df[ selected.features[ j ], ]
      
      # To be consistent with ssGSEA
      #NMI.matrix[i,j] <- signif(mutual.inf.P(x, y, n.grid=100)$NMI, 2)
      NMI.matrix[i,j] <- mutual.inf.P(x, y, n.grid=100)$NMI
      #dist.matrix[i, j] <- 1 - NMI.matrix[i, j]
    }
  }
  
  phen.NMI <- data.frame('Phenotype' = signif(NMI.matrix[,'Phenotype']/NMI.matrix['Phenotype','Phenotype'], 4))
  row.names(phen.NMI) <- selected.features
  
  NMI.matrix <- signif(NMI.matrix, 3)
  dist.matrix <- dist.matrix - NMI.matrix
  
  return(list('dist.matrix'=dist.matrix, 'phen.NMI'=phen.NMI))
} # End NMI.dist.mat

# This function will calculate the size of the gene set and normalize it to the range of 0.1 - 0.6
gs.size <- function(GSDB, genesets) {
  size.G <- GSDB$size.G
  gs.names <- GSDB$gs.names
  locs <- match(genesets, gs.names)
  size.G <- size.G[locs]	
  min.G <- min(size.G)
  max.G <- max(size.G)
  size.G <- signif((size.G - min.G) * 0.5 / (max.G - min.G) + 0.1)
  
  return(size.G)
} # End gs.size

# This function will calculate the connection between two genesets, based on Jaccard Coefficient
jaccard.coef <- function(GSDB, genesets){
#   write.table(GSDB,"GSDB.txt",append = F,quote = F,sep='\t')
#   write.table(genesets,"genesets.txt",append = F,quote = F,sep='\t')
  # Read gene set databases
  max.G <- 0
  size.G <- GSDB$size.G
  N.gs <- GSDB$N.gs
  max.G <- max(GSDB$size.G)
  gs.names <- GSDB$gs.names
  
  gs <- matrix("null", nrow = N.gs, ncol = max.G)
  gs[1:N.gs, 1:max.G] <- GSDB$gs[1:N.gs, 1:max.G]
  
  # Select desired genesets
  locs <- match(genesets, gs.names)
  gs <- gs[locs, ]
  gs.names <- gs.names[locs]
  size.G <- size.G[locs]
  
  N.gs <- sum(!is.na(locs))
  
  # Loop over genesets to generate the Jaccard coefficient matrix
  jac.coef <- matrix(1, nrow = N.gs, ncol = N.gs)
  for (i in 1:(N.gs - 1)){
    gene.set.i <- gs[i, 1:size.G[i]]
    for (j in (i + 1):N.gs){
      gene.set.j <- gs[j, 1:size.G[j]]
      overlap <- length(intersect(gene.set.i, gene.set.j))
      union <- length(union(gene.set.i, gene.set.j))
      coef <- signif(overlap / union, 4)
      jac.coef[i, j] <- coef
      jac.coef[j, i] <- coef
    }	
  }
  
  return(jac.coef)	
} # End jaccard.coef

########
# This function will calculate the jaccard index as well as report the number of genes in the gene set
# gene.set.stats <- function(gsdb, genesets){
#   max.G <- 0
#   GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
#   size.G <- GSDB$size.G
#   N.gs <- GSDB$N.gs
#   max.G <- max(GSDB$size.G)
#   gs.names <- GSDB$gs.names
#   
#   gs <- matrix("null", nrow = N.gs, ncol = max.G)
#   gs[1:N.gs, 1:max.G] <- GSDB$gs[1:N.gs, 1:max.G]
#   
#   # Select desired genesets
#   locs <- match(genesets, gs.names)
#   gs <- gs[locs, ]
#   gs.names <- gs.names[locs]
#   size.G <- size.G[locs]
#   
#   N.gs <- sum(!is.na(locs))
#   
#   # Loop over genesets to generate the Jaccard coefficient matrix
#   jac.coef <- matrix(1, nrow = N.gs, ncol = N.gs)
#   for (i in 1:(N.gs - 1)){
#     gene.set.i <- gs[i, 1:size.G[i]]
#     for (j in (i + 1):N.gs){
#       gene.set.j <- gs[j, 1:size.G[j]]
#       overlap <- length(intersect(gene.set.i, gene.set.j))
#       union <- length(union(gene.set.i, gene.set.j))
#       coef <- signif(overlap / union, 4)
#       jac.coef[i, j] <- coef
#       jac.coef[j, i] <- coef
#     }	
#   }
#   
#   return(list(coef = jac.coef, size = size.G))
# } # End gene.set.stats
########

extract.gene.set <- function(report.file, top.n, target){
  report <- read.table(report.file, header=T, sep="\t")
  Ns <- dim(report)[1]
  if(target == "UP" && top.n>sum(report$NMI>0)){
    exit.str1 <- paste0("top.n must not exceed number of gene sets with positive MI scores. Check that top.n is less than ",
                        sum(report$NMI>0), ".")
    stop(exit.str1)
  } else if(target == "DOWN" && top.n>sum(report$NMI<0)) {
    exit.str2 <- paste0("top.n must not exceed number of gene sets with negative MI scores. Check that top.n is less than ",
                        sum(report$NMI<0), ".")
    stop(exit.str2)
  }
  
  gene.sets <- report[1:top.n, 1]
  
  return(gene.sets)
} # End extract.gene.set

write.odf <- function(input.df, out.filename, header.list=NULL, convert.integers=FALSE)
{
  # Basic ODF file writer. Builds header.
  # Authors: Felix Wu
  # Args:
  #   -out.filename: name of final outputted ODF file.
  #   -input.df: input dataframe to be outputted as main body of data in ODF file.
  #   -header.list: list whose members are character vectors; each member denotes a header line.
  #     Names of list members are set as keywords in header lines.
  #     If no list is provided the basic header lines are created based on the column names of input.df.
  #   -convert.integers: logical telling the function whether or not to assign columns of integer-like doubles to java data type 'int'
  #
  # Returns:
  #   Writes ODF file to out.filename.
  
  
  # Basic error checking
  ######################################
  if(!is.character(out.filename)) {stop("out.filename must be object of type character")}
  
  out.filename.split <- strsplit(out.filename, split='\\.')
  out.filename.end <- tail(out.filename.split[[1]], n=1)
  if(out.filename.end!='odf' && out.filename.end!='ODF') {stop("out.filename must have extension .odf")}
  
  if(!is.data.frame(input.df)) {
    cat("input.df not data frame object, try to coerce to data frame.\n")
    input.df <- try(as.data.frame(input.df), silent=T)
    if(!is.data.frame(input.df)) {
      write(input.df[1], stderr())
      stop("Failed to coerce input.df to data frame.")
    }
    else {
      cat("input.df successfully coerced to data frame.\n")
      warning("input.df coerced to data frame. Data errors may have been introduced")
    }
  }#End data.frame check
  
  if(!is.list(header.list) && !is.null(header.list)) {
    warning("header.list not list object, ignoring header.list")
    header.list <- NULL
  }
  
  ######################################
  # Internal functions
  ######################################
  double.is.int <- function(x) {x %% 1 == 0} # function to check if number of type double is integer
  get.col.types <- function(df) {
    # gets the column types from data frame and returns character vector of columns in data frame order
    get.type <- function(clmn) {
      ch <- typeof(clmn)
      if(ch == "integer") {return("int")}
      else if(ch == "character") {return("String")}
      else if(ch == "double") {
        if(convert.integers && all(sapply(clmn, double.is.int))) {
          return("int")}
        else {return("double")}
      }
      else if(ch == "logical") {return("boolean")}
    }#End get.type()
    
    out.types <- lapply(df, get.type)
    #     names(out.types) <- NULL
    return(unlist(out.types))
  }#End get.col.types()
  
  ######################################
  # Make list of header line data if no header list provided; else make from provided header list
  ######################################
  header.fin <- list()
  numrow <- dim(input.df)[1]
  numcol <- dim(input.df)[2]
  if(is.null(header.list)) {
    cat("No header.list found. Creating basic header from provided data.\n")
    header.fin[["Model"]] <- "Dataset"
    header.fin[["DataLines"]] <- numrow
    header.fin[["COLUMN_TYPES"]] <- get.col.types(input.df)
    names(header.fin[["COLUMN_TYPES"]]) <- NULL
    if(!is.null(names(input.df))) {
      header.fin[["COLUMN_NAMES"]] <- names(input.df)
    }
    if(!all(row.names(input.df) == as.character(1:header.fin[["DataLines"]]))) {
      header.fin[["ROW_NAMES"]] <- row.names(input.df)
    }
  }
  else {
    # Get Model
    if(is.null(header.list[["Model"]])) {stop("Keyword \"Model\" required, not found in header list.")}
    else {header.list[["Model"]] <- header.list[["Model"]][1]}
    
    # Get DataLines
    if(is.null(header.list[["DataLines"]])) {
      cat("Keyword \"DataLines\" required, not found in header list. Calculating from data frame.\n")
      header.list[["DataLines"]] <- numrow
    }
    else {
      if(dim(input.df)[1] != header.list[["DataLines"]]) {stop("Number of DataLines specified in header.list does not match number of rows in input data frame.")}
      else {header.list[["DataLines"]] <- header.list[["DataLines"]][1]}
    }
    
    # Get COLUMN_TYPES
    if(is.null(header.list[["COLUMN_TYPES"]])) {
      cat("Keyword \"COLUMN_TYPES\" required, not found in header list. Parsing from data frame.\n")
      header.list[["COLUMN_TYPES"]] <- get.col.types(input.df)
    }
    else {
      if(length(header.list[["COLUMN_TYPES"]]) != numcol) {stop("Number of elements in \"COLUMN_TYPES\" does not match number of columns in input data frame.")}
      else {header.list[["COLUMN_TYPES"]] <- header.list[["COLUMN_TYPES"]]}
    }
    names(header.list[["COLUMN_TYPES"]]) <- NULL
    
    # Get COLUMN_NAMES
    if(is.null(header.list[["COLUMN_NAMES"]])) {
      cat("Keyworkd \"COLUMN_NAMES\" not found in header list. Parsing from data frame.\n")
      if(!is.null(names(input.df))) {header.list[["COLUMN_NAMES"]] <- names(input.df)}
      else(cat("No names found in data frame. COLUMN_NAMES will not be included in output header.\n"))
    }
    else {
      if(length(header.list[["COLUMN_NAMES"]]) != numcol) {stop("Number of elements in \"COLUMN_NAMES\" does not match number of columns in input data frame.")}
    }
    
    # Get ROW_NAMES
    if(!is.null(header.list[["ROW_NAMES"]])) {
      if(length(header.list[["ROW_NAMES"]]) != numrow) {stop("Number of elements in \"ROW_NAMES\" does not match number of rows in input data frame.")}
    }
    
    # Check that COLUMN_TYPES are valid java datatypes
    java.types <- c('byte','short','int','long','float','double','boolean','char','String')
    col.types.check <- sapply(header.list[["COLUMN_TYPES"]], function(x) x %in% java.types)
    if(!all(col.types.check)) {stop("Elements of \"COLUMN_TYPES\" must be valid java data types")}
    
    header.fin <- header.list
  }
  
  ######################################
  # Convert logicals to java boolean (e.g., TRUE > true)
  ######################################
  bool.check <- header.fin[["COLUMN_TYPES"]] %in% 'boolean'
  if(any(bool.check)) {
    input.df[,bool.check] <- sapply(input.df[,bool.check], tolower)
  }
  
  ######################################
  # Write ODF file
  ######################################
  con <- file(out.filename, 'w', blocking=F)
  line.1 <- "ODF 1.0"  # pre-header line 1
  line.2 <- paste0("HeaderLines=", length(header.fin))  # pre-header line 2
  
  write(c(line.1, line.2), file=con, append=F, sep='\n')
  
  keywords <- names(header.fin)  # get keywords for headers
  out.lines <- vector('list', length(header.fin))
  
  # Write header lines
  for(i in 1:length(keywords)) {
    temp.line <- NULL
    key <- keywords[[i]]
    vals <- header.fin[[i]]
    if(length(vals) > 1) {
      temp.line <- paste(vals, collapse='\t')
      temp.line <- paste0(key, ':', temp.line)
    }
    else {
      temp.line <- paste0(key, '=', vals)
    }
    write(temp.line, file=con, append=T)
  }
  
  # Write data (tab-delimited)
  write.table(input.df, file=con, append=T, quote=F, sep='\t', row.names=F, col.names=F)
  
  close(con)  # close connection
  return(0)
} # End write.odf

plot.data2odf.file <- function(
  x,  # x coordinates
  y,  # y coordinates
  gs,  # gene set names (indexed to match x, y)
  gs.prefix,  # gene set file/database name
  GSDB, # parsed gene set database
  jaccard.threshold,  # lower jaccard index bound for plotting connector lines
  target.class,  # phenotype/class compared against
  gct.filename,  # name of input gct file
  parsed.cls,  # parsed.cls: cls file pre-parsed by MSIG.ReadPhenFile.2()
  cls.filename,  # cls.filename: name of input cls file
  direction,  # UP or DOWN
  edgepairs, # paired indexes of nodes (gene sets) that will contain an edge
  jcoef # jaccard index matrix for all possible gene set pairs
  ) {
  # Writes constellation.map plot data to an ODF file
  # Authors: Felix Wu
  # Returns:
  #   -Writes plot data to ODF file

  # Build header list for nodes odf
  header.list <- list()
  header.list$COLUMN_NAMES <- c("Gene.Set.Name", "x", "y", "gene.set.size", "url")
  header.list$COLUMN_TYPES <- c("String", "double", "double", "int", "String")
  header.list$Model <- "ConstellationMap"
  header.list[["Plot Data"]] <- "Nodes"
  header.list[["Dataset File"]] <- gct.filename
  header.list[["Class File"]] <- cls.filename
  header.list[["Gene Set Database"]] <- gs.prefix
#   for(i in 1:length(parsed.cls$phen.list)){
#     key <- paste0("Class ", parsed.cls$phen[i])
#     header.list[[key]] <- parsed.cls$phen.list[i]
#   }
  header.list[["Target Class"]] <- target.class
  header.list[["Jaccard Threshold"]] <- jaccard.threshold
  header.list$Direction <- direction
  header.list[["Gene Set Names"]] <- gs
  ext <- max(c(abs(x), abs(y)))
  header.list[["Max Abs Extent"]] <- ext
  
  # Build output data.frame and rest of header file for nodes odf
  df <- data.frame(matrix(nrow=length(gs), ncol=5, NA))
  names(df) <- header.list$COLUMN_NAMES
  df$Gene.Set.Name <- gs
  df$x <- x
  df$y <- y
  all.sets <- GSDB
  gs.slice <- all.sets$gs.names %in% gs
  gs.genes <- all.sets$gs[gs.slice,]
  row.names(gs.genes) <- all.sets$gs.names[gs.slice]
  if(!is.matrix(gs.genes)) {gs.genes <- t(as.matrix(gs.genes))}
  gs.sizes <- all.sets$size.G[gs.slice]
  names(gs.sizes) <- all.sets$gs.names[gs.slice]

  for(i in 1:length(gs)){
    #### IMPORTANT: if uncommenting code below, must load "XML" libray
#     temp.url <- paste0("http://www.broadinstitute.org/gsea/msigdb/cards/", gs[i], ".html")
#     html.str <- toString.XMLNode(htmlTreeParse(temp.url))
#     if(!grepl("Gene Set Not Found", html.str)) {df$url[i] <- temp.url}
    ####
    temp.gs <- gs[i]
    df$url[i] <- paste0("http://www.broadinstitute.org/gsea/msigdb/cards/", temp.gs, ".html")
    temp.gene.vec <- gs.genes[temp.gs,]
    header.list[[gs[i]]] <- temp.gene.vec[!(temp.gene.vec %in% "null")]  # remove "null" values
    df$gene.set.size[match(temp.gs,df$Gene.Set.Name)] <- gs.sizes[temp.gs]
  }
  header.list$DataLines <- length(gs)
  
  # Specify odf output filename for nodes data
  out.filename <- "ConstellationMap.plot.data.nodes.odf"
  
  # Write nodes odf file
  write.odf(input.df=df, out.filename=out.filename, header.list=header.list, convert.integers=F)

  # Coerce edgepairs to matrix if not already matrix
#   if(!is.matrix(edgepairs)) {edgepairs <- t(as.matrix(edgepairs))}

  # Build header list for edges odf
  header.list2 <- header.list
  header.list2$COLUMN_NAMES <- c("Index1", "x1", "y1", "Index2", "x2", "y2", "Jaccard")
  header.list2$COLUMN_TYPES <- c("int", "double", "double", "int", "double", "double", "double")
  header.list2[["Plot Data"]] <- "Edges"
  if(is.null(edgepairs)) {
    header.list2$DataLines <- 0
    df2 <- data.frame(matrix(nrow=header.list2$DataLines, ncol=7, NA))
    names(df2) <- header.list2$COLUMN_NAMES
  } else {
    header.list2$DataLines <- as.numeric(dim(edgepairs)[1])
    # Build output data.frame for edges odf
    df2 <- data.frame(matrix(nrow=header.list2$DataLines, ncol=7, NA))
    names(df2) <- header.list2$COLUMN_NAMES
    for(i in 1:header.list2$DataLines) {
      temp.ind1 <- edgepairs[i,1]
      temp.ind2 <- edgepairs[i,2]
      
      df2$Index1[i] <- as.integer(temp.ind1)
      df2$x1[i] <- x[temp.ind1]
      df2$y1[i] <- y[temp.ind1]
      df2$Index2[i] <- as.integer(temp.ind2)
      df2$x2[i] <- x[temp.ind2]
      df2$y2[i] <- y[temp.ind2]
      df2$Jaccard[i] <- jcoef[temp.ind1, temp.ind2]
    }
  }
  
  

  # Specify odf output filename for edges data
  out.filename2 <- "ConstellationMap.plot.data.edges.odf"

  # Write edges odf file
  write.odf(input.df=df2, out.filename=out.filename2, header.list=header.list2, convert.integers=F)

  return(0)
} # End plot.data2odf.file

# Retrieves a gene set database file (GMT format)
# from GSEA MSigDB. FTP download done by java
# app, which writes db to dest.filename

Get.GeneSets.db <- function(
  javaexec,
  jardir,
  gene.sets.db,
  dest.filename
  ) {
  ConstellationMap.jar <- paste(jardir, "ConstellationMap.jar", sep="")
  edtftpj.jar <- paste(jardir, "edtftpj.jar", sep="")

  if (.Platform$OS.type == "windows") {classpath.sep <- ";"}
  else {classpath.sep <- ":"}

  classpath <- paste(".", ConstellationMap.jar, edtftpj.jar, sep=classpath.sep)

  command <- paste(javaexec, "-cp", classpath , "org.genepattern.modules.ConstellationMap.GeneSetsDownloader", gene.sets.db, dest.filename)

  if ((retval <- system(command)) != 0) {stop("failed to download gene sets db file from GSEA MSigDB; system retval: ", retval)}
} # End Get.GeneSets.db

copy.html.file <- function(
  libdir,
  html.filename
  ) {
  in.filepath <- paste0(libdir, html.filename)
  
  file.copy(from=in.filepath, to=html.filename, copy.mode=T)
} # End copy.html.file
