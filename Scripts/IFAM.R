#! /usr/bin/env Rscript

########################################################################
# Integrating Functional Annotation information by the genomic BLUP    #
# with Multiple random effects (FIAM)                                  #
# Copyright (C) 2024  Zhenshuang Tang, Lilin Yin and Xiaolei Liu       #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http://www.gnu.org/licenses/>. #
########################################################################

IFAM.version <- function()
{
cat(paste(paste(rep("-", 27), collapse=""), " Welcome to IFAM ", paste(rep("-", 26), collapse=""), sep=""), "\n")
cat("        Integrating Functional Annotation information by genomic        \n")
cat("              BLUP model with Multiple random effects                   \n")
cat("               ____  _____                                        \n")
cat("                ||   ||      ^       ^  ^                 \n")
cat("                ||   ||___  /|\\     /\\  /\\                  \n")
cat("                ||   ||    //_\\\\   //\\\\//\\\\                 \n")
cat("               _||_  ||   //   \\\\ //  \\/  \\\\  Version: 1.0.0\n")
cat("        Designed and Maintained by ZHneshuang Tang, Lilin Yin, and      \n") 
cat("        Xiaolei Liu                                                     \n")
cat("        Contact: xiaoleiliu@mail.hzau.edu.cn                            \n")
cat("        Website: https://github.com/Zhenshuang/IFAM/                    \n")
cat(paste(rep("-", 70), collapse=""), "\n")
}
IFAM.version()

cat(paste0("Start time: ", Sys.time(), "\n"))
options(stringsAsFactors=F)
library("optparse")

## making GRMs for some sets of SNPs lists
## args
## set_snplists: a matrix list for filenames of each annotation region
## bfile: the prefix of genotype file (plink binary format)
## weight: optional, the filename of SNP weight file (default: NULL)
## outPath: the output path
## output_prefix: the prefix of output
## thread: the number of threads (default: 1)
make_GRM <- function(set_snplists=NULL, bfile=NULL, weight=NULL, outPath=NULL, output_prefix=NULL, thread=1){
  GRM_list <- c()
  for(i in 1:nrow(set_snplists)){    
    snplist <- set_snplists[i,]
    GRM_name <- rownames(set_snplists)[i]
    GRM_name <- paste0(outPath, output_prefix, ".", GRM_name) 
    if (!is.null(weight)){
      makeGRM_cmd <- paste0("hiblup --bfile ", bfile, " --extract ", snplist, " --make-xrm --snp-weight ", weight, " --threads ", thread, " --out ", GRM_name, ".w")
      if (!file.exists(paste0(GRM_name, ".w.GA.bin"))){
        system(makeGRM_cmd, ignore.stdout=TRUE)
      }
      GRM_list <- c(GRM_list, paste0(GRM_name, ".w.GA"))
    } else{
      makeGRM_cmd <- paste0("hiblup --bfile ", bfile, " --extract ", snplist, " --make-xrm --threads ", thread, " --out ", GRM_name)
      if (!file.exists(paste0(GRM_name, ".GA.bin"))){
        system(makeGRM_cmd, ignore.stdout=TRUE)
      }
      GRM_list <- c(GRM_list, paste0(GRM_name, ".GA"))
    } 
  }
  GRM_list <- paste(GRM_list, collapse=",") 
  return(GRM_list)
}

## Optimizing the number of random effects based on the values of variance components
## args
## vc_list: the values of variance components for each annotations
## randomMax: the maximium number of random effects in the model (default:5)
opt_random <- function(vc_list=NULL, randomMax=5){
  var_max <- max(vc_list[,2])
  random_list <- list()
  for (j in 1:randomMax){
    index <- which(vc_list[,2] <= var_max/10^(j-1) & vc_list[,2] > var_max/10^j)
    if (length(index)==0){
        random_list[[j]] <- NULL
    } else {
        random_list[[j]] <- vc_list[index, 1]
    }
  }
  random_list <- random_list[vapply(random_list, Negate(is.null), NA)]
  if(length(unlist(random_list)) != nrow(vars)){
    cat(paste0("----- There are some annotations were romoved because the estimates were too small -----", "\n"))
    print(setdiff(vars[,1], unlist(random_list)))
  }
  return(random_list)
}

## Parameter setting
args_list <- list(
  make_option("--bfile", type = "character", default = NULL,
              help = "INPUT: the prefix of genotype file (plink binary format)", metavar = "character"),
  make_option("--pheno", type = "character", default = NULL,
              help = "INPUT: the filename of phenotype file", metavar = "character"),
  make_option("--anno_folder", type = "character", default = NULL,
              help = "INPUT: the folder where annotation files are stored", metavar = "character"),
  make_option("--anno_spec", type = "character", default = NULL,
              help = "INPUT: the filename of special annotation file, this file doesn't participate 
              in the optimization of random effects in the model", metavar = "character"),
  make_option("--GRMs_folder", type = "character", default = NULL,
              help = "INPUT: the folder where GRM of each annotation are stored", metavar = "character"),             
  make_option("--weight", type = "character", default = NULL,
              help = "INPUT: the filename of SNP weight file", metavar = "character"),
  make_option("--VCfile", type = "character", default = NULL,
              help = "INPUT: the filename of the results of variance component estimation 
              for all anotations", metavar = "character"),
  make_option("--outPath", type = "character", default = NULL,
              help = "INPUT: the path of output", 
              metavar = "character"), 
  make_option("--output_prefix", type = "character", default = "IFAM",
              help = "INPUT: the prefix of output (default:IFAM)", 
              metavar = "character"), 
  make_option("--pheno_pos", type = "integer", default = "2",
              help = "INPUT: the position of the analyzed phenotype in columns of phenotype file (default:2)", 
              metavar = "character"),
  make_option("--randomMax", type = "integer", default = "5",
              help = "INPUT: the maximium number of random effects in the model (default:5)", 
              metavar = "character"),
  make_option("--VCmethod", type = "character", default = "AI",
              help = "INPUT: the method of variance component estimation (default: AI method)", 
              metavar = "character"),
  make_option("--thread", type = "integer", default = "1",
              help = "INPUT: the number of threads (default: 1)", 
              metavar = "character"),
  make_option("--tmp_files", type = "logical", default = TRUE,
              help = "INPUT: Whether temporary files are stored (default: TRUE)")
)

opt_parser <- OptionParser(option_list=args_list)
opt <- parse_args(opt_parser)

## check the options
bfile_str <- paste0(opt$bfile, c(".bed", ".bim", ".fam"))
if (!file.exists(bfile_str[1])){
  cat(paste0("ERROR: ", opt$bfile, " does not exist! Please check!\n"))
  q()
} else {
  cat(paste0("The number of SNPs witnin genotype file: ", "\n"))
  cat(paste0(system(paste0("wc -l ", bfile_str[2]), intern=TRUE), "\n"))
  cat(paste0("The number of individuals within genotype file: ", "\n"))
  cat(paste0(system(paste0("wc -l ", bfile_str[3]), intern=TRUE),"\n"))
}

if (!file.exists(opt$pheno)){
  cat(paste0("ERROR: ", opt$pheno, " does not exist! Please check!\n"))
  q()
} else {
  cat(paste0("The number of individuals within phenotype file: ", "\n"))
  cat(paste0(system(paste0("wc -l ", opt$pheno), intern=TRUE), "\n"))
}

if (!file.exists(opt$anno_folder)){
  cat(paste0("ERROR: ", opt$anno_folder, " does not exist! Please check!\n"))
  q()
} else {
  anno_str <- list.files(path=opt$anno_folder, pattern = '*.txt', full.names=TRUE, recursive=FALSE)
  cat(paste0("Analysis ", length(anno_str), " annotations\n"))
  anno_str_names <- c()
  cat(paste0("The number of records within each annotation file: ", "\n"))
  for (i in 1:length(anno_str)){
    if (!file.exists(anno_str[i])){
      cat(paste0("ERROR: ", anno_str[i], " does not exist! Please check!\n"))
      q()
    } else {
      cat(paste0(system(paste0("wc -l ", anno_str[i]), intern=TRUE), "\n"))
      anno_str_str <- unlist(strsplit(anno_str[i], "/"))
      anno_str_names <- c(anno_str_names, gsub(".txt", "", anno_str_str[length(anno_str_str)]))
    }
  }
  anno_all <- matrix(anno_str, ncol=1, dimnames = list(anno_str_names, "Annotations"))
}

if (is.null(opt$outPath)){
  cat(paste0("ERROR: the path of output does not exist! Please check! \n"))
  q()
}

if (!is.null(opt$anno_spec)){
  cat(paste0("the annotation: ", opt$anno_spec,", which doesn't participate in the optimization of random effects in the model!\n"))
  anno_spec_str <- unlist(strsplit(opt$anno_spec, ","))
  anno_spec_names <- c()
  if (length(anno_spec_str) > 1){
    for (i in 1:length(anno_spec_str)){
      if (!file.exists(paste0(opt$anno_folder, "/", anno_spec_str[i], ".txt"))){
        cat(paste0("ERROR: ", anno_spec_str[i], " does not exist! Please check!\n"))
        q()
      }
    anno_spec_str_str <- unlist(strsplit(anno_spec_str[i], "/"))
    anno_spec_names <- c(anno_spec_names, gsub(".txt", "", anno_spec_str_str[length(anno_spec_str_str)]))
    }
  } else {
    anno_spec_file <- paste0(opt$anno_folder, "/", anno_spec_str, ".txt")
    if (!file.exists(anno_spec_file)){
      cat(paste0("ERROR: ", opt$anno_spec, " does not exist! Please check!\n"))
      q()
    }
    anno_spec_str_str <- unlist(strsplit(anno_spec_str, "/"))
    anno_spec_names <- c(anno_spec_names, gsub(".txt", "", anno_spec_str_str[length(anno_spec_str_str)]))
  }
  anno_spec_str2 <- matrix(anno_spec_file, ncol=1, dimnames = list(anno_spec_names, "Annotations"))
}

if (!is.null(opt$GRMs_folder)){
  GRM_bin <- list.files(path=opt$GRMs_folder, pattern = '.GA.bin', full.names=TRUE, recursive=FALSE)
  GRM_id <- list.files(path=opt$GRMs_folder, pattern = '.GA.id', full.names=TRUE, recursive=FALSE)
  if (length(anno_str) != length(GRM_bin) | length(anno_str) != length(GRM_id)) {
    cat(paste0("ERROR: the number of GRMs and all annotations must be the same! Please check!\n"))
    q()
  }
  GRM_str <- gsub(".bin", "", GRM_bin)
  GRM_names <- c()
  for (i in 1:length(GRM_bin)){
    GRM_str_str <- unlist(strsplit(GRM_str[i], "/"))
    GRM_names <- c(GRM_names, gsub(".GA", "", GRM_str_str[length(GRM_str_str)]))
  }
  GRM_str2 <- matrix(GRM_str, ncol=1, dimnames = list(GRM_names, "GRMs"))
}

if (!is.null(opt$weight)){
  if (!file.exists(opt$weight)){
    cat(paste0("ERROR: ", opt$weight, " does not exist! Please check!\n"))
    q()
  }
}

if (!is.null(opt$VCfile)){
  if (!file.exists(opt$VCfile)){
    cat(paste0("ERROR: ", opt$VCfile, " does not exist! Please check!\n"))
    q()
  }
}

## phenotype file 
phe_header <- unlist(strsplit(readLines(opt$pheno, n=1), "\t"))
trait_name <- phe_header[opt$pheno_pos]
cat(paste(paste(rep(" ", 5), collapse=""), " Analysis Trait: ", trait_name, sep=""), "\n")

# genome map file
map <- read.delim(bfile_str[2], head=FALSE)
map_SNP <- map[,2]
cat(paste(paste(rep(" ", 5), collapse=""), " The genome map file has ", nrow(map), " SNPs!", sep=""), "\n")

# summary of annotations
cat(paste(paste(rep("-", 27), collapse=""), " The summary information about annotations ", paste(rep("-", 26), collapse=""), sep=""), "\n")
anno_list <- list(Annotations=NULL, Number_annotaions=NULL, Number_overlap_anno_map=NULL, Anno_cover_percent=NULL)
SNP_list <- list()
for (i in 1:nrow(anno_all)){
  anno <- read.delim(anno_all[i,], head=FALSE)
  anno_list$Annotations[i] <- rownames(anno_all)[i]
  anno_list$Number_annotaions[i] <- nrow(anno)
  anno_list$Number_overlap_anno_map[i] <- length(intersect(anno[,1], map_SNP))
  anno_list$Anno_cover_percent[i] <- anno_list$Number_overlap_anno_map[i]/length(map_SNP)
  SNP_list[[i]] <- anno
}
anno_list$Anno_cover_percent <- paste(round(100*anno_list$Anno_cover_percent, 4), "%", sep="")
print(do.call(cbind, anno_list))
write.table(anno_list, paste0(opt$outPath, opt$output_prefix, ".annotaion.information.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

overlap_matrix <- matrix(NA, nrow=nrow(anno_all), ncol=nrow(anno_all), dimnames=list(rownames(anno_all), rownames(anno_all)))
for (i in 1:nrow(anno_all)){
  for (j in 1:nrow(anno_all)){
    if (i == j){
      overlap_matrix[i,j] <- "--"
    }else if (i < j){
      overlap_matrix[i,j] <- length(intersect(unlist(SNP_list[[i]]), unlist(SNP_list[[j]])))
    }
  }
}
print(overlap_matrix)
overlap_matrix[is.na(overlap_matrix)] <- " "
write.table(overlap_matrix, paste0(opt$outPath, opt$output_prefix, ".overlaped.SNPs.between.anntations.txt"), row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")


start <- proc.time()
cat(paste(paste(rep("-", 27), collapse=""), " Optimizing random effects ", paste(rep("-", 26), collapse=""), sep=""), "\n")

## Optimizing random effects in two scenarios
if (is.null(opt$anno_spec)){
  cat(paste(paste(rep(" ", 5), collapse=""), " There are no special annotations, so optimizing random effects using all annotations ", sep=""), "\n")
  cat(paste(paste(rep(" ", 5), collapse=""), " Estimating variance components using multiple random effects model ", sep=""), "\n")
  if (!is.null(opt$VCfile)){
    cat(paste(paste(rep(" ", 5), collapse=""), " The variance components of each annotation were provided ", sep=""), "\n")
    vars <- read.delim(opt$VCfile, head=TRUE)
    print(vars)
  } else{
    if (is.null(opt$GRMs_folder)){
      cat(paste(paste(rep(" ", 5), collapse=""), " Constructing GRMs for each annotation ", sep=""), "\n")
      GRMs <- make_GRM(set_snplists=anno_all, bfile=opt$bfile, weight=NULL, outPath=opt$outPath, output_prefix=opt$output_prefix, thread=opt$thread)
    } else{
      cat(paste(paste(rep(" ", 5), collapse=""), " The GRM of each annotation was provided ", sep=""), "\n")
      GRMs <- paste(GRM_str2, collapse=",") 
      print(GRM_str2)
    }
    cat(paste(paste(rep(" ", 5), collapse=""), " Estimating variance components ", sep=""), "\n")
    VC_cmd <- paste0("hiblup --single-trait --threads ", opt$thread, " --pheno ", opt$pheno, " --pheno-pos ", opt$pheno_pos, 
    " --xrm ", GRMs, " --vc-method ", opt$VCmethod, " --out ", opt$outPath, opt$output_prefix, "_", trait_name, "_vc")
    system(VC_cmd, ignore.stdout=TRUE)
    rm_cmd <- paste0("rm ", opt$outPath, opt$output_prefix, "_", trait_name, "_vc.beta; rm ", 
                      opt$outPath, opt$output_prefix, "_", trait_name, "_vc.rand")
    system(rm_cmd)
    vars <- read.delim(paste0(opt$outPath, opt$output_prefix, "_", trait_name, "_vc.vars"), head=TRUE)
    vars$Item <- gsub(paste0(opt$output_prefix, "."), "", vars$Item)
    vars$Item <- gsub(".GA", "", vars$Item)
    vars <- vars[vars$Item!="e", 1:2]
    print(vars)
  }      
  random_list <- opt_random(vc_list=vars, randomMax=opt$randomMax)
}else{
  cat(paste(paste(rep(" ", 5), collapse=""), " There have special annotations, which doesn't participate in the optimization of random effects ", sep=""), "\n")
  index <- match(setdiff(rownames(anno_all), rownames(anno_spec_str2)), rownames(anno_all))
  anno_remaining <- matrix(anno_all[index,], ncol=1)
  rownames(anno_remaining) <- rownames(anno_all)[index]
  if (!is.null(opt$VCfile)){
    cat(paste(paste(rep(" ", 5), collapse=""), " The variance components of each annotation were provided ", sep=""), "\n")
    vars <- read.delim(opt$VCfile, head=TRUE)
    if(nrow(vars) != nrow(anno_remaining)) {
      cat(paste0("ERROR: the number of variance components is wrong! Please check!\n"))
      q()
    }
    print(vars)
  }else{
    if (is.null(opt$GRMs_folder)){
      cat(paste(paste(rep(" ", 5), collapse=""), " Constructing GRMs for each annotation ", sep=""), "\n")
      GRMs <- make_GRM(set_snplists=anno_remaining, bfile=opt$bfile, weight=NULL, outPath=opt$outPath, output_prefix=opt$output_prefix, thread=opt$thread)
    } else{
      cat(paste(paste(rep(" ", 5), collapse=""), " The GRM of each annotation was provided ", sep=""), "\n")
      GRM_remaining <- matrix(GRM_str2[index,], ncol=1)
      rownames(GRM_remaining) <- rownames(GRM_str2)[index]
      GRMs <- paste(GRM_remaining, collapse=",") 
      print(GRM_remaining)
    }
    cat(paste(paste(rep(" ", 5), collapse=""), " Estimating variance components for remaining annotations ", sep=""), "\n")
    VC_cmd <- paste0("hiblup --single-trait --threads ", opt$thread, " --pheno ", opt$pheno, " --pheno-pos ", opt$pheno_pos, 
    " --xrm ", GRMs, " --vc-method ", opt$VCmethod, " --out ", opt$outPath, opt$output_prefix, "_", trait_name, "_vc")
    system(VC_cmd, ignore.stdout=TRUE)
    rm_cmd <- paste0("rm ", opt$outPath, opt$output_prefix, "_", trait_name, "_vc.beta; rm ", 
                      opt$outPath, opt$output_prefix, "_", trait_name, "_vc.rand")
    system(rm_cmd)
    vars <- read.delim(paste0(opt$outPath, opt$output_prefix, "_", trait_name, "_vc.vars"), head=TRUE)
    vars$Item <- gsub(paste0(opt$output_prefix, "."), "", vars$Item)
    vars$Item <- gsub(".GA", "", vars$Item)
    vars <- vars[vars$Item!="e", 1:2]
    print(vars)
  }
  random_list <- opt_random(vc_list=vars, randomMax=opt$randomMax - nrow(anno_spec_str2))
  for(i in 1:nrow(anno_spec_str2)){
    anno_spec_str2_s <- anno_spec_str2[i]
    anno_spec_str2_s <- unlist(strsplit(anno_spec_str2_s, "/"))
    anno_spec_str2_s <- gsub(".txt", "", anno_spec_str2_s[length(anno_spec_str2_s)])
    random_list[[length(random_list) + i]] <- anno_spec_str2_s
  }
}
cat(paste(paste(rep(" ", 5), collapse=""), " The number of random effects after optimization: ", length(random_list), sep=""), "\n")
end1 <- proc.time()
cat(paste(paste(rep(" ", 5), collapse=""), " Optimization time: ",  end1[3]-start[3], sep=""), "s \n")


cat(paste(paste(rep("-", 10), collapse=""), " Predicting additive genetic values using multiple random effects model ", paste(rep("-", 10), collapse=""), sep=""), "\n")
random_set_snplist_files <- list()
random_set_snplist_files_names <- c()
for (i in 1:length(random_list)){ 
  random_set <- random_list[[i]]
  random_set_snplist <- list()
  random_set_snplist_files[[i]] <- opt$outPath
  for(j in 1:length(random_set)){
    random_set_snplist_files[[i]] <- paste0(random_set_snplist_files[[i]], ".", random_set[j])
    random_set_file <- anno_all[rownames(anno_all)==random_set[j]]
    random_set_snplist[[j]] <- read.delim(random_set_file, head=FALSE)
  }
  random_set_snplist <- unique(unlist(random_set_snplist))
  cat(paste0("Random effects set ", i, " including ", length(random_set_snplist), " SNPs!\n"))
  random_set_snplist_files[[i]] <- paste0(random_set_snplist_files[[i]], ".txt")
  if (!file.exists(random_set_snplist_files[[i]])){
    write.table(random_set_snplist, random_set_snplist_files[[i]], row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  }
  temp <- unlist(strsplit(random_set_snplist_files[[i]], "/"))
  temp <- gsub(".txt", "", temp[length(temp)])
  temp <- unlist(strsplit(temp, "\\."))[-1]
  random_set_snplist_files_names <- c(random_set_snplist_files_names, paste(temp, collapse=".") )
}
random_set_snplist_files <- do.call(rbind, random_set_snplist_files)
rownames(random_set_snplist_files) <- random_set_snplist_files_names
cat(paste(paste(rep(" ", 5), collapse=""), " The list of random effects ", sep=""), "\n")
print(random_set_snplist_files)
cat(paste(paste(rep(" ", 5), collapse=""), " Constructiing GRM for each random effect ", sep=""), "\n")
if (is.null(opt$weight)){
  random_GRMs <- make_GRM(set_snplists=random_set_snplist_files, bfile=opt$bfile, weight=NULL, outPath=opt$outPath, output_prefix=opt$output_prefix, thread=opt$thread)
} else {
  random_GRMs <- make_GRM(set_snplists=random_set_snplist_files, bfile=opt$bfile, weight=opt$weight, outPath=opt$outPath, output_prefix=opt$output_prefix, thread=opt$thread)
}
cat(paste(paste(rep(" ", 5), collapse=""), " Running Genomic BLUP model with multiple random effects ", sep=""), "\n")
multipleBLUP_cmd <- paste0("hiblup --single-trait --threads ",  opt$thread, " --pheno ", opt$pheno, " --pheno-pos ", opt$pheno_pos, 
" --xrm ", random_GRMs, " --vc-method ", opt$VCmethod, " --out ", opt$outPath, output_prefix=opt$output_prefix, "_", trait_name)
system(multipleBLUP_cmd, ignore.stdout=TRUE)
rand <- read.delim(paste0(opt$outPath, output_prefix=opt$output_prefix, "_", trait_name, ".rand"), head=TRUE)
genetic_values <- apply(rand[2:(length(random_set_snplist_files)+1)], 1, sum)
output <- cbind(rand[, 1:(length(random_set_snplist_files)+1)], genetic_values, rand[,ncol(rand)])
colnames(output)[ncol(output)] <- "residuals"

## handling temporary files
if (!opt$tmp_files){
  rm_cmd <- paste0("rm ", opt$outPath, output_prefix=opt$output_prefix, "**")
  system(rm_cmd)
}

write.table(output, paste0(opt$outPath, output_prefix=opt$output_prefix, "_", trait_name, ".rand"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

cat(paste(paste(rep("-", 27), collapse=""), " IFAM is end! ", paste(rep("-", 26), collapse=""), sep=""), "\n")

end2 <- proc.time()
cat(paste(paste(rep(" ", 5), collapse=""), " Estimating genetic values time: ",  end2[3]-end1[3], sep=""), "s \n")
cat(paste0("End time: ", Sys.time(), "\n"))
