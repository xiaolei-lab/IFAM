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


## Parameter setting
args_list <- list(
  make_option("--bfile", type = "character", default = NULL,
              help = "INPUT: the prefix of genotype file (plink binary format)", metavar = "character"),
  make_option("--pheno", type = "character", default = NULL,
              help = "INPUT: the filename of phenotype file", metavar = "character"),
  make_option("--anno_folder", type = "character", default = NULL,
              help = "INPUT: the folder where annotation files are stored", metavar = "character"),
  make_option("--GRMs_folder", type = "character", default = NULL,
              help = "INPUT: the folder where GRM of each annotation are stored", metavar = "character"),          
  make_option("--outPath", type = "character", default = NULL,
              help = "INPUT: the path of output", metavar = "character"), 
  make_option("--output_prefix", type = "character", default = "IFAM",
              help = "INPUT: the prefix of output (default:IFAM)", metavar = "character"), 
  make_option("--Pruning", type = "logical", default = TRUE,
              help = "INPUT: TRUE represents to perform LD pruning (default: TRUE)"),
  make_option("--indep_pairwise", type = "character", default = NULL,
              help = "INPUT: the paremeters of LD pruning, please see more details about these files at PLINK user manual", metavar="character"),  
  make_option("--plink", type = "character", default = NULL,
              help = "INPUT: the perfix of Plink software", metavar = "character"),
  make_option("--pheno_pos", type = "integer", default = "2",
              help = "INPUT: the position of the analyzed phenotype in columns of phenotype file (default:2)", 
              metavar = "character"),
  make_option("--VCmethod", type = "character", default = "AI",
              help = "INPUT: the algorithms (AI, EM, EMAI, HE, HI) for variance component estimation (default: AI method)", 
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

if (is.null(opt$outPath)){
  cat(paste0("ERROR: the path of output does not exist! Please check! \n"))
  q()
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

if (!is.null(opt$indep_pairwise)){
  indep_pairwise_str <- unlist(strsplit(opt$indep_pairwise, ","))
  if (length(indep_pairwise_str) == 3){
    indep_pairwise <- paste(indep_pairwise_str[1], indep_pairwise_str[2], sep=" ")
    indep_pairwise <- paste(indep_pairwise, indep_pairwise_str[3], sep=" ")
  }else{
    cat(paste0("ERROR: ", opt$indep_pairwise, " is wrong! Please check!\n"))
    q()
  }
}


## phenotype file 
phe_header <- unlist(strsplit(readLines(opt$pheno, n=1), "\t"))
trait_name <- phe_header[opt$pheno_pos]
cat(paste("  ", "\n"))
cat(paste("Analysis Trait: ", trait_name, sep=""), "\n")

# genome map file
map <- read.delim(bfile_str[2], head=FALSE)
map_SNP <- map[,2]

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


## Estimating variance components 
cat(paste("  ", "\n"))
cat(paste(paste(rep("-", 27), collapse=""), " Estimating variance components ", paste(rep("-", 26), collapse=""), sep=""), "\n")

if (is.null(opt$GRMs_folder)){
  cat(paste("  ", "\n"))
  cat(paste(paste(rep("-", 11), collapse=""), " Constructing GRMs for each annotation ", paste(rep("-", 10), collapse=""), sep=""), "\n")
  GRMs <- make_GRM(set_snplists=anno_all, bfile=opt$bfile, weight=NULL, outPath=opt$outPath, output_prefix=opt$output_prefix, thread=opt$thread)
} else {
  cat(paste(paste(rep("-", 11), collapse=""), " The GRM of each annotation was provided!", paste(rep("-", 10), collapse=""), sep=""), "\n")
  GRMs <- paste(GRM_str2, collapse=",") 
  print(GRM_str2)
}

VC_cmd <- paste0("hiblup --single-trait --threads ", opt$thread, " --pheno ", opt$pheno, " --pheno-pos ", opt$pheno_pos, 
" --xrm ", GRMs, " --vc-method ", opt$VCmethod, " --out ", opt$outPath, opt$output_prefix, "_", trait_name, "_vc")
system(VC_cmd, ignore.stdout=TRUE)
rm_cmd <- paste0("rm ", opt$outPath, opt$output_prefix, "_", trait_name, "_vc.beta; rm ", 
                opt$outPath, opt$output_prefix, "_", trait_name, "_vc.rand")
system(rm_cmd)
vars <- read.delim(paste0(opt$outPath, opt$output_prefix, "_", trait_name, "_vc.vars"), head=TRUE)
vars$Item <- gsub(paste0(opt$output_prefix, "."), "", vars$Item)
vars$Item <- gsub(".GA", "", vars$Item)
vars <- vars[vars$Item!="e", c(1:2, 4)]
print(vars)

cat(paste("  ", "\n"))
cat(paste(paste(rep("-", 27), collapse=""), " Calcuating the per-SNP variance componts and per-SNP heritability ", paste(rep("-", 26), collapse=""), sep=""), "\n")
#bim <- read.delim(bfile_str[2], head=FALSE)
if (is.null(opt$GRMs_folder)){
  anno_name <- rownames(anno_all)
} else {
  anno_name <- rownames(GRM_str2)
}
if (opt$Pruning){
  cat(paste(paste(rep("-", 11), collapse=""), " Prunning SNPs within annotational region ", paste(rep("-", 10), collapse=""), sep=""), "\n")
  EV_results <- list(Annotations=NULL, Number_annotated_SNPs=NULL, Number_prunning_annotated_SNPs=NULL, Per_SNP_vc=NULL, Per_SNP_h2=NULL)
  for (i in 1:nrow(anno_all)){
    anno <- read.delim(anno_all[i, 1], head=FALSE)
    EV_results$Annotations[i] <- anno_name[i]
    EV_results$Number_annotated_SNPs[i] <- length(intersect(anno[,1], map_SNP))
    prunnig_cmd <- paste0(opt$plink, " --bfile ", opt$bfile, " --extract ", anno_all[i, 1], " --indep-pairwise ", 
                          indep_pairwise, " --out ", opt$outPath, opt$output_prefix, "_", EV_results$Annotations[i])
    system(prunnig_cmd, ignore.stdout=TRUE)
    prunning_SNPlist <- read.delim(paste0(opt$outPath, opt$output_prefix, "_", EV_results$Annotations[i], ".prune.in"), head=FALSE)
    EV_results$Number_prunning_annotated_SNPs[i] <- nrow(prunning_SNPlist)
    EV_results$Per_SNP_vc[i] <- vars[anno_name[i]==vars$Item, 2] / EV_results$Number_prunning_annotated_SNPs[i]
    EV_results$Per_SNP_h2[i] <- vars[anno_name[i]==vars$Item, 3] / EV_results$Number_prunning_annotated_SNPs[i]
  }
} else {
  EV_results <- list(Annotations=NULL, Number_annotated_SNPs=NULL, Per_SNP_vc=NULL, Per_SNP_h2=NULL)
  for (i in 1:nrow(anno_all)){
    anno <- read.delim(anno_all[i, 1], head=FALSE)
    EV_results$Annotations[i] <- anno_name[i]
    EV_results$Number_annotated_SNPs[i] <- length(intersect(anno[,1], map_SNP))
    EV_results$Per_SNP_vc[i] <- vars[anno_name[i]==vars$Item, 2] / EV_results$Number_annotated_SNPs[i]
    EV_results$Per_SNP_h2[i] <- vars[anno_name[i]==vars$Item, 3] / EV_results$Number_annotated_SNPs[i]
  }
}
print(do.call(cbind, EV_results))

## handling temporary files
if (!opt$tmp_files){
  rm_cmd <- paste0("rm ", opt$outPath, output_prefix=opt$output_prefix, "**")
  system(rm_cmd)
}
write.table(EV_results, paste0(opt$outPath, opt$output_prefix, ".EA.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

cat(paste0("End time: ", Sys.time(), "\n"))
