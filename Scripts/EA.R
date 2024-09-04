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

cat(paste0("Start time: ", Sys.time(), "\n"))
options(stringsAsFactors=F)
library("optparse")

## making GRMs for some sets of SNPs lists
## args
## set_snplists: sets of SNPs list for each annotation region
## bfile: the prefix of genotype file (plink binary format)
## weight: optional, the filename of SNP weight file (default: NULL)
## outPath: the output path
## thread: the number of threads (default: 1)
make_GRM <- function(set_snplists=NULL, bfile=NULL, weight=NULL, outPath=NULL, thread=1){
  set_snplists_str <- unlist(strsplit(set_snplists, ","))
  GRM_list <- c()
  for(i in 1:length(set_snplists_str)){    
    snplist <- set_snplists_str[i]
    GRM_name <- unlist(strsplit(snplist, "/"))
    GRM_name <- gsub(".txt", "", GRM_name[length(GRM_name)])
    GRM_name <- paste0(outPath, GRM_name) 
    if (!is.null(opt$weight)){
      makeGRM_cmd <- paste0("hiblup --bfile ", bfile, " --extract ", snplist, " --make-xrm --snp-weight ", weight, " --threads ", thread, " --out ", GRM_name, ".w")
      if (!file.exists(paste0(GRM_name, ".w.GA.bin"))){
        system(makeGRM_cmd)
      }
      GRM_list <- c(GRM_list, paste0(GRM_name, ".w.GA"))
    } else{
      makeGRM_cmd <- paste0("hiblup --bfile ", bfile, " --extract ", snplist, " --make-xrm --threads ", thread, " --out ", GRM_name)
      if (!file.exists(paste0(GRM_name, ".GA.bin"))){
        system(makeGRM_cmd)
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
  make_option("--anno", type = "character", default = NULL,
              help = "INPUT: the filename of annotation files", metavar = "character"),
  make_option("--anno_GRM", type = "character", default = NULL,
              help = "INPUT: the filename of GRM for each annotation", metavar = "character"),             
  make_option("--outPath", type="character", default=NULL,
              help="INPUT: the output path", metavar="character"),
  make_option("--output_prefix", type = "character", default = "IFAM",
              help = "INPUT: the prefix of output (default:IFAM)", 
              metavar = "character"), 
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
}
if (!file.exists(opt$pheno)){
  cat(paste0("ERROR: ", opt$pheno, " does not exist! Please check!\n"))
  q()
}
anno_str <- unlist(strsplit(opt$anno, ","))
for (i in 1:length(anno_str)){
  if (!file.exists(anno_str[i])){
    cat(paste0("ERROR: ", anno_str[i], " does not exist! Please check!\n"))
    q()
  }
}
if (!is.null(opt$anno_GRM)){
  anno_GRM_str <- unlist(strsplit(opt$anno_GRM, ","))
  for (i in 1:length(anno_GRM_str)){
    anno_GRM_str_str <- paste0(anno_GRM_str[i], c(".bin", ".id"))
    for(j in 1:length(anno_GRM_str_str)){
      if (!file.exists(anno_GRM_str_str[j])){
        cat(paste0("ERROR: ", anno_GRM_str_str[j], " does not exist! Please check!\n"))
        q()
      }
    }
  }
}

## phenotype file 
phe_header <- unlist(strsplit(readLines(opt$pheno, n=1), "\t"))
trait_name <- phe_header[opt$pheno_pos]
cat(paste0("Analysis Trait: ", trait_name, "\n"))

## annotation files
cat(paste0("Analysis ", length(anno_str), " annotations: \n"))
anno_str_names <- c()
for (i in 1:length(anno_str)){
  anno_str_str <- unlist(strsplit(anno_str[i], "/"))
  anno_str_names <- c(anno_str_names, gsub(".txt", "", anno_str_str[length(anno_str_str)]))
  cat(paste0(anno_str[i], "\n"))
}
anno_str2 <- matrix(anno_str, ncol=1, dimnames = list(anno_str_names, "Annotations"))

## GRMs files
if (!is.null(opt$anno_GRM)){
  anno_GRM_names <- c()
  for (i in 1:length(anno_GRM_str)){
    anno_GRM_str_str <- unlist(strsplit(anno_GRM_str[i], "/"))
    anno_GRM_names <- c(anno_GRM_names, gsub(".GA", "", anno_GRM_str_str[length(anno_GRM_str_str)]))
  }
  anno_GRM_str2 <- matrix(anno_GRM_str, ncol=1, dimnames = list(anno_GRM_names, "GRMs"))
}

## indep_pairwis
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

## making GRMs for each annotation for variance component estimation
if (is.null(opt$anno_GRM)){
  cat(paste0("------------ Constructing GRMs ------------\n"))
  GRMs <- make_GRM(set_snplists=opt$anno, bfile=opt$bfile, weight=NULL, outPath=opt$outPath, thread=opt$thread)
} else {
  GRMs <- opt$anno_GRM
}

## estimating variance components 
cat(paste0("------------ estimating variance components ------------\n"))
VC_cmd <- paste0("hiblup --single-trait --threads ", opt$thread, " --pheno ", opt$pheno, " --pheno-pos ", opt$pheno_pos, 
" --xrm ", GRMs, " --vc-method ", opt$VCmethod, " --out ", opt$outPath, opt$output_prefix, "_", trait_name, "_vc")
system(VC_cmd)
rm_cmd <- paste0("rm ", opt$outPath, opt$output_prefix, "_", trait_name, "_vc.beta; rm ", 
                opt$outPath, opt$output_prefix, "_", trait_name, "_vc.rand")
system(rm_cmd)
vars <- read.delim(paste0(opt$outPath, opt$output_prefix, "_", trait_name, "_vc.vars"), head=TRUE)
vars$Item <- gsub(".GA", "", vars$Item)
vars <- vars[vars$Item!="e", c(1:2,4)]

## calcuating the per-SNP variance componts and per-SNP heritability
bim <- read.delim(bfile_str[2], head=FALSE)
if (is.null(opt$anno_GRM)){
  anno_name <- rownames(anno_str2)
} else {
  anno_name <- rownames(anno_GRM_str2)
}
if (opt$Pruning){
  EV_results <- list(Annotations=NULL, Number_annotated_SNPs=NULL, Number_prunning_annotated_SNPs=NULL, Per_SNP_vc=NULL, Per_SNP_h2=NULL)
  for (i in 1:length(anno_str)){
    anno <- read.delim(anno_str[i], head=FALSE)
    EV_results$Annotations[i] <- anno_name[i]
    EV_results$Number_annotated_SNPs[i] <- length(intersect(anno[,1], bim[,2]))
    prunnig_cmd <- paste0(opt$plink, " --bfile ", opt$bfile, " --extract ", anno_str[i], " --indep-pairwise ", 
                          indep_pairwise, " --out ", opt$outPath, opt$output_prefix, "_", EV_results$Annotations[i])
    system(prunnig_cmd)
    prunning_SNPlist <- read.delim(paste0(opt$outPath, opt$output_prefix, "_", EV_results$Annotations[i], ".prune.in"), head=FALSE)
    EV_results$Number_prunning_annotated_SNPs[i] <- nrow(prunning_SNPlist)
    EV_results$Per_SNP_vc[i] <- vars[anno_name[i]==vars$Item, 2] / EV_results$Number_prunning_annotated_SNPs[i]
    EV_results$Per_SNP_h2[i] <- vars[anno_name[i]==vars$Item, 3] / EV_results$Number_prunning_annotated_SNPs[i]
  }
} else {
  EV_results <- list(Annotations=NULL, Number_annotated_SNPs=NULL, Per_SNP_vc=NULL, Per_SNP_h2=NULL)
  for (i in 1:length(anno_str)){
    anno <- read.delim(anno_str[i], head=FALSE)
    EV_results$Annotations[i] <- anno_name[i]
    EV_results$Number_annotated_SNPs[i] <- length(intersect(anno[,1], bim[,2]))
    EV_results$Per_SNP_vc[i] <- vars[anno_name[i]==vars$Item, 2] / EV_results$Number_annotated_SNPs[i]
    EV_results$Per_SNP_h2[i] <- vars[anno_name[i]==vars$Item, 3] / EV_results$Number_annotated_SNPs[i]
  }
}
print(do.call(cbind, EV_results))
write.table(EV_results, paste0(opt$outPath, opt$output_prefix, ".EA.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

## handling temporary files
if (!opt$tmp_files){
  if (is.null(opt$anno_GRM)){
    rm_cmd <- paste0("rm ", opt$outPath, "**.GA.**")
    system(rm_cmd)
  }
  if(opt$Pruning){
    rm_cmd <- paste0("rm ", opt$outPath, "**.prune.**")
    system(rm_cmd)
    rm_cmd <- paste0("rm ", opt$outPath, "**.nosex")
    system(rm_cmd)
    rm_cmd <- paste0("rm ", opt$outPath, "**.log")
    system(rm_cmd)
  }
}
cat(paste0("End time: ", Sys.time(), "\n"))
