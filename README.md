### IFAM: **I**ntegrating **F**unctional **A**nnotation information by genomic BLUP model with **M**ultiple random effects




### Repos:
**Github:** https://github.com/xiaolei-lab/IFAM  

### Authors:
**Design and Maintenance:** [Zhenshuang Tang](https://github.com/Zhenshuang), [Lilin Yin](https://github.com/YinLiLin), Shuhong Zhao, and [**Xiaolei Liu**](https://github.com/XiaoleiLiuBio).  

Questions, suggestions, and bug reports are welcome and appreciated: [xiaoleiliu@mail.hzau.edu.cn](mailto:xiaoleiliu@mail.hzau.edu.cn)




## About
IFAM extends the genomic best linear unbiased prediction (GBLUP) model with multiple random effects to accommodate massive types of functional annotations, these random effects include the groups of genetic markers from different functional annotations and an additional group of genetic markers that were associated significantly with objective traits, and the groups of genetic markers with similar contributions to phenotypic variance were automatically merged into one single random effect during variance assessment. For the aspect of speed, IFAM efficiently handles UKB scale datasets by using our previously proposed the phenotypic variance-covariance V matrix based "HE+PCG" strategy, which is pretty friendly to multiple random effect model since its computational complexity remains unchanged with the number of random effects increased. <br>

## Tutorial for IFAM
In this softare, we use `IFAM.R` script to make the usage of IFAM. WE STRONGLY RECOMMEND TO use this script in the Linux or Mac operating system.

### Input files and formats
* Genotype file: IFAM only accept the genotype in PLINK binary format, e.g. `./demo_data/geno/demo`, please see more details about these files at PLINK user manual. Users can convert any other format of genotype (e.g. VCF, HapMap, PED/MAP) to binary format by PLINK. 

* Phenotype file: e.g. `./demo_data/phe/phenotype.txt`, the first column must be the individual id, the second column is phenotypic records (optional), header should be included in the file.

* Annotation files: a list of the annotation files, the file name must be "annotation_names.txt", e.g. `./demo_data/annotations/A1.txt`. Each annotation file has only one column representing the annotationed SNPs, note that the genome version of the annotation file and genome file should be the same. If the annotation file is "*.bed" format, users can filtered all annotationed SNPs using bedtools or PLINK softwares.
```bash
# Using bedtools to filter annotationed SNPs
bedtools intersect -a ./demo_data/geno/demo.bed -b ./demo_data/annotations/A8.bed -wa -u > ./demo_data/test/A8.txt
# Using PLINK to filter annotationed SNPs
plink --bfile ./demo_data/geno/demo --extract ./demo_data/annotations/A8.plink.bed --range --make-bed --out ./demo_data/test/A8
awk '{print $2}' ./demo_data/test/A8.bim > ./demo_data/test/A8.txt
````


### Tutorial for running the IFAM model
Please install [HIBLUP v1.1.0](https://www.hiblup.com/tutorials#running-hiblup) software and the [optparse v1.7.5](https://cran.r-project.org/web/packages/optparse/index.html) R packages in advance

### Basic
```bash
# Set parameters
# The annotation files are stored in anno_folder
IFAM=./Scripts/IFAM.R
bfile=./demo_data/geno/demo
pheno=./demo_data/phe/phenotype.txt
anno_folder=./demo_data/annotations/
pheno_pos=2
outPath=./demo_data/test/
output_prefix=IFAM

# Run IFAM
Rscript ${IFAM} --bfile ${bfile} --pheno ${pheno} --anno_folder ${anno_folder} \
        --pheno_pos ${pheno_pos} --outPath ${outPath} --output_prefix ${output_prefix}
````

### Advanced
```bash
# Set parameters
# the pre-constructed GRM of each annotation is stored in GRMs_folder, and the number of GRMs must be the same as that of annotations
IFAM=./Scripts/IFAM.R
bfile=./demo_data/geno/demo
pheno=./demo_data/phe/phenotype.txt
anno_folder=./demo_data/annotations/
anno_spec=A8
GRMs_folder=./demo_data/GRMs/
weight=./demo_data/SNP_weight.txt
VCfile=./demo_data/trait.vars
outPath=./demo_data/test/
output_prefix=IFAM
pheno_pos=2
randomMax=5
thread=2
VCmethod=AI
tmp_files=FALSE

# Run IFAM
Rscript ${IFAM} --bfile ${bfile} --pheno ${pheno} --anno_folder ${anno_folder} --anno_spec ${anno_spec}\
        --pheno_pos ${pheno_pos} --GRMs_folder ${GRMs_folder} --weight ${weight} --randomMax ${randomMax}\
        --VCmethod ${VCmethod} --thread ${thread} --tmp_files ${tmp_files} --outPath ${outPath}\
        --VCfile ${VCfile} --output_prefix ${output_prefix}
````




### Tutorial for Evaluating the Annotations (EA)
Please install [HIBLUP v1.1.0](https://www.hiblup.com/tutorials#running-hiblup), [PLINK v1.90](https://zzz.bwh.harvard.edu/plink/), and the [optparse v1.7.5](https://cran.r-project.org/web/packages/optparse/index.html) R packages in advance
```bash
EA=./Scripts/EA.R
bfile=./demo_data/geno/demo
pheno=./demo_data/phe/phenotype.txt
anno_folder=./demo_data/annotations/
anno_GRM=./demo_data/GRMs/
Pruning=FALSE
indep_pairwise=1000,100,0.2
plink=plink
outPath=./demo_data/test/
output_prefix=IFAM_EA
pheno_pos=2
thread=2
VCmethod=AI
tmp_files=TRUE

# Run EA
Rscript ${EA} --bfile ${bfile} --pheno ${pheno} --anno_folder ${anno_folder} --anno_GRM ${anno_GRM} \
        --Pruning ${Pruning} --indep_pairwise ${indep_pairwise} --plink ${plink} --pheno_pos ${pheno_pos} \
        --VCmethod ${VCmethod} --thread ${thread} --tmp_files ${tmp_files} --outPath ${outPath} --output_prefix ${output_prefix}
````

 
# Citation
For HIBLUP software:
Lilin Yin, Haohao Zhang, Zhenshuang Tang, Dong Yin, Yuhua Fu, Xiaohui Yuan, Xinyun Li, Xiaolei Liu, Shuhong Zhao, HIBLUP: an integration of statistical models on the BLUP framework for efficient genetic evaluation using big genomic data, Nucleic Acids Research, Volume 51, Issue 8, 8 May 2023, Pages 3501–3512, https://doi.org/10.1093/nar/gkad074.


For IFAM software:
Zhenshuang Tang, Xiong Xiong, Haohao Zhang, Dong Yin, Yuhua Fu, Yunxia Zhao, Jingjin Li, Yuan Quan, Xiang Zhou, Xinyun Li, Lilin Yin, Shuhong Zhao, and Xiaolei Liu, IFAM:Improving genomic prediction accuracy of complex traits by integrating massive types of functional annotation information, https://github.com/xiaolei-lab/IFAM, DOI: 10.5281/zenodo.18802710, 2026.



