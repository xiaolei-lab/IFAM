# IFAM
**I**ntegrating **F**unctional **A**nnotation information by genomic BLUP model with **M**ultiple random effects

## About
IFAM adopts a multiple random effects model to  integrate massive types of genomic functional annotation information such as enhancer, promoter, and transcription factor binding site et al. to improve the accuracy of genomic prediction for complex traits. In IFAM, functional annotation markers with similar contributions to phenotype variance are automatically merged to construct random effects, the association between markers and trait is utilized as another random effect. <br>
Details of IFAM could be found in our [IFAM manuscript](https:****).

## Tutorial for IFAM
In this softare, we use `**.R` scripts to make the usage of IFAM. 

### Input files and formats
* Genotype file: IFAM only accept the genotype in PLINK binary format, e.g. demo.fam, demo.bim and demo.bed, please see more details about these files at PLINK user manual. Users can convert any other format of genotype (e.g. VCF, HapMap, PED/MAP) to binary format by PLINK2.

* Phenotype file (e.g. `demo_data/pheno.txt`): this file includes the phenotypic records, the environmental covariates, fixed and random effects. The first column must be the individual id, the second column is phenotypic records (optional), header should be included in the file.

* Annotation files: a list of the annotation files, the file name must be "annotation_names.txt" eg. "enhancer.txt". Each annotation file has only one column representing the position of annotationed SNPs, note that the genome version of the annotation file and genome file should be the same. If the annotation file is "*.bed" format, users can filtered all annotationed SNPs using bedtools or PLINK softwares.
```bash
    # Using bedtools to filter annotationed SNPs
    bedtools intersect -a demo.bed -b enhancer.bed -wa -u > enhancer.txt
    # Using PLINK to filter annotationed SNPs
    plink --bfile demo --extract enhancer.bed --range --make-bed --out enhancer
    awk '{print $2}' enhancer.bim > enhancer.txt
````

For other parapeters, please use the commond "--help"

### Tutorial for running the IFAM model
Please install [HIBLUP](https://www.hiblup.com/tutorials#running-hiblup) software in advance
```bash
# Set parameters
IFAM=./Scripts/IFAM.R
bfile=./demo_data/geno/demo
pheno=./demo_data/phe/phenotype.txt
anno_folder=./demo_data/annotations
anno_spec=A8
GRMs_folder=./demo_data/GRMs
weight=./demo_data/SNP_weight.txt
VCfile=./demo_data/trait.vars
outPath=./test/
output_prefix=IFAM
pheno_pos=2
randomMax=5
thread=2
VCmethod=AI
tmp_files=TRUE

# Run IFAM
Rscript ${IFAM} --bfile ${bfile} --pheno ${pheno} --anno_folder ${anno_folder} --anno_spec ${anno_spec}\
        --pheno_pos ${pheno_pos} --GRMs_folder ${GRMs_folder} --weight ${weight} --randomMax ${randomMax}\
        --VCmethod ${VCmethod} --thread ${thread} --tmp_files ${tmp_files} --outPath ${outPath}\
        --VCfile ${VCfile} --output_prefix ${output_prefix}
````

### Tutorial for Evaluating the Annotations (EA)
Please install [HIBLUP](https://www.hiblup.com/tutorials#running-hiblup) and [PLINK](https://zzz.bwh.harvard.edu/plink/) softwares in advance
```bash
# Set parameters
EA=/Scripts/EA.R
bfile=/demo_data/demo
pheno=/demo_data/phenotype.txt
anno=/demo_data/A1.txt,/demo_data/A2.txt,/demo_data/A3.txt,/demo_data/A4.txt,/demo_data/A5.txt,/demo_data/A6.txt,/demo_data/A7.txt
anno_GRM=/demo_data/A1.GA,/demo_data/A2.GA,/demo_data/A3.GA,/demo_data/A4.GA,/demo_data/A5.GA,/demo_data/A6.GA,/demo_data/A7.GA
Pruning=FALSE
indep_pairwise=1000,100,0.2
plink=plink
outPath=/output_path/test/
output_prefix=test
pheno_pos=2
thread=2
VCmethod=AI
tmp_files=FALSE

# Run EA
Rscript ${EA} --bfile ${bfile} --pheno ${pheno} --anno ${anno} --anno_GRM ${anno_GRM}\
        --Pruning ${Pruning} --indep_pairwise ${indep_pairwise} --plink ${plink} --pheno_pos ${pheno_pos} \
        --VCmethod ${VCmethod} --thread ${thread} --tmp_files ${tmp_files} --outPath ${outPath}\
        --output_prefix ${output_prefix}
````
 
# Citation
For IFAM:
...........   <br>
For HIBLUP software:
Lilin Yin, Haohao Zhang, Zhenshuang Tang, Dong Yin, Yuhua Fu, Xiaohui Yuan, Xinyun Li, Xiaolei Liu, Shuhong Zhao, HIBLUP: an integration of statistical models on the BLUP framework for efficient genetic evaluation using big genomic data, Nucleic Acids Research, Volume 51, Issue 8, 8 May 2023, Pages 3501–3512, https://doi.org/10.1093/nar/gkad074.
