### Tutorial for running the IFAM model
## Basic
# Set parameters
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

## Advanced
# Set parameters
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


## Tutorial for evaluating the annotations
# Set parameters
EA=./Scripts/EA.R
bfile=./demo_data/geno/demo
pheno=./demo_data/phe/phenotype.txt
anno_folder=./demo_data/annotations/
GRMs_folder=./demo_data/GRMs/
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
Rscript ${EA} --bfile ${bfile} --pheno ${pheno} --anno_folder ${anno_folder} --GRMs_folder ${GRMs_folder} \
        --Pruning ${Pruning} --indep_pairwise ${indep_pairwise} --plink ${plink} --pheno_pos ${pheno_pos} \
        --VCmethod ${VCmethod} --thread ${thread} --tmp_files ${tmp_files} --outPath ${outPath} --output_prefix ${output_prefix}


# Convert annotations using bedtools
bedtools intersect -a ./demo_data/geno/demo.bed -b ./demo_data/annotations/A8.bed -wa -u > ./demo_data/test/A8.txt
# Convert annotations using plink
plink --bfile ./demo_data/geno/demo --extract ./demo_data/annotations/A8.plink.bed --range --make-bed --out ./demo_data/test/A8
awk '{print $2}' ./demo_data/test/A8.bim > ./demo_data/test/A8.txt
