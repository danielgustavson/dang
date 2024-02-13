#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=0-15:00:00     # 1 day - 1 hrs
#SBATCH --output=CATSLife_PRScs_EA4_Log.txt
#SBATCH --job-name=CATSLife_PRScs_EA4.txt

#########################
####   User Input    ####
#########################
# Name of the phenotype
PHENO=EA4                                          
# Your directory where this script is run
PGS_PATH=/pl/active/IBG/dang/CATSLife_PGSs/PRScs/EA4/   
# N for this GWAS
N_GWAS=765283                                        
# Path for Sumstat file 
SUM_PATH=/pl/active/IBG/dang/CATSLife_PGSs/PRScs/EA4/EA4_cs_format.txt

#### Probably leave these as is ####
# path to directory for genetic data
BIM_PATH=/pl/active/IBG/dang/CATSLife_PGSs/CATSLife_hrc_filtered_Oct_2023
# path to reference directory (EUR)
REF_PATH=/pl/active/IBG/dang/PRScs/ld_files/ldblk_1kg_eur

########################
####   Run PRScs    ####
########################

ml anaconda 

python /pl/active/IBG/dang/PRScs/PRScs.py \
--ref_dir=$REF_PATH \
--bim_prefix=$BIM_PATH \
--n_gwas=$N_GWAS \
--sst_file=$SUM_PATH \
--out=$PGS_PATH/${PHENO}

cat $PGS_PATH/${PHENO}_pst_eff_a1_b0.5_phiauto_chr* > $PGS_PATH/${PHENO}_prs_cs_gwas.txt

plink --bfile $BIM_PATH \
--score $PGS_PATH/${PHENO}_prs_cs_gwas.txt 2 4 6 \
--out $PGS_PATH/${PHENO}_scores
