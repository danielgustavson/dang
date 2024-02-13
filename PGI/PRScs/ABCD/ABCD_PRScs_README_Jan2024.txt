
###################################################################
##############   Generating PGS using PRSCS       #################
###   Daniel Gustavson - daniel.gustavson@colorado.edu - 7/5/23 ###
###################################################################


###################################################################
###            STEP 1: REFORMAT SUMMARY STATISTICS              ###
###################################################################

# Use an awk command to keep only the columns you need for prs-cs
# Replace the column numbers in the following command with the columns from you sst file
# COLUMNS NEEDED: SNP A1 A2 BETA P
# NOTE: You want to use summary statistics in Build 37 format!

awk -F"\t" 'NR==1{print "SNP","A1","A2","BETA","P";next}{print $5,$4,$3,$9,$11}' file > newfile



# EXAMPLE:
# Current columns:
# ROW SNP     A1      A2      EAF     beta    SE      p

# Needed columns:
# SNP A1 A2 BETA P

# awk -F"\t" 'NR==1{print "SNP","A1","A2","BETA","P";next}{print $2,$3,$4,$6,$8}' EF_summstats.txt > EF_cs_format.txt



###################################################################
###            STEP 2: RUN PRS-CS                               ###
###################################################################

# Submit a SLURM script to run PRS-CS. 
# An example can be found here: /pl/active/IBG/dang/ABCD_PGSs/PRScs_Neuroticism_ABCD_Example.slurm

### information about the settings (Note: the "#" sign is actually needed here - it's not just a comment!!)

#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=0-15:00:00     # 1 day - 1 hrs
#SBATCH --output=ABCD_PRScs_Neuroticism_Log.txt
#SBATCH --job-name=ABCD_PRScs_Neuroticism.txt

#########################
####   User Input    ####
#########################
# NOTE: this section should be all you need to change to update to your GWAS and your folder!

# Name of the phenotype
PHENO=Neuro                                          
# Your directory where this script is run
PGS_PATH=/pl/active/IBG/dang/ABCD_PGSs/Neuroticism/   
# N for this GWAS
N_GWAS=523783                                        
# Path for Sumstat file 
SUM_PATH=/pl/active/IBG/dang/CATSLife_PGSs/PRScs/Neuroticism/MA_NEU_cs_format.txt 

#### Probably leave these as is ####
# path to directory for genetic data
BIM_PATH=/pl/active/IBG/dang/ABCD_imputed_build37/ABCD_imputed_b37
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


###################################################################
###            STEP 3: MERGE WITH ANCESTRY & PC FILE            ###
###################################################################

## IF everything worked out, you should have a XXXXX.profile file in your output directory
# Download this file which has the ancestry and PC information for ABCD:
# /pl/active/IBG/dang/ABCD/ABCD_Ancestry_1KGPCs_Jan2024.txt 

# Note: PCs are plotted against the 1KG Phase 3 PCs
#   The Ancestry groups are determined based on +/- 5 SD from the population mean from 1KG
#   (except for AMR which is based on 3 SD and do not use people already classified in other groups.
#    This is because the AMR group has a really wide distribution, especially for the first 1 or 2 PCs)

# Here is an example of how you would merge you PRScs output with the PC/Ancestry file

# Read PC/Ancestry and PGS score files & merge
PCs <- read.table("ABCD_Ancestry_1KGPCs_Jan2024.txt",header=T)
PGS <- read.table("GCA_scores.profile",header=T)
PGS2 <- select(PGS2, IID, SCORE)
PGS_PC <- merge(PGS2, PCs)

# Fix Typo for one observation: `NDAR_INVF3FYXH1G (the ` should not be there)
PGS_PC$FID<-gsub("`","",as.character(PGS_PC$FID))
PGS_PC$IID<-gsub("`","",as.character(PGS_PC$IID))


## You should be able to merge this file with additional phenotypic information using the 'ID' variable. Good luck!

