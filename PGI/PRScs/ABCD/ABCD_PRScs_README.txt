
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
# An example can be found here: /pl/active/IBG/dang/ABCD_PGSs/GCA/prs_cs_GCA_ABCD.slurm

# Your script should look like this:
## Basic commands for SLURM (edit as necessary) & load anaconda

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=0-12:00:00     # 12 hrs
#SBATCH --output=ABCD_GCA.txt
#SBATCH --job-name=PGS_GCA_ABCD

ml anaconda 

## Call PRS-CS. 
# You may need to update the reference directory (if your GWAS is not from Europeans) and the bim prefix (if you want to apply this outside of ABCD). Otherwise keep as is!
# You definitely need to update the n_gwas (based on the GWAS you are using), the sst_file location (to your directory), and the output directory 

python /pl/active/IBG/dang/PRScs/PRScs.py --ref_dir=/pl/active/IBG/dang/PRScs/ld_files/ldblk_1kg_eur --bim_prefix=/pl/active/IBG/dang/ABCD_imputed_build37/ABCD_imputed_b37 --n_gwas=331679 --sst_file=/pl/active/IBG/dang/ABCD_PGSs/GCA/GCA_cs_format.txt --out=/pl/active/IBG/dang/ABCD_PGSs/GCA/

## Concatinate PRS-CS output
# Update the directory names and output filename. PRS-CS generates "_pst_eff_a1_b0.5_phiauto_chr" files for each chromosome, so this just concatinates them together.

cat /pl/active/IBG/dang/ABCD_PGSs/GCA/_pst_eff_a1_b0.5_phiauto_chr* > /pl/active/IBG/dang/ABCD_PGSs/GCA/GCA_prs_cs_gwas.txt

## Generate scores in PLINK
# update the --score statement to the file from the previous line and update the output directory. The --bfile needs to stay the same as what you used above

plink --bfile /pl/active/IBG/dang/ABCD_imputed_build37/ABCD_imputed_b37 --score /pl/active/IBG/dang/ABCD_PGSs/GCA/GCA_prs_cs_gwas.txt 2 4 6 --out /pl/active/IBG/dang/ABCD_PGSs/GCA/GCA_scores



###################################################################
###            STEP 3: MERGE WITH ANCESTRY & PC FILE            ###
###################################################################

## IF everything worked out, you should have a XXXXX.profile file in your output directory
# Download this file which has the ancestry and PC information for ABCD:
# /pl/active/IBG/dang/ABCD/ABCD_ancestries_PRIMUS_PCs 

# You will need to strip part of the subject ID in order to merge with other phenotypic data.
# Here is some example R code (this may be updated in the future)

# Read files
PCs <- read.csv("ABCD_ancestries_PRIMUS_PCs.csv",header=T)
PGS <- read.table("GCA_scores.profile",header=T)


# Split the subject key from the full ID
PGS2 <- separate_wider_delim(as.data.frame(PGS), delim="_", cols=IID, names=c(NA,"B","C"), too_few="align_start")
PGS2$IID <- paste(PGS2$B, PGS2$C, sep="_")
head(PGS2)

# Select just the new IID & PGS (SCORE) and merge with ancestry file 
PGS3 <- select(PGS2, IID, SCORE)
PGS_PC <- merge(PGS3, PCs)


## You should be able to merge this file with additional phenotypic information. Good luck!

