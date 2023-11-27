##################################################################################################
# PART 1. preparing the GWAS summary statistics (weights)
# Chandra Reynolds / Adapted by Daniel Gustavson
# Date: 11 Nov 2023
##################################################################################################

###############################################
####################	1. Prep Work ##########
###############################################
### 1.1 Assign phenotype, date, author of discovery GWAS, directories
PHENO=EF
DATE=Oct2023
AUTHOR=Friedman
YEAR=2021
USER=gustavsd

# Work Areas
HOME_DIR=/pl/active/IBG/dang/CATSLife_PGSs/SbayesR/EF
HM3=/pl/active/IBG/dang/CATSLife_PGSs
SS_DIR=/pl/active/IBG/dang/ABCD_PGSs/$PHENO

#mkdir /Users/chandrar/chrs/catslife/sbayesr/scratch
#mkdir /Users/chandrar/chrs/catslife/sbayesr/scratch/$PHENO
cd $HOME_DIR/SbayesR/EF/

###############################################
### 2. Preparing GWAS summary stats (weights)
###############################################

# 2.1 Get rid of duplicates
# (base) gpvpn-general-172-18-4-115:EF_IBG chandrar$ zless $SS_DIR/EF_BoltImputeFinal.bgen.stats.gz
# SNP     CHR     BP      GENPOS  ALLELE1 ALLELE0 A1FREQ  INFO    CHISQ_LINREG    P_LINREG        BETA    SE      CHISQ_BOLT_LMM_INF      P_BOLT_LMM_INF
# rs28544273      1       751343  0.00487391      T       A       0.877114        0.98065 0.0101028       9.2E-01 8.57115e-05     0.00173069      0.00245267      9.6E-01
# rs28527770      1       751756  0.00487859      T       C       0.877007        0.982639        0.00554277      9.4E-01 3.77144e-05     0.00172833      0.000476168     9.8E-01
# rs3115860       1       753405  0.00489273      C       A       0.12913 0.993454        0.513794        4.7E-01 -0.00099343     0.00168405      0.347987        5.6E-01
# rs3131970       1       753425  0.00489285      T       C       0.124491        0.990894        0.0174722       8.9E-01 -4.54247e-05    0.00171218      0.000703857     9.8E-01


gunzip -c $HOME_DIR/EF_summstats.txt |\
awk '{seen[$1]++; if(seen[$1]==1){ print}}' |\
gzip - > $PHENO.$AUTHOR.$YEAR.nodup.gz

wc -l EF_summstats.txt
#(base) gpvpn-general-172-18-6-47:EF_IBG chandrar$ gunzip -c $SS_DIR/EF_BoltImputeFinal.bgen.stats.gz | wc -l
# 7382983

gunzip -c $PHENO.$AUTHOR.$YEAR.nodup.gz | wc -l
#(base) gpvpn-general-172-18-6-47:EF_IBG chandrar$ gunzip -c $PHENO.$AUTHOR.$YEAR.nodup.gz | wc -l
# 7382983

# 2.1.2 check if dups by CHR:POS
gunzip -c $HOME_DIR/$PHENO.$AUTHOR.$YEAR.nodup.gz |\
awk '{seen[$2, $3]++; if(seen[$2, $3]==1){ print}}' |\
gzip - > $PHENO.$AUTHOR.$YEAR.nodup23.gz

gunzip -c $PHENO.$AUTHOR.$YEAR.nodup23.gz | wc -l

#(base) gpvpn-general-172-18-6-47:EF_IBG chandrar$ gunzip -c $PHENO.$AUTHOR.$YEAR.nodup23.gz | wc -l
# 7382983
# It's the same!

# 2.2 Get rid of low MAF

gunzip -c $PHENO.$AUTHOR.$YEAR.nodup.gz |\
awk 'NR==1 || ($5>0.01 && $5<0.99) {print}' |\
gzip  > $PHENO.$AUTHOR.$YEAR.MAF.gz

gunzip -c $PHENO.$AUTHOR.$YEAR.MAF.gz | wc -l

#(base) gpvpn-general-172-18-6-47:EF_IBG chandrar$ gunzip -c $PHENO.$AUTHOR.$YEAR.MAF.gz | wc -l
# 7382484


# 2.3 Get rid of ambiguous
gunzip -c $PHENO.$AUTHOR.$YEAR.MAF.gz |\
awk '!( ($3=="A" && $4=="T") || \
        ($3=="T" && $4=="A") || \
        ($3=="G" && $4=="C") || \
        ($3=="C" && $4=="G")) {print}' |\
gzip > $PHENO.$AUTHOR.$YEAR.QC.gz

gunzip -c $PHENO.$AUTHOR.$YEAR.QC.gz | wc -l

#(base) gpvpn-general-172-18-6-47:EF_IBG chandrar$ gunzip -c $PHENO.$AUTHOR.$YEAR.QC.gz | wc -l
# 6251403
 
 #dropped 7382983-6251403 = 1,131,580 ambiguous SNPs

SUMSTATS=/pl/active/IBG/dang/CATSLife_PGSs/SbayesR/EF/$PHENO.$AUTHOR.$YEAR.QC.gz
# To view the sumstat file:
zless -S $SUMSTATS 
# To leave view:
q

#(base) gpvpn-general-172-18-4-115:EF_IBG chandrar$ zless -S $SUMSTATS 
#SNP     CHR     BP      GENPOS  ALLELE1 ALLELE0 A1FREQ  INFO    CHISQ_LINREG    P_LINREG        BETA    SE      CHISQ_BOLT_LMM_INF      P_BOLT_LMM_INF
#rs28527770      1       751756  0.00487859      T       C       0.877007        0.982639        0.00554277      9.4E-01 3.77144e-05     0.00172833      0.000476168     9.8E-01
#rs3115860       1       753405  0.00489273      C       A       0.12913 0.993454        0.513794        4.7E-01 -0.00099343     0.00168405      0.347987        5.6E-01
#rs3131970       1       753425  0.00489285      T       C       0.124491        0.990894        0.0174722       8.9E-01 -4.54247e-05    0.00171218      0.000703857     9.8E-01


## drop the row names, rename columns, and add N column (if it doesn't exist already)
R
library(tidyverse)
dat <- read.table("EF.Friedman.2021.QC",row.names=1)
colnames(dat) <- c("RSID","A1","A2","freq","b","se","p")
dat$N <- 427037
write.table(dat,"EF.Friedman.2021.QC.formatted",sep="\t",row.names=F,quote=F)


gzip EF.Friedman.2021.QC.formatted


########## 2.4 Specify the correct columns in the sumstat file 

### MODIFY
##### 2.7.1 Assign the columns
#CHR=2
#POS=3
SNP=1
EA=2
NEA=3
FRQEA=4 # A1FREQ
BETA=5
SE=6
PVAL=7		#P_BOLT_LMM_INF = p value for BETA?
#N=15
# NOTE: CR - Will need CHR:POS later for merging with catslife as Target, will get that from SNPs_IGEMS_selected_20210228.txt

### Select one of the below depending on if N is reported in the sumstats (affects the SBayesR code):
# If yes:
#NREP=Nrep
# If no:
NREP=noN

### MODIFY (depending on sumstat format)
##### 2.7.2 Select columns from the the sumstats
### NOTE! One of these options are selected, depending on whether MAF and N are provided in the sumstats

##### A) If Freq Tested Allele is available:
### A1) If N is also available, use this code: 

## Start here if QC work already completed on target set, and have Freq Tested Allele already

#gunzip -c ${SUMSTATS} | awk 'BEGIN {print "RSID" "\t" "A1" "\t" "A2" "\t" "freq" "\t" "b" "\t" "se" "\t" "p" "\t" "N"} \
#NR!=1 {print $'"$SNP"' "\t" $'"$EA"' "\t" $'"$NEA"' "\t" $'"$FRQEA"' "\t" $'"$BETA"' "\t" $'"$SE"' "\t" $'"$PVAL"' "\t" $'"$N"'}' > sumstatfile.txt

### A2) If N is not available, use this code
# Add the N from the GWAS paper here:
#EF (https://www.biorxiv.org/content/10.1101/674515v2):  N=427,037
Ngwas=427037
# Then run this
#zcat ${SUMSTATS} | awk 'BEGIN {print "RSID" "\t" "A1" "\t" "A2" "\t" "freq" "\t" "b" "\t" "se" "\t" "p" "\t" "N"} \
#NR!=1 {print $'"$SNP"' "\t" $'"$EA"' "\t" $'"$NEA"' "\t" $'"$MAF"' "\t" $'"$BETA"' "\t" $'"$SE"' "\t" $'"$PVAL"' "\t" '"$Ngwas"'}' > sumstatfile.txt
gunzip -c ${SUMSTATS} | awk 'BEGIN {print "RSID" "\t" "A1" "\t" "A2" "\t" "freq" "\t" "b" "\t" "se" "\t" "p" "\t" "N"} \
NR!=1 {print $'"$SNP"' "\t" $'"$EA"' "\t" $'"$NEA"' "\t" $'"$FRQEA"' "\t" $'"$BETA"' "\t" $'"$SE"' "\t" $'"$PVAL"' "\t" '"$Ngwas"'}' > sumstatfile.txt

#RSID    A1      A2      freq    b       se      p       N
#rs28527770      T       C       0.877007        3.77144e-05     0.00172833      9.8E-01 427037
#rs3115860       C       A       0.12913 -0.00099343     0.00168405      5.6E-01 427037
#rs3131970       T       C       0.124491        -4.54247e-05    0.00171218      9.8E-01 427037
#rs2073813       G       A       0.871339        0.00090427      0.00168808      5.9E-01 427037
#rs3131969       A       G       0.12955 -0.000984466    0.00168215      5.6E-01 427037
#rs3131968       A       G       0.129453        -0.00106321     0.00168283      5.3E-01 427037
#rs3131967       T       C       0.129551        -0.000991024    0.00168212      5.6E-01 427037
#rs3131962       A       G       0.1299  -0.00109204     0.00167511      5.1E-01 427037
#rs3115853       G       A       0.13062 -0.00106686     0.00167357      5.2E-01 427037



########## B) If Freq Tested Allele is NOT available:
### Get MAF from the 1kg reference
#plink --bfile /nfs/AGE/twins.psychchip.data/GRS/Data/DerivedData/ldref/1kgEUR --freq --out 1kG_freq
#awk 'BEGIN {print "RSID" "\t" "MA" "\t" "freq"} NR!=1 {print $2 "\t" $3 "\t" $5}' 1kG_freq.frq > 1kg_MAF.txt

### Then select one of the following, depending on if N is reported in the sumstats
### B1) If N is available, use this code: 
#zcat ${SUMSTATS} | awk 'BEGIN {print "RSID" "\t" "A1" "\t" "A2" "\t" "b" "\t" "se" "\t" "p" "\t" "N"} \
#NR!=1 {print $'"$SNP"' "\t" $'"$EA"' "\t" $'"$NEA"' "\t" $'"$BETA"' "\t" $'"$SE"' "\t" $'"$PVAL"' "\t" '"$N"'}' > sumstatfile_pre.txt
#awk 'NR==FNR {a[$1] = $0; next} ($1) in a {print $0 "\t" a[$1]}' 1kg_MAF.txt sumstatfile_pre.txt > sumstatfile_pre2.txt
#awk 'NR==1 {print $1 "\t" $2 "\t" $3 "\t" $10 "\t" $4 "\t" $5 "\t" $6 "\t" $7} \
#NR!=1 {if ($2==$9) print $1 "\t" $2 "\t" $3 "\t" $10 "\t" $4 "\t" $5 "\t" $6 "\t" $7} \
##{if ($3==$9) print $1 "\t" $2 "\t" $3 "\t" 1-$10 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' sumstatfile_pre2.txt > sumstatfile.txt

# Report number of SNPs in full GWAS
#echo 'N SNPs in full sumstats, prior to merging with 1kg for MAF' >> report.log
#wc -l < sumstatfile_pre.txt >> report.log

### NOTE! If e.g lower case letters are used for A1 and A2, or if OR is reported instead of beta, see below for modifications

### NOTE! Depending on the format of the summary statistics, some of the options below may be needed 
##### NOTE! If OR is reported rather than beta, add the following:
# mv sumstatfile.txt sumstat_pre.txt
# awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4, log($5), $6, $7, $8}' sumstat_pre.txt > sumstatfile.txt
# head sumstatfile.txt

##### NOTE! If alleles are not in capital letters, add the following:
# Convert alleles to capital letters
# mv sumstatfile.txt sumstat_pre.txt
# awk 'BEGIN{OFS="\t"}{print $1, toupper($2), toupper($3), $4, $5, $6, $7, $8}' sumstat_pre.txt > sumstatfile.txt
# head sumstatfile.txt

################################################################################
############### FORMAT MODIFICATIONS DONE! Make no changes to the below
##### Remove code not used above (left commented out, optional) and save as specified in the header 
##### Then can run as is with the nohup bash command as specified in the header
################################################################################

#echo 'N SNPs in SUMSTATS at start' >> report.log
#gunzip -c $SS_DIR/EF_BoltImputeFinal.bgen.stats.gz | wc -l >> report.log
#echo 'N SNPs after QC (sumstatfile.txt) ' >> report.log
#wc -l < sumstatfile.txt >> report.log

############ 3. Keep only those that are in HM3 SNPs##########
# Notes on IGEMS list, with header hg19

#echo 'N SNPs in CATSLife HM3 target list' >> report.log
#wc -l CATSLife_HRC_HM3_SNPs_12Feb22.nodup.txt >> report.log

# 3.1 Filter out SNPs with alleles different from those in the target data 
#awk '{print $0 "\t" $2$3 "\t" $3$2}' sumstatfile.txt > sumstatfile2.txt
#awk 'NR==FNR {a[$4] = $3; next} ($1) in a {print $0 "\t" a[$1]}' CATSLife_HRC_HM3_SNPs_12Feb22.nodup.txt sumstatfile2.txt > sumstatfile3.txt
#awk 'NR==1 {print $0} NR!=1 {if ($9==$11 || $10==$11) print $0}' sumstatfile3.txt > sumstatfile4.txt

#awk '{print $9}' sumstatfile3.txt >temp9

# Number of SNPs left after merging with IGEMS SNPs
#echo 'N SNPs after merging with IGEMS SNPs and dropping mismatches' >> report.log
#wc -l < sumstatfile4.txt >> report.log


#(base) gpvpn-general-172-18-6-47:EF_IBG chandrar$ head sumstatfile4.txt 
#RSID	A1	A2	freq	b	se	p	N	A1A2	A2A1	ALLELES
#rs3934834	C	T	0.857183	0.000449849	0.00160874	7.8E-01	427037	CT	TC	CT
#rs3766192	C	T	0.438264	-0.000970202	0.00114073	4.0E-01	427037	CT	TC	CT
#rs3766191	C	T	0.864001	0.000616973	0.00164824	7.1E-01	427037	CT	TC	CT
#rs9442372	A	G	0.425688	-0.000824158	0.00113724	4.7E-01	427037	AG	GA	AG
#rs10907177	A	G	0.856094	0.000384987	0.00161853	8.1E-01	427037	AG	GA	AG
#rs3737728	A	G	0.283103	-0.000328164	0.00126135	7.9E-01	427037	AG	GA	AG
#rs9442398	A	G	0.280574	-0.000487751	0.00126068	7.0E-01	427037	AG	GA	AG
#rs6689308	A	G	0.842895	-0.000326473	0.001546	8.3E-01	427037	AG	GA	AG
#rs6687776	C	T	0.843193	-0.00042625	0.00154311	7.8E-01	427037	CT	TC	CT

#(base) gpvpn-general-172-18-6-47:EF_IBG chandrar$ wc -l sumstatfile2.txt 
# 6251403 sumstatfile2.txt
#(base) gpvpn-general-172-18-6-47:EF_IBG chandrar$ wc -l sumstatfile3.txt 
#  931263 sumstatfile3.txt
#(base) gpvpn-general-172-18-6-47:EF_IBG chandrar$ wc -l sumstatfile4.txt 
#  931263 sumstatfile4.txt
  
# 3.1.1Drop the added columns
#awk '{$9="";$10="";$11="";print $0}' OFS='\t' sumstatfile4.txt > QCd_SUMSTATS_$PHENO.txt


#(base) gpvpn-general-172-18-6-47:EF_IBG chandrar$ head QCd_SUMSTATS_$PHENO.txt
#RSID	A1	A2	freq	b	se	p	N			
#rs3934834	C	T	0.857183	0.000449849	0.00160874	7.8E-01	427037			
#rs3766192	C	T	0.438264	-0.000970202	0.00114073	4.0E-01	427037			
#rs3766191	C	T	0.864001	0.000616973	0.00164824	7.1E-01	427037			
#rs9442372	A	G	0.425688	-0.000824158	0.00113724	4.7E-01	427037			
#rs10907177	A	G	0.856094	0.000384987	0.00161853	8.1E-01	427037			
#rs3737728	A	G	0.283103	-0.000328164	0.00126135	7.9E-01	427037			
#rs9442398	A	G	0.280574	-0.000487751	0.00126068	7.0E-01	427037			
#rs6689308	A	G	0.842895	-0.000326473	0.001546	8.3E-01	427037			
#rs6687776	C	T	0.843193	-0.00042625	0.00154311	7.8E-01	427037			

#N SNPs check after removing extra columns (should be same as prior step)-- Yes!
#wc -l < QCd_SUMSTATS_$PHENO.txt

#(base) gpvpn-general-172-18-6-47:EF_IBG chandrar$ wc -l < QCd_SUMSTATS_$PHENO.txt
#  931263

############ 4.	Rescale summary statistics using GCTB ##########
#a.	Get sparse LD matrix generated using 50k random UKB samples by SBayesR authors (22GB file downloadable via the following link: 10.5281/zenodo.3350914)

#b.	Calculate the re-weighted effect size estimates by chromosome of your summary statistics using SBayesR with the sparse LD matrix [use and adapt example script \SBayesR_steps\gctb_sbayesr_BMI_chr1to22.sh.]

mkdir OUTPUT

for i in {22..22}
do
#grun -n ${PHENO}_sbr_${i} -q hugemem.q 
gctb --sbayes R  \
     --ldm /pl/active/IBG/dang/SbayesR_ref/ukb_50k_bigset_2.8M/ukb50k_shrunk_chr${i}_mafpt01.ldm.sparse  \
	 --pi 0.95,0.02,0.02,0.01  \
	 --gamma 0.0,0.01,0.1,1  \
	 --gwas-summary /pl/active/IBG/dang/CATSLife_PGSs/SbayesR/EF/EF.Friedman.2021.QC.formatted  \
	 --seed 123  \
	 --chain-length 10000  \
	 --burn-in 2000  \
	 --robust  \
	 --impute-n \
	 --out-freq 10  \
	 --out /pl/active/IBG/dang/CATSLife_PGSs/SbayesR/EF/OUTPUT/EF_Friedman_2021_${i}
done

### used the  GCTB sparse shrunk LD matrices from 2.8M common variants from the UKB
### https://zenodo.org/records/3375373

#for i in {1..22}
#do
#grun -n ${PHENO}_sbr_${i} -q hugemem.q 
#gctb --sbayes R  \
#     --ldm /Users/chandrar/chrs/vetsa/sbayesr/SparseLD/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_chr${i}_v3_50k.ldm.sparse  \
#	 --pi 0.95,0.02,0.02,0.01  \
#	 --gamma 0.0,0.01,0.1,1  \
#	 --gwas-summary /Users/chandrar/chrs/catslife/sbayesr/scratch/${PHENO}/QCd_SUMSTATS_${PHENO}.txt  \
#	 --seed 123  \
#	 --chain-length 10000  \
#	 --burn-in 2000  \
#	 --robust  \
#	 --impute-n \
#	 --out-freq 10  \
#	 --out /Users/chandrar/chrs/catslife/sbayesr/scratch/${PHENO}/OUTPUT/${PHENO}_${AUTHOR}_${DATE}_${i}
#done


# c.	Concatenate the output .snpRes files from each chromosome to one data set [use example script \SBayesR_steps\concatenate_snpRes.sh]

cd OUTPUT

grep 'Name' ${PHENO}_${AUTHOR}_${DATE}_22.snpRes > header.txt
cat ${PHENO}_${AUTHOR}_${DATE}_*.snpRes | grep -v -w 'Name' > tmp.txt 
cat header.txt tmp.txt > ../SBayesR_sumstats_${PHENO}_${AUTHOR}_${DATE}.snpRes 
rm tmp.txt

cd ..

echo 'N SNPs after SBayesR with impute-n' >> report.log
wc -l < SBayesR_sumstats_${PHENO}_${AUTHOR}_${DATE}.snpRes >> report.log

# d.	Extract the columns "Name", "A1" and "A1Effect" from the file needed to calculate individual sum scores [use and adapt example script \SBayesR_steps\get_betas_BMI.sh]

awk '{print $2,$5,$8}' SBayesR_sumstats_EF_Friedman_2021.snpRes > EF_effects_scores.txt
wc -l < EF_effects_scores.txt

awk '{print $2,$5,$8}' SBayesR_sumstats_${PHENO}_${AUTHOR}_${DATE}.snpRes > ${PHENO}_effects_scores.txt
wc -l < ${PHENO}_effects_scores.txt




#(base) gpvpn-general-172-18-6-47:EF_IBG chandrar$ awk '{print $2,$5,$8}' SBayesR_sumstats_${PHENO}_${AUTHOR}_${DATE}.snpRes > ${PHENO}_effects_scores.txt
#(base) gpvpn-general-172-18-6-47:EF_IBG chandrar$ wc -l < ${PHENO}_effects_scores.txt
#  763703
  
#(base) gpvpn-general-172-18-6-47:EF_IBG chandrar$ head ${PHENO}_effects_scores.txt
#Name A1 A1Effect
#rs3934834 C 0.000003
#rs3766192 C -0.000017
#rs3766191 C 0.000001
#rs9442372 A -0.000024
#rs9442398 A -0.000016
#rs6689308 A -0.000009
#rs6687776 C -0.000017
#rs6678318 G -0.000018
#rs6671356 T -0.000011

############ STOPPED HERE - RESOLVE PGEN FILES and RSID ############
############ 5. replace SNP column for score generation ############
########### Make Column to match with CATSLife ID CHROM:POS:REF:ALT ##########

head $iHM3/SNPs_IGEMS_selected_20210228.txt

#(base) MacBook-Pro-4-449:EF_IBG chandrar$ head $iHM3/SNPs_IGEMS_selected_20210228.txt
#CHR	POS	ALLELES	RSID
#1	100000827	CT	rs6678176
#1	100005477	AG	rs12069019
#1	100006117	AG	rs6686057
#1	100008607	AC	rs11166268
#1	100017701	AG	rs6656134
#1	100018587	AC	rs11807493
#1	100019269	AC	rs1339866
#1	100021743	CT	rs1339865
#1	100023385	AG	rs7530721


awk '{s=substr($3,1,1)}{g=substr($3,2,length($3))}{print $1,$2, s, g, $4}' $iHM3/SNPs_IGEMS_selected_20210228.txt > SNPs_A1A2_IGEMS_selected_20210228.txt

#(base) MacBook-Pro-4-449:EF_IBG chandrar$ head SNPs_A1A2_IGEMS_selected_20210228.txt
#CHR POS A LLELES RSID
#1 100000827 C T rs6678176
#1 100005477 A G rs12069019
#1 100006117 A G rs6686057
#1 100008607 A C rs11166268
#1 100017701 A G rs6656134
#1 100018587 A C rs11807493
#1 100019269 A C rs1339866
#1 100021743 C T rs1339865
#1 100023385 A G rs7530721


awk 'BEGIN {print "SNP" "\t" "CHR" "\t" "POS" "\t" "A1" "\t" "A2" "\t" "RSID"} NR!=1 \
{print $1":"$2":"$3":"$4 "\t" $1"\t" $2"\t" $3"\t" $4"\t" $5"\t" $6}' SNPs_A1A2_IGEMS_selected_20210228.txt > SNPs_IGEMS_selected_A.txt
 
awk 'NR==FNR{a[$6]=$0;next} ($1 in a){print (a[$1]) $0}' SNPs_IGEMS_selected_A.txt ${PHENO}_effects_scores.txt > temp_A.txt
 
(base) MacBook-Pro-4-449:EF_IBG chandrar$ head temp_A.txt 
1:1005806:C:T	1	1005806	C	T	rs3934834	rs3934834 C 0.000003
1:1017197:C:T	1	1017197	C	T	rs3766192	rs3766192 C -0.000017
1:1017587:C:T	1	1017587	C	T	rs3766191	rs3766191 C 0.000001
1:1018704:A:G	1	1018704	A	G	rs9442372	rs9442372 A -0.000024
1:1021695:A:G	1	1021695	A	G	rs9442398	rs9442398 A -0.000016
1:1029805:A:G	1	1029805	A	G	rs6689308	rs6689308 A -0.000009
1:1030565:C:T	1	1030565	C	T	rs6687776	rs6687776 C -0.000017
1:1030633:A:G	1	1030633	A	G	rs6678318	rs6678318 G -0.000018
1:1040026:C:T	1	1040026	C	T	rs6671356	rs6671356 T -0.000011
1:1041700:A:G	1	1041700	A	G	rs6604968	rs6604968 A -0.000015


awk 'BEGIN {print "SNP" "\t" "CHR" "\t" "POS" "\t" "A1" "\t" "A2" "\t" "RSID"} NR!=1 \
{print $1":"$2":"$4":"$3 "\t" $1"\t" $2"\t" $3"\t" $4"\t" $5"\t" $6}' SNPs_A1A2_IGEMS_selected_20210228.txt > SNPs_IGEMS_selected_B.txt
 
awk 'NR==FNR{a[$6]=$0;next} ($1 in a){print (a[$1]) $0}' SNPs_IGEMS_selected_B.txt ${PHENO}_effects_scores.txt > temp_B.txt

(base) MacBook-Pro-4-449:EF_IBG chandrar$ head temp_B.txt 
1:1005806:T:C	1	1005806	C	T	rs3934834	rs3934834 C 0.000003
1:1017197:T:C	1	1017197	C	T	rs3766192	rs3766192 C -0.000017
1:1017587:T:C	1	1017587	C	T	rs3766191	rs3766191 C 0.000001
1:1018704:G:A	1	1018704	A	G	rs9442372	rs9442372 A -0.000024
1:1021695:G:A	1	1021695	A	G	rs9442398	rs9442398 A -0.000016
1:1029805:G:A	1	1029805	A	G	rs6689308	rs6689308 A -0.000009
1:1030565:T:C	1	1030565	C	T	rs6687776	rs6687776 C -0.000017
1:1030633:G:A	1	1030633	A	G	rs6678318	rs6678318 G -0.000018
1:1040026:T:C	1	1040026	C	T	rs6671356	rs6671356 T -0.000011
1:1041700:G:A	1	1041700	A	G	rs6604968	rs6604968 A -0.000015

(base) MacBook-Pro-4-449:EF_IBG chandrar$ wc -l temp_A.txt 
  763702 temp_A.txt
(base) MacBook-Pro-4-449:EF_IBG chandrar$ wc -l temp_B.txt 
  763702 temp_B.txt 

cat temp_A.txt temp_B.txt > temp.txt

(base) MacBook-Pro-4-449:EF_IBG chandrar$ cat temp_A.txt temp_B.txt > temp.txt
(base) MacBook-Pro-4-449:EF_IBG chandrar$ wc -l temp.txt
 1527404 temp.txt

# 5.1 Doing it the long way, and removing columns not needed:
awk '{print $1"\t" $8"\t" $9}' temp.txt > temp2.txt
 
(base) MacBook-Pro-4-449:EF_IBG chandrar$ head temp2.txt 
1:1005806:C:T	C	0.000003
1:1017197:C:T	C	-0.000017
1:1017587:C:T	C	0.000001
1:1018704:A:G	A	-0.000024
1:1021695:A:G	A	-0.000016
1:1029805:A:G	A	-0.000009
1:1030565:C:T	C	-0.000017
1:1030633:A:G	G	-0.000018
1:1040026:C:T	T	-0.000011
1:1041700:A:G	A	-0.000015

## 5.1.1 Add Headers

grep 'Name' ${PHENO}_effects_scores.txt  > header2.txt
cat header2.txt temp2.txt > ${PHENO}_effects_scores_CHRPOS.txt

Name A1 A1Effect
1:1005806:C:T	C	0.000003
1:1017197:C:T	C	-0.000017
1:1017587:C:T	C	0.000001
1:1018704:A:G	A	-0.000024
1:1021695:A:G	A	-0.000016
1:1029805:A:G	A	-0.000009
1:1030565:C:T	C	-0.000017
1:1030633:A:G	G	-0.000018
1:1040026:C:T	T	-0.000011


echo 'N SNPs after SBayesR x 2 (to be sure to get a match with Alt/Ref SNP order), check change in Col 1 to CHR:POS ' >> report.log
wc -l < ${PHENO}_effects_scores_CHRPOS.txt >> report.log

(base) MacBook-Pro-4-449:EF_IBG chandrar$ wc -l < ${PHENO}_effects_scores_CHRPOS.txt
 1527405
 
rm temp*
rm head*

################################################################################
############################### Go To Part II ####################################
################################################################################

##############################################################################
################################## Tasks ##################################### 
#1.	Prepare the target genotype data for polygenic score calculation [several ways to do these steps so no tool specific details here]:
#a.	Restrict SNPs to HapMap3 SNPs available across all IGEMS samples using the provided list \SBayesR_steps\ SNPs_IGEMS_selected_20210228.txt. Note that this list already has filtered out any out any rare (MAF<1%) and poorly imputed (INFO<0.8) variants from each chromosome across the IGEMS consortium and only includes snps available in all samples!!]
#b.	Filter out any variant with duplicate position (the final scoring will not work if there are duplicate variants in the data)
#c.	Merge all chromosomes to one dataset
##############################################################################

# Do all work in the "scratch" folder on statgen
### 1.1 Assign phenotype, date, author of discovery GWAS, directories
PHENO=EF_IBG
DATE=Feb2022
AUTHOR=Friedman
YEAR=2021
USER=chandrar

cd /home/reynoldc/PRS_work/SBayesR/scratch

mkdir $PHENO

cd $PHENO

################## 1. Generate PRS

### plink2 code run 08NOV21, rerunning the following locations
#Already made scratch

pgen=/home/reynoldc/PRS_work/Data/pgen
plink2 \
--pfile $pgen/CUdata_08Feb2022 \
--score ${PHENO}_effects_scores_CHRPOS.txt list-variants \
--rm-dup exclude-mismatch \
--out prs2_CUdata_08Feb2022_${PHENO}

awk 'FNR==NR {a[$1]; next} $1 in a' temp_A.txt prs2_CUdata_08Feb2022_EF_IBG.sscore.vars >dupsA
(base) MacBook-Pro-4-449:EF_IBG chandrar$ wc -l dupsA
  381273 dupsA

awk 'FNR==NR {a[$1]; next} $1 in a' temp_B.txt prs2_CUdata_08Feb2022_EF_IBG.sscore.vars >dupsB
(base) MacBook-Pro-4-449:EF_IBG chandrar$ wc -l dupsB
  382181 dupsB

pgen=/home/reynoldc/PRS_work/Data/pgen
plink2 \
--pfile $pgen/CUdata_08Feb2022 \
--score ${PHENO}_effects_scores_CHRPOS.txt cols=scoresums \
--rm-dup exclude-mismatch \
--out prs2_SS_CUdata_08Feb2022_${PHENO}


## my version
plink2 \
--pfile /pl/active/IBG/dang/CATSLife_PGSs/CATSLife_hrc_filtered_Oct_2023 \
--score /pl/active/IBG/dang/CATSLife_PGSs/SbayesR/EF/EF_effects_scores.txt 1 2 3 \
--rm-dup exclude-mismatch \
--out /pl/active/IBG/dang/CATSLife_PGSs/SbayesR/EF/EF_SbayesR_scores



## Clean up
#cp prs_VETSA_${PHENO}.profile $HOME_DIR/IGEMS_pheno_PRS/Data/Derived/PRS_files/${PHENO}_${AUTHOR}_PRS_SBayesR_VETSA_${DATE}.txt
#cp report.log $HOME_DIR/IGEMS_pheno_PRS/Logs/${PHENO}_${AUTHOR}_PRS_SBayesR_VETSA_report_${DATE}.log
#cp  prs_VETSA_BMI.log $HOME_DIR/IGEMS_pheno_PRS/Logs/${PHENO}_${AUTHOR}_PRS_SBayesR_VETSA_plink_${DATE}.log



