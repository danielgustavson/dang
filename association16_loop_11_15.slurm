#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=0-12:00:00     # 1 day - 1 hrs
#SBATCH --output=cats_gwas16_11_15.txt
#SBATCH --job-name=cats_gwas16_11_15


for i in {11..15}
do
	rvtest --inVcf /pl/active/IBG/dang/CATSLife_PGSs/origVCFs/chr$i.MZincluded.vcf.gz --pheno /pl/active/IBG/dang/CATSLife_PGSs/GWAS/musicGWAS_Nov2023.pheno --pheno-name z16talent --covar musicGWAS_Nov2023.cov --covar-name sex,age16,age16_sq,sex_age16,byear,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --out /pl/active/IBG/dang/CATSLife_PGSs/GWAS/output16/output16.chr$i --kinship CATSLife_kinship_Oct2023.kinship --meta score
done