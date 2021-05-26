#!/bin/sh

# Adapted from https://github.com/gusevlab/fusion_twas/blob/master/examples/GEUV.compute_weights.sh

VCF=$1 # VCF WITH GENOTYPES OF ALL SAMPLES
COVAR=$2 # COVARIATE MATRIX2
RPKM=$3 # PATH TO QTL PHENOTYPES
LDREF=$4 # THIS IS USED TO RESTRICT INPUT SNPS TO REFERENCE IDS ONLY
CHR=$5
WIN=$6
OUT_BASE=$7
FUSION_PATH=$8 # path to folder with FUSION.compute_weights.R


# PATH TO OUTPUT DIRECTORY 
OUT_DIR="$OUT_BASE/WEIGHTS"

# --- BEGIN SCRIPT:

#reformat covariate matrix:
Rscript $FUSION_PATH/reformat_covar.R $COVAR $OUT_BASE/covar.format.mat  #should move this so it's not done with every batch

#create file with ids in expected format
cut -f1-2 $OUT_BASE/covar.format.mat  -d ' ' | sed '1d' > $OUT_BASE/ids.txt #should move this so it's not done with every batch

mkdir --parents $OUT_BASE/tmp/$CHR #should be able to remove this
mkdir --parents $OUT_BASE/hsq/$CHR #should be able to remove this
mkdir --parents $OUT_BASE/out/$CHR

# THIS IS DIRECTORY WHERE THE OUTPUT WILL GO:
mkdir -p $OUT_DIR

# Loop through each gene expression phenotype in the batch
zcat $RPKM | awk -vc=$CHR '$1==c' |  while read PARAM; do

# Get the gene positions +/- win
CHR=`echo $PARAM | awk '{ print $1 }'`
P0=`echo $PARAM | awk -v win=$WIN '{ print $2 - win }'`
P1=`echo $PARAM | awk -v win=$WIN '{ print $3 + win }'`
GNAME=`echo $PARAM | awk '{ print $4 }'`

if [ $P0 -lt 0 ]
then
	P0=0
fi

OUT="$OUT_BASE/tmp/$CHR/$GNAME"

echo $GNAME $CHR $P0 $P1
echo

# Pull out the current peak phenotype
echo $PARAM | tr ' ' '\n' |  tail -n+7 | paste $OUT_BASE/ids.txt - > $OUT.pheno  

# Get the locus genotypes for all samples and set current peak as the phenotype
plink2 --vcf $VCF --const-fid --pheno $OUT.pheno --make-bed --out $OUT --chr $CHR --from-bp $P0 --to-bp $P1 --extract $LDREF/1000G.EUR.$CHR.bim --allow-no-vars #allow-no vars is tmp

# Process all samples together 
mkdir $OUT_DIR/ALL
FINAL_OUT="$OUT_DIR/ALL/ALL.$GNAME"

Rscript $FUSION_PATH/FUSION.compute_weights.R --bfile $OUT --tmp $OUT --out $FINAL_OUT --verbose 0 --save_hsq --PATH_gcta gcta64 --PATH_gemma gemma --models lasso,top1 --covar $OUT_BASE/covar.format.mat --PATH_plink plink2

# Append heritability output to hsq file
if [ -f $FINAL_OUT.hsq ]; then
	cat $FINAL_OUT.hsq >> $OUT_BASE/hsq/$CHR.hsq
fi

rm $OUT.*

# GO TO THE NEXT GENE
done
