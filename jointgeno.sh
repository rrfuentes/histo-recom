#!/usr/bin/bash

#-----------------------------Output files-----------------------------
#SBATCH --output=output_%A.%a.txt
#SBATCH --error=error_%A.%a.txt
#-----------------------------Required resources-----------------------
#SBATCH --time=0-200:00:00
#SBATCH --array=1-12
#SBATCH --cpus-per-task=20
#SBATCH --mem=20000
#----------------------------------------------------------------------

threads=20
module load java/jre/1.8.0/144
module load gcc/7.1.0
module load bedtools/gcc/64/2.28.0
GATK=Tools/gatk-4.1.8.1/gatk

samplelist=FINAL_PimpSamples_filtered
chromfile=S_lycopersicum_chromosomes.4.00.chrom
reference=S_lycopersicum_chromosomes.4.00.fa
repeats=ITAG4.0_RepeatModeler_repeats_light.gff
outdir=Pimp_Genomes/Merged

chrom=$(sed -n "$SLURM_ARRAY_TASK_ID"p "$chromfile")
mkdir -p $outdir/$chrom
#rm -r $dbname
outfile=$chrom/${chrom}_merged
dbname=$outdir/$chrom/${chrom}_DB

$GATK --java-options "-Xmx50g" GenomicsDBImport --genomicsdb-workspace-path $dbname --intervals $chrom --sample-name-map $samplelist --tmp-dir /path/Tmp

$GATK --java-options "-Xmx50g" GenotypeGVCFs -R $reference -V gendb://${dbname} -O $outdir/${outfile}.vcf.gz

$GATK VariantFiltration -V $outdir/${outfile}.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O $outdir/${outfile}_tmp1.vcf.gz

#More filtering plus setting het to missing
bedtools intersect -header -v -a $outdir/${outfile}_tmp1.vcf.gz -b $repeats | bcftools +fill-tags -Ov -- -t AC_Het,MAF | bcftools view -m2 -M2 -v snps -i '%FILTER=="PASS" && AC_Het<=1 && F_MISSING<0.1 && MAF>0.05' -Ov | bcftools +setGT -Oz -o $outdir/${outfile}_filtered.vcf.gz -- -t q -i 'GT="het"' -n .

#Impute missing genotypes
#default ne=1000000; longest distance between consecutive pair of SNPs across all chrom is below 1Mb
#5cM or 5Mb allows more SNPs per window given varying SNP density across the chromosome 
java -Xmx20g -jar Beagle/beagle.18May20.d20.jar gt=$outdir/${outfile}_filtered.vcf.gz out=$outdir/${outfile}_imputed nthreads=$threads window=5 overlap=2 iterations=30 burnin=10 gp=true err=0.001 

vcftools --gzvcf $outdir/${outfile}_imputed.vcf.gz --chr $chrom --ldhat --out $outdir/${outfile}

rm $outdir/${outfile}_tmp1.vcf.gz*
