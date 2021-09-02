#!/usr/bin/bash

#-----------------------------Output files-----------------------------
#SBATCH --output=output_%A.%a.txt
#SBATCH --error=error_%A.%a.txt
#-----------------------------Required resources-----------------------
#SBATCH --time=0-200:00:00
#SBATCH --array=0-80
#SBATCH --cpus-per-task=2 #20
#SBATCH --mem-per-cpu=2000 #3000
#----------------------------------------------------------------------
  
module load WUR/BIOINF/bwa/0.7.17
module load bamtools/gcc/64/2.4.0
module load java/jre/1.8.0/144

configfile=PIPE_PARAM.cfg
if [ -f $configfile ];then
   echo "Loading config file..." >&2
   #Check for incorrectly formatted input
   CONFIG_SYNTAX="(^\s*#|^\s*$|^\s*[a-z_][^[:space:]]*=[^;&\(\`]*$)"
   if egrep -q -iv "$CONFIG_SYNTAX" "$configfile"; then
      echo "Config file is unclean, Please  cleaning it..." >&2
      exit 1
   fi
   source "$configfile"
else
   echo "Configuration file ${configfile} is missing."
   exit 1
fi

threads=2

SV_dir=$dest/StructuralVariant
dellyout=$SV_dir/Delly
pindelout=$SV_dir/Pindel
gromout=$SV_dir/GROM
lumpyout=$SV_dir/Lumpy

if [ ! -d $SV_dir ]; then
    mkdir $SV_dir
    mkdir $dellyout
    mkdir $pindelout
    mkdir $gromout
    mkdir $lumpyout
fi

BAM=$dest/BAM
mkdir -p $BAM

intermediate=$dest/intermediate
mkdir -p $intermediate

VCF=$dest/VCF
mkdir -p $VCF

samplelist=()
while read lib
do
    samplelist+=($lib)
done < $input

## Align to reference
sam_info=${samplelist[${SLURM_ARRAY_TASK_ID}]}
sample=$(echo $sam_info | cut -f1 -d",")
   
readfile=$(echo $sam_info | cut -f2 -d",")
   
java -jar $Trimmomatic PE -threads $threads $input_path/${readfile}_1.fastq.gz $input_path/${readfile}_2.fastq.gz $intermediate/${sample}_1.fastq.gz $intermediate/${sample}_1u.fastq.gz $intermediate/${sample}_2.fastq.gz $intermediate/${sample}_2u.fastq.gz ILLUMINACLIP:$adapter:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:3
   
bwa mem -t $threads -o $intermediate/${sample}.sam $ref $intermediate/${sample}_1.fastq.gz $intermediate/${sample}_2.fastq.gz

${samtools_dir}/samtools sort -@ $threads -n -O BAM $intermediate/${sample}.sam | ${samtools_dir}/samtools fixmate -@ $threads -m -O BAM - - | ${samtools_dir}/samtools sort -@ $threads -O BAM - | ${samtools_dir}/samtools markdup -@ $threads -O BAM - - | ${samtools_dir}/samtools addreplacerg -@ $threads -r ID:${sample} -r PL:Illumina -r SM:${sample} -r LB:${sample} -O BAM - > $BAM/${sample}.bam

samtools index $BAM/${sample}.bam
rm $intermediate/${sample}*
    
/home/WUR/fuent011/Tools/gatk-4.1.0.0/gatk --java-options "-Xmx60g" HaplotypeCaller --native-pair-hmm-threads $threads -R $ref -I $BAM/${sample}.bam -O $VCF/${sample}_HC.vcf.gz -ERC GVCF

bamfile=$BAM/${sample}.bam

bamtools stats -in $bamfile -insert > $BAM/${sample}.bam.stats
bedtools genomecov -ibam $bamfile -bga | awk '$1!="SL4.0ch00"{if($4>0) sum+=$3-$2;} END{print "Total covered bases: "sum;}' - >> $BAM/${sample}.bam.stats
samtools view -F 4 $bamfile | awk '{if(NR==10000000) exit; sum+=length($10);} END{printf("Read length: %d\n",sum/10000000);}' - >> $BAM/${sample}.bam.stats

echo "Done: $sample"
