#!/usr/bin/bash

#-----------------------------Output files-----------------------------
#SBATCH --output=seqLD_output_%A.%a.txt
#SBATCH --error=seqLD_error_%A.%a.txt
#-----------------------------Required resources-----------------------
#SBATCH --time=0-200:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=30000
#----------------------------------------------------------------------

while getopts c:i:r:h:l:w:p: option
do
  case "${option}"
  in
    c) chrom=${OPTARG};;
    i) inputsnp=${OPTARG};;
    r) inputrate=${OPTARG};;
    #h) haplotypes=${OPTARG};;
    #l) loci=${OPTARG};;
    w) windowsize=${OPTARG};; #smoothing window
    p) prefix=${OPTARG};;
  esac
done

idx=$SLURM_ARRAY_TASK_ID
locfile=${inputsnp}_${idx}.locs
genfile=${inputsnp}_${idx}.sites
recfile=${inputrate}_${idx}.res.txt
recback=${inputrate}_${idx}.bgres
tmpfile=${prefix}_${idx}_input

#PREPARE input for sequenceLDhot
#Compute md5dum hash values for each haplotype, then group and count similar haplotype
awk -F"\t" -v outfile=$tmpfile 'BEGIN{seqname="";} NR==1{chr=$1; loci=$2;} NR>1{if($1~/^>/){ if(seqname!=""){ cmd="echo -n "sq[seqname]" | md5sum"; cmd|getline md5; close(cmd); gsub(/ .*/,"",md5); haplotypes[md5]++; hapseq[md5]=sq[seqname];} seqname=$1;} else{ sq[seqname]=sq[seqname]""$0;}} END{cmd="echo -n "sq[seqname]" | md5sum"; cmd|getline md5; close(cmd); gsub(/ .*/,"",md5); haplotypes[md5]++; hapseq[md5]=sq[seqname]; for(i in haplotypes){ print hapseq[i]" "haplotypes[i]; distinct++;} print "#"; printf("Distinct = %d\nGenes = %d\nLoci = %d\nI = 1\nK = -2\nPositions of loci:\n",distinct,chr,loci) > outfile;}' $genfile >> ${tmpfile}_tmp

awk -v ORS="" 'NR>1{if(NR==2) printf("%.0f",$1*1000); else printf(" %.0f",$1*1000);} END{print "\nHaplotypes\n"}' $locfile >> $tmpfile
less ${tmpfile}_tmp >> $tmpfile
rm ${tmpfile}_tmp 

#REFORMAT LDhat output and use as background recombination
lastsnp=$(tail -n1 $locfile | awk '{printf("%d",$1*1000);}')

awk -v lastp=$lastsnp -v chr=$chrom 'NR==3{start=$1*1000; rate=$2;} 
		NR>3{ win=1000; end=$1*1000; i=start; 
			do{if(i+win>=end) win=end-i; 
				printf("%s\t%.0f\t%.0f\t%f\t%d\n",chr,i,i+win-1,rate,NR-2); 
				i+=win;
			}while(i+1<end) 
			start=end; rate=$2;
		 }END{win=1000; end=lastp; i=start;
                        do{if(i+win>=end) win=end-i+1;
                                printf("%s\t%.0f\t%.0f\t%f\t%d\n",chr,i,i+win-1,rate,NR-2);
                                i+=win;
                        }while(i<end)
                        start=end; rate=$2;
		}' $recfile > ${recfile}_1kbwindows

#Smooth the recombination rate using 50-kb window and the median rho
echo -e "Lastsnp: $lastsnp"
awk -v lastp=$lastsnp -v win=$windowsize -v chr=$chrom 'BEGIN{window=win/2;} 
		NR==3{start=$1*1000; exit;} 
		END{lastp+=2000; #include extra bases to allow inclusion of last position
			for(i=start;i<=lastp;i+=1000){
				mid=i+500; 
				if(mid-window<1) x=1; 
				else x=mid-window; 
				printf("%s\t%.0f\t%.0f\t%d\n",chr,x,mid+window,mid);
			}
		}'  $recfile | bedtools map -c 4 -o median -a - -b ${recfile}_1kbwindows -null 0.1 | awk '{printf("%s\t%d\t%d\t%d\t%.10f\n", $1,$2,$3,$4,$5);}' - > $recback

rm ${recfile}_1kbwindows
awk '{print $5;}' $recback > ${recback}_rate

#RUN sequenceLDhot
sequenceLDhot /path/Infile ${tmpfile} -V ${recback}_rate


