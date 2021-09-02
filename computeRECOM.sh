#!/usr/bin/bash

#-----------------------------Output files-----------------------------
#SBATCH --output=output_%A.%a.txt
#SBATCH --error=error_%A.%a.txt
#-----------------------------Required resources-----------------------
#SBATCH --time=0-200:00:00
#SBATCH --array=1-12
#SBATCH --cpus-per-task=2 #10
#SBATCH --mem=5000 #50000
#----------------------------------------------------------------------

chromfile=S_lycopersicum_chromosomes.4.00.chrom
chrom=$(sed -n "$SLURM_ARRAY_TASK_ID"p "$chromfile")
inputdir=Pimp_Genomes/Merged/$chrom
outdir=Pimp_Genomes
sites=$inputdir/${chrom}_merged.ldhat.sites
locs=$inputdir/${chrom}_merged.ldhat.locs
window=5000
overlap=500
pos=$(awk 'NR==1{print $1; exit;}' $locs)
splitdir=$inputdir/Windows
LDout=$inputdir/LDhat
seqLDout=$inputdir/seqLDhot
rhoout=$outdir/LDhat
seqout=$outdir/seqLD
mkdir -p $splitdir
mkdir -p $LDout
mkdir -p $seqLDout
mkdir -p $rhoout
mkdir -p $seqout

numhap=$(awk 'NR==1{print $1; exit;}' $sites) #retrieve the number of haplotypes

#GENERATE lookup table separately
#if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
#       lkgen -lk /home/WUR/fuent011/Tools/LDhat/lk_files/lk_n192_t0.001 -nseq $numhap -prefix $outdir/lk_n${numhap}_t0.001
#fi

lkfile=Pimp_Genomes/Merged/lk_n${numhap}_t0.001new_lk.txt

lim=$(((pos/(window-overlap))-1)) #0-based file/window naming
remainder=$((pos%(window-overlap)))
if [ $remainder -ne 0 ];then
	if [ $remainder -gt $overlap ];then (( lim++ )); fi
	#if remainder is less than overlap, adjust the size of the last window
fi

if [ $pos -lt $window ];then #Number of SNPs is less than window size
	window=$pos
	overlap=0
	lim=0
fi

echo -e "$chrom\t$((lim+1))"

skip=1
if [ $skip -eq 0 ]; then
for i in `seq 0 $lim`
do
	strpos=$((i*window - overlap*$i + 1)) #1-based
	endpos=$((strpos+window - 1))
	cur_window=$window
	if [ $remainder -ne 0 ] && [ $i -eq $lim ]; then
		#Last window contains the remaining SNPs but spans the same size as other window#
		if [ $remainder -gt $overlap ];then 
			strpos=$((pos-window+1)) #take the last (window)-SNPs from the end
			endpos=$pos #The very last SNP position
		else #otherwise, use same start but shorter window size
			cur_window=$((window-overlap+remainder))
			endpos=$((strpos+window-overlap+remainder-1))
		fi
	fi
	echo -e "Window:$i;\t$strpos-$endpos Size:$cur_window"
        
	#Prepare locs file per window
	awk -v winidx=$i -v win=$cur_window -v sp=$strpos -v ep=$endpos 'BEGIN{cnt=0;}
		{if(NR-1>=sp && NR-1<=ep){
			row[++cnt]=$0; 
			if(NR-1==sp) x=$1; if(NR-1==ep) y=$1; #Get genomic position of the bounding SNPs
		} 
	  	} END{
			#Add header
			if(winidx) printf("%d\t%.3f\tL\n",win,y-x+1); #Add 1 to enable LDhat to consider the last SNP per window
			else printf("%d\t%.3f\tL\n",win,y+1);
			
			#Extract blocks of SNP positions
			for(i=1;i<=win;i++){
				if(winidx) printf("%.4f\n",row[i]-x+0.001); #Second window to last window
				else  printf("%.4f\n",row[i]); #First window
			}
		}' $locs > $splitdir/${chrom}_${i}.locs

	awk -v winidx=$i -v win=$cur_window -v sp=$strpos -v ep=$endpos 'NR==1{print $1"\t"win"\t1";} 
		NR>1{
			if($0~/^>/){
				if(length($0)>30){
					print "ERROR: Sequence name is too long\n";
					exit 1;
				}
				print $0; #Print sequence name 
			}else{
                                lines=int(win/2000); #only 2k characters per line allowed in LDhat
                                for(i=0;i<lines;i++){
                                        print substr($0,sp+i*2000,2000);
                                        if(i+1==lines && win%2000) #Print the remaining characters
						print substr($0,sp+((i+1)*2000),win%2000);
                                }
			}
	     	}' $sites > $splitdir/${chrom}_${i}.sites
done
fi

#RUN LDhat using a job array
sbatch --array=0-$lim runLDhat.sh -i $splitdir/${chrom} -l $lkfile -x 20000000 -y 2000 -b 5 -p $LDout/${chrom}
wait

#RUN LDhat stat on each window
if [ $skip -eq 0 ]; then
cd $LDout
for i in `seq 0 $lim`
do
        stat -input $LDout/${chrom}_${i}.rates.txt -loc $splitdir/${chrom}_${i}.locs -burn 2000
        mv res.txt $LDout/${chrom}_${i}.res.txt
done
fi

#for i in `less /lustre/nobackup/WUR/BIOINF/fuent011/References/Tomato/S_lycopersicum_chromosomes.4.00.chrom`; do inputdir="/lustre/nobackup/WUR/BIOINF/fuent011/Project2.0/Pimp_Genomes/Merged"; /home/WUR/fuent011/Tools/LDhat/convert -seq $inputdir/$i/${i}_merged.ldhat.sites -loc $inputdir/$i/${i}_merged.ldhat.locs | awk -v chr=$i '$1~/Segregating/{n=split($0,a," "); if(n==4){numsnp=a[4];}} $1~/^Watterson/{split($0,a," "); printf("%s\t%d\t%f\t%.6f\n",chr,numsnp,a[4],(a[4]/numsnp)/1e-8);}' -; done

#Merge windows
if [ $skip -eq 0 ]; then
printf "" > $rhoout/${chrom}_rho
offset=0
for i in `seq 0 $lim`
do
        strpos=$((i*window - overlap*$i + 1)) #1-based
        endpos=$((strpos+window - 1))
	cur_window=$window
	if [ $remainder -ne 0 ] && [ $i -eq $lim ]; then
		#Last window contains the remaining SNPs but spans the same size as other window#
                if [ $remainder -gt $overlap ];then
                        strpos=$((pos-window+1)) #take the last (window)-SNPs from the end
                        endpos=$pos #The very last SNP position
                else #otherwise, use same start but shorter window size
                        cur_window=$((window-overlap+remainder))
                        endpos=$((strpos+window-overlap+remainder-1))
                fi
        fi
	if [ $i -ne 0 ]; then 
		offset=$(awk -v sp=$strpos '{if(NR-1==sp) print $1;}' $locs ) 
	fi
        echo -e "Offset:$offset; Window=$i; Startpos:$strpos"
	awk -v os=$offset -v winidx=$i -v win=$cur_window -v overl=$overlap -v lm=$lim -v rem=$remainder 'NR>2{ #offset of 2 rows for res.txt
		if(winidx==0){		
			#first window or single-window chromosome
                        if(NR<=win+2-(overl/2)) printf("%.3f\t%.5f\n",$1,$2);
		}else if(winidx==lm){
			#if(NR==win-rem+(overl/2)+2) print "Last window";
			#only include the remainder positions; exclude additional overlap/2 because the first window has size window-(overlap/2)
			if(rem>overl){ #last block of size equal to other windows
				if(NR>win-rem+(overl/2)+2) printf("%.3f\t%.5f\n",$1+os-0.001,$2);
			}else{ #last block with shorter size
				if(NR>(overl/2)+2) printf("%.3f\t%.5f\n",$1+os-0.001,$2);
			}	
		}else{ #if(winidx)
			#exclude overlap/2 positions from both ends of the window
			if(NR>(overl/2)+2 && NR<=win+2-(overl/2)) printf("%.3f\t%.5f\n",$1+os-0.001,$2);
		} 

	}' $LDout/${chrom}_${i}.res.txt >> $rhoout/${chrom}_rho 
done
fi 

#printf "" > FINAL_results/Pimp/Pimp_ALLCHR_rho; for i in `ls Pimp_Genomes/LDhat/*_rho`; do awk '{tmp=FILENAME; gsub(/^.*\//,"",tmp); gsub("_rho","",tmp); printf("%s\t%d\t%.10f\n", tmp,$1*1000,$2);}'  $i >> FINAL_results/Pimp/Pimp_ALLCHR_rho; done

#printf "" > FINAL_results/Pimp/Pimp_ALLCHR_rho.bed; for i in `ls Pimp_Genomes/LDhat/*_rho`; do awk 'NR==1{tmp=FILENAME; gsub(/^.*\//,"",tmp); gsub("_rho","",tmp); startp=$1*1000; rho=$2;} NR>1{endp=$1*1000; printf("%s\t%d\t%d\t%.10f\n",tmp,startp,endp,rho); startp=$1*1000; rho=$2;}' $i >> FINAL_results/Pimp/Pimp_ALLCHR_rho.bed; done

#printf "" > FINAL_results/Pimp_ALLCHR_rho_100kbwindow; for i in `ls Pimp_Genomes/LDhat/SL4.0ch*_rho`; do printf "" > tmp; awk 'NR==1{tmp=FILENAME; gsub(/^.*\//,"",tmp); gsub("_rho","",tmp); pos=$1; rho=$2;} NR>1{startp=pos*1000; endp=$1*1000; rate=($1-pos)*rho; intervallen=($1-pos)*1000; for(i=startp;i<endp;i++){printf("%s\t%d\t%d\t%.10f\n",tmp,i,i,rate/intervallen);} pos=$1; rho=$2;}' $i >> tmp; awk -v chr=$i 'BEGIN{gsub(/^.*\//,"",chr); gsub("_rho","",chr);} NR>1{if($1==chr) print $0;}' ../References/Tomato/S_lycopersicum_chromosomes.4.00.len | bedtools makewindows -g - -w 100000 -s 10000 > tmpwin; bedtools map -a tmpwin -b tmp -c 4 -o sum >> FINAL_results/Pimp_ALLCHR_rho_100kbwindow; done

#awk '{print $0"\t"$2;}' FINAL_results/Pimp/Pimp_SNPpositions | bedtools intersect -v -a - -b ../References/Tomato/ITAG4.0_RepeatModeler_repeats_light.gff | bedtools intersect -a FINAL_results/Window100kb_step50kb -b - -c > FINAL_results/Pimp/Pimp_SNPs_100kb

#RUN sequenceLDhot for using a job array
sbatch --array=0-$lim runseqLDhot.sh -c $chrom -i $splitdir/${chrom} -r $LDout/${chrom} -w 50000 -p $seqLDout/${chrom}
wait

if [ $skip -eq 1 ]; then
printf "" > $seqout/${chrom}_sum
printf "" > $seqout/${chrom}_bgres
for i in `seq 0 $lim`
do
        strpos=$((i*window - overlap*$i + 1)) #1-based
        endpos=$((strpos+window - 1))
	if [ $remainder -ne 0 ] && [ $i -eq $lim ]; then
		#Last window contains the remaining SNPs but spans the same size as other window#
		if [ $remainder -gt $overlap ];then
			strpos=$((pos-window+1)) #take the last (window)-SNPs from the end
			endpos=$pos #The very last SNP position
		else #otherwise, use same start but shorter window size
			cur_window=$((window-overlap+remainder))
			endpos=$((strpos+window-overlap+remainder-1))
		fi
	fi
        if [ $i -ne 0 ]; then
		offset=$(awk -v sp=$strpos '{if(NR-1==sp) print $1*1000;}' $locs )
		if [ $i -ne $lim ]; then
			if [ $remainder -gt $overlap ];then
                		pos1=$(awk -v sp=$((strpos+overlap/2)) '{if(NR-1==sp) print $1*1000;}' $locs )
				pos2=$(awk -v ep=$((endpos-overlap/2)) '{if(NR-1==ep) print $1*1000;}' $locs )
			else #otherwise, use same start but shorter window size
				pos1=$(awk -v sp=$((strpos+overlap/2)) '{if(NR-1==sp) print $1*1000;}' $locs )
				pos2=$(awk -v ep=$pos '{if(NR-1==ep) print $1*1000;}' $locs )
			fi
		else
			pos1=$(awk -v sp=$((pos-remainder+1+overlap/2)) '{if(NR-1==sp) print $1*1000;}' $locs )
                        pos2=$(awk -v ep=$pos '{if(NR-1==ep) print $1*1000;}' $locs )
		fi
	else
		offset=0
		pos1=$(awk -v sp=$strpos '{if(NR-1==sp) print $1*1000;}' $locs )
                pos2=$(awk -v ep=$((endpos-overlap/2)) '{if(NR-1==ep) print $1*1000;}' $locs )
        fi
        echo -e "$strpos; $endpos; $pos1; $pos2; $offset"
        awk -F" " -v p1=$pos1 -v p2=$pos2 -v os=$offset '{if($1+os>=p1 && $2+os<=p2) print $1+os"\t"$2+os"\t"$3"\t"$4;}' $seqLDout/${chrom}_${i}_input.sum >> $seqout/${chrom}_sum
	awk -F" " -v p1=$pos1 -v p2=$pos2 -v os=$offset '{if($3+os>=p1 && $2+os<=p2) printf("%s\t%d\t%d\t%.10f\n",$1,$2+os,$3+os,$5);}' $LDout/${chrom}_${i}.bgres >> $seqout/${chrom}_bgres
done
fi

exit
#Combine all chrom background recombination rates
cat Pimp_Genomes/seqLD/SL4.0ch*_bgres | sort -k1,1 -k2,2n > FINAL_results/Pimp/Pimp_ALLCHR_bgres

printf "" > FINAL_results/Pimp/Pimp_seqLD; for i in `ls Pimp_Genomes/seqLD/SL4.0ch*_sum`; do awk 'NR==1{LR=$3; RHO=$4; tmp=FILENAME; gsub(/.*\//,"",tmp); gsub("_sum","",tmp);} NR>1{if($3==""){print tmp"\t"$1"\t"$2"\t"LR"\t"RHO;}else{print tmp"\t"$0; LR=$3; RHO=$4;}}' $i >> FINAL_results/Pimp/Pimp_seqLD; done

LR_lim=$(awk '$4>0{print $4;}' FINAL_results/Pimp/Pimp_seqLD | sort -k1,1n | awk '{row[NR]=$0;}END{perc=NR*0.95; for(i=1;i<=NR;i++) if(i>perc){print row[i]; exit;}}' - )

awk -v lim=$LR_lim '$4>lim{print $0;}' FINAL_results/Pimp/Pimp_seqLD | sort -k1,1 -k2,2n | bedtools merge -i - -d 500 -c 4,5 -o max | bedtools map -a - -b FINAL_results/Pimp/Pimp_ALLCHR_rho.bed -c 4 -o max | bedtools map -a - -b FINAL_results/Pimp/Pimp_ALLCHR_bgres -c 4 -o mean | awk '$7>0 && $6/$7>10{print $0;}' - > FINAL_results/Pimp/Pimp_seqLD_hotspots

#Convert p/kb to pb; Sum p per map interval, then divide by interval width (Mb); compute p/Mb
printf "" > Compare2Map/SLVintage_vs_EXPIM2012; for i in `ls SLVintage_Genomes/LDhat/SL4.0ch*_rho`; do printf "" > tempfile; awk 'NR==1{tmp=FILENAME; gsub(/^.*\//,"",tmp); gsub("_rho","",tmp); startp=$1*1000; rho=$2;} NR>1{endp=$1*1000; intervallen=endp-startp; rate=(rho/1000); for(i=startp;i<endp;i++){printf("%s\t%d\t%d\t%.10f\n",tmp,i,i,rate);} startp=$1*1000; rho=$2;}' $i >> tempfile; awk -v chr=$i 'BEGIN{gsub(/^.*\//,"",chr); gsub("_rho","",chr);} $1==chr{if(NR==1){pos=$2;cM=$3;}else{if(cM==$3) next; printf("%s\t%d\t%d\t%.10f\n", $1,pos,$2,(($3-cM)/($2-pos))*1000000); pos=$2; cM=$3;}}' GeneticMap/EXPIM2012_physicalpos.txt > tempfile2; bedtools map -a tempfile2 -b tempfile -c 4 -o sum | awk '{printf("%s\t%d\t%d\t%.10f\n",$1,$2,$3,($4/($3-$2))*1000000); }' >> Compare2Map/SLVintage_vs_EXPIM2012; done

