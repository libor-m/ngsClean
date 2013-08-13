#!/bin/bash

ID=$1
N_IND=$2
CHR=$3
BAM_FILES=$4
REF_SEQ=$5
ANC_SEQ=$6


BASE_DIR=/space/tmp1/fgvieira
DATA_DIR=$BASE_DIR

SCRIPTS_DIR=$BASE_DIR/scripts
APPZ_DIR=$BASE_DIR/appz
ANGSD_DIR=$APPZ_DIR/angsd0.204_ORI
SAMTOOLS_DIR=$APPZ_DIR/samtools-0.1.18


if [[ -z $ID ]]; then
    echo "ERROR: second field (ID) must be a non-null string!"
    exit
elif [[ ! $N_IND =~ ^[0-9]+$ ]]; then
    echo "ERROR: third field (number of individuals) must be an integer!"
    exit
elif [[ -z $CHR ]]; then
    echo "ERROR: forth field (CHR) must be a non-null string!"
    exit
elif [[ -z $BAM_FILES ]] || [ ! -s $BAM_FILES ]; then
    echo "ERROR: fifth field (file with path to all BAM files) must be an existing non-zero file!"
    exit
elif [[ -z $REF_SEQ ]] || [ ! -s $REF_SEQ ]; then
    echo "ERROR: sixth field (file with REF sequence) must be an existing non-zero file!"
    exit
elif [[ -z $ANC_SEQ ]] || [ ! -s $ANC_SEQ ]; then
    echo "ERROR: seventh field (file with ANCESTRAL sequence) must be an existing non-zero file!"
    exit
fi



# SAMTools options
MPILEUP_OPTS='-AIDS -q 0 -Q 20 -C 50'
# SNPcleaner options
VARFILTER_OPTS="-v -a 0 -k $((N_IND/2)) -u 2 -h 0 -H 1e-4"

if [[ $CHR != "allchr" ]]; then
    MPILEUP_OPTS=$MPILEUP_OPTS" -r "$CHR
fi



######################
### Filtering Step ###
######################
if [[ ! -s $ID.bed ]]; then
    echo "==> Left alignemnt reads (TODO)"

    echo "==> Local realignment with GATK (TODO)"
    #./gatk -I $i  -T RealignerTargetCreator -R hg19/hg19Chr.fa -o $ONAME.intervals -known $KIND -known $KIND2 -nt 3
    #./gatk -T IndelRealigner -R hg19/hg19Chr.fa -I $i -o $ONAME -targetIntervals $ONAME.intervals -known $KIND -known $KIND2  -LOD 0.4 -model KNOWNS_ONLY -compress 0 --disable_bam_indexing


    echo "==> Filter out individuals with extreme low/high depths (TODO)"
    

    echo "==> Parse BAM files into BCF (apply read quality filters)"
    time $SAMTOOLS_DIR/samtools mpileup -u $MPILEUP_OPTS -f $REF_SEQ -b $BAM_FILES | $SAMTOOLS_DIR/bcftools/bcftools view -gI - | pbzip2 -p5 -cz --best > $ID.raw.vcf.bz2
    
    
    echo "==> Set depth limits based on empirical distribution"
    time pbzip2 -p5 -dc $ID.raw.vcf.bz2 | perl $SCRIPTS_DIR/NGS/SNPcleaner.pl $VARFILTER_OPTS -A $ANC_SEQ -p $ID.filter_out1.vcf.bz2 | pbzip2 -p5 -cz --best > /tmp/$ID.TMP.vcf.bz2
    time pbzip2 -p5 -dc /tmp/$ID.TMP.vcf.bz2 | tr ";" "\t" | awk 'BEGIN{print "chr_pos\tdepth"} !/#/{sub("DP=","",$8); if(rand()<=0.05) print $1"_"$2"\t"$8}' > $ID.sdepth
    time Rscript --vanilla --slave $SCRIPTS_DIR/get_depth_thresh.R --in_file $ID.sdepth --out_file $ID.sdepth.fit.pdf --rnd_sample 20000 > /tmp/$ID.depth_limits.tsv
    MIN_DEPTH=`tail -n 1 /tmp/$ID.depth_limits.tsv | cut -f 1 -d " "`
    if [[ $MIN_DEPTH -lt 3 ]]; then
        MIN_DEPTH=3
    fi
    MAX_DEPTH=`tail -n 1 /tmp/$ID.depth_limits.tsv | cut -f 2 -d " "`
    echo -e "\t=> using $MIN_DEPTH and $MAX_DEPTH as depth cut-offs"
    
    
#    echo "==> Estimate siteF and do LRT between lkl(F=0) and lkl(F=1) (TODO)"
#    time nice $ANGSD_DIR/angsd.g++ -nThreads 15 -ref $REF_SEQ -anc $ANC_SEQ -outfiles $ID.NOFILTER -doHWE 1 -doMaf 2 -doMajorMinor 1 mpileup -g $MPILEUP_OPTS -f $REF_SEQ -b $BAM_FILES > /dev/null
    
    
    echo "==> Filter positions (based on quality, depth, bias, HWE, ...)"
    time bzcat /tmp/$ID.TMP.vcf.bz2 | perl $SCRIPTS_DIR/NGS/SNPcleaner.pl $VARFILTER_OPTS -d $MIN_DEPTH -D $MAX_DEPTH -p $ID.filter_out2.vcf.bz2 | pbzip2 -p5 -cz --best > $ID.vcf.bz2
    
    
    echo "==> Make BED file with \"good\" positions"
    time pbzip2 -p5 -dc $ID.vcf.bz2 | awk '!/#/ {print $1"\t"$2}' > $ID.keep
    time pbzip2 -p5 -dc $ID.vcf.bz2 | awk '!/#/ {print $1"\t"$2-1"\t"$2}' | $APPZ_DIR/bedtools/bin/mergeBed -i stdin > $ID.bed
else
    echo 'BED file '$ID.bed' already exists. Skipping filtering step...'
fi


############################
### Preliminary Analysis ###
############################
if [[ ! -s $ID.sfs.pdf ]]; then
    ## Do fold SFS?
    FOLD=""
    if [[ $REF_SEQ == $ANC_SEQ ]]; then
	N_IND=$(echo "scale=1; $N_IND/2" | bc)
        FOLD="-dofold 1"
    fi


    echo "==> Estimating \"realSFS 1\" on selected positions (BED)"
    time $ANGSD_DIR/angsd.g++ -chunkSize 100000 -ref $REF_SEQ -anc $ANC_SEQ -outfiles $ID -realSFS 1 $FOLD mpileup -g $MPILEUP_OPTS -f $REF_SEQ -b $BAM_FILES -l $ID.bed > /dev/null


    echo "==> Estimating SFS prior (optimSFS)"
    rm -f $ID.sfs.ml*
    N_CHR=$(echo "scale=0;($N_IND*2)/1" | bc)
    time $ANGSD_DIR/misc/optimSFS.gcc -nThreads 10 -binput $ID.sfs -nChr $N_CHR -outnames $ID.sfs -nSites 30000000

    
    echo "==> Parse realSFS (assume highest posterior prob categ to be true)"
    N_SFS_CATEG=$(echo "(2*$N_IND+1)/1" | bc)
    time hexdump -v -e "$N_SFS_CATEG/8 \"%.10g\t\"\"\n\"" $ID.sfs | perl -n -a -e 'next LINE if($F[0] =~ m/nan/i); $sum=0; $sum+=exp($_) for @F; for($i=0; $i<=$#F; $i++){ $tot[$i]+=exp($F[$i])/$sum }; END{$sum_tot+=$_ for @tot; $_/=$sum_tot for @tot; print join("\t",@tot)."\n"}' > $ID.sfs.sum


    echo "==> Plot SFS"
    seq -s $'\t' 0 $N_CHR | cat - $ID.sfs.sum $ID.sfs.ml | R --slave --vanilla -e 'library(ggplot2); d<-t(read.table("stdin",header=T,check.names=F)[-1]); colnames(d)<-c("realSFS","optimSFS"); d<-rename(melt(d), c(X1="categ",X2="SFS")); ggplot(d) + geom_bar(aes(x=categ, y=value, fill=SFS), stat = "identity", position = "dodge")'; mv Rplots.pdf $ID.sfs.pdf
elif [[ -s $ID.sfs.pdf ]]; then
    echo 'SFS plot '$ID.sfs.pdf' already exists. Skipping preliminary analysis step...'
fi

echo "==> Clean up"
rm -f $ID.arg $ID.sfs.ml.local* $ID.sfs.ml.log $ID.NOFILTER.*
