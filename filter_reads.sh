#!/bin/bash
################################################################
# a script for read QC (adaptor/duplicate)                     #
# external dependencies: trimmomatic, cutadapt, flash, picard  #
# Based on script by S. Singhal, sonal.singhal1 [at] gmail.com #
################################################################


#################
# General Paths #
#################

BASE_DIR=/space/F4/fgvieira/appz/ngsClean
EXTERN_DIR=$BASE_DIR/external_progs
ADAPT_DIR=$BASE_DIR/suppl_files

SCYTHE_DIR=$EXTERN_DIR/scythe
SICKLE_DIR=$EXTERN_DIR/sickle

CUTADAPT_DIR=$EXTERN_DIR/cutadapt/bin

FLASH_DIR=$EXTERN_DIR/FLASH-1.2.6
PICARD_DIR=$EXTERN_DIR/picard-tools-1.90


###########################################################
########## Don't change anything pass this point ##########
###########################################################


#################
# Input Options #
#################
RUN=$1
FASTQ1=$2
FASTQ2=$3

if [[ ! -s $FASTQ1 ]]; then
    echo "ERROR: Cannot find FASTQ1 file: "$FASTQ1
    exit
fi



#############
# Filenames #
#############
ID=`basename $FASTQ1`
ID=${ID%\.gz}
ID=${ID%\.*}
ID=${ID%_*}

R1p=$ID"_1.TRIM.fq"
R2p=$ID"_2.TRIM.fq"
Ru=$ID"_U.TRIM.fq"

OUT1p=$ID"_1.FINAL.fq"
OUT2p=$ID"_2.FINAL.fq"
OUTu=$ID"_U.FINAL.fq"



###########
# Options #
###########
N_THREADS=10
PHRED=33
MIN_READ_LEN=$((`zcat $FASTQ1 | head -n 2 | tail -n 1 | wc -c`/3))
MIN_QUAL=10


#################
# Read Cleaning #
#################

if [[ -s $FASTQ2 ]]; then
    echo "========== Paired-End reads =========="
    echo "===> Removing adaptors..."
    $SCYTHE_DIR/scythe -a $ADAPT_DIR/TruSeq2-PE.fa -M 0 -p 0.1 $FASTQ1 | pigz --best > $ID"_1.ADAPT.fq.gz"
    $SCYTHE_DIR/scythe -a $ADAPT_DIR/TruSeq2-PE.fa -M 0 -p 0.1 $FASTQ2 | pigz --best > $ID"_2.ADAPT.fq.gz"
    
    echo "===> Trimming..."
    $SICKLE_DIR/sickle pe -t sanger --pe-file1 $ID"_1.ADAPT.fq.gz" --pe-file2 $ID"_2.ADAPT.fq" -q $MIN_QUAL -l $MIN_READ_LEN --output-pe1 $R1p --output-pe2 $R2p --output-single $Ru
else
    echo "========== Single-End reads =========="
    echo "===> Removing adaptors..."
    $SCYTHE_DIR/scythe -a $ADAPT_DIR/TruSeq2-SE.fa -M 0 -p 0.1 $FASTQ1 | pigz --best > $ID"_U.ADAPT.fq.gz"

    echo "===> Trimming..."
    $SICKLE_DIR/sickle se -t sanger --fastq-file $ID"_U.ADAPT.fq.gz" -q $MIN_QUAL -l $MIN_READ_LEN --output-file $Ru
fi



if [[ $RUN == "FAST" ]]; then
    if [[ -s $FASTQ2 ]]; then
	mv $R1p $OUT1p
	mv $R2p $OUT2p
    fi
    mv $Ru $OUTu
else
    echo "====> Performing exhaustive adaptor search... (NOT IMPLEMENTED YET)"
    #for FILE in $R1p $R2p $Ru
    #do
         #$CUTADAPT_DIR/cutadapt --quality-base $PHRED -f fastq -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -e 0.1 -O 5 -n 5 -m 0 --match-read-wildcards $FILE > /tmp/$FILE
    #done


    echo "===> Converting FASTQ to SAM..."
    java -jar $PICARD_DIR/FastqToSam.jar QUALITY_FORMAT=Standard SAMPLE_NAME=$ID FASTQ=$Ru OUTPUT=/tmp/$ID.u.sam
    cp /tmp/$ID.u.sam /tmp/$ID.sam
    if [[ -s $FASTQ2 ]]; then
	java -jar $PICARD_DIR/FastqToSam.jar QUALITY_FORMAT=Standard SAMPLE_NAME=$ID FASTQ=$R1p FASTQ2=$R2p OUTPUT=/tmp/$ID.p.sam
	java -jar $PICARD_DIR/MergeSamFiles.jar INPUT=/tmp/$ID.p.sam INPUT=/tmp/$ID.u.sam SORT_ORDER=queryname ASSUME_SORTED=true OUTPUT=/tmp/$ID.sam
    fi
    rm -f /tmp/$ID.u.sam /tmp/$ID.p.sam


    echo "===> Removing low complexity regions..."
    # Remove low complexity regions (stretches of same base >= 30% read length)
    perl -an -e '$tot = int(0.3*length($F[9])); print if($F[9] =~ m/A{$tot}|T{$tot}|C{$tot}|G{$tot}/gi)' /tmp/$ID.sam > $ID.REM.mono
    # Remove low complexity regions (read base composition: any base > 50%; N > 30%)
    perl -an -e 'my %n; $l=length($F[9]); map{$n{$_}++} split("",uc($F[9])); print if($n{A}>=0.5*$l || $n{T}>=0.5*$l || $n{C}>=0.5*$l || $n{G}>=0.5*$l || $n{N}>=0.3*$l)' /tmp/$ID.sam > $ID.REM.compos
    # Join all filter results
    cut -f 1 $ID.REM.mono $ID.REM.compos | sort -u | fgrep -v '@' > $ID.REM.read_id
    
    # Remove reads that failed filters
    java -jar $PICARD_DIR/FilterSamReads.jar INPUT=/tmp/$ID.sam FILTER=excludeReadList READ_LIST_FILE=$ID.REM.read_id WRITE_READS_FILES=false OUTPUT=/tmp/$ID.filtered.sam
    # Convert SE back to FASTQ
    awk '/^@/ || $2 == 4' /tmp/$ID.filtered.sam | java -jar $PICARD_DIR/SamToFastq.jar INPUT=/dev/stdin FASTQ=$OUTu
    
    if [[ -s $FASTQ2 ]]; then
	# Convert PE back to FASTQ
	awk '/^@/ || $2 != 4' /tmp/$ID.filtered.sam | java -jar $PICARD_DIR/SamToFastq.jar INPUT=/dev/stdin FASTQ=$OUT1p SECOND_END_FASTQ=$OUT2p

	echo "===> Merging reads (FLASH)..."
	$FLASH_DIR/flash $OUT1p $OUT2p -t $N_THREADS -M 70 -x 0.1 -o $ID
	cat $ID.extendedFrags.fastq >> $OUTu
	mv $ID.notCombined_1.fastq $OUT1p
	mv $ID.notCombined_2.fastq $OUT2p
	rm -f $ID.extendedFrags.fastq $ID.hist $ID.histogram
    fi
fi



############
# Clean-up #
############
echo "====> Clean up..."
if [[ -s $FASTQ2 ]]; then
    pigz --best $OUT1p
    pigz --best $OUT2p
fi
pigz --best $OUTu
rm -f /tmp/$ID.sam /tmp/$ID.filtered.sam $ID*.trim.fastq $ID.REM.*
