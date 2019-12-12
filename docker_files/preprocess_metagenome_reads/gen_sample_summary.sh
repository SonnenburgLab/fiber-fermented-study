#!/usr/bin/env bash
# This script gathers info from log files and makes a sample summary table
# It requires doing "export SampleName={sample} in the shell snakemake rule.
RAW_READS=$((`grep "Input:" /mnt/OUTPUT/00_LOGS/TRIM_MARKED_${SampleName}.err.txt | awk '{print $2}'`/2)) # Divide by 2
RAW_BASES=`grep "Input:" /mnt/OUTPUT/00_LOGS/TRIM_MARKED_${SampleName}.err.txt | awk '{print $4}'`
TRIM_READS=$((`grep "Result:" /mnt/OUTPUT/00_LOGS/TRIM_MARKED_${SampleName}.err.txt | awk '{print $2}'`/2)) # Divide by 2
TRIM_BASES=`grep "Result:" /mnt/OUTPUT/00_LOGS/TRIM_MARKED_${SampleName}.err.txt | awk '{print $5}'`
UNIQUE_READS=$((`grep "Reads Out:" /mnt/OUTPUT/00_LOGS/NODUP_HMN_UNMAPPED_TRIM_MARKED_${SampleName}.err.txt | awk '{print $3}'`/2)) # Divide by 2
UNIQUE_BASES=`grep "Bases Out:" /mnt/OUTPUT/00_LOGS/NODUP_HMN_UNMAPPED_TRIM_MARKED_${SampleName}.err.txt | awk '{print $3}'`
DUP_READS=$((`grep "Reads Out:" /mnt/OUTPUT/00_LOGS/DUP_HMN_UNMAPPED_TRIM_MARKED_${SampleName}.err.txt | awk '{print $3}'`/2)) # Divide by 2
DUP_BASES=`grep "Bases Out:" /mnt/OUTPUT/00_LOGS/DUP_HMN_UNMAPPED_TRIM_MARKED_${SampleName}.err.txt | awk '{print $3}'`
HOST_READS=$((TRIM_READS - UNIQUE_READS - DUP_READS))
K1_SIZE=`grep "Estimated genome size:" /mnt/OUTPUT/00_LOGS/MASH_k1_NODUP_HMN_UNMAPPED_TRIM_MARKED_${SampleName}.err.txt | awk '{print $4}'`
K1_COV=`grep "Estimated coverage:" /mnt/OUTPUT/00_LOGS/MASH_k1_NODUP_HMN_UNMAPPED_TRIM_MARKED_${SampleName}.err.txt | awk '{print $3}'`
K2_SIZE=`grep "Estimated genome size:" /mnt/OUTPUT/00_LOGS/MASH_k2_NODUP_HMN_UNMAPPED_TRIM_MARKED_${SampleName}.err.txt | awk '{print $4}'`
K2_COV=`grep "Estimated coverage:" /mnt/OUTPUT/00_LOGS/MASH_k2_NODUP_HMN_UNMAPPED_TRIM_MARKED_${SampleName}.err.txt | awk '{print $3}'`

printf "SampleName\traw_reads\traw_bases\ttrimmed_reads\ttrimmed_bases\tunique_reads\tunique_bases\tduplicate_reads\tduplicate_bases\thost_reads\tmash_k1_genome_size\tmash_k1_genome_coverage\tmash_k2_genome_size\tmash_k2_genome_coverage\n" > /mnt/OUTPUT/00_LOGS/QC_SUMMARY_${SampleName}.txt
printf "$SampleName\t$RAW_READS\t$RAW_BASES\t$TRIM_READS\t$TRIM_BASES\t$UNIQUE_READS\t$UNIQUE_BASES\t$DUP_READS\t$DUP_BASES\t$HOST_READS\t$K1_SIZE\t$K1_COV\t$K2_SIZE\t$K2_COV\n" >> /mnt/OUTPUT/00_LOGS/QC_SUMMARY_${SampleName}.txt
