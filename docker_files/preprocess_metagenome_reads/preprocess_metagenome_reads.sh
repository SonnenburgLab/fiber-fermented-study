#!/bin/bash -x

set -euo pipefail

START_TIME=$SECONDS

LOCAL=$(pwd)

# The next two lines are necessary because there are two different /mnt locations
# Without this odd copy step, Snakemake fails (other things do too).
cp -pr * $LOCAL/
cd $LOCAL

export PATH="/opt/conda/bin:${PATH}"

# Setting up bbtools dependency to speed things up
# which bbmap.sh
# cd `which bbmap.sh | sed 's/bbmap.sh//'`../opt/bbmap*/jni/ && make -f makefile.linux -isk && cd $LOCAL
# export JAVA_HOME=/opt/conda && make -f /opt/conda/opt/bbmap*/jni/makefile.linux -C /opt/conda/opt/bbmap*/jni/

coreNum=${coreNum:-16}
snakefile=${SNAKEFILE:-Snakefile}
CONFIG=${CONFIG:-config.yaml}
mem_gb=$(( $mem_mb / 1024 ))
READS_PREFIX=${READS_PREFIX:-$SAMPLE}

echo $SAMPLE
echo ${READS_PREFIX}
echo ${READPATHS[@]}


echo "Getting the specified reference index for host read removal ..."
human_index=s3://czbiohub-microbiome/Sonnenburg_Lab/Fiber_Study_Reference_Data/human_index/
host_index_path=${HOST_INDEX_PATH:-$human_index}
echo ${host_index_path}
mkdir -p $LOCAL/host_index/
aws s3 sync ${host_index_path} $LOCAL/host_index/
echo "Got the host index. Getting all the reads now..."

# Downloading all of the reads
mkdir -p $LOCAL/OUTPUT
x=1
for RUN in ${READPATHS[@]}; do
RUN=${RUN%/}
echo "Getting reads from $RUN ..."
mkdir -p run$x
aws s3 sync ${RUN}/ $LOCAL/run${x}/ --exclude "*" --include "*${READS_PREFIX}*"
echo "Got reads ${READS_PREFIX} from ${RUN}/" | tee -a $LOCAL/OUTPUT/FASTQ_READS_USED.txt
x=$(( $x + 1 ))
done

echo "Combining FASTQ files for ${SAMPLE}..."
cat $LOCAL/run*/${READS_PREFIX}*R1*.gz > $LOCAL/${SAMPLE}_R1.fastq.gz
cat $LOCAL/run*/${READS_PREFIX}*R2*.gz > $LOCAL/${SAMPLE}_R2.fastq.gz
rm -r $LOCAL/run*/*.gz

echo "Running snakemake pipeline ..."
# Snakemake needs a minimum of two processors
time snakemake -p -j $coreNum --resources mem_mb=$mem_mb -s $snakefile --config SAMPLE=$SAMPLE mem_gb=$mem_gb all_mem_mb=$mem_mb all_cpu=$coreNum --configfile $CONFIG


# Cleaning up snakemake output
#mv $LOCAL/OUTPUT/03_NODUP/*fastqc* $LOCAL/OUTPUT/01_FASTQC
cp LATEST $LOCAL/OUTPUT/00_LOGS/VERSION__preprocess_metagenome_reads


#echo "Generating sample summary ..." # maybe make this a Snakemake rule???
#export SampleName=$SAMPLE
#source gen_sample_summary.sh

echo "Syncing data back to S3 ..."
aws s3 sync $LOCAL/OUTPUT/ ${DATA_DEST}/${SAMPLE}/

echo "File transfer complete. Pipeline has finished."
