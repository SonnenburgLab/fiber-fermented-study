# Variables to set before running this script
# NO TRAILING SLASHES!!!
coreNum=16
mem_mb=64000
SNAKEFILE=Snakefile
HOST_INDEX_PATH=s3://czbiohub-microbiome/Sonnenburg_Lab/Fiber_Study_Reference_Data/human_index/
SAMPLE=Hadza_MoBio_hadza-I-L_N_4_1315
STUDY=Hadza
READPATHS=(s3://czbiohub-microbiome/Original_Sequencing_Data/180615_A00111_0161_BH5LMWDSXX/Bryan_Merrill s3://czbiohub-microbiome/Original_Sequencing_Data/180815_A00111_0190_BH7J22DSXX/${STUDY} s3://czbiohub-microbiome/Original_Sequencing_Data/181002_A00111_0217_BHCC2KDSXX s3://czbiohub-microbiome/Original_Sequencing_Data/181203_A00111_0236_BHCWTNDSXX/Bryan_Merrill s3://czbiohub-microbiome/Original_Sequencing_Data/181203_A00111_0235_AHFHTHDSXX/Bryan_Merrill s3://czb-seqbot/fastqs/190122_A00111_0259_AHHC2MDSXX_Redemux/190122_A00111_0259_AHHC2MDSXX/rawdata/Bryan_Merrill_DualIndex s3://czb-seqbot/fastqs/190517_A00111_0312_AHKJW3DSXX/rawdata/Bryan_Merrill s3://czb-seqbot/fastqs/190517_A00111_0313_BHKL2LDSXX/rawdata/Bryan_Merrill)
DATA_DEST=s3://czbiohub-microbiome/Sonnenburg_Lab/Hadza_Nepal_Metagenomics/FINAL_PROJECT_DATA/${STUDY}
CONFIG=config.yaml

# Submit script
aegea batch submit --queue microbiome-lowPriority \
--image bmerrill9/preprocess_metagenome_reads:latest \
--storage /mnt=500 \
--vcpus $coreNum \
--memory $mem_mb \
--command="\
#Set time zone \
ln -fs /usr/share/zoneinfo/America/Los_Angeles /etc/localtime/; \
dpkg-reconfigure --frontend noninteractive tzdata; \
export coreNum=$coreNum; \
export mem_mb=${mem_mb}; \
export SNAKEFILE=$SNAKEFILE; \
export CONFIG=$CONFIG; \
export export HOST_INDEX_PATH=$HOST_INDEX_PATH; \
export READS_PREFIX=${READS_PREFIX}; \
export SAMPLE=$SAMPLE; \
export READPATHS=$READPATHS; \
export DATA_DEST=${DATA_DEST}; \
source preprocess_metagenome_reads.sh"

# Delete host reads? y/n
# Add config json variable to specify which file you want to use

# LATEST
aegea batch submit --queue microbiome-highPriority         --image bmerrill9/preprocess_metagenome_reads:latest         --storage /mnt=500         --vcpus 16         --memory 64000         --command="        ln -fs /usr/share/zoneinfo/America/Los_Angeles /etc/localtime;         dpkg-reconfigure --frontend noninteractive tzdata;         export coreNum=16;         export mem_mb=64000;         export SNAKEFILE=Snakefile;         export CONFIG=config.yaml;         export HOST_INDEX_PATH=s3://czbiohub-microbiome/Sonnenburg_Lab/Fiber_Study_Reference_Data/human_index/;         export SAMPLE=Hadza_MoBio_hadza-I-L_N_4_1315;         export READPATHS=(s3://czbiohub-microbiome/Original_Sequencing_Data/180615_A00111_0161_BH5LMWDSXX/Bryan_Merrill s3://czbiohub-microbiome/Original_Sequencing_Data/180815_A00111_0190_BH7J22DSXX/Hadza s3://czbiohub-microbiome/Original_Sequencing_Data/181002_A00111_0217_BHCC2KDSXX s3://czbiohub-microbiome/Original_Sequencing_Data/181203_A00111_0236_BHCWTNDSXX/Bryan_Merrill s3://czb-seqbot/fastqs/190122_A00111_0259_AHHC2MDSXX_Redemux/190122_A00111_0259_AHHC2MDSXX/rawdata/Bryan_Merrill_DualIndex s3://czb-seqbot/fastqs/190517_A00111_0312_AHKJW3DSXX/rawdata/Bryan_Merrill s3://czb-seqbot/fastqs/190517_A00111_0313_BHKL2LDSXX/rawdata/Bryan_Merrill);         export DATA_DEST=s3://czbiohub-microbiome/Sonnenburg_Lab/Hadza_Nepal_Metagenomics/FINAL_PROJECT_DATA/Hadza;  source preprocess_metagenome_reads.sh" 2>> logs/20190525_233221__preprocess_metagenome_reads/AEGEA.err 1>> logs/20190525_233221__preprocess_metagenome_reads/AEGEA.out


# MUST USE source preprocess_metagenome_reads.sh and NOT ./preprocess_metagenome_reads.sh or else the bash array ${READPATHS} will not be accessible inside preprocess_metagenome_reads.sh !!!!!!!!!!


# qtrim_reads.sh
aegea batch submit --queue microbiome-highPriority         --image bmerrill9/preprocess_metagenome_reads:latest         --storage /mnt=500         --vcpus 16         --memory 64000         --command="\
export coreNum=16;         export mem_mb=64000; \
export TRIMQ=30; export MINLEN=55; \
export READS_PATH=s3://czbiohub-microbiome/Sonnenburg_Lab/Hadza_Nepal_Metagenomics/FINAL_PROJECT_DATA/Hadza/Hadza_MoBio_hadza-I-L_N_4_1315/03_MERGED/; \
./qtrim_reads.sh"


# Must set:
# mem_mb, coreNum, TRIMQ, MINLEN, 
# READS_PATH=