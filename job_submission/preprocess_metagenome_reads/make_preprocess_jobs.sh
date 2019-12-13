touch 191212_aegea_batch_submit_script.sh
rm 191212_aegea_batch_submit_script.sh
for SAMPLE in `cat sample_list.txt`; do 
cat >> 191212_aegea_batch_submit_script.sh << EOL
echo "Submitting $SAMPLE for preprocessing ... " 2>&1 | tee -a 191212_preprocess_log_file.txt
aegea batch submit --queue microbiome-lowPriority  --retry-attempts 3  --image hwastyk/preprocess_metagenome_reads:latest         --storage /mnt=500         --vcpus 8         --memory 32000         --command="export coreNum=8;         export mem_mb=32000;         export SNAKEFILE=Snakefile;         export CONFIG=config.yaml;         export HOST_INDEX_PATH=s3://czbiohub-microbiome/Sonnenburg_Lab/Fiber_Study_Reference_Data/human_index/;         export SAMPLE=$SAMPLE;         export READPATHS=(s3://czb-seqbot/fastqs/191125_A00111_0398_AHWM7CDSXX/Pilot);         export DATA_DEST=s3://czbiohub-microbiome/Sonnenburg_Lab/hwastyk/201912-fiber-fermented-study;  source preprocess_metagenome_reads.sh" &>> 191212_preprocess_log_file.txt
sleep 10

EOL
done
