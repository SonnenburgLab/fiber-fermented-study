parallel 'aws s3 cp s3://czbiohub-microbiome/Sonnenburg_Lab/hwastyk/201912-fiber-fermented-study/{}/00_LOGS/ QC_SUMMARIES/ --exclude "*" --include "QC*" --recursive' ::: `aws s3 ls s3://czbiohub-microbiome/Sonnenburg_Lab/hwastyk/201912-fiber-fermented-study/ | awk '{print $2}' | sed 's|/||'`

head -n 1 `ls QC_SUMMARIES/*.txt | head -n 1` | cut -f 1-10 > ALL_QC_SUMMARIES.txt
for i in QC_SUMMARIES/*.txt; do sed 1d ${i} | cut -f 1-10 >> ALL_QC_SUMMARIES.txt; done
