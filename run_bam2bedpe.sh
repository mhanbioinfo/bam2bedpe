#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=ming.han@uhn.ca
#SBATCH -D ./logs_slurm/
#SBATCH -t 1-00:00:00
#SBATCH --mem=16G
#SBATCH -J bam2bedpe
#SBATCH -p himem
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -o %j-%x.out
#SBATCH -e %j-%x.err

if [ $# != 1 ]; then
    echo "Please specify a config.yml file."
    exit 1
fi

## PARSE CONFIG ######################################################

## https://github.com/ash-shell/yaml-parse/blob/master/lib/yaml_parse.sh
## yml file must be 2 spaces indented
function parse_yaml {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         gsub(/\s*#.*$/, "", $3);
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}

CONFIG_F=$1
params=($(parse_yaml $CONFIG_F))

for param in "${params[@]}"; do
   #echo $param
   lhs=${param%=*}
   rhs=${param#*=}
   rhs="${rhs%%\#*}"    # Del in line right comments
   rhs="${rhs%%*( )}"   # Del trailing spaces
   rhs="${rhs%\"*}"     # Del opening string quotes 
   rhs="${rhs#\"*}"     # Del closing string quotes 
   declare $lhs=$rhs
done


## Main program ####################################################

source ${CONDA_ACTIVATE} ${CONDA_ENV}

echo "Processing bam2bedpe..."
echo "output path:      $OUT_DIR"
echo "number of chunks: $NUM_OF_CHUNKS"
echo "processing on:    $SLURMLOCAL"
echo ""

SRC_DIR="$(pwd)"

OLDIFS=$IFS
IFS=','
[ ! -f $SAMPLESHEET_PATH ] && { echo "$SAMPLESHEET_PATH file not found"; exit 99; }
sed 1d $SAMPLESHEET_PATH | \
while read -r SAMPLE_NAME INPUT_BAM_PATH; do


    #SAMPLE_NAME="toy01"
    #INPUT_BAM_PATH="/cluster/home/t110409uhn/git/cfmedip_medremix_bedpe_git/toy_files/toy01.bam"
    
    echo "Processing sample: ${SAMPLE_NAME}"
    echo "Input bam path:    ${INPUT_BAM_PATH}"
    echo "Job started at "$(date) 
    time1=$(date +%s)
    
    PY_SCRIPT_DIR=${SRC_DIR}
    PY_SCRIPT_PATH="${PY_SCRIPT_DIR}/bam2bedpe_pysam.py"
    INPUT_DIR="${INPUT_BAM_PATH%/*}"
    INPUT_BAM="${INPUT_BAM_PATH##*/}"
    
    TMP_DIR="${OUT_DIR}/tmp_dir/${INPUT_BAM%.*}"
    mkdir -p $TMP_DIR
    
    OUT_FRAG_NAMES="${INPUT_BAM%.*}.fragNames"
    OUT_MERGED_SORTD_BEDPE="${INPUT_BAM%.*}_coordSortd.bedpe"
    
    ## -------------------------------------- ##
    ## get fragNames
    samtools view ${INPUT_DIR}/${INPUT_BAM} \
        | cut -f1 \
        | awk '{ a[$1]++ } END { for (b in a) { print b } }' \
        > ${TMP_DIR}/${OUT_FRAG_NAMES}
    
    ## -------------------------------------- ##
    ## split bam into chunks
    NUM_ROWS=$(cat ${TMP_DIR}/${OUT_FRAG_NAMES} | wc -l)
    CHUNK_SIZE=$(( $NUM_ROWS / $NUM_OF_CHUNKS ))   ## won't have decimal, sometimes will be short, last chunk has to go to last row
    echo "bam file has $NUM_ROWS fragments."
    echo "Number of chunks was set to ${NUM_OF_CHUNKS}."
    echo "Each chunk will be $CHUNK_SIZE rows."
    
    for ((i=0; i<=$(( $NUM_OF_CHUNKS - 2 )); i++)); do
        CHUNK=$(( i + 1 ))
        ROW_START=$(( i*CHUNK_SIZE + 1))
        ROW_END=$(( CHUNK*CHUNK_SIZE ))
        echo "Processing CHUNK${CHUNK}... starting with row ${ROW_START}, ending in row ${ROW_END}."
        sed -n "$ROW_START,$ROW_END p" ${TMP_DIR}/${OUT_FRAG_NAMES} \
            > ${TMP_DIR}/${OUT_FRAG_NAMES}.CHUNK${CHUNK}
    
    done
    
    ## since chunk might not be evenly divided, last chunk will just be to last row
    ith=$(( $NUM_OF_CHUNKS - 1 ))
    ROW_START=$(( ith*CHUNK_SIZE + 1 ))
    echo "Processing CHUNK$NUM_OF_CHUNKS... starting with row ${ROW_START}, ending in row ${NUM_ROWS}."
    sed -n "$ROW_START,$NUM_ROWS p" ${TMP_DIR}/${OUT_FRAG_NAMES} \
        > ${TMP_DIR}/${OUT_FRAG_NAMES}.CHUNK${NUM_OF_CHUNKS}
    
    ## make log dir
    mkdir -p ${OUT_DIR}/logs_slurm
    
    ## -------------------------------------- ##
    ## picard FilterSamReads on chunks
    for CHUNK in ${TMP_DIR}/${OUT_FRAG_NAMES}.CHUNK*; do
        echo "Picard FilterSamReads for ${CHUNK##*/}..."
        echo "Output is ${INPUT_BAM%.*}_${CHUNK##*.}.bam"
        
        echo "Writing out ${INPUT_BAM%.*}_bam2bedpe_${CHUNK##*.}.sh"
        ## write out sample sbatch script
        cat <<- EOF > "${TMP_DIR}/${INPUT_BAM%.*}_bam2bedpe_${CHUNK##*.}.sh" 
	#!/bin/bash
	#SBATCH -t 3-00:00:00
	#SBATCH -J bam2bedpe_chunks_${CHUNK##*.}
	#SBATCH -D ${OUT_DIR}/logs_slurm
	#SBATCH --mail-type=ALL
	#SBATCH --mail-user=ming.han@uhn.ca
	#SBATCH -p himem
	#SBATCH -c 1
	#SBATCH --mem=16G
	#SBATCH -o ./%j-%x.out
	#SBATCH -e ./%j-%x.err
	
	source ${CONDA_ACTIVATE} ${CONDA_ENV}
	
	echo "Job started at "\$(date) 
	time1=\$(date +%s)
	
	picard FilterSamReads \
	    -I ${INPUT_DIR}/${INPUT_BAM} \
	    -O ${TMP_DIR}/"${INPUT_BAM%.*}_${CHUNK##*.}.bam" \
	    -READ_LIST_FILE ${CHUNK} \
	    -FILTER includeReadList \
	    -WRITE_READS_FILES false \
	    -USE_JDK_DEFLATER true \
	    -USE_JDK_INFLATER true \
	
	python ${PY_SCRIPT_PATH} \
	    --sort_bam_by_qname \
	    --bam_input ${TMP_DIR}/"${INPUT_BAM%.*}_${CHUNK##*.}.bam" \
	    --bedpe_output ${TMP_DIR}/"${INPUT_BAM%.*}_${CHUNK##*.}.bedpe"
	
	rm ${TMP_DIR}/"${INPUT_BAM%.*}_${CHUNK##*.}.bam"
	
	time2=\$(date +%s)
	echo "Job ended at "\$(date) 
	echo "Job took \$(((time2-time1)/3600)) hours \$((((time2-time1)%3600)/60)) minutes \$(((time2-time1)%60)) seconds"
	EOF
    
        if [ $SLURMLOCAL == "slurm" ]; then
            sbatch "${TMP_DIR}/${INPUT_BAM%.*}_bam2bedpe_${CHUNK##*.}.sh"
        elif [ $SLURMLOCAL == "local" ]; then
            nohup bash "${TMP_DIR}/${INPUT_BAM%.*}_bam2bedpe_${CHUNK##*.}.sh" &> "${TMP_DIR}/${INPUT_BAM%.*}_bam2bedpe_${CHUNK##*.}.log" &
        fi
    
    done
    
    echo $TMP_DIR
    echo $OUT_DIR
    echo $OUT_MERGED_SORTD_BEDPE
    
    
    ## periodically check if chunks have been written to completely
    while true; do
        if [ $(find ${TMP_DIR} -maxdepth 1 -mmin +3 -type f -regex ".*_CHUNK[1-9][0-9]*\.bedpe" | wc -l) -eq $NUM_OF_CHUNKS ] 
        then
            echo "All bedpe chunks have been written to, merging bedpe chunks..."
            find ${TMP_DIR} -maxdepth 1 -mmin +3 -type f -regex ".*_CHUNK[1-9][0-9]*\.bedpe" -exec cat {} + \
                | sort -k1,1V -k2,2n -k3,3n -k5,5n -k6,6n \
                | gzip -c \
                > ${OUT_DIR}/${OUT_MERGED_SORTD_BEDPE}.gz
            break 1
        else
            echo "Sleeping for 5 minutes, then recheck if bedpe chunks were written to completely."
            sleep 5m
        fi
    done
    
    ## remove TMP_DIR
    if "$KEEP_TMP"; then
        echo "Keeping temp dir"
    else
        echo "Removed temp dir after completion"
        rm -rf ${TMP_DIR}
    fi
    
    
    time2=$(date +%s)
    echo "Job ended at "$(date) 
    echo "Job took $(((time2-time1)/3600)) hours $((((time2-time1)%3600)/60)) minutes $(((time2-time1)%60)) seconds"
    echo ""


done
IFS=$OLDIFS


echo "Finished processing bam2bedpe for all samples."
