#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=ming.han@uhn.ca
#SBATCH -t 1-00:00:00
#SBATCH -D ./logs_slurm/
#SBATCH --mem=16G
#SBATCH -J bedpe7plus12validation
#SBATCH -p himem
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -o %j-%x.out
#SBATCH -e %j-%x.err

# getopts ###################################################
usage(){
    echo 
    echo "Usage: bash bedpe7plus12validation.sh -i BEDPE_INPUT"
    echo 
}
no_args="true"

## Help 
Help()
{
    # Display Help
    echo 
    echo "Checks if .bedpe.gz follows BEDPE7+12 format (including correct number of columns,tab-delimited,appropriate values for each column)."
    echo
    echo "Usage: bash bedpe7plus12validation.sh -i BEDPE_INPUT"
    echo "options:"
    echo "-h   [HELP]      print help"
    echo "-i   [REQUIRED]  input .bedpe.gz file (full path)"
    echo
}

## Get the options
while getopts ":hi:" option; do
    case "${option}" in
        h) Help
           exit;;
        i) INPUT_BEDPE=${OPTARG};;
       \?) echo "Error: Invalid option"
           exit;;
    esac
    no_args="false"
done

[[ "$no_args" == "true" ]] && { usage; exit 1; }

echo "Input .bedpe.gz: $INPUT_BEDPE"

# Main program ##############################################

echo "Job started at "$(date) 
time1=$(date +%s)

NUM_COLS=19
NUM_LINES=$(zcat ${INPUT_BEDPE} | wc -l)
#echo $NUM_LINES

LINE_NUM_NF_neq19=$(zcat ${INPUT_BEDPE} | awk -F'\t' 'BEGIN { ORS="," } NF != 19 {print NR}')
if [ ${#LINE_NUM_NF_neq19} -eq 0 ]; then
    echo "All rows have required 19 fields, tab-delimited."
else
    echo "These rows do not have the correct number of fields, and/or is not tab-delimited: ${LINE_NUM_NF_neq19%,*}."
fi

BAD_CHR_ROWS=$(zcat ${INPUT_BEDPE} | awk -F'\t' 'BEGIN { ORS="," } !(($1 ~/[A-Za-z]+[1-9]+[0-9]*|\./) && ($4 ~/[A-Za-z]+[1-9]+[0-9]*|\./)) {print NR}')
if [ ${#BAD_CHR_ROWS} -eq 0 ]; then
    echo "All rows have correct format for CHROMOSOME fields."
else
    echo "These rows have INCORRECT format for CHROMOSOME fields: ${BAD_CHR_ROWS%,*}."
fi

BAD_START_END_ROWS=$(zcat ${INPUT_BEDPE} | awk -F'\t' 'BEGIN { ORS="," } !(($2 ~/^([0-9]+)$|\-1/) && ($3 ~/^([0-9]+)$|\-1/) && ($5 ~/^([0-9]+)$|\-1/) && ($6 ~/^([0-9]+)$|\-1/)) {print NR}')
if [ ${#BAD_START_END_ROWS} -eq 0 ]; then
    echo "All rows have correct format for START and END fields."
else
    echo "These rows have INCORRECT format for START and END fields: ${BAD_START_END_ROWS%,*}."
fi

BAD_MAPQ_STRAND_ROWS=$(zcat ${INPUT_BEDPE} | awk -F'\t' 'BEGIN { ORS="," } !(($8 >= 0 && $8 < 2^8) && ($9 >= 0 && $9 < 2^8) && ($10 == "-" || $10 == "+" || $10 == ".") && ($11 == "-" || $11 == "+" || $11 == ".")) {print NR}')
if [ ${#BAD_MAPQ_STRAND_ROWS} -eq 0 ]; then
    echo "All rows have correct format for MAPQ and STRAND fields."
else
    echo "These rows have INCORRECT format for MAPQ and STRAND fields: ${BAD_MAPQ_STRAND_ROWS%,*}."
fi

BAD_CIGAR_ROWS=$(zcat ${INPUT_BEDPE} | awk -F'\t' 'BEGIN { ORS="," } !(($12 ~ /*|([0-9]+[MIDNSHPX=])+/) && ($13 ~ /*|([0-9]+[MIDNSHPX=])+/)) {print NR}')
if [ ${#BAD_CIGAR_ROWS} -eq 0 ]; then
    echo "All rows have correct format for CIGAR fields."
else
    echo "These rows have INCORRECT format for CIGAR fields: ${BAD_CIGAR_ROWS%,*}."
fi

BAD_FLAG_ROWS=$(zcat ${INPUT_BEDPE} | awk -F'\t' 'BEGIN { ORS="," } !(($14 >= 0 && $14 < 2^16) && ($15 >= 0 && $15 < 2^16)) {print NR}')
if [ ${#BAD_FLAG_ROWS} -eq 0 ]; then
    echo "All rows have correct format for FLAG field."
else
    echo "These rows have INCORRECT format for FLAG field: ${BAD_FLAG_ROWS%,*}."
fi

BAD_TLEN_ROWS=$(zcat ${INPUT_BEDPE} | awk -F'\t' 'BEGIN { ORS="," } !(($16 ~ /^(-*[0-9]+)$/) && ($17 ~ /^(-*[0-9]+)$/)) {print NR}')
if [ ${#BAD_TLEN_ROWS} -eq 0 ]; then
    echo "All rows have correct format for TLEN field."
else
    echo "These rows have INCORRECT format for TLEN field: ${BAD_TLEN_ROWS%,*}."
fi

BAD_NM_ROWS=$(zcat ${INPUT_BEDPE} | awk -F'\t' 'BEGIN { ORS="," } !(($18 ~ /NA|^([0-9]+)$/) && ($19 ~ /NA|^([0-9]+)$/)) {print NR}')
if [ ${#BAD_NM_ROWS} -eq 0 ]; then
    echo "All rows have correct format for NM field."
else
    echo "These rows have INCORRECT format for NM field: ${BAD_NM_ROWS%,*}."
fi


EMPTY_FIELDS=($(zcat ${INPUT_BEDPE} | awk -F'\t' 'BEGIN { ORS=" " } {for(N=1; N<=NF; N++) if($N=="" || $N ~/\s+/) print NR}'))
EMPTY_FIELDS_UNIQ_ROWS=$(printf "%s\n" "${EMPTY_FIELDS[@]}" | sort -V -u | tr '\n' ',')
EMPTY_FIELDS_UNIQ_ROWS2=${EMPTY_FIELDS_UNIQ_ROWS%,*}
if [ ${#EMPTY_FIELDS_UNIQ_ROWS2} -eq 0 ]; then
    echo "All rows have NO EMPTY or field with only WHITESPACES."
else
    echo "These rows have cells that are EMPTY or contain only WHITESPACES: ${EMPTY_FIELDS_UNIQ_ROWS%,*}."
fi


time2=$(date +%s)
echo "Job ended at "$(date) 
echo "Job took $(((time2-time1)/3600)) hours $((((time2-time1)%3600)/60)) minutes $(((time2-time1)%60)) seconds"

## EOF
