#!/bin/bash

while getopts ":n:f:s:e:u" opt; do
  case $opt in
    n) run_name="$OPTARG"
    ;;
    f) flowcell="$OPTARG"
    ;;
    s) start_lane="$OPTARG"
    ;;
    e) end_lane="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

if [[ $# -ne 8 ]]; then
    echo "USAGE: -n NAME -f FLOWCELL -s START_LANE -e END_LANE"
    exit 1
fi

if [ $(echo $HOSTNAME | cut -f1 -d"." | cut -c 1-8) != "transfer" ]; then
    echo "Must be on a transfer node to access the Genetics filesystem!"
    exit 1
fi

echo "name="${run_name}
echo "flowcell="$flowcell
echo "start_lane="${start_lane}
echo "end_lane="${end_lane}

SN_dir=$(find /n/files/Genetics/reichseq/reich/broad/get.broadinstitute.org/pkgs -name "*${flowcell}*" -execdir pwd \; -quit)

if [ ${start_lane} -gt ${end_lane} ]; then
    echo "start_lane should be less than end_lane"
    exit 1
elif [ ${start_lane} -eq ${end_lane} ]; then
    lanes=${start_lane}
elif [ ${start_lane} -lt ${end_lane} ]; then
    lanes=$(seq -s , ${start_lane} ${end_lane})
else
    echo "something is messed up with your lanes..."
    exit 1
fi

date_string=$(echo 20$(less ${SN_dir}/${start_lane}_${flowcell}.${start_lane}.1.fastq.gz | head -1 | cut -f2 -d":" | cut -c 10-15))
scratch_dir=/n/scratch3/users/$(echo $USER | cut -c 1)/$USER/automated_pipeline/${date_string}_${run_name}/fastq

mkdir -p ${scratch_dir}
eval rsync --ignore-existing --progress ${SN_dir}/{${lanes}}_*.fastq.gz ${scratch_dir} || rsync --ignore-existing --progress ${SN_dir}/${start_lane}_*.fastq.gz ${scratch_dir}

if [ ${start_lane} -eq ${end_lane} ]; then
    echo "--name ${run_name} --date_string ${date_string} --illumina_directory ${date_string}_${start_lane}_${flowcell}"
else
    echo "--name ${run_name} --date_string ${date_string} --illumina_directory ${date_string}_${start_lane}-${end_lane}_${flowcell}"
fi
