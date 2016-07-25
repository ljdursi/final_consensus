#!/bin/bash

readonly METADATA=/oicr/data/pancanxfer/consensus/final_consensus/metadata/pcawg_summary_internal_freeze_160323.tsv
readonly ID=$1

if [[ -z "$ID" ]]
then
    >&2 echo "Usage: ${0} id "
    >&2 echo "        Returns 'tcga' for samples that are in tcga, 'icgc' otherwise"
    exit 1
fi

countrycode=$( grep $ID ${METADATA} \
                    | cut -f 2 -d $'\t' \
                    | cut -f 2 -d - )

if [[ "$countrycode" == "US" ]]
then
    echo "tcga" 
else
    echo "icgc" 
fi
