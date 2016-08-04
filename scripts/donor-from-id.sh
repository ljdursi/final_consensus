#!/bin/bash

readonly METADATA=/oicr/data/pancanxfer/consensus/final_consensus/metadata/release.tsv
readonly ID=$1

if [[ -z "$ID" ]]
then
    >&2 echo "Usage: ${0} id "
    >&2 echo "        Returns gender for given donor id"
    exit 1
fi

donor=$( grep $ID ${METADATA} \
                | cut -f 1 -d $'\t' )

echo "${donor}" 
