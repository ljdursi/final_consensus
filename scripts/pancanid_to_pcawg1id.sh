#!/bin/bash

readonly METADATA=/oicr/data/pancanxfer/consensus/final_consensus/metadata/sample_metadata.csv

if [ $# -eq 0 ] || [ -z "$1" ] 
then
    echo "$0 - converts global pancancer donor id to tumour analysis id, used by PCAWG-1"
    echo "Usage: $0 pancan_id"
    echo "invocation was: $0 $1 $2"
    exit 
fi

readonly ID=$1

grep ${ID} ${METADATA} | awk -F, '{print $16}'

