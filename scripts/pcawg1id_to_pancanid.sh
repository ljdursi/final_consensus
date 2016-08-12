#!/bin/bash

readonly METADATA=/oicr/data/pancanxfer/validation/metadata/sample_metadata.csv

if [ $# -eq 0 ] || [ -z "$1" ] 
then
    echo "$0 - converts pcawg1d (tumour analysis id) to global pancancer donor id (tumour analyzed sample UUID)"
    echo "Usage: $0 pcawg1_id"
    echo "invocation was: $0 $1 $2"
    exit 
fi

readonly ID=$1

grep ${ID} ${METADATA} | awk -F, '{print $14}'

