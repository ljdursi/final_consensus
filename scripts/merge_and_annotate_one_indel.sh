#!/bin/bash -l

module purge

function usage {
    >&2 echo "usage: $0 sample-id "
    >&2 echo "       merges and annotates indel VCFs for a single tumour"
    exit 1
}

readonly ID=$1
if [ -z $ID ]
then
    >&2 echo "argument missing"
    usage
fi

./scripts/merge-one-tumour-indel.sh $ID
./scripts/dbsnp_annotate_one.sh $ID "indel"
./scripts/vaf_annotate_one_indel.sh $ID
