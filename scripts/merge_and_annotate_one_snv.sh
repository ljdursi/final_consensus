#!/bin/bash -l

module purge

function usage {
    >&2 echo "usage: $0 sample-id "
    >&2 echo "       merges and annotates SNV VCFs for a single tumour"
    exit 1
}

readonly VARIANT=snv_mnv

readonly ID=$1
if [ -z $ID ]
then
    >&2 echo "argument missing"
    usage
fi

./scripts/merge-one-tumour-snv.sh $ID
./scripts/dbsnp_annotate_one.sh $ID "snv_mnv"
./scripts/vaf_oxog_annotate_one_snv.sh $ID
