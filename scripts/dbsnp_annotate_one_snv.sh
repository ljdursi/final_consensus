#!/bin/bash
module load vcfanno

function usage {
    >&2 echo "usage: $0 sample-id "
    >&2 echo "       annotates one merged snv tumor"
    exit 1
}

readonly VARIANT=snv_mnv
readonly INDIR=merged/${VARIANT}

readonly ID=$1
if [ -z $ID ]
then
    >&2 echo "argument missing"
    usage
fi

readonly input_file=${INDIR}/${ID}.merged.${VARIANT}.vcf.gz

if [ ! -f $input_file ] 
then
    >&2 echo "file missing: ${input_file} not found"
    usage
fi

readonly OUTDIR=dbsnp_annotated/${VARIANT}
mkdir -p $OUTDIR
readonly output_file=${OUTDIR}/${variant}/${ID}.annotated.${VARIANT}.vcf

if [[ -f $output_file ]] && [[ $outfile -nt $input_file ]]
then
    >&2 echo "$0: ${output_file} exists and is newer than inputs; cowardly refusing to overwrite."
    exit 1
fi

vcfanno -p 1 annotation/dbsnp.annotations.conf ${input_file} > ${output_file}
