#!/bin/bash
module load vcfanno

function usage {
    >&2 echo "usage: $0 sample-id {snv_mnv|indel}"
    >&2 echo "       annotates one merged tumor"
    exit 1
}

readonly ID=$1
if [ -z $ID ]
then
    >&2 echo "argument missing"
    usage
fi

readonly VARIANT=${2:-snv_mnv}
if [[ "${VARIANT}" != "snv_mnv" ]] && [[ "${VARIANT}" != "indel" ]]
then
    >&2 echo "Invalid variant type"
    usage
fi

readonly INDIR=merged/${VARIANT}
readonly input_file=${INDIR}/${ID}.merged.${VARIANT}.vcf.gz

if [ ! -f $input_file ] 
then
    >&2 echo "file missing: ${input_file} not found"
    usage
fi

readonly OUTDIR=dbsnp_annotated/${VARIANT}
mkdir -p $OUTDIR
readonly output_file=${OUTDIR}/${ID}.annotated.${VARIANT}.vcf

if [[ -f $output_file ]] && [[ $outfile -nt $input_file ]]
then
    >&2 echo "$0: ${output_file} exists and is newer than inputs; cowardly refusing to overwrite."
    exit 1
fi

function fix_vcfanno_header {
    sed \
        -e 's/^##INFO=<ID=1000genomes_AF.*$/##INFO=<ID=1000genomes_AF,Number=1,Type=Float,Description="Thousand Genomes phase 3 occurance fraction if found: ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz">/' \
        -e 's/^##INFO=<ID=1000genomes_ID.*$/##INFO=<ID=1000genomes_ID,Number=1,Type=String,Description="Thousand Genomes phase 3 ID if found: ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz">/' \
        -e 's/^##INFO=<ID=cosmic,.*$/##INFO=<ID=cosmic,Number=1,Type=String,Description="(first) cosmic ID if found, COSMICv76">/' \
        -e 's/^##INFO=<ID=dbsnp,.*$/##INFO=<ID=dbsnp,Number=1,Type=String,Description="(first) dbSNP ID if found, build 147, All_20160408.vcf.gz">/' \
        -e 's/^##INFO=<ID=repeat_masker,.*$/##INFO=<ID=repeat_masker,Number=1,Type=String,Description="Repeat masker region if in one">/' 
}

vcfanno -p 1 annotation/dbsnp.annotations.conf ${input_file} \
    | fix_vcfanno_header \
    > ${output_file}
