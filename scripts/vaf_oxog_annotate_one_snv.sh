#!/bin/bash -l
module load vcfanno
module load tabix/0.2.6
module load python-packages/2

function usage {
    >&2 echo "usage: $0 sample-id "
    >&2 echo "       annotates one merged snv tumor vcf"
    exit 1
}

readonly VARIANT=snv_mnv
readonly INDIR=dbsnp_annotated/${VARIANT}
readonly DKFZBIASDIR=/oicr/data/pancanxfer/preliminary_final_release_dkfz_bias/${VARIANT}

readonly ID=$1
if [ -z $ID ]
then
    >&2 echo "argument missing"
    usage
fi

readonly input_file=${INDIR}/${ID}.annotated.${VARIANT}.vcf

if [ ! -f $input_file ] 
then
    >&2 echo "file missing: ${input_file} not found"
    usage
fi

readonly OUTDIR=annotated/${VARIANT}
mkdir -p $OUTDIR
readonly tmp_output_file=${OUTDIR}/${variant}/tmp.${ID}.annotated.${VARIANT}.vcf
readonly output_file=${OUTDIR}/${variant}/${ID}.annotated.${VARIANT}.vcf
sed -e "s/@@SAMPLE@@/${ID}/" annotation/vaf_oxog.annotations.conf.template > annotation/vaf.${ID}.conf

if [[ -f $output_file ]] && [[ $outfile -nt $input_file ]]
then
    >&2 echo "$0: ${output_file} exists and is newer than inputs; cowardly refusing to overwrite."
    exit 1
fi

vcfanno -p 1 annotation/vaf.${ID}.conf ${input_file} > ${tmp_output_file}
rm annotation/vaf.${ID}.conf
python ./scripts/clean_snv_calls.py ${tmp_output_file} -o ${output_file}
rm ${tmp_output_file}

readonly dkfz_bias_file=${DKFZBIASDIR}/${ID}.annotated.dkfz_bias.${VARIANT}.vcf.gz
readonly germline_maf=/oicr/data/pancanxfer/consensus/filters/germline/Broad_germline_site_filter_failed_mutations.tsv
readonly classification_maf=/oicr/data/pancanxfer/consensus/annotations/classifications/${VARIANT}.maf

if [ ! -f ${dkfz_bias_file} ] || [ ! -f ${germline_maf} ] || [ ! -f ${classification_maf} ]
then
    >&2 echo "$0: Warning: ${dkfz_bias_file} or ${germline_maf} or ${classification_maf} not found, not applying filters/annotations"
else
    mv ${output_file} ${tmp_output_file}
    python ./scripts/apply_bias_filters.py ${dkfz_bias_file} -i ${tmp_output_file} \
        | python ./scripts/omit_germline_from_maf.py ${germline_maf} ${ID} \
        | python ./scripts/variant_classification_from_maf.py ${classification_maf} \
            -o ${output_file}
    rm ${tmp_output_file}
fi

bgzip -f ${output_file}
tabix -p vcf ${output_file}.gz
