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

readonly dkfz_bias_file=/oicr/data/pancanxfer/consensus/filters/dkfz_bias/${VARIANT}/${ID}.annotated.dkfz_bias.${VARIANT}.vcf.gz
readonly classification_maf=/oicr/data/pancanxfer/consensus/annotations/classifications/${VARIANT}/${ID}.${VARIANT}.maf
readonly remapfilter=/oicr/data/pancanxfer/consensus/filters/realignment/${VARIANT}/${ID}.${VARIANT}.maf
readonly normalpanel=/oicr/data/pancanxfer/consensus/filters/panel_of_normals/${VARIANT}/${ID}.${VARIANT}.maf
readonly GERMLINE_MAF=/oicr/data/pancanxfer/consensus/filters/germline/Broad_germline_site_filter_failed_mutations.tsv
readonly GERMLINE_OVLP_MAF=/oicr/data/pancanxfer/consensus/filters/germline_overlap/somatic_germline_overlap_by_patient.tsv
readonly SIGR1=/oicr/data/pancanxfer/consensus/annotations/pcawg7_artifacts/R1.tsv
readonly SIGR2=/oicr/data/pancanxfer/consensus/annotations/pcawg7_artifacts/R2.tsv
readonly SIGN3=/oicr/data/pancanxfer/consensus/annotations/pcawg7_artifacts/N3.tsv
readonly TIN=/oicr/data/pancanxfer/consensus/annotations/TiN/release_may2016.v1.1.TiN__donor.TiNsorted.20Jul2016.tsv
readonly STARS=/oicr/data/pancanxfer/consensus/annotations/star/PAWG_QC_Summary_of_Measures.tsv

PCAWG1ID=$( ./scripts/pancanid_to_pcawg1id.sh ${ID} )
readonly VALIDATION_FILE=/oicr/data/pancanxfer/validation/vcfs/quality-filtered/${PCAWG1ID}.${VARIANT}.vcf

if [ ! -f ${dkfz_bias_file} ] || [ ! -f ${germline_maf} ] || [ ! -f ${classification_maf} ]
then
    >&2 echo "$0: Warning: ${dkfz_bias_file} or ${germline_maf} or ${classification_maf} not found, not applying filters/annotations"
else
    mv ${output_file} ${tmp_output_file}
    python ./scripts/apply_bias_filters.py ${dkfz_bias_file} -i ${tmp_output_file} \
        | python ./scripts/apply_validation_calls.py ${VALIDATION_FILE} \
        | python ./scripts/filter_by_presence_in_maf.py -d "Evidence of Germline Call Based on 1000Genomes" ${GERMLINE_MAF} ${ID} 1000GENOMES \
        | python ./scripts/filter_by_presence_in_maf.py -d "Overlaps Germline Call" ${GERMLINE_OVLP_MAF} ${ID} GERMLINEOVLP \
        | python ./scripts/filter_by_presence_in_maf.py -d "Presence in Panel of Normals" ${normalpanel} ${ID} NORMALPANEL \
        | python ./scripts/filter_by_presence_in_maf.py --info -d "Sanger Tower: Possible Artifact" ${SIGR1} ${ID} signature_R1 \
        | python ./scripts/filter_by_presence_in_maf.py --info -d "Suspected C>A oxo-guanine signature" ${SIGR2} ${ID} signature_R2 \
        | python ./scripts/filter_by_presence_in_maf.py --info -d "T>A mutation often in bleed through context" ${SIGN3} ${ID} signature_N3 \
        | python ./scripts/info_or_filter_from_MAF.py -a info -c Variant_Classification -d "Variant Classification" ${classification_maf} Variant_Classification \
        | python ./scripts/info_or_filter_from_MAF.py -a filter -d "Variant no longer seen under remapping" ${remapfilter} REMAPFAIL \
        | ./scripts/annotate_vcf_from_tsv_column.sh -c 11 -i ${ID} -n TumourInNormalEstimate -t ${TIN} \
        | ./scripts/annotate_vcf_from_tsv_column.sh -c 5 -i ${ID} -n BAMQCStars -t ${STARS} \
            > ${output_file}
    rm ${tmp_output_file}
fi

bgzip -f ${output_file}
tabix -p vcf ${output_file}.gz