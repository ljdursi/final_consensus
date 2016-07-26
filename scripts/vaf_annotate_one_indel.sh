#!/bin/bash -l
module load vcfanno
module load tabix/0.2.6

function usage {
    >&2 echo "usage: $0 sample-id "
    >&2 echo "       annotates one merged indel tumor vcf"
    exit 1
}

readonly VARIANT=indel
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
readonly output_file=${OUTDIR}/${variant}/${ID}.annotated.${VARIANT}.vcf
readonly tmp_output_file=${OUTDIR}/${variant}/tmp.${ID}.annotated.${VARIANT}.vcf

template_file=annotation/vaf.indel.annotations.single.template
rm -f annotation/vaf.indel.${ID}.conf
for caller in broad dkfz sanger smufin
do
    if grep ${caller} $input_file > /dev/null
    then
        sed -e "s/@@SAMPLE@@/${ID}/" -e "s/@@CALLER@@/${caller}/" $template_file >> annotation/vaf.indel.${ID}.conf
        echo " " >> annotation/vaf.indel.${ID}.conf
    fi
done

if [[ -f $output_file ]] && [[ $outfile -nt $input_file ]]
then
    >&2 echo "$0: ${output_file} exists and is newer than inputs; cowardly refusing to overwrite."
    exit 1
fi

vcfanno -p 1 annotation/vaf.indel.${ID}.conf ${input_file} \
    | sed -e 's/\tFORMAT.*$//' \
    | sed -e 's/^##INFO=<ID=TumorVAF,.*$/##INFO=<ID=TumorVAF,Number=1,Type=Float,Description="VAF of variant in tumor from sga">/' \
    | sed -e 's/^##INFO=<ID=TumorVarDepth,.*$/##INFO=<ID=TumorVarDepth,Number=1,Type=Integer,Description="Tumor alt count from sga">/' \
    | sed -e 's/^##INFO=<ID=TumorTotalDepth,.*$/##INFO=<ID=TumorTotalDepth,Number=1,Type=Integer,Description="Tumor total read depth from sga">/' \
    | sed -e 's/^##INFO=<ID=NormalVAF,.*$/##INFO=<ID=NormalVAF,Number=1,Type=Float,Description="VAF of variant in normal from sga">/' \
    | sed -e 's/^##INFO=<ID=NormalVarDepth,.*$/##INFO=<ID=NormalVarDepth,Number=1,Type=Integer,Description="Normal alt count from sga">/' \
    | sed -e 's/^##INFO=<ID=NormalTotalDepth,.*$/##INFO=<ID=NormalTotalDepth,Number=1,Type=Integer,Description="Normal total read depth from sga">/' \
    > ${tmp_output_file}
rm annotation/vaf.indel.${ID}.conf

readonly classification_maf=/oicr/data/pancanxfer/consensus/annotations/classifications/${VARIANT}/${ID}.${VARIANT}.maf
readonly normalpanel=/oicr/data/pancanxfer/consensus/filters/panel_of_normals/${VARIANT}/${ID}.${VARIANT}.maf
readonly TIN=/oicr/data/pancanxfer/consensus/annotations/TiN/release_may2016.v1.1.TiN__donor.TiNsorted.20Jul2016.tsv
readonly STARS=/oicr/data/pancanxfer/consensus/annotations/star/PAWG_QC_Summary_of_Measures.tsv
PCAWG1ID=$( ./scripts/pancanid_to_pcawg1id.sh ${ID} )
readonly VALIDATION_FILE=/oicr/data/pancanxfer/validation/vcfs/quality-filtered/${PCAWG1ID}.${VARIANT}.vcf

python ./scripts/clean_indel_calls.py ${tmp_output_file} \
    | python ./scripts/apply_validation_calls.py ${VALIDATION_FILE} \
    | python ./scripts/filter_by_presence_in_maf.py -d "Presence in Panel of Normals" ${normalpanel} ${ID} NORMALPANEL \
    | python ./scripts/info_or_filter_from_MAF.py -a info -c Variant_Classification -d "Variant Classification" ${classification_maf} Variant_Classification \
    | ./scripts/annotate_vcf_from_tsv_column.sh -c 11 -i ${ID} -n TumourInNormalEstimate -t ${TIN} \
    | ./scripts/annotate_vcf_from_tsv_column.sh -c 5 -i ${ID} -n BAMQCStars -t ${STARS} \
        > ${output_file}
bgzip -f ${output_file}
tabix -p vcf ${output_file}.gz
rm ${tmp_output_file} 
