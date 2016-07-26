#!/bin/bash -l
module load vcfanno
module load tabix/0.2.6
module load python-packages/2

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

readonly OUTDIR=annotated_pre_model/${VARIANT}
mkdir -p $OUTDIR
readonly output_file=${OUTDIR}/${ID}.annotated.${VARIANT}.vcf

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
    > ${output_file}
rm annotation/vaf.indel.${ID}.conf
