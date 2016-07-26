#!/bin/bash -l

export PATH=/u/jdursi/.local/bin:${PATH}
module load python-packages/2
module load tabix/0.2.6

function usage {
    >&2 echo "usage: $0 sample-id "
    >&2 echo "       merges snv VCFs for sample $0"
    exit 1
}

readonly VARIANT=snv_mnv

readonly ID=$1
if [ -z $ID ]
then
    >&2 echo "argument missing"
    usage
fi

readonly INDIR=input_vcfs/processed/${VARIANT}

readonly musefile=${INDIR}/${ID}.muse.snv_mnv.vcf.gz  
readonly broadfile=${INDIR}/${ID}.broad.snv_mnv.vcf.gz
readonly dkfzfile=${INDIR}/${ID}.dkfz.snv_mnv.vcf.gz
readonly sangerfile=${INDIR}/${ID}.sanger.snv_mnv.vcf.gz

if [ ! -f $musefile ] || [ ! -f $broadfile ] || [ ! -f $dkfzfile ] || [ ! -f $sangerfile ]
then
    >&2 echo "files missing: one of ${musefile} ${broadfile} ${dkfzfile} ${sangerfile} not found"
    usage
fi

readonly OUTDIR=merged/${VARIANT}
mkdir -p ${OUTDIR}

newest=$musefile
if [[ $broadfile -nt $newest ]]; then newest=$broadfile; fi
if [[ $dkfzfile -nt $newest ]]; then newest=$dkfzfile; fi
if [[ $sangerfile -nt $newest ]]; then newest=$sangerfile; fi

readonly outfile=${OUTDIR}/${ID}.merged.${VARIANT}.vcf

if [[ -f "${outfile}.gz" ]] && [[ "${outfile}.gz" -nt "${newest}" ]]
then
    >&2 echo "$0: ${outfile} exists and is newer than inputs; cowardly refusing to overwrite."
    exit 1
fi

mergevcf -l broad,dkfz,muse,sanger \
    ${broadfile} ${dkfzfile} ${musefile} ${sangerfile} \
    --ncallers --mincallers 2 \
    -o ${outfile}.tmp

grep "^##" ${outfile}.tmp > ${outfile}
echo '##FILTER=<ID=LOWSUPPORT,Description="Not called by enough callers in ensemble">' >> ${outfile}
grep "^#C" ${outfile}.tmp >> ${outfile}
grep -v "^#" ${outfile}.tmp \
    | grep -v "Callers=muse;" \
    | sort -k1,1d -k2,2n \
    >> ${outfile}
rm ${outfile}.tmp

bgzip -f ${outfile}
tabix -p vcf ${outfile}.gz
