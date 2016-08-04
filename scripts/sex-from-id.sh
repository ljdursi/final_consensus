#!/bin/bash

readonly METADATA=/oicr/data/pancanxfer/consensus/final_consensus/metadata/release.tsv
readonly DONORS=/oicr/data/pancanxfer/consensus/final_consensus/metadata/pcawg_donor.csv
readonly ID=$1

if [[ -z "$ID" ]]
then
    >&2 echo "Usage: ${0} id "
    >&2 echo "        Returns gender for given donor id"
    exit 1
fi

readonly donor=$( ./scripts/donor-from-id.sh ${ID} )
sex=$( grep ${donor} ${DONORS} \
        | cut -f 4 -d, \
        | head -n1 )

echo "${sex}"
