#!/usr/bin/env python
from __future__ import print_function
import argparse
import vcf  
import sys
import numpy

def passed(record):
    """ Did this variant pass all of its filters? """
    return record.FILTER is None or len(record.FILTER) == 0 or \
            record.FILTER == "PASS" or record.FILTER[0] == "PASS" or \
            record.FILTER[0][:4] == "Tier"

def round_three(num):
    scale=1000
    return int(num*scale)*1./scale

def get_vaf(record, broad=False, sanger=False, muse=False, dkfz=False, SNV=True):
    """Returns the VAF corresponding to the record and the caller. """
    if SNV:
        if broad:
            return float(record.INFO['tumor_f'][0])
        if dkfz:
            return record.INFO['AF'][1]
        if sanger:
            return record.samples[1]['PM']
        if muse:
            return 1.0*record.samples[0]['AD'][1]/record.samples[0]['DP']
        return None
    else:  #indel
        if broad:
            try:
                value = float(record.INFO['TFRAC'])
            except:
                value = None
            return value
        if dkfz:
            return float(record.INFO['FR'][0])
        if sanger:
            return None  
        return None

def get_depth(record, broad=False, sanger=False, muse=False, dkfz=False, SNV=True):
    """ Returns the depth corresponding to the record and the caller.
        Indels not yet implemented."""
    if SNV:
        if broad:
            return int(record.samples[0]['alt_count']) + int(record.samples[0]['ref_count'])
        if dkfz:
            return int(record.samples[1]['DP'])
        if sanger:
            keys = ['FAZ', 'RAZ', 'FCZ', 'RCZ', 'FGZ', 'RGZ', 'FTZ', 'RTZ']
            return sum([int(record.samples[1][key]) for key in keys])
        if muse:
            return int(record.samples[0]['DP'])
    else:
        return None

def get_readcount(record, broad=False, sanger=False, muse=False, dkfz=False, SNV=True):
    """ Returns the alt-supporting read counts corresponding to the record and the caller.
        Indels not yet implemented."""
    if SNV:
        if broad:
            return int(record.samples[0]['alt_count'])
        if dkfz:
            return int(record.samples[1]['DP4'][2]) + int(record.samples[1]['DP4'][3])
        if sanger:
            altbase = str(record.ALT[0])
            keys = ['F'+altbase+'Z', 'R'+altbase+'Z']
            return sum([int(record.samples[1][key]) for key in keys])
        if muse:
            return int(record.samples[0]['AD'][1])
    else:
        return None

def populate_dict(filename, *args, **kwargs):
    """
    Read each variant in filename, creating a dictionary
    with key (chrom, pos, ref, alt) and value VAF
    """
    vaf_dict = {}
    readcount_dict = {}
    depth_dict = {}
    if filename is None:
        return {}
    vcf_reader = vcf.Reader(open(filename, 'r'))
    for variant in vcf_reader:
        if not passed(variant):
            continue
        vaf = get_vaf(variant, *args, **kwargs)
        rc  = get_readcount(variant, *args, **kwargs)
        dp  = get_depth(variant, *args, **kwargs)
        for alt in variant.ALT:
            key = (variant.CHROM, variant.POS, variant.REF, str(alt))
            vaf_dict[key] = vaf
            readcount_dict[key] = rc
            depth_dict[key] = dp
    return vaf_dict, readcount_dict, depth_dict

def pick_best(vafs, dps, rcs):
    idxs = numpy.argsort(vafs)
    assert len(vafs) == len(dps) == len(rcs)
    nvafs = len(vafs)

    # odd number: choose the median
    if nvafs % 2 == 1:
        best = idxs[(nvafs-1)//2]
        return vafs[best], dps[best], rcs[best]

    # even number: choose the best
    case1 = nvafs // 2
    case2 = case1 - 1

    # in case of tie in depths, go by lowest read count:
    if dps[case2] == dps[case1]:
        if rcs[case2] < rcs[case1]:
            best = case2
        else:
            best = case1
    # otherwise go by lowest depth:
    elif dps[case2] < dps[case1]:
        best = case2
    else:
        best = case1

    return vafs[best], dps[best], rcs[best]
    

def main():
    parser = argparse.ArgumentParser(description='Annotate merged vcf with VAF information where available')
    parser.add_argument('mergedvcf', type=argparse.FileType('r'), default=sys.stdin, help="Merged VCF file")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout, help="Specify output file (default:stdout)")
    parser.add_argument('-b', '--broad', type=str, help="Broad file")
    parser.add_argument('-d', '--dkfz', type=str, help="DKFZ file")
    parser.add_argument('-s', '--sanger', type=str, help="Sanger file")
    parser.add_argument('-m', '--muse', type=str, help="Muse file")
    parser.add_argument('-i', '--indel', action='store_true',
                        help="Variant type == indel (default:snv_mnv)")
    args = parser.parse_args()
    snvs = not args.indel

    vaf_dicts, readcount_dicts, depth_dicts = \
            tuple(zip(*[populate_dict(args.broad, broad=True, SNV=snvs),
                        populate_dict(args.dkfz, dkfz=True, SNV=snvs),
                        populate_dict(args.sanger, sanger=True, SNV=snvs),
                        populate_dict(args.muse, muse=True, SNV=snvs)]))

    vcf_reader = vcf.Reader(args.mergedvcf)
    vcf_writer = vcf.Writer(args.output, vcf_reader)
    for variant in vcf_reader:
        key = variant.CHROM, variant.POS, variant.REF, str(variant.ALT[0])
        vafs = [vaf_dict[key]
                for vaf_dict in vaf_dicts
                if key in vaf_dict
                if vaf_dict[key] is not None]
        dps = [depth_dict[key]
                for depth_dict in depth_dicts
                if key in depth_dict
                if depth_dict[key] is not None]
        rcs = [rc_dict[key]
                for rc_dict in readcount_dicts
                if key in rc_dict
                if rc_dict[key] is not None]

        for oldkey in ['t_alt_count', 't_ref_count']:
            if oldkey in variant.INFO:
                del variant.INFO[oldkey]

        if len(dps)==len(rcs) and len(rcs)==len(vafs) and len(vafs) > 0:
            variant.INFO['caller_depths'] = dps
            variant.INFO['caller_variant_readcounts'] = rcs
            roundvafs = [round_three(vaf) for vaf in vafs]
            variant.INFO['caller_VAFs'] = roundvafs

            variant.INFO['VAF'], variant.INFO['t_depth'], variant.INFO['t_alt_count'] = pick_best(vafs, dps, rcs)

        vcf_writer.write_record(variant)

    return 0

if __name__ == "__main__":
    sys.exit(main())
