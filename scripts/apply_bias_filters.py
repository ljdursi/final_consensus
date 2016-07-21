#!/usr/bin/env python
from __future__ import print_function
import argparse
import vcf  
import sys

def variant_tuple(record, alt):
    return (record.CHROM, record.POS, record.REF, str(alt))

def get_bias_failures(biasfile, filtername=None):
    failures = {}
    reader = vcf.Reader(biasfile)
    for record in reader:
        if not record.FILTER:
            continue
        for alt in record.ALT:
            variant=variant_tuple(record, alt)
            failures[variant] = filtername if filtername is not None else record.FILTER
    return failures

def main():
    parser = argparse.ArgumentParser(description='Fix dbsnp VP calls and add OXOG filter')
    parser.add_argument('filtervcf', type=argparse.FileType('r'), help="Filter vcf file")
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), default=sys.stdin, help="Merged and annotated VCF file (default: stdin)")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout, help="Specify output file (default:stdout)")
    parser.add_argument('-f', '--filtername', help="Specify filter name to use (default: use filter field from VCF file)")
    parser.add_argument('-d', '--filterdesc', default="", help="Specify description of filter")
    args = parser.parse_args()

    reader = vcf.Reader(args.input)
    if args.filtername is not None:
        reader.filters[args.filtername] = vcf.parser._Filter(id=args.filtername, desc=args.filterdesc)
    writer = vcf.Writer(args.output, reader)

    failuredict = get_bias_failures(args.filtervcf, args.filtername)
    for record in reader:
        assert len(record.ALT) == 1
        variant = variant_tuple(record, record.ALT[0])
        if variant in failuredict:
            if not record.FILTER:
                record.FILTER = failuredict[variant]
            else:
                record.FILTER = record.FILTER + failuredict[variant]

        writer.write_record(record)

    return 0

if __name__ == "__main__":
    sys.exit(main())
