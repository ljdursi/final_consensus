#!/usr/bin/env python
from __future__ import print_function
import argparse
import vcf  
import vcf.parser
import csv
import sys

def variant_tuple(sample, record):
    return (sample, record.CHROM, record.POS)

def get_entries_from_MAF(maf):
    entries = set()
    mafreader = csv.DictReader(maf, delimiter='\t')
    for record in mafreader:
        items = ['tumor_aliquot_id', 'Chromosome', 'Start_position' ]
        if not all([item in record for item in items]):
            continue
        variant=(record['tumor_aliquot_id'], record['Chromosome'], int(record['Start_position']))
        entries.add(variant)
    return entries

def main():
    parser = argparse.ArgumentParser(description='Fix dbsnp VP calls and add OXOG filter')
    parser.add_argument('MAF', type=argparse.FileType('r'), help="MAF file for filtering")
    parser.add_argument('sample', type=str, help="tumour aliquot id")
    parser.add_argument('filtername', type=str, help="Filter name to apply")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout, help="Specify output file (default:stdout)")
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), default=sys.stdin, help="Merged and annotated VCF file (default: stdin)")
    parser.add_argument('-d', '--desc', type=str, default="", help="Description of filter")
    args = parser.parse_args()

    reader = vcf.Reader(args.input)
    reader.filters[args.filtername] = vcf.parser._Filter(id=args.filtername, desc=args.desc)
    writer = vcf.Writer(args.output, reader)

    entries = get_entries_from_MAF(args.MAF)
    for record in reader:
        variant = variant_tuple(args.sample, record)
        if variant in entries:
            if not record.FILTER:
                record.FILTER = [args.filtername]
            else:
                record.FILTER = record.FILTER + [args.filtername]

        writer.write_record(record)

    return 0

if __name__ == "__main__":
    sys.exit(main())
