#!/usr/bin/env python
from __future__ import print_function
import argparse
import vcf  
import sys

def variant_tuple(record, alt):
    chrom = record.CHROM
    if record.CHROM[:3] == "chr":
        chrom = record.CHROM[3:]
    return (chrom, record.POS, record.REF, str(alt))

def get_info_field(inputfile, fieldname):
    infos = {}
    header = None
    try:
        reader = vcf.Reader(filename=inputfile)
        header = reader.infos[fieldname]
        for record in reader:
            if not fieldname in record.INFO:
                continue
            for alt in record.ALT:
                variant=variant_tuple(record, alt)
                infos[variant] = record.INFO[fieldname]
    except:
        pass
    return header, infos

def main():
    parser = argparse.ArgumentParser(description='Fix dbsnp VP calls and add OXOG filter')
    parser.add_argument('validationvcf', help="Validation vcf file")
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), default=sys.stdin, help="Merged and annotated VCF file (default: stdin)")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout, help="Specify output file (default:stdout)")
    parser.add_argument('-f', '--fieldname', default="Validation_status", help="Specify INFO field name to use (default: Validation_status")
    parser.add_argument('-s', '--skipif', default="LOWDEPTH", help="Comma-delimted list of items which won't get carried over (default: LOWDEPTH)")
    args = parser.parse_args()

    skips = args.skipif.split(',')
    header, infos = get_info_field(args.validationvcf, args.fieldname)
    reader = vcf.Reader(args.input)
    if len(infos) > 0:
        reader.infos[args.fieldname] = header
    writer = vcf.Writer(args.output, reader)

    for record in reader:
        assert len(record.ALT) == 1
        variant = variant_tuple(record, record.ALT[0])
        if variant in infos:
            items = infos[variant]
            apply_items = [item for item in items if item not in skips]
            if apply_items:
                record.INFO[args.fieldname] = ','.join(apply_items)

        writer.write_record(record)
    return 0

if __name__ == "__main__":
    sys.exit(main())
