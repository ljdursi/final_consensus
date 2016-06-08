#!/usr/bin/env python
from __future__ import print_function
import argparse
import vcf  
import sys

def main():
    parser = argparse.ArgumentParser(description='Fix dbsnp VP calls and add OXOG filter')
    parser.add_argument('inputvcf', type=argparse.FileType('r'), default=sys.stdin, help="Merged and annotated VCF file")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout, help="Specify output file (default:stdout)")
    args = parser.parse_args()

    reader = vcf.Reader(args.inputvcf)
    writer = vcf.Writer(args.output, reader)

    for record in reader:
        new_info = {}

        # copy some records over directly
        for item in ['dbsnp', 'cosmic', 'Callers', 'NumCallers', 'repeat_masker', '1000genomes_AF', '1000genomes_ID', 'VAF', 't_alt_count', 't_ref_count']:
            if item in record.INFO:
                new_info[item] = record.INFO[item]

        if 'dbsnp_VP' in record.INFO:
            qualbyte = int(record.INFO['dbsnp_VP'],16) & 255
            somatic = (qualbyte & 2**5) > 0
            if somatic:
                new_info['dbsnp_somatic'] = True

        if 'OXOG_Fail' in record.INFO:
            if record.INFO['OXOG_Fail'] == 'True':
                if record.FILTER is None:
                    record.FILTER = ['OXOGFAIL']
                else:
                    record.FILTER.append('OXOGFAIL')

        record.INFO = new_info
        writer.write_record(record)

    return 0

if __name__ == "__main__":
    sys.exit(main())
