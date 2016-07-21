#!/usr/bin/env python
from __future__ import print_function
import argparse
import vcf  
import vcf.parser
import csv
import sys

def variant_tuple(record, alt):
    return (record.CHROM, record.POS, record.REF, str(alt))

def get_classification_dict_from_MAF(maf):
    classifications = {}
    mafreader = csv.DictReader(maf, delimiter='\t')
    for record in mafreader:
        items = ['Chromosome', 'Start_position', 'Variant_Classification', 'Reference_Allele', 'Tumor_Seq_Allele2' ]
        if not all([item in record for item in items]):
            continue
        variant=(record['Chromosome'], int(record['Start_position']), record['Reference_Allele'], record['Tumor_Seq_Allele2'])
        classifications[variant] = record['Variant_Classification']
    return classifications

def main():
    parser = argparse.ArgumentParser(description='Fix dbsnp VP calls and add OXOG filter')
    parser.add_argument('MAF', type=argparse.FileType('r'), help="MAF file for variant classification annotation")
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), default=sys.stdin, help="Merged and annotated VCF file (default: stdin)")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout, help="Specify output file (default:stdout)")
    args = parser.parse_args()

    infofield = 'Classification'

    reader = vcf.Reader(args.input)
    reader.infos[infofield] = vcf.parser._Info(id=infofield, num='1', type='String', desc='Variant Classification', source=None, version=None)
    writer = vcf.Writer(args.output, reader)
    classification_dict = get_classification_dict_from_MAF(args.MAF)

    for record in reader:
        assert len(record.ALT) == 1
        variant = variant_tuple(record, record.ALT[0])
        if variant in classification_dict:
            record.INFO[infofield] = classification_dict[variant]

        writer.write_record(record)

    return 0

if __name__ == "__main__":
    sys.exit(main())
