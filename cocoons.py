#!/bin/python

import argparse
import re
import subprocess

import vcf_tools

parser = argparse.ArgumentParser()

parser.add_argument('--vcf', type=str,
        help='The input vcf file')

parser.add_argument('--maf', type=str,
        help='The input maf file')

parser.add_argument('--distance', type=int, default=5,
        help='Maximum distance between two mutations to be considered [default: 5]')

parser.add_argument('--mapq', type=int, default=10,
        help='Minimum read mapping quality; used if given bam [default: 10]')

parser.add_argument('--snv-only', action='store_true',
        help='Only SNVs are considered [default: no]')

parser.add_argument('--min-cooccuring-freq', type=float, default=.8,
        help='Minimum cooccuring frequency of the adjacent mutations; used if given bam [default: 0.8]')

parser.add_argument('--max-vaf-diff', type=int, default=10,
        help='Maximum variant allele frequency (VAF) difference to define DNP; used if not given bam [default:10]')

parser.add_argument('--bam', type=str,
        help='File containing bam file locations for the samples used in MAF')

parser.add_argument('--threads', type=int, default=4,
        help='Maximum number of processes')



parser.add_argument("--out", type=str, default='cocoons',
        help='output file prefix')

args = parser.parse_args()


# def filter_input_vcf(input_fp, compressed_input, position_filter_fp, variant_filter_fp=None, out_prefix=''):
#     tool_args = ['--positions', args.position_filter]
# 
#     if compressed_input:
#         tool_args += ['--gzvcf', input_fp]
#     else:
#         tool_args += ['--vcf', input_fp]
# 
#     tool_args += ['--out', out_prefix + '.filtered']
#     tool_args += ['--recode', '--recode-INFO-all']
# 
#     print('position filter')
#     subprocess.check_output(['vcftools'] + tool_args)
# 
#     if variant_filter_fp:
#         tool_args = [out_prefix + '.filtered.recode.vcf', variant_filter_fp]
#         tool_args += ['--out', out_prefix + '.variant-filtered.vcf']
# 
#         print('variant filter')
#         subprocess.check_output(['python', 'filtering_tool.py'] + tool_args)
#         
#         return out_prefix + '.variant-filtered.vcf'
# 
#     return  out_prefix + '.filtered.recode.vcf'
def get_input_file(cli_args):
    """Returns input file and boolean indicating whether input file is .vcf

    Returns:
        - input_fp
            filepath to input file
        - is_vcf
            bool indicating whether input file is vcf or not
    """
    if args.vcf:
        return args.vcf, True
    elif args.maf:
        return args.maf, False

    raise ValueError('Invalid inputs - must provide a .vcf or .maf file as input')
    




def main():
    input_fp, is_vcf = get_input_file(args)

    possible_cooccuring = vcf_tools.get_cooccuring_variants(input_fp, 5, threads=4)

    print(len(list(possible_cooccuring.values())[0]))


if __name__ == '__main__':
    main()
