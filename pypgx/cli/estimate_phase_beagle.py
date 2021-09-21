import sys

from ..api import utils

import fuc
import pysam

description = f"""
##########################################################################
# Estimate haplotype phase of observed variants with the Beagle program. #
##########################################################################

If your input data is GRCh37, I recommend using the 1000 Genomes Project phase 3 reference panel. You can easily download it thanks to the authors of Beagle:

$ wget -r --no-parent http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/

Usage examples:
  $ pypgx {fuc.api.common._script_name()} imported-variants.zip ref.vcf phased-variants.zip
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Estimate haplotype phase of observed variants with the Beagle program.',
        description=description,
    )
    parser.add_argument(
        'imported_variants',
        metavar='imported-variants',
        help='Archive file with the semantic type VcfFrame[Imported].'
    )
    parser.add_argument(
        'panel',
        help='Reference haplotype panel.'
    )
    parser.add_argument(
        'phased_variants',
        metavar='phased-variants',
        help='Archive file with the semantic type VcfFrame[Phased].'
    )
    parser.add_argument(
        '--impute',
        action='store_true',
        help='Whether to perform imputation of missing genotypes.'
    )

def main(args):
    result = utils.estimate_phase_beagle(
        args.imported_variants, args.panel, impute=args.impute
    )
    result.to_file(args.phased_variants)
