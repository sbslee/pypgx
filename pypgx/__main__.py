import sys

from .common import get_parser, logging
from .pgkb import pgkb
from .report import report
from .sdf2gdf import sdf2gdf
from .bam2sdf import bam2sdf
from .bam2gdf import bam2gdf
from .minivcf import minivcf
from .merge import merge

logger = logging.getLogger(__name__)

def _pgkb(args):
    result = pgkb(args.t)
    if args.o:
        with open(args.o, "w") as f:
            f.write(result)
    else:
        sys.stdout.write(result)

def _report(args):
    f = open(args.gt)
    result = report(f)
    f.close()
    if args.o:
        with open(args.o, "w") as f:
            f.write(result)
    else:
        sys.stdout.write(result)

def _bam2gdf(args):
    result = bam2gdf(args.tg, args.cg, args.bam)
    if args.o:
        with open(args.o, "w") as f:
            f.write(result)
    else:
        sys.stdout.write(result)

def _minivcf(args):
    result = minivcf(args.vcf, args.region)
    if args.o:
        with open(args.o, "w") as f:
            f.write(result)
    else:
        sys.stdout.write(result)

def _merge(args):
    result = merge(args.vcf, args.r)
    if args.o:
        with open(args.o, "w") as f:
            f.write(result)
    else:
        sys.stdout.write(result)

def main():
    parser = get_parser()
    args = parser.parse_args()

    logger.info(f"""Command: '{" ".join(sys.argv)}'""")
    
    if args.tool == "pgkb":
        _pgkb(args)

    elif args.tool == "report":
        _report(args)

    elif args.tool == "sdf2gdf":
        sdf2gdf(args)

    elif args.tool == "bam2sdf":
        bam2sdf(args)

    elif args.tool == "bam2gdf":
        _bam2gdf(args)

    elif args.tool == "minivcf":
        _minivcf(args)

    elif args.tool == "merge":
        _merge(args)

if __name__ == "__main__":
    main()
