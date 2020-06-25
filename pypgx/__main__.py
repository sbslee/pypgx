import logging

from .common import get_parser

from .pgkb import pgkb
from .report import report
from .sdf2gdf import sdf2gdf
from .bam2sdf import bam2sdf
from .bam2gdf import bam2gdf

def main():
    parser = get_parser()
    args = parser.parse_args()

    if args.tool == "pgkb":
        pgkb(args)

    elif args.tool == "report":
        report(args)

    elif args.tool == "sdf2gdf":
        sdf2gdf(args)

    elif args.tool == "bam2sdf":
        bam2sdf(args)

    elif args.tool == "bam2gdf":
        bam2gdf(args)

if __name__ == "__main__":
    main()
