import argparse
import logging

from .pgkb import pgkb
from .report import report

def main():
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(
        dest='tool',
        metavar='tool',
        help='name of tool',
    )
    subparsers.required = True

    pgkb_parser = subparsers.add_parser(
        "pgkb",
        help="extract CPIC guidelines using PharmGKB API",
    )
    pgkb_parser.add_argument(
        "-o",
        metavar="FILE",
        help="output to FILE [stdout]",
    )

    report_parser = subparsers.add_parser(
        "report",
        help="create HTML report using data from Stargazer",
    )
    report_parser.add_argument(
        "gt",
        help="Stargazer genotype file",
    )
    report_parser.add_argument(
        "-o",
        metavar="FILE",
        help="output to FILE [stdout]",
    )

    args = parser.parse_args()

    if args.tool == "pgkb":
        pgkb(args)

    elif args.tool == "report":
        report(args)

if __name__ == "__main__":
    main()
