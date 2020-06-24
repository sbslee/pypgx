import argparse
import logging

LINE_BREAK1 = "-" * 70
LINE_BREAK2 = "*" * 70

logging.basicConfig(level=logging.INFO)

def get_parser():
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
    return parser

def is_namespace(x):
    return isinstance(x, argparse.Namespace)