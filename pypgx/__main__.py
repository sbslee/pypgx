import argparse
from pypgx import *
from .version import __version__

def get_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="Show the version and exit."
    )

    subparsers = parser.add_subparsers(
        dest="command",
        metavar="COMMAND",
        help="Name of the command."
    )

    subparsers.required = True

    compare_genotype_calls_parser = subparsers.add_parser(
        "compare-stargazer-calls",
        help=("Compute the concordance between two 'genotype-calls.tsv' "
              "files created by Stargazer."),
        description=("Compute the concordance between two 'genotype-calls.tsv' "
              "files created by Stargazer.")
    )

    compare_genotype_calls_parser.add_argument(
        "-r",
        "--ref-file",
        metavar="PATH",
        required=True,
        help=("Path to the reference or truth 'genotype-calls.tsv' file "
              "created by Stargazer. [required]")
    )

    compare_genotype_calls_parser.add_argument(
        "-t",
        "--test-file",
        metavar="PATH",
        required=True,
        help=("Path to the test 'genotype-calls.tsv' file "
              "created by Stargazer. [required]")
    )

    compare_genotype_calls_parser.add_argument(
        "-o",
        "--output-file",
        metavar="PATH",
        required=True,
        help=("Path to the output file. [required]")
    )

    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()
    command = args.command
    delattr(args, "command")
    commands[command](**vars(args))

if __name__ == "__main__":
    main()
