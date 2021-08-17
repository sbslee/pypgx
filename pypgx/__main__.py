import argparse

import pypgx
from .cli import commands
from .version import __version__

def _func(l):
    """Convert a list to a pretty text."""
    return '{' + ", ".join([f"'{x}'" for x in l]) + '}'

def _get_parser():
    parser = argparse.ArgumentParser(add_help=False)

    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="Show the version and exit."
    )

    parser.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit."
    )

    subparsers = parser.add_subparsers(
        dest="command",
        metavar="COMMAND",
        help="Name of the command."
    )

    subparsers.required = True

    compare_stargazer_calls_parser = subparsers.add_parser(
        "compare-stargazer-calls",
        add_help=False,
        help=("Compute the concordance between two 'genotype-calls.tsv' "
              "files created by Stargazer."),
        description=("Compute the concordance between two "
                     "'genotype-calls.tsv' files created by Stargazer.")
    )

    compare_stargazer_calls_parser._optionals.title = "Arguments"

    compare_stargazer_calls_parser.add_argument(
        "-r",
        "--ref-file",
        metavar="PATH",
        required=True,
        help=("Path to the reference or truth 'genotype-calls.tsv' file "
              "created by Stargazer. [required]")
    )

    compare_stargazer_calls_parser.add_argument(
        "-t",
        "--test-file",
        metavar="PATH",
        required=True,
        help=("Path to the test 'genotype-calls.tsv' file "
              "created by Stargazer. [required]")
    )

    compare_stargazer_calls_parser.add_argument(
        "-o",
        "--output-file",
        metavar="PATH",
        required=True,
        help=("Path to the output file. [required]")
    )

    compare_stargazer_calls_parser.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit."
    )
    
    return parser

def main():
    parser = _get_parser()
    args = parser.parse_args()
    command = args.command
    delattr(args, "command")
    commands[command](**vars(args))

if __name__ == "__main__":
    main()
