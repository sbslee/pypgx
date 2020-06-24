import logging

from .common import get_parser

from .pgkb import pgkb
from .report import report

def main():
    parser = get_parser()
    args = parser.parse_args()

    if args.tool == "pgkb":
        pgkb(args)

    elif args.tool == "report":
        report(args)

if __name__ == "__main__":
    main()
