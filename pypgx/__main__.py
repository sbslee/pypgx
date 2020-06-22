import argparse
import logging

from pgkb import pgkb

def main():
    logging.basicConfig(level=logging.INFO)
    logging.info("Who called me?")

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(
        dest='tool',
        metavar='tool',
        help='name of tool',
    )
    pgkb_parser = subparsers.add_parser("pgkb")
    args = parser.parse_args()

    if args.tool == "pgkb":
        pgkb(args)

if __name__ == "__main__":
    main()
