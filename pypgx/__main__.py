import argparse

from .version import __version__
from .cli import commands

def main():
    parser = argparse.ArgumentParser(
        add_help=False,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        '-h',
        '--help',
        action='help',
        default=argparse.SUPPRESS,
        help='Show this help message and exit.',
    )
    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version=f'%(prog)s {__version__}',
        help='Show the version number and exit.'
    )
    subparsers = parser.add_subparsers(
        dest='command',
        metavar='COMMAND',
        required=True,
    )
    for name, command in commands.items():
        command.create_parser(subparsers)
    args = parser.parse_args()
    commands[args.command].main(args)

if __name__ == '__main__':
    main()
