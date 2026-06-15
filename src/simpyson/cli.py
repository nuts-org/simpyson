from __future__ import annotations

import argparse

from simpyson.gui import main as gui_main


def main() -> None:
    """
    Entry point for the ``simpyson`` command-line interface.

    Subcommands
    -----------
    gui
        Launch the simpyson GUI. Optionally pass file paths to open on start.
    """
    parser = argparse.ArgumentParser(prog="simpyson")
    subparsers = parser.add_subparsers(dest="command")

    parser_gui = subparsers.add_parser("gui", help="Launch the GUI")
    parser_gui.add_argument('files', nargs='*', help='Optional file(s) to open.')
    parser_gui.set_defaults(func=lambda args: gui_main(args.files))

    args = parser.parse_args()
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()
