"""Command-line entry point for VarTrace."""

import argparse

from vartrace import __version__


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="vartrace",
        description="Estimate whether variants seen in single-cell RNA-seq "
                    "are present in the genome.",
    )
    parser.add_argument("--version", action="version",
                        version=f"vartrace {__version__}")
    parser.parse_args(argv)
    print(f"vartrace {__version__}: under development, "
          "no analysis is available yet.")
    return 0
