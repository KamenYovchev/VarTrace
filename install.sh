#!/bin/sh
# Installs VarTrace on macOS and Linux.
#   curl -LsSf https://raw.githubusercontent.com/KamenYovchev/VarTrace/main/install.sh | sh
set -eu

say() { printf '%s\n' "$*"; }
die() { printf '%s\n' "$*" >&2; exit 1; }

have() { command -v "$1" >/dev/null 2>&1; }

python_ok() {
    have "$1" || return 1
    "$1" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 12) else 1)' 2>/dev/null
}

SHELL_PATH="$PATH"

say "Installing VarTrace."
say ""

INSTALLER=""
for py in python3.14 python3.13 python3.12 python3; do
    if python_ok "$py"; then
        if have pipx; then
            INSTALLER="pipx"
        fi
        break
    fi
done

if [ -z "$INSTALLER" ] && have uv; then
    INSTALLER="uv"
fi

if [ -z "$INSTALLER" ]; then
    say "Getting uv, which brings a Python of its own."
    curl -LsSf https://astral.sh/uv/install.sh | sh >/dev/null 2>&1 \
        || die "Could not install uv. See https://docs.astral.sh/uv/"
    PATH="$HOME/.local/bin:$PATH"
    export PATH
    have uv || die "uv installed but is not on the path. Open a new terminal and run this again."
    INSTALLER="uv"
fi

case "$INSTALLER" in
    pipx) pipx install --force vartrace >/dev/null 2>&1 || die "pipx could not install vartrace." ;;
    uv)   uv tool install --force vartrace >/dev/null 2>&1 || die "uv could not install vartrace." ;;
esac

BIN="$HOME/.local/bin"
# against the path the shell had, not the one this script widened
if ! PATH="$SHELL_PATH" command -v vartrace >/dev/null 2>&1; then
    [ -x "$BIN/vartrace" ] || die "vartrace was installed but cannot be found. Open a new terminal and try again."
    say "One thing first. vartrace is in $BIN, which your shell does not look in."
    say "Add this line to your ~/.zshrc or ~/.bashrc, then open a new terminal:"
    say ""
    say "    export PATH=\"\$HOME/.local/bin:\$PATH\""
    say ""
fi

VERSION=$("$BIN/vartrace" --version 2>/dev/null || vartrace --version 2>/dev/null || echo "vartrace")

say "$VERSION is installed."
say ""
say "  vartrace              open the interface in your browser"
say "  vartrace --help       what it can do from the command line"
say ""
say "Your files need to be human variants on GRCh38, from single-cell RNA"
say "sequencing, called with GATK, carrying AD or VAF, covered by 5 reads or"
say "more. VarTrace checks all of that and says so before it scores anything."
say ""
say "Upload one file for each cell, or one VCF with a column per cell barcode."
say "Nothing leaves your machine except, if you allow it, the positions sent to"
say "Ensembl for dbSNP and gnomAD."
