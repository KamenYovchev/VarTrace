"""Command-line entry point for VarTrace."""

import argparse
import os
import socket
import subprocess
import sys
import threading
import time
import webbrowser

from vartrace import __version__, annotate, features, io, model as models, report


def score(variants, n_cells, models_dir, use_ensembl, log=print):
    """Everything between reading the files and having a score per variant."""
    features.compute(variants, n_cells)

    blind = features.without_allele_counts(variants)
    if blind > 0.2 * len(variants):
        log(f"warning: {blind} of {len(variants)} variants have no allele "
            "fraction, because the input carries neither AD nor VAF. The "
            "scores will not mean much.")

    available = features.available(variants)

    annotation = None
    if use_ensembl and not (available & features.DB_GROUPS):
        log(f"no database annotation in the input, asking Ensembl about "
            f"{len(variants)} variants")
        try:
            annotation = annotate.annotate(variants)
            for v in variants:
                v.features.update(features.database_features(v))
            available = features.available(variants)
            log(f"Ensembl had a record for {annotation['with_a_record']} of "
                f"{annotation['queried']}")
        except annotate.EnsemblUnavailable as problem:
            log(f"Ensembl did not answer ({problem}), carrying on without it")

    model, why = models.choose(models.load_all(models_dir), available)
    return model, why, model.predict(variants), annotation


def open_when_ready(port, url, timeout=30):
    deadline = time.time() + timeout
    while time.time() < deadline:
        with socket.socket() as probe:
            probe.settimeout(0.5)
            if probe.connect_ex(("127.0.0.1", port)) == 0:
                webbrowser.open(url)
                return
        time.sleep(0.3)


def serve(argv):
    parser = argparse.ArgumentParser(prog="vartrace serve",
                                     description="Open the web interface.")
    parser.add_argument("-p", "--port", type=int, default=8501)
    parser.add_argument("--expose", action="store_true",
                        help="listen on every network interface instead of "
                             "this machine only, for serving others")
    parser.add_argument("--no-browser", action="store_true")
    args = parser.parse_args(argv)
    try:
        import streamlit  # noqa: F401
    except ImportError:
        sys.exit('the web interface needs extra packages, install with '
                 'pipx install "vartrace[web]"')

    # people upload genomic data here, so by default nothing outside this
    # machine can reach it
    address = "0.0.0.0" if args.expose else "127.0.0.1"
    url = f"http://localhost:{args.port}"
    if not args.no_browser:
        threading.Thread(target=open_when_ready, args=(args.port, url),
                         daemon=True).start()

    print(f"open {url} in a browser, stop with Ctrl-C")
    if args.expose:
        print("warning: reachable from the network, anyone who can reach this "
              "port can upload files")
    app = os.path.join(os.path.dirname(__file__), "web", "app.py")
    # headless stops Streamlit asking for an email address on first run, so the
    # browser is opened above instead; statistics are off because of the data
    return subprocess.call([sys.executable, "-m", "streamlit", "run", app,
                            "--server.port", str(args.port),
                            "--server.address", address,
                            "--server.headless=true",
                            "--browser.gatherUsageStats=false"])


def run(argv):
    parser = argparse.ArgumentParser(
        prog="vartrace",
        description="Estimate whether variants seen in single-cell RNA-seq "
                    "are present in the genome.",
    )
    parser.add_argument("files", nargs="+",
                        help="one VCF or ANNOVAR table per cell")
    parser.add_argument("-o", "--output", default="vartrace_results.tsv")
    parser.add_argument("-t", "--threshold", type=float, default=0.5,
                        help="0.5 by default, 0.867 for the high confidence "
                             "mode where precision is about 0.99")
    parser.add_argument("--no-annotate", action="store_true",
                        help="do not ask Ensembl for dbSNP and gnomAD when the "
                             "input has no annotation of its own")
    parser.add_argument("--models", default=models.MODELS_DIR)
    parser.add_argument("--version", action="version",
                        version=f"vartrace {__version__}")
    args = parser.parse_args(argv)

    variants, n_cells, skipped = io.read_all(args.files)
    if not variants:
        sys.exit("no usable variants in the input")

    model, why, scores, annotation = score(
        variants, n_cells, args.models, not args.no_annotate)

    report.write(args.output, report.rows(variants, scores, model, args.threshold))
    info = report.provenance(model, why, args.threshold, annotation)
    with open(args.output + ".info.txt", "w") as fh:
        fh.write(info + "\n")

    print(report.summary(scores, args.threshold, model, why, n_cells, skipped))
    print(f"written to {args.output} and {args.output}.info.txt")
    return 0


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] == "serve":
        return serve(argv[1:])
    return run(argv)
