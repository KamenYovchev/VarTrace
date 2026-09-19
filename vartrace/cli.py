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


MOST = 0.2  # a fifth of the input is enough to make the run a different thing


class Refused(Exception):
    """The input is outside what the models were trained on."""

    def __init__(self, problems):
        self.problems = problems
        super().__init__("\n\n".join(p["what"] for p in problems))


def check(paths, variants):
    """The conditions that can be read off the files themselves.

    Whether the reads came from single cells or from bulk, which chemistry was
    used and what tissue it was leave no trace in a VCF and cannot be checked
    here. They belong in the instructions. The four below do leave a trace, and
    each of them changes what the scores mean rather than merely weakening
    them, which is why the answer is to stop rather than to warn.
    """
    n = len(variants)
    problems = []

    if "GRCh37" in {io.genome_build(path) for path in paths}:
        problems.append({
            "what": "The files say GRCh37. The models and the Ensembl lookup "
                    "both work on GRCh38, so every position here means "
                    "something else.",
            "do": "Lift the variants over to GRCh38, with CrossMap or the "
                  "Ensembl Assembly Converter, and upload them again.",
        })

    plain = features.without_read_level(variants)
    if plain > MOST * n:
        problems.append({
            "what": f"{plain} of {n} variants carry none of the read-level "
                    "annotations, so 7 of the 35 features would be empty. FS, "
                    "SOR, MQ, MQRankSum, ReadPosRankSum, BaseQRankSum and QD "
                    "are written by GATK, and unlike the databases no model "
                    "was trained without them.",
            "do": "Call the variants with GATK HaplotypeCaller, or with a "
                  "caller that writes the same INFO fields. If you already "
                  "have BAM files, only the calling step has to be redone.",
        })

    blind = features.without_allele_counts(variants)
    if blind > MOST * n:
        problems.append({
            "what": f"{blind} of {n} variants carry neither AD nor VAF, so "
                    "they have no allele fraction. It is the feature the "
                    "models lean on most after the databases.",
            "do": "Keep the AD field when you filter or convert the VCF. Some "
                  "tools drop it; bcftools annotate -x never touches AD "
                  "unless it is asked to.",
        })

    shallow = sum(1 for v in variants
                  if v.features["dp"] is not None and v.features["dp"] < 5)
    if shallow > MOST * n:
        problems.append({
            "what": f"{shallow} of {n} variants are covered by fewer than 5 "
                    "reads. Everything below 5 was removed from the training "
                    "data, so these sit outside what the models have seen.",
            "do": "Filter them out first: bcftools view -e 'FMT/DP<5'.",
        })

    return problems


def score(variants, n_cells, models_dir, use_ensembl, log=print, paths=(),
          anyway=False):
    """Everything between reading the files and having a score per variant."""
    features.compute(variants, n_cells)

    # before the Ensembl call, so a refused run costs nobody anything
    problems = check(paths, variants)
    if problems and not anyway:
        raise Refused(problems)
    for problem in problems:
        log("warning: " + problem["what"])

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


# One accent colour and nothing else. Set here rather than in a config file,
# because Streamlit looks for one in the working directory, which is wherever
# the user happens to be standing when they type vartrace.
THEME = {
    "base": "light",
    "primaryColor": "#2F6B45",
    "backgroundColor": "#FFFFFF",
    "secondaryBackgroundColor": "#F4F6F5",
    "textColor": "#1F2933",
    "linkColor": "#2F6B45",
    "borderColor": "#D4DAD6",
    "dataframeBorderColor": "#D4DAD6",
    "baseRadius": "medium",
    "buttonRadius": "medium",
}


def serve(argv):
    parser = argparse.ArgumentParser(prog="vartrace serve",
                                     description="Open the web interface.")
    parser.add_argument("-p", "--port", type=int, default=8501)
    parser.add_argument("--expose", action="store_true",
                        help="listen on every network interface instead of "
                             "this machine only, for serving others")
    parser.add_argument("--no-browser", action="store_true")
    args = parser.parse_args(argv)

    # people upload genomic data here, so by default nothing outside this
    # machine can reach it
    address = "0.0.0.0" if args.expose else "127.0.0.1"
    url = f"http://localhost:{args.port}"
    if not args.no_browser:
        threading.Thread(target=open_when_ready, args=(args.port, url),
                         daemon=True).start()

    print(f"open {url} in a browser, stop with Ctrl-C")
    print("to work from the command line instead, run: vartrace --help")
    if args.expose:
        print("warning: reachable from the network, anyone who can reach this "
              "port can upload files")
    app = os.path.join(os.path.dirname(__file__), "web", "app.py")
    # headless stops Streamlit asking for an email address on first run, so the
    # browser is opened above instead; statistics are off because of the data.
    # The minimal toolbar removes Streamlit's own Deploy button, which offers
    # to publish the app to their cloud and has no business here. A VCF holding
    # a few hundred cells passes the default upload limit easily, hence 4 GB.
    command = [sys.executable, "-m", "streamlit", "run", app,
               "--server.port", str(args.port),
               "--server.address", address,
               "--server.headless=true",
               "--server.maxUploadSize=4096",
               "--client.toolbarMode=minimal",
               "--browser.gatherUsageStats=false"]
    command += [f"--theme.{name}={value}" for name, value in THEME.items()]
    return subprocess.call(command)


def run(argv):
    parser = argparse.ArgumentParser(
        prog="vartrace",
        description="Estimate whether variants seen in single-cell RNA-seq "
                    "are present in the genome.",
        epilog="run vartrace with no file at all and the web interface opens\n"
               "in the browser. Use serve instead of a file to give it a port\n"
               "of its own:\n"
               "\n"
               "    vartrace serve --port 8600\n",
        formatter_class=argparse.RawDescriptionHelpFormatter,
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
    parser.add_argument("--anyway", action="store_true",
                        help="score even when the input is outside what the "
                             "models were trained on, having read why")
    parser.add_argument("--models", default=models.MODELS_DIR)
    parser.add_argument("--version", action="version",
                        version=f"vartrace {__version__}")
    args = parser.parse_args(argv)

    try:
        variants, n_cells, skipped = io.read_all(args.files)
    except io.UnknownFormat as problem:
        sys.exit(str(problem))
    if not variants:
        sys.exit("no usable variants in the input")

    try:
        model, why, scores, annotation = score(
            variants, n_cells, args.models, not args.no_annotate,
            paths=args.files, anyway=args.anyway)
    except Refused as refusal:
        lines = ["These files are not what the models were trained on."]
        for problem in refusal.problems:
            lines += ["", problem["what"], "  -> " + problem["do"]]
        lines += ["", "If you understand this and want the numbers anyway, "
                      "run again with --anyway."]
        sys.exit("\n".join(lines))

    report.write(args.output, report.rows(variants, scores, model, args.threshold))
    info = report.provenance(model, why, args.threshold, annotation)
    with open(args.output + ".info.txt", "w") as fh:
        fh.write(info + "\n")

    print(report.summary(scores, args.threshold, model, why, n_cells, skipped))
    print(f"written to {args.output} and {args.output}.info.txt")
    return 0


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    # plain "vartrace" opens the interface, because that is what most people
    # who install this will want and the command line is the harder path
    if not argv:
        return serve([])
    if argv[0] == "serve":
        return serve(argv[1:])
    return run(argv)
