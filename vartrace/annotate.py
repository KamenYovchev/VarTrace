"""Fetching dbSNP and gnomAD from Ensembl when the input has no annotation.

Only these two are fetched. Of the five database features they carry nearly all
of the weight, they are free, and they come from one place. REDIportal and
COSMIC are used when the input already has them but are never fetched, because
COSMIC is licensed and REDIportal adds nothing measurable.

The public Ensembl service is shared and often busy, so temporary errors are
retried and a request that keeps failing is reported rather than hidden. The
caller can then carry on without annotation instead of stopping.

The models were trained on gnomAD as ANNOVAR supplies it. Ensembl may serve a
different release, so which one answered is written into the output.
"""

import time

import requests

SERVER = "https://rest.ensembl.org"
BATCH = 100
TIMEOUT = 60
RETRY_ON = {429, 500, 502, 503, 504}
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}


class EnsemblUnavailable(RuntimeError):
    pass


def region(v):
    chrom = v.chrom[3:] if v.chrom.startswith("chr") else v.chrom
    return f"{chrom} {v.pos} . {v.ref} {v.alt} . . ."


def post(session, regions, attempts=5):
    # af_gnomade alone is accepted; asking for the genome frequencies as well
    # makes the service answer with an error, so only the exomes are requested
    body = {"variants": regions, "af_gnomade": 1}
    wait = 2
    for attempt in range(attempts):
        try:
            r = session.post(f"{SERVER}/vep/human/region", headers=HEADERS,
                             json=body, timeout=TIMEOUT)
        except requests.RequestException:
            time.sleep(wait)
            wait *= 2
            continue
        if r.status_code in RETRY_ON:
            time.sleep(float(r.headers.get("Retry-After", wait)))
            wait *= 2
            continue
        if r.status_code != 200:
            raise EnsemblUnavailable(f"Ensembl answered {r.status_code}")
        return r.json()
    raise EnsemblUnavailable("Ensembl did not answer after several attempts")


def entry_keys(entry):
    """Every way of recognising which variant an answer belongs to.

    The service echoes the input line, but not in every release, so the
    position and alleles are used as well and the order of the answers is
    never relied on.
    """
    keys = []
    if entry.get("input"):
        keys.append(entry["input"])
    chrom, start = entry.get("seq_region_name"), entry.get("start")
    alleles = entry.get("allele_string", "")
    if chrom is not None and start is not None and "/" in alleles:
        ref, _, alt = alleles.partition("/")
        keys.append(f"{chrom}:{start}:{ref}:{alt}")
    return keys


def variant_keys(v):
    chrom = v.chrom[3:] if v.chrom.startswith("chr") else v.chrom
    return [region(v), f"{chrom}:{v.pos}:{v.ref}:{v.alt}"]


def gnomad_af(entry, alt):
    best = 0.0
    for col in entry.get("colocated_variants") or []:
        freqs = (col.get("frequencies") or {}).get(alt) or {}
        for key in ("gnomade", "gnomadg", "gnomad"):
            value = freqs.get(key)
            if isinstance(value, (int, float)):
                best = max(best, float(value))
    return best


def rs_id(entry):
    for col in entry.get("colocated_variants") or []:
        name = col.get("id", "")
        if name.startswith("rs"):
            return name
    return "."


def release(session):
    try:
        r = session.get(f"{SERVER}/info/software", headers=HEADERS, timeout=15)
        if r.status_code == 200:
            return str(r.json().get("release", "unknown"))
    except requests.RequestException:
        pass
    return "unknown"


def annotate(variants, batch=BATCH, progress=None):
    """Adds dbSNP and gnomAD to every variant and says what answered."""
    session = requests.Session()
    answers = {}
    for start in range(0, len(variants), batch):
        chunk = variants[start:start + batch]
        for entry in post(session, [region(v) for v in chunk]):
            for key in entry_keys(entry):
                answers[key] = entry
        if progress:
            progress(min(start + batch, len(variants)), len(variants))

    found = 0
    for v in variants:
        entry = None
        for key in variant_keys(v):
            if key in answers:
                entry = answers[key]
                break
        rs = rs_id(entry) if entry else "."
        af = gnomad_af(entry, v.alt) if entry else 0.0
        v.dbsnp = rs
        v.gnomad_af = [str(af)]
        v.databases |= {"dbsnp", "gnomad"}
        if rs != "." or af > 0:
            found += 1

    return {
        "source": "Ensembl VEP",
        "server": SERVER,
        "release": release(session),
        "queried": len(variants),
        "with_a_record": found,
    }
