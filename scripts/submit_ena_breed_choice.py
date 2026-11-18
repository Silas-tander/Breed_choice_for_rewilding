#!/usr/bin/env python3
"""Interactive ENA submission of the mols*/modoc* reads: registers the study and samples, writes the
Webin-CLI manifests and runs webin-cli -validate. Run on an interactive node after
`snakemake fastq_formatting`:
    ./scripts/submit_ena_breed_choice.py [--production]
Credentials: getpass prompt or ENA_USERNAME / ENA_PASSWORD. webin-cli only accepts the password as a
command-line flag, so it is visible in `ps` while each validate call runs.
"""
import argparse
import csv
import getpass
import os
import subprocess
import sys
import xml.etree.ElementTree as ET

import requests
from requests.auth import HTTPBasicAuth

FASTQ_META = "metadata/sample_sheet.tsv"
FASTQ_IND_PREFIXES = ("mols", "modoc")
BREED_CHOICE_META = "metadata/samples.csv"
WEBIN_CLI_JAR = "scripts/webin-cli-9.0.3.jar"

ENA_SUBMIT_TEST_URL = "https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit/"
ENA_SUBMIT_PROD_URL = "https://www.ebi.ac.uk/ena/submit/drop-box/submit/"

ENA_PROJECT_ALIAS = "Breed choice for rewilding"
ENA_PROJECT_TITLE = "Informing rewilding with genomics: choosing horse breeds for long-term sustainable populations"
ENA_PROJECT_DESCRIPTION = ("Whole-genome sequencing of Exmoor ponies from Mols laboratory Denmark and feral "
                            "mustangs from Modoc National Forest California.")

ENA_TAXON_ID = "9796"  # Equus caballus
ENA_SCIENTIFIC_NAME = "Equus caballus"
ENA_INSERT_SIZE = "280"  # TODO: confirm against the sequencing run

ACCDIR = "data/ENA/accessions"
MANIFESTDIR = "data/ENA/manifests"
OUTDIR = "data/ENA/ena_output"

SUBMISSION_ADD_XML = """<SUBMISSION>
  <ACTIONS>
    <ACTION>
      <ADD/>
    </ACTION>
  </ACTIONS>
</SUBMISSION>
"""


def load_line_fastq_pairs():
    """ind -> sorted [basename, ...] (e.g. modoc1 -> [149978x_1.fastq,
    149978x_2.fastq]) for every mols*/modoc* row in FASTQ_META -- a file
    that lives under this repo (Silas/Horses), read purely for its ind/
    raw_file columns. raw_file's own path may point outside this repo (into
    a sibling project's raw-data tree), but this script never opens or
    depends on anything at that path -- only os.path.basename() of the
    string matters, to line up with the already-staged data/ENA/*.gz
    files."""
    pairs = {}
    with open(FASTQ_META) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["ind"].startswith(FASTQ_IND_PREFIXES):
                pairs.setdefault(row["ind"], []).append(os.path.basename(row["raw_file"]))
    return {ind: sorted(basenames) for ind, basenames in pairs.items()}


def load_line_ena_sample_meta(line_fastq_pairs):
    """ind -> {instrument, library_strategy, library_source,
    library_selection, breed, gender}, read from metadata/samples.csv (the
    project's ENA-ready sample sheet, keyed by sample_alias)."""
    meta = {}
    with open(BREED_CHOICE_META, newline="") as fh:
        for row in csv.DictReader(fh):
            alias = row["sample_alias"]
            if alias in line_fastq_pairs:
                meta[alias] = {
                    "instrument": row["instrument_platform"],
                    "library_strategy": row["library_strategy"],
                    "library_source": row["library_source"],
                    "library_selection": row["library_selection"],
                    "breed": row["breed"],
                    "gender": row["sex"],
                }
    return meta


def ena_submit_xml(url, user, password, submission_path, payload_field, payload_path, label):
    """POST a SUBMISSION + one payload (PROJECT or SAMPLE) XML to an ENA
    drop-box URL and return the accession ENA assigned. Uses `requests` with
    HTTP Basic auth so the credentials travel in-process (request headers),
    never as a visible command-line argument the way a `curl -u user:pass`
    call would expose them via `ps`. ENA reports failure inside a 200
    response (success="false" + <ERROR> messages), not via HTTP status, so
    that's checked explicitly rather than trusting the HTTP status code."""
    with open(submission_path, "rb") as submission_fh, open(payload_path, "rb") as payload_fh:
        response = requests.post(
            url,
            auth=HTTPBasicAuth(user, password),
            files={
                "SUBMISSION": (os.path.basename(submission_path), submission_fh, "text/xml"),
                payload_field: (os.path.basename(payload_path), payload_fh, "text/xml"),
            },
        )
    root = ET.fromstring(response.text)
    if root.get("success") != "true":
        errors = [e.text for e in root.iter("ERROR")] or [response.text]
        raise RuntimeError(f"ENA {label} registration failed:\n" + "\n".join(errors))
    accession = root.find(payload_field).get("accession")
    if not accession:
        raise RuntimeError(f"ENA {label} response had no accession:\n{response.text}")
    return accession, response.text


def register_study(url, user, password):
    os.makedirs(ACCDIR, exist_ok=True)
    project_xml = os.path.join(ACCDIR, "project.xml")
    submission_xml = os.path.join(ACCDIR, "project_submission.xml")
    with open(project_xml, "w") as fh:
        fh.write(f"""<PROJECT_SET>
  <PROJECT alias="{ENA_PROJECT_ALIAS}">
    <TITLE>{ENA_PROJECT_TITLE}</TITLE>
    <DESCRIPTION>{ENA_PROJECT_DESCRIPTION}</DESCRIPTION>
    <SUBMISSION_PROJECT>
      <SEQUENCING_PROJECT/>
    </SUBMISSION_PROJECT>
  </PROJECT>
</PROJECT_SET>
""")
    with open(submission_xml, "w") as fh:
        fh.write(SUBMISSION_ADD_XML)
    accession, receipt = ena_submit_xml(url, user, password, submission_xml, "PROJECT", project_xml, "study")
    with open(os.path.join(ACCDIR, "project_receipt.xml"), "w") as fh:
        fh.write(receipt)
    with open(os.path.join(ACCDIR, "study_accession.txt"), "w") as fh:
        fh.write(accession + "\n")
    return accession


# ENA requires collection date and country; unknown here, so 'NA' (use an INSDC missing-value
# term if ENA rejects it)
def register_sample(url, user, password, ind, meta):
    sample_xml = os.path.join(ACCDIR, f"{ind}_sample.xml")
    submission_xml = os.path.join(ACCDIR, f"{ind}_sample_submission.xml")
    with open(sample_xml, "w") as fh:
        fh.write(f"""<SAMPLE_SET>
  <SAMPLE alias="{ind}">
    <TITLE>{ind}</TITLE>
    <SAMPLE_NAME>
      <TAXON_ID>{ENA_TAXON_ID}</TAXON_ID>
      <SCIENTIFIC_NAME>{ENA_SCIENTIFIC_NAME}</SCIENTIFIC_NAME>
    </SAMPLE_NAME>
    <SAMPLE_ATTRIBUTES>
      <SAMPLE_ATTRIBUTE>
        <TAG>sex</TAG>
        <VALUE>{meta["gender"]}</VALUE>
      </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE>
        <TAG>breed</TAG>
        <VALUE>{meta["breed"]}</VALUE>
      </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE>
        <TAG>collection date</TAG>
        <VALUE>NA</VALUE>
      </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE>
        <TAG>geographic location (country and/or sea)</TAG>
        <VALUE>NA</VALUE>
      </SAMPLE_ATTRIBUTE>
    </SAMPLE_ATTRIBUTES>
  </SAMPLE>
</SAMPLE_SET>
""")
    with open(submission_xml, "w") as fh:
        fh.write(SUBMISSION_ADD_XML)
    accession, receipt = ena_submit_xml(url, user, password, submission_xml, "SAMPLE", sample_xml,
                                         f"sample ({ind})")
    with open(os.path.join(ACCDIR, f"{ind}_sample_receipt.xml"), "w") as fh:
        fh.write(receipt)
    with open(os.path.join(ACCDIR, f"{ind}_sample_accession.txt"), "w") as fh:
        fh.write(accession + "\n")
    return accession


def write_manifest(ind, meta, study_accession, sample_accession, basenames):
    os.makedirs(MANIFESTDIR, exist_ok=True)
    manifest_path = os.path.join(MANIFESTDIR, f"{ind}_manifest.txt")
    with open(manifest_path, "w") as fh:
        fh.write(f"STUDY\t{study_accession}\n")
        fh.write(f"SAMPLE\t{sample_accession}\n")
        fh.write(f"NAME\t{ind}\n")
        fh.write(f"INSTRUMENT\t{meta['instrument']}\n")
        fh.write(f"LIBRARY_STRATEGY\t{meta['library_strategy']}\n")
        fh.write(f"LIBRARY_SOURCE\t{meta['library_source']}\n")
        fh.write(f"LIBRARY_SELECTION\t{meta['library_selection']}\n")
        fh.write(f"INSERT_SIZE\t{ENA_INSERT_SIZE}\n")
        for basename in basenames:
            fh.write(f"FASTQ\t{os.path.abspath(os.path.join('data/ENA', basename + '.gz'))}\n")
    return manifest_path


def webin_validate(ind, manifest_path, user, password):
    """Run `webin-cli -validate` for one individual. NOTE: ENA's Webin-CLI
    only accepts credentials via the -userName/-password command-line
    flags -- there is no supported environment-variable, stdin, or .netrc
    alternative (the webin-cli-username/webin-cli-password env vars in its
    own repo are for testing webin-cli's build, not for submitters). That
    means -password is visible in `ps`/`/proc/<pid>/cmdline` to other users
    on the same node for as long as this call runs. This is a limitation
    of webin-cli itself and cannot be avoided from this script; the window
    is limited to the duration of each individual validate call."""
    outdir = os.path.join(OUTDIR, ind)
    os.makedirs(outdir, exist_ok=True)
    subprocess.run(
        ["java", "-jar", WEBIN_CLI_JAR,
         "-context", "reads",
         "-userName", user,
         "-password", password,
         "-manifest", manifest_path,
         "-outputDir", outdir,
         "-validate"],
        check=True
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--production", action="store_true",
                         help="submit to ENA's real drop-box instead of the wwwdev sandbox")
    args = parser.parse_args()
    url = ENA_SUBMIT_PROD_URL if args.production else ENA_SUBMIT_TEST_URL

    line_fastq_pairs = load_line_fastq_pairs()
    line_ena_sample_meta = load_line_ena_sample_meta(line_fastq_pairs)
    missing_meta = sorted(set(line_fastq_pairs) - set(line_ena_sample_meta))
    if missing_meta:
        sys.exit(f"No metadata/samples.csv row for: {', '.join(missing_meta)}")

    missing_gz = [os.path.join("data/ENA", b + ".gz")
                  for basenames in line_fastq_pairs.values() for b in basenames
                  if not os.path.exists(os.path.join("data/ENA", b + ".gz"))]
    if missing_gz:
        sys.exit("Missing staged fastq(s), run the Snakemake staging rule first "
                  "(snakemake fastq_formatting):\n  " + "\n  ".join(missing_gz))
    if not os.path.exists(WEBIN_CLI_JAR):
        sys.exit(f"{WEBIN_CLI_JAR} not found")

    user = os.environ.get("ENA_USERNAME") or input("ENA Webin username: ").strip()
    password = os.environ.get("ENA_PASSWORD") or getpass.getpass("ENA Webin password: ")
    if not user or not password:
        sys.exit("ENA username/password cannot be empty")

    print(f"Registering study against {url} ...", file=sys.stderr)
    study_accession = register_study(url, user, password)
    print(f"  study accession: {study_accession}", file=sys.stderr)

    results = {}
    for ind in sorted(line_fastq_pairs):
        meta = line_ena_sample_meta[ind]
        try:
            print(f"[{ind}] registering sample ...", file=sys.stderr)
            sample_accession = register_sample(url, user, password, ind, meta)
            print(f"[{ind}]   sample accession: {sample_accession}", file=sys.stderr)
            manifest_path = write_manifest(ind, meta, study_accession, sample_accession, line_fastq_pairs[ind])
            print(f"[{ind}] running webin-cli -validate ...", file=sys.stderr)
            webin_validate(ind, manifest_path, user, password)
            results[ind] = "OK"
        except (RuntimeError, subprocess.CalledProcessError) as exc:
            print(f"[{ind}] FAILED: {exc}", file=sys.stderr)
            results[ind] = "FAILED"

    print("\nSummary:", file=sys.stderr)
    for ind, status in results.items():
        print(f"  {ind}: {status}", file=sys.stderr)
    if any(status == "FAILED" for status in results.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
