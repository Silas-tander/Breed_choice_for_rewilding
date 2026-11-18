# Stage the in-house mols*/modoc* fastqs (outside this project) gzipped into data/ENA/ for ENA upload.
# Not in rule all: snakemake fastq_formatting

import csv
import os
import re


FASTQ_META = "metadata/sample_sheet.tsv"
FASTQ_SRC_PREFIX = "data/raw_fastq/"
FASTQ_IND_PREFIXES = ("mols", "modoc")  # UoN* are already on ENA

def _load_line_fastqs():
    """basename (e.g. 149978x_1.fastq) -> source path relative to the repo
    root, for every raw_file row under the Line/ tree whose sample (ind) is
    an Exmoor (mols*) or Modoc (modoc*) individual."""
    mapping = {}
    with open(FASTQ_META) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            raw = row["raw_file"]
            if raw.startswith(FASTQ_SRC_PREFIX) and row["ind"].startswith(FASTQ_IND_PREFIXES):
                basename = os.path.basename(raw)
                mapping[basename] = os.path.normpath(raw)  # raw_file is relative to the repository root
    return mapping

LINE_FASTQS = _load_line_fastqs()
LINE_FASTQ_ORDER = sorted(LINE_FASTQS.keys())

# First rule, so it is the default target when this file is run on its own
rule fastq_formatting:
    input:
        expand("data/ENA/{basename}.gz", basename=LINE_FASTQ_ORDER),
        "data/ENA/checksums.md5"

def _prev_ena_output(wildcards):
    """The previous file's .gz output in LINE_FASTQ_ORDER, as an extra input
    -- not read by the shell command, just an edge in the DAG. Forces each
    stage_line_fastq_to_ena_breed_choice job to wait for the one before it,
    so only one of these 64 GB gzips ever runs at a time regardless of
    --cores. The first file has no predecessor, so it gets none."""
    idx = LINE_FASTQ_ORDER.index(wildcards.basename)
    if idx == 0:
        return []
    return [f"data/ENA/{LINE_FASTQ_ORDER[idx - 1]}.gz"]

rule stage_line_fastq_to_ena_breed_choice:
    input:
        src = lambda wc: LINE_FASTQS[wc.basename],
        prev = _prev_ena_output
    output:
        "data/ENA/{basename}.gz"
    wildcard_constraints:
        basename = "|".join(re.escape(b) for b in LINE_FASTQS)
    threads: 1
    resources:
        mem_mb = 2000,
        runtime = 1440
    shell:
        r"""
        gzip -c {input.src} > {output}
        """

# Webin-CLI wants md5s keyed by bare file name, hence the cd
rule checksums_ena_breed_choice:
    input:
        expand("data/ENA/{basename}.gz", basename=LINE_FASTQ_ORDER)
    params:
        basenames = " ".join(f"{b}.gz" for b in LINE_FASTQ_ORDER)
    output:
        "data/ENA/checksums.md5"
    threads: 1
    resources:
        mem_mb = 1000,
        runtime = 120
    shell:
        r"""
        cd data/ENA && md5sum {params.basenames} > checksums.md5
        """
