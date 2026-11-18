"""Per-autosome length, N bp, mask bp and their union (bp that can never carry a called SNP).
Usage: python autosome_uncallable_bp.py <ref.fa> <ref.fa.fai> <mask.bed> <x_contig> <out.tsv>
"""
import sys

import numpy as np

fa, fai, bed, x_contig, out = sys.argv[1:6]

lengths = {}
for line in open(fai):
    name, length = line.split("\t")[:2]
    if name != x_contig:
        lengths[name] = int(length)

mask = {c: [] for c in lengths}
for line in open(bed):
    c, s, e = line.split("\t")[:3]
    if c in mask:
        mask[c].append((int(s), int(e)))


def report(name, seq_chunks, fh):
    seq = np.frombuffer(b"".join(seq_chunks), dtype=np.uint8)
    assert len(seq) == lengths[name], f"{name}: fasta length {len(seq)} != fai {lengths[name]}"
    is_n = (seq == ord("N")) | (seq == ord("n"))
    masked = np.zeros(len(seq), dtype=bool)
    for s, e in mask[name]:
        masked[s:e] = True
    fh.write(f"{name}\t{len(seq)}\t{is_n.sum()}\t{masked.sum()}\t{(is_n | masked).sum()}\n")


with open(out, "w") as fh:
    fh.write("contig\tlength\tN_bp\tmask_bp\tN_or_mask_bp\n")
    name, chunks = None, []
    for line in open(fa, "rb"):
        if line.startswith(b">"):
            if name in lengths:
                report(name, chunks, fh)
            name, chunks = line[1:].split()[0].decode(), []
        elif name in lengths:
            chunks.append(line.rstrip())
    if name in lengths:
        report(name, chunks, fh)
