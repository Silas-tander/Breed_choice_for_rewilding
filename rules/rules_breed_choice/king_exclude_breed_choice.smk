# Greedy KING exclusion (PO/FS/Dup-MZ/2nd pairs), alternative to the ngsRelate keep list:
# A) re-creates the previous hand-curated keep list, B) same on KING of all --mind 0.20 samples.

KING_EXCL_DIR = "results/breed_choice/relatedness_king"

KING_SEG = "data/king/king_ibs.seg"
KING_EXCL_TYPES = ["PO", "FS", "Dup/MZ", "2nd"]
KING_MANUAL_EXCLUDE = {
    "TG1111D2629": "Mongolian; extreme heterozygosity outlier, removed by hand in the previous keep list",
}


def _king_greedy_exclude(input, output, params, log):
    import pandas as pd

    seg = pd.read_csv(input.seg, sep="\t", dtype={"ID1": str, "ID2": str})
    pairs = seg[seg["InfType"].isin(params.types)]

    # igraph::graph_from_data_frame vertex order, for its which.max tie-breaking
    order = list(dict.fromkeys(list(pairs["ID1"]) + list(pairs["ID2"])))
    rank = {s: i for i, s in enumerate(order)}
    adj = {s: set() for s in order}
    for id1, id2 in zip(pairs["ID1"], pairs["ID2"]):
        adj[id1].add(id2)
        adj[id2].add(id1)

    drop = []
    while any(adj.values()):
        v = max(order, key=lambda n: (len(adj[n]), -rank[n]))
        drop.append(v)
        for nbr in adj[v]:
            adj[nbr].discard(v)
        adj[v] = set()

    names = [l.strip() for l in open(input.samples) if l.strip()]
    reasons = {}
    n_premind = 0
    if "king_fam" in input.keys():  # B: samples removed by --mind before KING
        in_king = {line.split()[1] for line in open(input.king_fam) if line.strip()}
        reasons = {s: "F_MISS > 0.2 (removed before KING)" for s in names if s not in in_king}
        n_premind = len(reasons)
    reasons.update({s: "related (KING " + "/".join(params.types) + ")" for s in drop if s in names})
    for s, why in params.manual.items():
        reasons.setdefault(s, f"manual: {why}")
    keep = sorted(set(names) - set(reasons))

    os.makedirs(os.path.dirname(output.keep), exist_ok=True)
    with open(output.drop, "w") as fh:
        fh.writelines(s + "\n" for s in drop)
    with open(output.keep, "w") as fh:
        fh.writelines(s + "\n" for s in keep)
    pd.DataFrame(sorted(reasons.items()), columns=["sample_alias", "reason"]).to_csv(output.excluded, index=False)
    pairs.assign(kinship=pairs["PropIBD"] / 2,
                 ID1_dropped=pairs["ID1"].isin(drop),
                 ID2_dropped=pairs["ID2"].isin(drop),
                 both_in_breed_choice=pairs["ID1"].isin(names) & pairs["ID2"].isin(names)) \
         [["ID1", "ID2", "InfType", "IBD1Seg", "IBD2Seg", "PropIBD", "kinship",
           "ID1_dropped", "ID2_dropped", "both_in_breed_choice"]] \
         .to_csv(output.pairs, index=False)

    hand_drop = [l.strip() for l in open(input.hand_drop) if l.strip()]
    hand_keep = sorted(l.strip() for l in open(input.hand_keep) if l.strip())
    in_seg = set(seg["ID1"]) | set(seg["ID2"])
    with open(log[0], "w") as lf:
        lf.write(
            f"{input.seg}: {len(pairs)} KING {'/'.join(params.types)} pairs involving {len(order)} samples; "
            f"{len(drop)} dropped by greedy removal ({len(reasons) - len(params.manual) - n_premind} of them in breed_choice).\n"
            f"{n_premind} samples removed for missingness before KING.\n"
            f"{len(names)} breed_choice samples, {len(set(names) & in_seg)} of them in the KING .seg; "
            f"{len(params.manual)} manual exclusions; {len(keep)} kept.\n"
            f"drop order == first {len(drop)} lines of {input.hand_drop}: {drop == hand_drop[:len(drop)]}\n"
            f"keep list == {input.hand_keep}: {keep == hand_keep}\n"
            f"  only in this keep list: {sorted(set(keep) - set(hand_keep))}\n"
            f"  only in the hand-curated keep list: {sorted(set(hand_keep) - set(keep))}\n"
        )


# A) from a previous KING run (data/king/king_ibs.seg; 67 of the 146 samples)
rule king_relatedness_exclude_breed_choice:
    input:
        seg = KING_SEG,
        samples = f"{REL_ANALYSIS_DIR}/relatedness/all_GTmiss.samples.txt",
        hand_drop = "metadata/previous_drop_order.txt",  # only for the comparison in the log
        hand_keep = "metadata/previous_keep_list.txt"
    output:
        drop = f"{KING_EXCL_DIR}/king_greedy_drop_order.txt",
        keep = f"{KING_EXCL_DIR}/unrelated_samples_king_breed_choice.txt",
        excluded = f"{KING_EXCL_DIR}/excluded_samples_king_breed_choice.csv",
        pairs = f"{KING_EXCL_DIR}/related_pairs_king.csv"
    params:
        types = KING_EXCL_TYPES,
        manual = KING_MANUAL_EXCLUDE
    log:
        "results/breed_choice/logs/relatedness/king_relatedness_exclude.log"
    resources:
        mem_mb = 4000,
        runtime = 30
    run:
        _king_greedy_exclude(input, output, params, log)


# B) from king_all_breed_choice.smk
rule king_mind_relatedness_exclude_breed_choice:
    input:
        seg = f"{KING_ALL_DIR}/king_ibs.seg",
        king_fam = f"{KING_ALL_DIR}/all_GTmiss.mind.fam",
        samples = f"{REL_ANALYSIS_DIR}/relatedness/all_GTmiss.samples.txt",
        hand_drop = "metadata/previous_drop_order.txt",
        hand_keep = "metadata/previous_keep_list.txt"
    output:
        drop = f"{KING_ALL_DIR}/king_greedy_drop_order.txt",
        keep = f"{KING_ALL_DIR}/unrelated_samples_king_mind_breed_choice.txt",
        excluded = f"{KING_ALL_DIR}/excluded_samples_king_mind_breed_choice.csv",
        pairs = f"{KING_ALL_DIR}/related_pairs_king_mind.csv"
    params:
        types = KING_EXCL_TYPES,
        manual = KING_MANUAL_EXCLUDE
    log:
        "results/breed_choice/logs/relatedness/king_mind_relatedness_exclude.log"
    resources:
        mem_mb = 4000,
        runtime = 30
    run:
        _king_greedy_exclude(input, output, params, log)
