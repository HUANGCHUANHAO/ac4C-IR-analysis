# Genomic-region distribution of differential peaks, peak-sequence extraction
# and 5-mer enrichment, and CCW (C-C-A/T) motif density against a
# mononucleotide-shuffled background.
# Outputs are written to figure2/results and rendered by Fig2_plot.py.

import bisect
import gzip
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.stats import fisher_exact, ranksums
import warnings
warnings.filterwarnings("ignore")

PROJ = Path(__file__).resolve().parents[2]
FIG1 = PROJ / "figure1"
FIG2 = PROJ / "figure2"
OUT = FIG2 / "results"
OUT.mkdir(parents=True, exist_ok=True)

DIFF = FIG1 / "results" / "all_diff_ac4c_sites.csv"
ALLX = FIG1 / "raw_data" / "Supplementary_Table_S1_differential_ac4C_sites.xlsx"
GTF = FIG2 / "raw_data/annotation/Mus_musculus_MM10_forRNAseq3875_20170608.gtf.gz"
FAI = FIG2 / "raw_data/MM10/BOWTIE2_MM10_Base.fa.fai"
GENOME_FA = FIG2 / "raw_data/MM10/BOWTIE2_MM10_Base.fa"
SEQ_CACHE = OUT / "peak_sequences.csv"


def merge(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def subtract_intervals(base, remove):
    if not base:
        return []
    if not remove:
        return [list(x) for x in base]
    res = []
    i = j = 0
    cur_s, cur_e = base[0]
    rem_s, rem_e = remove[0]
    while True:
        if cur_e <= rem_s:
            res.append([cur_s, cur_e])
            i += 1
            if i >= len(base):
                break
            cur_s, cur_e = base[i]
        elif cur_s >= rem_e:
            j += 1
            if j >= len(remove):
                res.append([cur_s, cur_e])
                i += 1
                while i < len(base):
                    res.append(list(base[i]))
                    i += 1
                break
            rem_s, rem_e = remove[j]
        else:
            if cur_s < rem_s:
                res.append([cur_s, rem_s])
            if cur_e > rem_e:
                cur_s = rem_e
                j += 1
                if j >= len(remove):
                    res.append([cur_s, cur_e])
                    i += 1
                    while i < len(base):
                        res.append(list(base[i]))
                        i += 1
                    break
                rem_s, rem_e = remove[j]
            else:
                i += 1
                if i >= len(base):
                    break
                cur_s, cur_e = base[i]
    return res


def step1_genomic_distribution():
    txs = {}
    with gzip.open(GTF, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            feat, chrom, s, e, strand, attrs = p[2], p[0], int(p[3]), int(p[4]), p[6], p[8]
            m = re.search(r'transcript_id "([^"]+)"', attrs)
            if not m:
                continue
            tid = m.group(1)
            t = txs.setdefault(tid, {
                "chrom": chrom, "strand": strand,
                "exon": [], "CDS": [], "UTR": [], "span": [10**18, 0]})
            t["span"][0] = min(t["span"][0], s)
            t["span"][1] = max(t["span"][1], e)
            if feat in ("exon", "CDS", "UTR"):
                t[feat].append((s, e))

    CLASSES = ["CDS", "5'UTR", "3'UTR", "lncRNA_exon", "intron"]
    CATS = ["CDS", "5'UTR", "3'UTR", "lncRNA_exon", "intron", "intergenic"]
    regions = defaultdict(lambda: defaultdict(list))
    for tid, t in txs.items():
        chrom, strand = t["chrom"], t["strand"]
        exons = merge(t["exon"])
        if t["CDS"]:
            cds_s = min(s for s, _ in t["CDS"])
            cds_e = max(e for _, e in t["CDS"])
            for s, e in t["UTR"]:
                mid = (s + e) / 2
                cls = None
                if mid <= cds_s:
                    cls = "5'UTR" if strand == "+" else "3'UTR"
                elif mid >= cds_e:
                    cls = "3'UTR" if strand == "+" else "5'UTR"
                if cls:
                    regions[chrom][cls].append((s, e))
            regions[chrom]["CDS"].extend(t["CDS"])
        else:
            regions[chrom]["lncRNA_exon"].extend(t["exon"])
        prev_e = None
        for s, e in exons:
            if prev_e is not None and s - prev_e > 1:
                regions[chrom]["intron"].append((prev_e + 1, s - 1))
            prev_e = e

    regions = {c: {k: sorted(iv) for k, iv in d.items()} for c, d in regions.items()}
    starts = {c: {k: [x[0] for x in iv] for k, iv in d.items()} for c, d in regions.items()}

    MM10_CHROM_LEN = {}
    with open(FAI, "rt") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            MM10_CHROM_LEN[parts[0]] = int(parts[1])

    bg_len = defaultdict(int)
    for c, d in regions.items():
        for k, iv in d.items():
            for s, e in merge(iv):
                bg_len[k] += e - s + 1

    def classify(chrom, s, e):
        if chrom not in regions:
            return "intergenic"
        d = regions[chrom]
        for k in CLASSES:
            ivs = d.get(k)
            if not ivs:
                continue
            st = starts[chrom][k]
            i = bisect.bisect_right(st, e) - 1
            if i >= 0:
                for j in range(i, max(i - 3, -1), -1):
                    a, b = ivs[j]
                    if b < s:
                        break
                    if a <= e and b >= s:
                        return k
        return "intergenic"

    def load_peaks_diff():
        df = pd.read_csv(DIFF)
        return df[["chrom", "txStart", "txEnd", "direction"]].dropna()

    def load_peaks_all():
        frames = []
        xl = pd.ExcelFile(ALLX)
        for sh in xl.sheet_names:
            raw = pd.read_excel(ALLX, sheet_name=sh, header=None)
            hdr = None
            for i in range(min(40, len(raw))):
                if str(raw.iloc[i, 0]).strip() == "chrom":
                    hdr = i
                    break
            if hdr is None:
                continue
            df = pd.read_excel(ALLX, sheet_name=sh, header=hdr)
            df = df.dropna(subset=["chrom"])
            frames.append(df[["chrom", "txStart", "txEnd"]])
        return pd.concat(frames, ignore_index=True)

    sets = {"up": None, "down": None, "all": None}
    diff = load_peaks_diff()
    sets["up"] = diff[diff["direction"] == "up"]
    sets["down"] = diff[diff["direction"] == "down"]
    sets["all"] = load_peaks_all()

    counts = pd.DataFrame(0, index=CATS, columns=list(sets))
    for name, df in sets.items():
        for chrom, s, e in zip(df["chrom"], df["txStart"], df["txEnd"]):
            counts.loc[classify(str(chrom), int(s), int(e)), name] += 1

    props = counts.div(counts.sum(axis=0), axis=1)
    counts.to_csv(OUT / "Fig2A_region_counts.csv")
    (props * 100).round(2).to_csv(OUT / "Fig2A_region_proportions.csv")

    priority = ["CDS", "5'UTR", "3'UTR", "lncRNA_exon", "intron"]
    bg_exclusive = {}
    occupied = {c: [] for c in MM10_CHROM_LEN}
    for cls in priority:
        total_len = 0
        for chrom in MM10_CHROM_LEN:
            ivs = merge(regions.get(chrom, {}).get(cls, []))
            exclusive = subtract_intervals(ivs, occupied[chrom])
            total_len += sum(e - s + 1 for s, e in exclusive)
            occupied[chrom] = merge(occupied[chrom] + exclusive)
        bg_exclusive[cls] = total_len
    bg_exclusive["intergenic"] = sum(MM10_CHROM_LEN.values()) - sum(bg_exclusive.values())
    bg_total_v2 = sum(bg_exclusive.values())
    bg_props_v2 = {k: v / bg_total_v2 for k, v in bg_exclusive.items()}

    rows = []
    for region in CATS:
        rows.append({"Region": region, "group": "up", "proportion": props.loc[region, "up"]})
        rows.append({"Region": region, "group": "down", "proportion": props.loc[region, "down"]})
        rows.append({"Region": region, "group": "background", "proportion": bg_props_v2.get(region, 0)})
    pd.DataFrame(rows).to_csv(OUT / "Fig2A_region_prop_long.csv", index=False)

    with open(OUT / "Fig2A_region_summary.txt", "w", encoding="utf-8") as f:
        f.write("Counts\n")
        f.write(counts.to_string() + "\n\n")
        f.write("Proportions (%)\n")
        f.write((props * 100).round(2).to_string() + "\n\n")
        f.write("Background proportions (%)\n")
        for k, v in sorted(bg_props_v2.items(), key=lambda x: -x[1]):
            f.write(f"{k}: {100*v:.2f}\n")
    print(counts)
    print((props * 100).round(2))
    print("step1 genomic distribution done")


def step2_extract_5mer():
    genome_dict = SeqIO.to_dict(SeqIO.parse(str(GENOME_FA), "fasta"))
    df_diff = pd.read_csv(DIFF)
    up_df = df_diff.loc[df_diff["direction"] == "up", ["chrom", "txStart", "txEnd"]].copy()
    down_df = df_diff.loc[df_diff["direction"] == "down", ["chrom", "txStart", "txEnd"]].copy()

    def fetch_seq(chrom, start, end):
        chrom = str(chrom)
        if chrom not in genome_dict:
            return None
        rec = genome_dict[chrom]
        s, e = int(start), int(end)
        if s < 0:
            return None
        return str(rec.seq[s:e]).upper()

    def get_seq_table(input_df, label):
        res = []
        for _, row in input_df.iterrows():
            dna_seq = fetch_seq(row["chrom"], row["txStart"], row["txEnd"])
            res.append({"chrom": row["chrom"], "txStart": row["txStart"], "txEnd": row["txEnd"],
                        "seq": dna_seq, "group": label})
        return pd.DataFrame(res)

    seq_cache = pd.concat([get_seq_table(up_df, "up"), get_seq_table(down_df, "down")],
                          ignore_index=True)
    seq_cache.to_csv(SEQ_CACHE, index=False)
    seq_cache = seq_cache[~seq_cache["seq"].isna()].copy()
    seq_cache = seq_cache[~seq_cache["seq"].str.contains("N")].copy()

    up_seqs = seq_cache.loc[seq_cache["group"] == "up", "seq"].tolist()
    down_seqs = seq_cache.loc[seq_cache["group"] == "down", "seq"].tolist()
    k = 5

    def count_kmer(seq_list, k_size):
        kmer_counts = {}
        for s in seq_list:
            for i in range(len(s) - k_size + 1):
                mer = s[i:i + k_size]
                kmer_counts[mer] = kmer_counts.get(mer, 0) + 1
        return kmer_counts

    up_counts = count_kmer(up_seqs, k)
    down_counts = count_kmer(down_seqs, k)
    all_mers = set(list(up_counts.keys()) + list(down_counts.keys()))
    out_rows = []
    pseudo = 0.5
    total_up = sum(up_counts.values())
    total_down = sum(down_counts.values())
    for mer in all_mers:
        n_up = up_counts.get(mer, 0)
        n_down = down_counts.get(mer, 0)
        table = [[n_up + pseudo, total_up - n_up + pseudo],
                 [n_down + pseudo, total_down - n_down + pseudo]]
        odds, p = fisher_exact(table)
        out_rows.append({"mer": mer, "up_n": n_up, "down_n": n_down,
                         "log2OR": np.log2(odds), "pvalue": p})
    result_df = pd.DataFrame(out_rows)
    n_test = result_df.shape[0]
    result_df["padj"] = result_df["pvalue"] * n_test
    result_df.loc[result_df["padj"] > 1, "padj"] = 1.0
    result_df = result_df.sort_values("log2OR", ascending=False).reset_index(drop=True)
    result_df.to_csv(OUT / "Fig2B_5mer_enrich.csv", index=False)
    print("step2 5-mer enrichment done:", len(result_df), "motifs")


def step3_ccw_density():
    RNG = np.random.default_rng(2026)
    df = pd.read_csv(DIFF)
    cache = pd.read_csv(SEQ_CACHE).drop_duplicates(subset=["chrom", "txStart", "txEnd"])
    cache["key"] = list(zip(cache["chrom"].astype(str),
                            cache["txStart"].astype(int), cache["txEnd"].astype(int)))
    seq_map = {k: str(v).upper() for k, v in zip(cache["key"], cache["seq"])}

    def get_seq(r):
        return seq_map.get((str(r["chrom"]), int(r["txStart"]), int(r["txEnd"])), "")

    def ccw_density(seq):
        if not seq:
            return 0.0
        matches = re.findall(r"(?=(CC[AT]))", seq)
        return len(matches) / (len(seq) / 100.0)

    def shuffle_preserve(seq):
        chars = list(seq)
        RNG.shuffle(chars)
        return "".join(chars)

    df["seq"] = df.apply(get_seq, axis=1)
    df = df[df["seq"] != ""].copy()
    rows = []
    for idx, r in df.iterrows():
        seq = r["seq"]
        shuf = shuffle_preserve(seq)
        rows.append({
            "PeakID": r.get("PeakID", f"peak_{idx}"),
            "direction": r["direction"],
            "ccw_density_real": ccw_density(seq),
            "ccw_density_shuffled": ccw_density(shuf),
            "seq_length": len(seq)})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "Fig2C_ccw_density.csv", index=False)
    with open(OUT / "Fig2C_ccw_summary.txt", "w", encoding="utf-8") as f:
        f.write("CCW (C-C-A/T) density: real peaks vs mononucleotide shuffled background\n\n")
        f.write(f"- Total peaks: {len(out)}\n")
        for dire, label in [("up", "Up"), ("down", "Down")]:
            sub = out[out["direction"] == dire]
            real = sub["ccw_density_real"]
            shuffled = sub["ccw_density_shuffled"]
            stat, p = ranksums(real, shuffled)
            line = (f"- {label}: n={len(sub)}, real={real.mean():.3f}\u00b1{real.std():.3f}, "
                    f"shuffled={shuffled.mean():.3f}\u00b1{shuffled.std():.3f}, Wilcoxon P={p:.3e}\n")
            print(line.strip())
            f.write(line)
    print("step3 CCW density done")


if __name__ == "__main__":
    step1_genomic_distribution()
    step2_extract_5mer()
    step3_ccw_density()
    print("Figure 2 data complete")
