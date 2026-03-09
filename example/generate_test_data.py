#!/usr/bin/env python3
"""Generate synthetic test data for a 6-chromosome genome (~180 Mb)."""

import random
import math
import os

random.seed(42)

CHROMS = [
    ('chr1', 45_000_000), ('chr2', 38_000_000), ('chr3', 32_000_000),
    ('chr4', 28_000_000), ('chr5', 22_000_000), ('chr6', 15_000_000),
]

OUTDIR = os.path.dirname(os.path.abspath(__file__))


def write_chrom_sizes():
    with open(os.path.join(OUTDIR, 'chrom.sizes'), 'w') as f:
        for name, size in CHROMS:
            f.write(f"{name}\t{size}\n")


def write_centromeres():
    with open(os.path.join(OUTDIR, 'centromeres.bed'), 'w') as f:
        for name, size in CHROMS:
            mid = size // 2 + random.randint(-size // 8, size // 8)
            w = random.randint(500_000, 2_000_000)
            f.write(f"{name}\t{mid - w // 2}\t{mid + w // 2}\n")


def write_tidk():
    with open(os.path.join(OUTDIR, 'tidk.tsv'), 'w') as f:
        f.write("id\twindow\tforward_repeat_number\treverse_repeat_number\n")
        window_size = 10_000
        for name, size in CHROMS:
            n_windows = size // window_size
            has_5 = random.random() > 0.2
            has_3 = random.random() > 0.2
            for i in range(n_windows):
                pos = i * window_size
                fwd, rev = 0, 0
                if has_5 and i < 3:
                    fwd = random.randint(20, 80)
                    rev = random.randint(10, 40)
                elif has_3 and i > n_windows - 4:
                    fwd = random.randint(20, 80)
                    rev = random.randint(10, 40)
                f.write(f"{name}\t{pos}\t{fwd}\t{rev}\n")


def write_gaps():
    with open(os.path.join(OUTDIR, 'gaps.bed'), 'w') as f:
        for name, size in CHROMS:
            n_gaps = random.randint(0, 3)
            for _ in range(n_gaps):
                pos = random.randint(1_000_000, size - 1_000_000)
                gap_len = random.randint(500, 5000)
                f.write(f"{name}\t{pos}\t{pos + gap_len}\n")


def write_repeatmasker():
    classes = [
        ('LINE/L1', 0.025), ('LINE/L2', 0.012),
        ('SINE/Alu', 0.018), ('SINE/MIR', 0.008),
        ('LTR/ERV1', 0.010), ('LTR/ERVL', 0.006),
        ('DNA/hAT', 0.008), ('DNA/TcMar', 0.005),
        ('Simple_repeat', 0.012), ('Satellite/centr', 0.004),
    ]
    with open(os.path.join(OUTDIR, 'repeats.out'), 'w') as f:
        f.write("   SW   perc perc perc  query     position in query    matching       repeat         position in repeat\n")
        f.write("score   div. del. ins.  sequence  begin end  (left)    repeat         class/family   begin end  (left)  ID\n")
        f.write("\n")
        for name, size in CHROMS:
            for class_family, base_rate in classes:
                # position-dependent rate
                n_features = int(size / 1000 * base_rate)
                for _ in range(n_features):
                    pos = random.randint(0, size - 1000)
                    frac = pos / size
                    local_rate = 1.0 + 0.8 * math.sin(frac * math.pi * 4)
                    if random.random() > local_rate * 0.5:
                        continue
                    length = random.randint(100, 800)
                    end = min(pos + length, size)
                    left = size - end
                    strand = '+' if random.random() > 0.3 else 'C'
                    rep_name = class_family.split('/')[1] if '/' in class_family else class_family
                    f.write(f"  {random.randint(200,2000)}  {random.uniform(5,30):.1f}  "
                            f"0.0  0.0  {name}  {pos}  {end}  ({left})  "
                            f"{strand}  {rep_name}  {class_family}  1  {length}  (0)  {random.randint(1,9999)}\n")


def write_segdups():
    with open(os.path.join(OUTDIR, 'segdups.bedpe'), 'w') as f:
        # Intra-chromosomal
        for name, size in CHROMS:
            n_sd = random.randint(3, 10)
            for _ in range(n_sd):
                s1 = random.randint(0, size - 5_000_000)
                len1 = random.randint(2000, 50_000)
                s2 = s1 + random.randint(1_000_000, 5_000_000)
                len2 = len1 + random.randint(-500, 500)
                if s2 + len2 > size:
                    continue
                error = random.uniform(0.5, 8.0)
                f.write(f"{name}\t{s1}\t{s1 + len1}\t{name}\t{s2}\t{s2 + len2}\t.\t{error:.2f}\n")
        # Inter-chromosomal
        for i in range(len(CHROMS)):
            for j in range(i + 1, len(CHROMS)):
                if random.random() > 0.4:
                    continue
                c1, s1 = CHROMS[i]
                c2, s2 = CHROMS[j]
                n_links = random.randint(1, 4)
                for _ in range(n_links):
                    p1 = random.randint(0, s1 - 50_000)
                    p2 = random.randint(0, s2 - 50_000)
                    length = random.randint(2000, 30_000)
                    error = random.uniform(0.5, 7.0)
                    f.write(f"{c1}\t{p1}\t{p1 + length}\t{c2}\t{p2}\t{p2 + length}\t.\t{error:.2f}\n")


def write_genes_gff():
    with open(os.path.join(OUTDIR, 'genes.gff3'), 'w') as f:
        f.write("##gff-version 3\n")
        gene_id = 0
        for name, size in CHROMS:
            # ~20-40 genes per Mb (spread across chromosome)
            n_genes = int(size / 1e6 * random.randint(20, 40))
            positions = sorted(random.sample(range(1000, size - 10_000), n_genes))
            for pos in positions:
                gene_id += 1
                length = random.randint(1000, 50_000)
                end = min(pos + length, size)
                strand = '+' if random.random() > 0.5 else '-'
                f.write(f"{name}\ttest\tgene\t{pos}\t{end}\t.\t{strand}\t.\t"
                        f"ID=gene{gene_id:05d};Name=GENE{gene_id:05d}\n")


def write_gc():
    window = 100_000
    with open(os.path.join(OUTDIR, 'gc.bed'), 'w') as f:
        for name, size in CHROMS:
            n_windows = size // window
            for i in range(n_windows):
                start = i * window
                end = start + window
                frac = i / n_windows
                # GC varies 0.35-0.55 with some spatial pattern
                gc = 0.40 + 0.08 * math.sin(frac * math.pi * 6) + random.uniform(-0.03, 0.03)
                gc = max(0.30, min(0.60, gc))
                f.write(f"{name}\t{start}\t{end}\t{gc:.4f}\n")


def write_t2t():
    statuses = ['T2T', 'telomere_only', 'incomplete']
    with open(os.path.join(OUTDIR, 't2t.tsv'), 'w') as f:
        f.write("chromosome\tstatus\n")
        for name, _ in CHROMS:
            # weight toward T2T for realism
            s = random.choices(statuses, weights=[5, 2, 1])[0]
            f.write(f"{name}\t{s}\n")


if __name__ == '__main__':
    write_chrom_sizes()
    write_centromeres()
    write_tidk()
    write_gaps()
    write_repeatmasker()
    write_segdups()
    write_genes_gff()
    write_gc()
    write_t2t()
    print("Generated test data in", OUTDIR)
