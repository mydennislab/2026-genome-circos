#!/usr/bin/env python3
"""
Circular chromosome ideogram with density tracks and segmental duplication links.
Requires: numpy, matplotlib, pycirclize
"""

import bisect
import math
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from pycirclize import Circos
import argparse
from pathlib import Path

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8, 'figure.dpi': 300, 'savefig.dpi': 300,
})

# --- Colors ---

COL_T2T_COMPLETE   = '#2E7D32'
COL_T2T_PARTIAL    = '#FFA000'
COL_T2T_INCOMPLETE = '#BDBDBD'
COL_CENTROMERE     = '#424242'
COL_GAP            = '#D32F2F'
COL_TELO_YES       = '#2E7D32'
COL_TELO_NO        = '#EF5350'
COL_GENE           = '#1976D2'
COL_REPEAT         = '#FB8C00'
COL_GC             = '#43A047'
COL_SD             = '#7B1FA2'
COL_INTRA_LINK      = '#AB47BC'
COL_INTRA_LINK_GENE = '#FF8F00'
COL_LINK             = '#00897B'
COL_LINK_GENE        = '#E65100'


# --- Parsers ---

def parse_chrom_sizes(filepath):
    """Parse .fai or 2-column TSV (chrom<TAB>size)."""
    sizes = {}
    with open(filepath) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                try:
                    sizes[parts[0]] = int(parts[1])
                except ValueError:
                    continue
    return sizes


def parse_bed(filepath):
    if filepath is None or not Path(filepath).exists():
        return {}
    regions = {}
    with open(filepath) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                try:
                    regions.setdefault(parts[0], []).append(
                        (int(parts[1]), int(parts[2])))
                except ValueError:
                    continue
    return regions


def parse_tidk(filepath):
    """Parse tidk search output — first and last window per chrom."""
    if filepath is None or not Path(filepath).exists():
        return {}
    telomeres = {}
    header = None
    chrom_rows = {}
    with open(filepath) as f:
        for line in f:
            if header is None:
                header = line.strip().split('\t')
                continue
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            chrom = parts[0]
            total = int(parts[2]) + int(parts[3])
            window = int(parts[1])
            if chrom not in chrom_rows:
                chrom_rows[chrom] = {'first_win': window, 'first_val': total,
                                     'last_win': window, 'last_val': total}
            else:
                if window < chrom_rows[chrom]['first_win']:
                    chrom_rows[chrom]['first_win'] = window
                    chrom_rows[chrom]['first_val'] = total
                if window > chrom_rows[chrom]['last_win']:
                    chrom_rows[chrom]['last_win'] = window
                    chrom_rows[chrom]['last_val'] = total
    for chrom, d in chrom_rows.items():
        telomeres[chrom] = {'5prime': d['first_val'], '3prime': d['last_val']}
    return telomeres


def parse_repeatmasker(filepath):
    """Parse RepeatMasker .out → {chrom: [(start, end), ...]}."""
    repeats = {}
    if filepath is None or not Path(filepath).exists():
        return repeats
    with open(filepath) as f:
        for i, line in enumerate(f):
            if i < 3:
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            try:
                chrom = parts[4]
                start, end = int(parts[5]), int(parts[6])
                repeats.setdefault(chrom, []).append((start, end))
            except (ValueError, IndexError):
                continue
    return repeats


def parse_gff_genes(filepath):
    """Parse GFF3 → gene positions {chrom: [(start, end), ...]}."""
    genes = {}
    if filepath is None or not Path(filepath).exists():
        return genes
    with open(filepath) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 9 and parts[2] == 'gene':
                try:
                    genes.setdefault(parts[0], []).append(
                        (int(parts[3]), int(parts[4])))
                except ValueError:
                    continue
    return genes


def parse_gff_genes_named(filepath):
    """Parse GFF3 → {chrom: [(start, end, name), ...]} sorted by start."""
    genes = {}
    if filepath is None or not Path(filepath).exists():
        return genes
    with open(filepath) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 9 and parts[2] == 'gene':
                try:
                    start, end = int(parts[3]), int(parts[4])
                except ValueError:
                    continue
                name = ''
                for token in parts[8].split(';'):
                    if token.startswith('Name='):
                        name = token[5:]
                        break
                    if token.startswith('ID=') and not name:
                        name = token[3:]
                genes.setdefault(parts[0], []).append((start, end, name))
    for chrom in genes:
        genes[chrom].sort(key=lambda x: x[0])
    return genes


def parse_segdups(filepath, chrom_prefix=""):
    """Parse BED or BEDPE segmental duplications.
    BEDPE (BISER/SEDEF format) with divergence in col 8: keeps <10% divergence and >=1 kb.
    Returns (intra_regions, intra_links, inter_regions, inter_links)."""
    intra_regions, intra_links = {}, []
    inter_regions, inter_links = {}, []
    if filepath is None or not Path(filepath).exists():
        return intra_regions, intra_links, inter_regions, inter_links
    is_bedpe = None
    with open(filepath) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            # detect format from first data line
            if is_bedpe is None:
                is_bedpe = (len(parts) >= 6
                            and not parts[3].isdigit()
                            and bool(re.match(r'^[A-Za-z_]', parts[3])))
            if is_bedpe:
                if len(parts) < 6:
                    continue
                try:
                    chrom1, s1, e1 = parts[0], int(parts[1]), int(parts[2])
                    chrom2, s2, e2 = parts[3], int(parts[4]), int(parts[5])
                except ValueError:
                    continue
                avg_len = ((e1 - s1) + (e2 - s2)) / 2.0
                if len(parts) > 7:
                    try:
                        divergence = float(parts[7])
                        if divergence >= 10.0 or avg_len < 1000:
                            continue
                    except ValueError:
                        pass
                # filter by prefix if given
                c1_ok = not chrom_prefix or chrom1.startswith(chrom_prefix)
                c2_ok = not chrom_prefix or chrom2.startswith(chrom_prefix)
                if chrom1 == chrom2:
                    for c, s, e in [(chrom1, s1, e1), (chrom2, s2, e2)]:
                        intra_regions.setdefault(c, []).append((s, e))
                    if c1_ok:
                        intra_links.append((chrom1, s1, e1, s2, e2, avg_len))
                else:
                    for c, s, e in [(chrom1, s1, e1), (chrom2, s2, e2)]:
                        inter_regions.setdefault(c, []).append((s, e))
                    if c1_ok and c2_ok:
                        inter_links.append((chrom1, s1, e1, chrom2, s2, e2, avg_len))
            else:
                try:
                    chrom = parts[0]
                    s, e = int(parts[1]), int(parts[2])
                    intra_regions.setdefault(chrom, []).append((s, e))
                except ValueError:
                    continue
    return intra_regions, intra_links, inter_regions, inter_links


def parse_t2t_status(filepath):
    """Parse T2T summary TSV → {chrom: 'complete'|'partial'|'incomplete'}."""
    status = {}
    if filepath is None or not Path(filepath).exists():
        return status
    with open(filepath) as f:
        header = None
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            if header is None:
                header = parts
                continue
            if len(parts) < 2:
                continue
            chrom = parts[0]
            s = parts[1].strip()
            if 'T2T' in s:
                status[chrom] = 'complete'
            elif 'only' in s:
                status[chrom] = 'partial'
            else:
                status[chrom] = 'incomplete'
    return status


def parse_gc(filepath):
    """Parse GC BED → {chrom: [(start, end, gc_frac), ...]}."""
    gc = {}
    if filepath is None or not Path(filepath).exists():
        return gc
    with open(filepath) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            try:
                gc.setdefault(parts[0], []).append(
                    (int(parts[1]), int(parts[2]), float(parts[3])))
            except ValueError:
                continue
    return gc


# --- Helpers ---

def strip_prefix(name, prefix):
    if prefix and name.startswith(prefix):
        return name[len(prefix):]
    return name


def chrom_sort_key(chrom, prefix=""):
    name = strip_prefix(chrom, prefix)
    if name.startswith('chr'):
        name = name[3:]
    if name == 'X':  return (1, 100, chrom)
    if name == 'Y':  return (1, 101, chrom)
    if name in ('M', 'MT'): return (1, 102, chrom)
    try:
        return (0, int(name), chrom)
    except ValueError:
        m = re.match(r'(\d+)', name)
        return (0, int(m.group(1)), chrom) if m else (2, 0, chrom)


def compute_density(features, chrom_size, window_size):
    n_win = int(np.ceil(chrom_size / window_size))
    positions = np.minimum(
        np.arange(n_win) * window_size + window_size / 2, chrom_size - 1)
    starts = np.array([s for s, e in features])
    ends = np.array([e for s, e in features])
    lengths = ends - starts
    bins = np.clip(starts // window_size, 0, n_win - 1).astype(int)
    density = np.zeros(n_win)
    np.add.at(density, bins, lengths)
    density /= window_size
    np.minimum(density, 1.0, out=density)
    return positions, density


def intersect_sd_with_genes(genes_named, intra_links, inter_links, output_tsv):
    """Intersect SD link endpoints with gene intervals using bisect search.
    Returns (intra_gene_flags, inter_gene_flags) as parallel bool lists."""
    chrom_starts = {
        chrom: [g[0] for g in ivs]
        for chrom, ivs in genes_named.items()
    }

    def overlapping_genes(chrom, s, e):
        ivs = genes_named.get(chrom, [])
        if not ivs:
            return []
        starts = chrom_starts[chrom]
        idx = bisect.bisect_left(starts, e)
        return [ivs[i][2] for i in range(idx) if ivs[i][1] > s]

    intra_gene_flags, inter_gene_flags = [], []
    rows = []

    for chrom, s1, e1, s2, e2, avg_len in intra_links:
        g1, g2 = overlapping_genes(chrom, s1, e1), overlapping_genes(chrom, s2, e2)
        has_gene = bool(g1 or g2)
        intra_gene_flags.append(has_gene)
        rows.append((chrom, s1, e1, chrom, s2, e2, 'intra', has_gene,
                     ';'.join(g1) or '.', ';'.join(g2) or '.'))

    for c1, s1, e1, c2, s2, e2, avg_len in inter_links:
        g1, g2 = overlapping_genes(c1, s1, e1), overlapping_genes(c2, s2, e2)
        has_gene = bool(g1 or g2)
        inter_gene_flags.append(has_gene)
        rows.append((c1, s1, e1, c2, s2, e2, 'inter', has_gene,
                     ';'.join(g1) or '.', ';'.join(g2) or '.'))

    with open(output_tsv, 'w') as fh:
        fh.write('chrom1\tstart1\tend1\tchrom2\tstart2\tend2\ttype\thas_gene\tgenes_region1\tgenes_region2\n')
        for row in rows:
            fh.write('\t'.join(str(x) for x in row) + '\n')

    n_gi, n_ge = sum(intra_gene_flags), sum(inter_gene_flags)
    print(f"  Gene-containing SD links: {n_gi} intra, {n_ge} inter")
    print(f"  Saved: {output_tsv}")
    return intra_gene_flags, inter_gene_flags


# --- Main ---

def plot_circos(chrom_sizes, centromeres=None, telomeres=None, gaps=None,
                genes=None, repeats=None, segdups_data=None,
                gc_data=None, t2t_status=None,
                output_prefix='assembly', title='Circos Ideogram',
                window_size=1_000_000, chrom_prefix="", min_size=0,
                show_links=True, show_sd_intra=True, show_sd_inter=True,
                sd_gene_links_data=None):

    centromeres = centromeres or {}
    telomeres = telomeres or {}
    gaps = gaps or {}
    genes = genes or {}
    repeats = repeats or {}
    gc_data = gc_data or {}
    t2t_status = t2t_status or {}

    if segdups_data is None:
        segdups_data = ({}, [], {}, [])
    sd_intra_regions, sd_intra_links, sd_inter_regions, sd_inter_links = segdups_data

    if sd_gene_links_data is None:
        intra_gene_flags = [False] * len(sd_intra_links)
        inter_gene_flags = [False] * len(sd_inter_links)
        has_gene_links = False
    else:
        intra_gene_flags, inter_gene_flags = sd_gene_links_data
        has_gene_links = True

    # Filter and sort chromosomes
    chromosomes = sorted(
        [c for c in chrom_sizes
         if (not chrom_prefix or c.startswith(chrom_prefix))
         and chrom_sizes[c] >= min_size],
        key=lambda c: chrom_sort_key(c, chrom_prefix))

    if not chromosomes:
        print("ERROR: No chromosomes matched. Check --prefix and --min-size.")
        return

    n_chrom = len(chromosomes)
    max_size = max(chrom_sizes[c] for c in chromosomes)
    print(f"  {n_chrom} chromosomes, largest = {max_size / 1e6:.1f} Mb")

    # Build sector dict with display names
    sectors = {}
    name_map = {}      # original → display
    reverse_map = {}   # display → original
    for c in chromosomes:
        display = strip_prefix(c, chrom_prefix)
        sectors[display] = chrom_sizes[c]
        name_map[c] = display
        reverse_map[display] = c

    # Precompute densities with genome-wide normalization
    has_genes = any(genes.get(c) for c in chromosomes)
    has_repeats = any(repeats.get(c) for c in chromosomes)
    has_gc = any(gc_data.get(c) for c in chromosomes)
    has_telo = any(telomeres.get(c) for c in chromosomes)
    has_sd = any(sd_intra_regions.get(c) or sd_inter_regions.get(c) for c in chromosomes)

    gene_densities, repeat_densities = {}, {}
    for chrom in chromosomes:
        cs = chrom_sizes[chrom]
        g = genes.get(chrom, [])
        if g:
            gene_densities[chrom] = compute_density(g, cs, window_size)
        r = repeats.get(chrom, [])
        if r:
            repeat_densities[chrom] = compute_density(r, cs, window_size)

    gene_gmax = max((d.max() for _, d in gene_densities.values()), default=1)
    repeat_gmax = max((d.max() for _, d in repeat_densities.values()), default=1)
    gc_gmax = max((gc for c in chromosomes for _, _, gc in gc_data.get(c, [])),
                  default=0.5)

    # Determine which tracks to draw
    track_list = []
    track_list.append('ideogram')
    if has_telo:
        track_list.append('telomere')
    if has_gc:
        track_list.append('gc')
    if has_genes:
        track_list.append('gene')
    if has_repeats:
        track_list.append('repeat')
    if has_sd:
        track_list.append('sd')

    # Compute radius ranges dynamically based on active tracks
    # Ideogram always at 94-100, others stack inward
    track_height = {
        'ideogram': 6, 'telomere': 3, 'gc': 7,
        'gene': 10, 'repeat': 10, 'sd': 8,
    }
    track_gap = 1
    radii = {}
    r_top = 100
    for t in track_list:
        h = track_height[t]
        radii[t] = (r_top - h, r_top)
        r_top = r_top - h - track_gap

    link_radius = radii[track_list[-1]][0] if track_list else 48

    # Initialize Circos
    circos = Circos(sectors, space=3)
    print("  Drawing tracks...")

    for sector in circos.sectors:
        orig = reverse_map[sector.name]

        # Track 1: Ideogram — colored by T2T status
        r_lo, r_hi = radii['ideogram']
        ideo_track = sector.add_track((r_lo, r_hi))
        status = t2t_status.get(orig, 'incomplete')
        if t2t_status:
            if status == 'complete':
                ideo_fc = COL_T2T_COMPLETE
            elif status == 'partial':
                ideo_fc = COL_T2T_PARTIAL
            else:
                ideo_fc = COL_T2T_INCOMPLETE
        else:
            ideo_fc = '#E8E8E8'
        ideo_track.axis(fc=ideo_fc, ec='#757575', lw=0.5, alpha=0.35)

        # Centromere
        for c_start, c_end in centromeres.get(orig, []):
            ideo_track.rect(c_start, c_end, fc=COL_CENTROMERE, ec='none')

        # Assembly gaps
        if orig in gaps:
            cs = chrom_sizes[orig]
            for gs, ge in gaps[orig]:
                raw_w = ge - gs
                min_w = sector.size * 0.003
                ge_c = min(gs + max(raw_w, min_w), cs - 1)
                if ge_c > gs:
                    ideo_track.rect(gs, ge_c, fc=COL_GAP, ec='none')

        sector.text(sector.name, r=radii['ideogram'][1] + 8, size=7, fontweight='bold')

        # Track 2: Telomere indicators
        if 'telomere' in radii:
            r_lo, r_hi = radii['telomere']
            telo_track = sector.add_track((r_lo, r_hi))
            telo_track.axis(fc='#FAFAFA', ec='none')
            td = telomeres.get(orig, {})
            tw = sector.size * 0.03
            telo_track.rect(0, tw,
                            fc=COL_TELO_YES if td.get('5prime', 0) >= 10 else COL_TELO_NO,
                            ec='none')
            telo_track.rect(sector.size - tw, sector.size,
                            fc=COL_TELO_YES if td.get('3prime', 0) >= 10 else COL_TELO_NO,
                            ec='none')

        # Track 3: GC content
        if 'gc' in radii:
            r_lo, r_hi = radii['gc']
            gc_track = sector.add_track((r_lo, r_hi))
            gc_track.axis(fc='#FAFAFA', ec='#E0E0E0', lw=0.3)
            gc_windows = gc_data.get(orig, [])
            if gc_windows:
                gc_pos = np.array([(s + e) / 2 for s, e, _ in gc_windows])
                gc_vals = np.array([g for _, _, g in gc_windows])
                gc_pos = np.minimum(gc_pos, sector.size - 1)
                gc_norm = gc_vals / gc_gmax
                gc_track.fill_between(gc_pos, gc_norm * (r_hi - r_lo), 0,
                                      fc=COL_GC, ec='none', alpha=0.8)

        # Track 4: Gene density
        if 'gene' in radii:
            r_lo, r_hi = radii['gene']
            gene_track = sector.add_track((r_lo, r_hi))
            gene_track.axis(fc='#FAFAFA', ec='#E0E0E0', lw=0.3)
            if orig in gene_densities:
                pos, dens = gene_densities[orig]
                dens_norm = dens / gene_gmax
                gene_track.fill_between(pos, dens_norm * (r_hi - r_lo), 0,
                                        fc=COL_GENE, ec='none', alpha=0.8)

        # Track 5: Repeat density
        if 'repeat' in radii:
            r_lo, r_hi = radii['repeat']
            repeat_track = sector.add_track((r_lo, r_hi))
            repeat_track.axis(fc='#FAFAFA', ec='#E0E0E0', lw=0.3)
            if orig in repeat_densities:
                pos, dens = repeat_densities[orig]
                dens_norm = dens / repeat_gmax
                repeat_track.fill_between(pos, dens_norm * (r_hi - r_lo), 0,
                                          fc=COL_REPEAT, ec='none', alpha=0.8)

        # Track 6: SD regions
        if 'sd' in radii:
            r_lo, r_hi = radii['sd']
            sd_track = sector.add_track((r_lo, r_hi))
            sd_track.axis(fc='#FAFAFA', ec='#E0E0E0', lw=0.3)
            sd_rects = []
            if show_sd_intra:
                sd_rects.extend(sd_intra_regions.get(orig, []))
            if show_sd_inter:
                sd_rects.extend(sd_inter_regions.get(orig, []))
            for sd_s, sd_e in sd_rects:
                sd_track.rect(sd_s, sd_e, fc=COL_SD, ec='none', alpha=0.7)

    # SD alpha scaled by log(length)
    _log_min = math.log10(1_000)
    _log_max = math.log10(500_000)

    def _sd_alpha(avg_len):
        t = (math.log10(max(avg_len, 1_000)) - _log_min) / (_log_max - _log_min)
        return 0.12 + 0.53 * min(max(t, 0.0), 1.0)

    # Intra-chromosomal arcs
    if show_sd_intra and sd_intra_links:
        print(f"  Drawing {len(sd_intra_links)} intra-chr SD arcs...")
        for i, (chrom, s1, e1, s2, e2, avg_len) in enumerate(sd_intra_links):
            name = name_map.get(chrom)
            if name:
                alpha = _sd_alpha(avg_len)
                has_gene = has_gene_links and intra_gene_flags[i]
                fc = COL_INTRA_LINK_GENE if has_gene else COL_INTRA_LINK
                circos.link((name, s1, e1), (name, s2, e2),
                            r1=link_radius, r2=link_radius,
                            fc=fc, ec=fc, lw=0.1, alpha=alpha, height_ratio=0.2)

    # Inter-chromosomal links
    if show_links and show_sd_inter and sd_inter_links:
        print(f"  Drawing {len(sd_inter_links)} inter-chr SD links...")
        for i, (c1, s1, e1, c2, s2, e2, avg_len) in enumerate(sd_inter_links):
            n1, n2 = name_map.get(c1), name_map.get(c2)
            if n1 and n2:
                alpha = _sd_alpha(avg_len)
                has_gene = has_gene_links and inter_gene_flags[i]
                fc = COL_LINK_GENE if has_gene else COL_LINK
                circos.link((n1, s1, e1), (n2, s2, e2),
                            fc=fc, ec=fc, lw=0.1, alpha=alpha)

    # Render figure with gridspec: left = circos (polar), right = legend
    fig = plt.figure(figsize=(18, 12))
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1.2],
                           left=0.02, right=0.98, top=0.93, bottom=0.03,
                           wspace=0.05)

    ax_circos = fig.add_subplot(gs[0, 0], projection='polar')
    circos.plotfig(ax=ax_circos)
    fig.suptitle(title, fontsize=12, fontweight='bold', y=0.97)

    # Legend
    ax_legend = fig.add_subplot(gs[0, 1])
    ax_legend.set_xlim(0, 1)
    ax_legend.set_ylim(0, 1)
    ax_legend.axis('off')

    track_legend = []
    if t2t_status:
        track_legend.append(("Ideogram (T2T status)", [
            (COL_T2T_COMPLETE, 0.45, "T2T complete"),
            (COL_T2T_PARTIAL,  0.45, "Partial (one telomere)"),
            (COL_T2T_INCOMPLETE, 0.45, "Incomplete"),
            (COL_CENTROMERE, 1.0, "Centromere"),
            (COL_GAP, 1.0, "Assembly gap"),
        ]))
    elif centromeres or gaps:
        items = []
        if centromeres:
            items.append((COL_CENTROMERE, 1.0, "Centromere"))
        if gaps:
            items.append((COL_GAP, 1.0, "Assembly gap"))
        track_legend.append(("Ideogram", items))

    if has_telo:
        track_legend.append(("Telomere status", [
            (COL_TELO_YES, 1.0, "Present"),
            (COL_TELO_NO, 1.0, "Absent"),
        ]))
    if has_gc:
        track_legend.append(("GC content", [
            (COL_GC, 0.8, f"GC fraction per window"),
        ]))
    if has_genes:
        track_legend.append(("Gene density", [
            (COL_GENE, 0.8, "Genes per window"),
        ]))
    if has_repeats:
        track_legend.append(("Repeat density", [
            (COL_REPEAT, 0.8, "Repeats per window"),
        ]))
    if has_sd:
        sd_items = []
        if show_sd_intra:
            sd_items.append((COL_SD, 0.7, "Intra-chr SD region"))
            sd_items.append((COL_INTRA_LINK, 0.65, "Intra-chr arc"))
            if has_gene_links:
                sd_items.append((COL_INTRA_LINK_GENE, 0.65, "Intra-chr arc (gene)"))
        if show_sd_inter:
            sd_items.append((COL_SD, 0.4, "Inter-chr SD region"))
            if show_links:
                sd_items.append((COL_LINK, 0.65, "Inter-chr link"))
                if has_gene_links:
                    sd_items.append((COL_LINK_GENE, 0.65, "Inter-chr link (gene)"))
        if sd_items:
            track_legend.append(("Segmental duplications", sd_items))

    y = 0.98
    line_h = 0.035
    swatch_w, swatch_h = 0.05, 0.02
    for track_label, items in track_legend:
        ax_legend.text(0.0, y, track_label, fontsize=8, fontweight='bold',
                       va='top', transform=ax_legend.transAxes)
        y -= line_h
        for color, alpha, label in items:
            ax_legend.add_patch(plt.Rectangle(
                (0.04, y - swatch_h * 0.5), swatch_w, swatch_h,
                fc=color, alpha=alpha, ec='none',
                transform=ax_legend.transAxes, clip_on=False))
            ax_legend.text(0.04 + swatch_w + 0.02, y, label, fontsize=7,
                           va='center', transform=ax_legend.transAxes)
            y -= line_h * 0.82
        y -= line_h * 0.35

    for ext in ['png', 'pdf']:
        outpath = f"{output_prefix}_circos.{ext}"
        fig.savefig(outpath, dpi=300, facecolor='white')
        print(f"Saved: {outpath}")
    plt.close()


if __name__ == '__main__':
    p = argparse.ArgumentParser(
        description='Circular chromosome ideogram with density tracks and SD links')
    p.add_argument('--chrom-sizes', required=True,
                   help='.fai or 2-column TSV (chrom, size)')
    p.add_argument('--centromere', default=None, help='Centromere BED')
    p.add_argument('--tidk', default=None, help='tidk search TSV')
    p.add_argument('--gaps', default=None, help='Assembly gaps BED')
    p.add_argument('--gff', default=None, help='Gene annotation GFF3')
    p.add_argument('--repeatmasker', default=None, help='RepeatMasker .out')
    p.add_argument('--segdups', default=None, help='Segmental duplications BED/BEDPE')
    p.add_argument('--gc', default=None,
                   help='GC content BED (chrom, start, end, gc_frac)')
    p.add_argument('--t2t', default=None, help='T2T status TSV (chromosome, status)')
    p.add_argument('-o', '--output', default='assembly', help='Output prefix')
    p.add_argument('-t', '--title', default='Circos Ideogram')
    p.add_argument('-w', '--window', type=int, default=1_000_000,
                   help='Density window size in bp (default: 1000000)')
    p.add_argument('--prefix', default='', help='Chromosome name prefix filter')
    p.add_argument('--min-size', type=int, default=0, help='Min chromosome size (bp)')
    p.add_argument('--no-links', action='store_true',
                   help='Suppress inter-chr SD links')
    p.add_argument('--sd-intra', action='store_true', default=None,
                   help='Show intra-chr SD arcs only')
    p.add_argument('--sd-inter', action='store_true', default=None,
                   help='Show inter-chr SD links only')
    p.add_argument('--sd-gene-links', action='store_true',
                   help='Recolor SD links by gene overlap; save TSV')
    args = p.parse_args()

    # If neither --sd-intra / --sd-inter given, default to both
    if not args.sd_intra and not args.sd_inter:
        show_sd_intra = True
        show_sd_inter = True
    else:
        show_sd_intra = bool(args.sd_intra)
        show_sd_inter = bool(args.sd_inter)

    print(f"Loading data for {args.title}...")
    chrom_sizes = parse_chrom_sizes(args.chrom_sizes)
    centromeres = parse_bed(args.centromere)
    gaps = parse_bed(args.gaps)
    telomeres = parse_tidk(args.tidk)
    repeats = parse_repeatmasker(args.repeatmasker)
    genes = parse_gff_genes(args.gff)
    gc_data = parse_gc(args.gc)
    t2t_status = parse_t2t_status(args.t2t)
    segdups_data = parse_segdups(args.segdups, chrom_prefix=args.prefix)

    sd_intra_links = segdups_data[1]
    sd_inter_links = segdups_data[3]
    print(f"  SD: {len(sd_intra_links)} intra-chr links, {len(sd_inter_links)} inter-chr links")

    # SD-gene intersection
    sd_gene_links_data = None
    if args.sd_gene_links and args.gff:
        genes_named = parse_gff_genes_named(args.gff)
        tsv_path = f"{args.output}_sd_gene_links.tsv"
        sd_gene_links_data = intersect_sd_with_genes(
            genes_named, sd_intra_links, sd_inter_links, tsv_path)

    if args.repeatmasker:
        print("  Parsing RepeatMasker (this may take a moment)...")

    plot_circos(
        chrom_sizes, centromeres=centromeres, telomeres=telomeres,
        gaps=gaps, genes=genes, repeats=repeats, segdups_data=segdups_data,
        gc_data=gc_data, t2t_status=t2t_status,
        output_prefix=args.output, title=args.title, window_size=args.window,
        chrom_prefix=args.prefix, min_size=args.min_size,
        show_links=not args.no_links,
        show_sd_intra=show_sd_intra, show_sd_inter=show_sd_inter,
        sd_gene_links_data=sd_gene_links_data)
