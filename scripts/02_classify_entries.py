#!/usr/bin/env python3
"""Phase 2: Classify 23,976 plant protein identifiers into buckets A-E."""

import csv
import re
import sys
from pathlib import Path

INPUT_CSV  = Path("/local/storage/jedric/1_GeneOptimization/data/genes_dedup.csv")
OUTPUT_CSV = Path("/local/storage/jedric/1_GeneOptimization/outputs/classified_entries.csv")

# ── Bucket A patterns ─────────────────────────────────────────────────────────
# UniProt: 6-char old format
RE_UNIPROT_6  = re.compile(r'^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9])$')
# UniProt: 10-char new format (e.g. A0A178WN35)
RE_UNIPROT_10 = re.compile(r'^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9]$')

# Locus tags — case-insensitive where biologically appropriate
RE_ARABIDOPSIS = re.compile(r'^AT[1-5CMG]G\d{5}(\.\d+)?$', re.IGNORECASE)
RE_RICE_RAP    = re.compile(r'^Os\d{2}g\d{7}$', re.IGNORECASE)
RE_RICE_MSU    = re.compile(r'^LOC_Os\d{2}g\d{5}$', re.IGNORECASE)
RE_TOMATO      = re.compile(r'^Solyc\d{2}g\d{6}(\.\d+)?(\.\d+)?$', re.IGNORECASE)
RE_MAIZE_GRMZM = re.compile(r'^GRMZM\d+G\d+$', re.IGNORECASE)
RE_MAIZE_ZM    = re.compile(r'^Zm\d+d\d+$', re.IGNORECASE)
RE_SOYBEAN     = re.compile(r'^Glyma\.\d+G\d+(\.\d+)?$', re.IGNORECASE)
RE_SORGHUM     = re.compile(r'^Sobic\.\d+G\d+(\.\d+)?$', re.IGNORECASE)
RE_BRACHY      = re.compile(r'^Bradi\d+g\d+(\.\d+)?$', re.IGNORECASE)
RE_MEDICAGO    = re.compile(r'^Medtr\d+\w+\d+(\.\d+)?$', re.IGNORECASE)
RE_REFSEQ      = re.compile(r'^(NP|XP|NM|XM|NR)_\d+(\.\d+)?$')
RE_GI          = re.compile(r'^gi\|\d+$')

BUCKET_A_CHECKS = [
    (RE_UNIPROT_6,   "uniprot_accession"),
    (RE_UNIPROT_10,  "uniprot_accession"),
    (RE_ARABIDOPSIS, "arabidopsis_locus"),
    (RE_RICE_RAP,    "rice_locus_rap"),
    (RE_RICE_MSU,    "rice_locus_msu"),
    (RE_TOMATO,      "tomato_locus"),
    (RE_MAIZE_GRMZM, "maize_locus"),
    (RE_MAIZE_ZM,    "maize_locus"),
    (RE_SOYBEAN,     "soybean_locus"),
    (RE_SORGHUM,     "sorghum_locus"),
    (RE_BRACHY,      "brachypodium_locus"),
    (RE_MEDICAGO,    "medicago_locus"),
    (RE_REFSEQ,      "refseq_id"),
    (RE_GI,          "ncbi_gi"),
]

# ── Bucket B: species-prefixed gene symbols ───────────────────────────────────
# 2-3 lowercase after first capital, then uppercase gene symbol
# Must not contain spaces, must be at least ~6 chars total
RE_SPECIES_PREFIX = re.compile(
    r'^([A-Z][a-z]{1,2})([A-Z][A-Za-z0-9][-A-Za-z0-9]*)$'
)

KNOWN_PREFIXES = {
    'At', 'Os', 'Sl', 'Zm', 'Gm', 'Mt', 'Ta', 'Gh', 'Cs', 'Hv',
    'Br', 'Nt', 'Vv', 'Sb', 'Si', 'Pp', 'Mp', 'Pt', 'Rc', 'Eg',
    'Dc', 'Cp', 'Md', 'Nt', 'Nb', 'Le', 'Ps', 'Pv', 'Fv', 'Aa',
    'Ac', 'Ah', 'Al', 'Am', 'Ap', 'Ar', 'Av', 'Bp', 'Bs', 'Ca',
    'Cb', 'Cc', 'Cd', 'Ce', 'Cf', 'Cg', 'Ch', 'Ci', 'Cj', 'Ck',
    'Cl', 'Cm', 'Cn', 'Co', 'Cr', 'Ct', 'Cu', 'Cv', 'Cx', 'Cy',
    'Dl', 'Dm', 'Dn', 'Ds', 'Dt', 'Du', 'Dv', 'Fa', 'Fb', 'Fc',
    'Fd', 'Ff', 'Fg', 'Fh', 'Fi', 'Fj', 'Fl', 'Fn', 'Fo', 'Fp',
    'Gs', 'Gt', 'Ha', 'Hb', 'Hc', 'Hi', 'Hj', 'Hk', 'Hl', 'Hm',
    'Hn', 'Ho', 'Hp', 'Hr', 'Hs', 'Ht', 'Hu', 'Ia', 'Ib', 'Jc',
    'Jn', 'Kf', 'La', 'Lb', 'Lc', 'Ld', 'Lf', 'Lg', 'Lh', 'Li',
    'Lj', 'Lk', 'Ll', 'Lm', 'Ln', 'Lo', 'Lp', 'Lr', 'Ls', 'Lt',
    'Lu', 'Lv', 'Ma', 'Mb', 'Mc', 'Mf', 'Mg', 'Mh', 'Mi', 'Mj',
    'Mk', 'Ml', 'Mm', 'Mn', 'Mo', 'Mq', 'Mr', 'Ms', 'Mu', 'Mv',
    'Mw', 'Mx', 'My', 'Na', 'Nc', 'Nd', 'Ne', 'Nf', 'Ng', 'Nh',
    'Ni', 'Nj', 'Nk', 'Nl', 'Nm', 'Nn', 'No', 'Np', 'Nr', 'Ns',
    'Nu', 'Nv', 'Ob', 'Oc', 'Od', 'Of', 'Og', 'Oh', 'Oi', 'Oj',
    'Ok', 'Ol', 'Om', 'On', 'Oo', 'Op', 'Oq', 'Or', 'Ot', 'Ou',
    'Ov', 'Ow', 'Ox', 'Oy', 'Pa', 'Pb', 'Pc', 'Pd', 'Pe', 'Pf',
    'Pg', 'Ph', 'Pi', 'Pj', 'Pk', 'Pl', 'Pm', 'Pn', 'Po', 'Pr',
    'Pu', 'Pv', 'Pw', 'Px', 'Py', 'Pz', 'Qa', 'Qb', 'Ra', 'Rb',
    'Rf', 'Rg', 'Rh', 'Ri', 'Rj', 'Rk', 'Rl', 'Rm', 'Rn', 'Ro',
    'Rp', 'Rr', 'Rs', 'Rt', 'Ru', 'Rv', 'Sa', 'Sc', 'Sd', 'Se',
    'Sf', 'Sg', 'Sh', 'Sk', 'Sm', 'Sn', 'So', 'Sp', 'Sq', 'Sr',
    'Ss', 'St', 'Su', 'Sv', 'Sw', 'Sx', 'Sy', 'Sz', 'Tc', 'Td',
    'Tf', 'Tg', 'Th', 'Ti', 'Tj', 'Tk', 'Tl', 'Tm', 'Tn', 'To',
    'Tp', 'Tq', 'Tr', 'Ts', 'Tt', 'Tu', 'Tv', 'Tw', 'Tx', 'Ty',
    'Tz', 'Ua', 'Ub', 'Va', 'Vb', 'Vc', 'Vd', 'Ve', 'Vf', 'Vg',
    'Vh', 'Vi', 'Vj', 'Vk', 'Vl', 'Vm', 'Vn', 'Vo', 'Vp', 'Vq',
    'Vr', 'Vs', 'Vt', 'Vu', 'Vw', 'Vx', 'Vy', 'Vz', 'Wa', 'Wb',
    'Wc', 'Wd', 'We', 'Wf', 'Wg', 'Wh', 'Wi', 'Wj', 'Wk', 'Wl',
    'Wm', 'Wn', 'Wp', 'Wq', 'Wr', 'Ws', 'Wt', 'Wu', 'Wv', 'Ww',
    'Wx', 'Wy', 'Wz', 'Xa', 'Xb', 'Xc', 'Xd', 'Xe', 'Xf', 'Xg',
    'Xh', 'Xi', 'Xj', 'Xk', 'Xl', 'Xm', 'Xn', 'Xo', 'Xp', 'Xr',
    'Xs', 'Xt', 'Xu', 'Xv', 'Xw', 'Xx', 'Xy', 'Xz', 'Ya', 'Yb',
    'Yc', 'Yd', 'Ze', 'Zf', 'Zg', 'Zh', 'Zi', 'Zj', 'Zk', 'Zl',
    'Zm', 'Zo', 'Zp', 'Zq', 'Zr', 'Zs', 'Zt', 'Zu', 'Zv', 'Zw',
    # 3-letter prefixes common in plant literature
    'Ath', 'Osl', 'Sol', 'Zma', 'Gma', 'Mtr', 'Hvl', 'Tae',
    'Bdi', 'Sbi', 'Set', 'Ptr', 'Rco', 'Egr', 'Vvi', 'Ppe',
    'Csi', 'Cpa', 'Mes', 'Ppa', 'Mpo', 'Smo',
}

# ── Bucket E: noise detectors ─────────────────────────────────────────────────
RE_NOISE_NUMBER   = re.compile(r'^\d+$')
RE_NOISE_SEQ      = re.compile(r'^\([A-Z]/[A-Z]\)')         # e.g. (C/T)ACGTGTC
RE_NOISE_GREEK    = re.compile(r'^[ωαβγδεζηθιλμνξπρστυφχψΩΑΒΓΔΕΖΗΘΙΛΜΝΞΠΡΣΤΥΦΧΨ]')

NOISE_GENERIC_TERMS = {
    'atp synthase complex', 'electron transport chain', 'electron transport chains',
    'molecular chaperone', 'transaminase', 'protein complex', 'ribosome',
    'photosystem i', 'photosystem ii', 'proteasome', 'nucleosome',
}

# ── Bucket C: bare gene symbols ───────────────────────────────────────────────
# Short all-caps (or mixed-caps) gene symbol, no spaces, no species prefix
RE_BARE_SYMBOL = re.compile(r'^[A-Z][A-Z0-9]{1,9}(\.\d+)?$')
RE_BARE_MIXED  = re.compile(r'^[A-Z][A-Za-z0-9]{1,9}(\d)$')


def classify(entity: str):
    """Return (bucket, sub_type, species_prefix, gene_symbol, notes)."""
    raw = entity.strip()
    s   = raw.strip('"').strip()   # strip surrounding quotes if any

    if not s:
        return ('E', 'noise_empty', '', '', '')

    # ── Bucket A ──────────────────────────────────────────────────────────────
    for pat, sub in BUCKET_A_CHECKS:
        if pat.match(s):
            return ('A', sub, '', '', '')

    # ── Bucket B ──────────────────────────────────────────────────────────────
    m = RE_SPECIES_PREFIX.match(s)
    if m:
        prefix = m.group(1)
        symbol = m.group(2)
        # Require at least 4 chars total and symbol has at least one digit or ≥4 chars
        # to avoid matching random CamelCase descriptive names
        if len(s) >= 5 and (any(c.isdigit() for c in symbol) or len(symbol) >= 3):
            return ('B', f'species_prefix_{prefix}', prefix, symbol, '')

    # ── Bucket E (noise) ──────────────────────────────────────────────────────
    if RE_NOISE_NUMBER.match(s):
        return ('E', 'noise_number', '', '', '')
    if RE_NOISE_SEQ.match(s):
        return ('E', 'noise_sequence_motif', '', '', '')
    if RE_NOISE_GREEK.match(s):
        return ('E', 'noise_greek_symbol', '', '', '')
    if s.lower() in NOISE_GENERIC_TERMS:
        return ('E', 'noise_generic_term', '', '', '')
    # Single letters or very short (<= 2 char)
    if len(s) <= 2:
        return ('E', 'noise_too_short', '', '', '')

    # ── Bucket C ──────────────────────────────────────────────────────────────
    if RE_BARE_SYMBOL.match(s) and ' ' not in s:
        return ('C', 'bare_gene_symbol', '', '', '')
    if RE_BARE_MIXED.match(s) and ' ' not in s:
        return ('C', 'bare_gene_symbol', '', '', '')

    # ── Bucket D ──────────────────────────────────────────────────────────────
    return ('D', 'descriptive_name', '', '', '')


# ── Main ──────────────────────────────────────────────────────────────────────
entries = []
with open(INPUT_CSV) as fh:
    reader = csv.DictReader(fh)
    for row in reader:
        entries.append(row['entity'])

print(f"Loaded {len(entries)} entries")

results = []
for entity in entries:
    bucket, sub_type, prefix, symbol, notes = classify(entity)
    results.append({
        'entity':                   entity,
        'bucket':                   bucket,
        'sub_type':                 sub_type,
        'extracted_species_prefix': prefix,
        'extracted_gene_symbol':    symbol,
        'notes':                    notes,
    })

OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
with open(OUTPUT_CSV, 'w', newline='') as fh:
    writer = csv.DictWriter(
        fh,
        fieldnames=['entity','bucket','sub_type','extracted_species_prefix','extracted_gene_symbol','notes']
    )
    writer.writeheader()
    writer.writerows(results)

print(f"Saved: {OUTPUT_CSV}")

# ── Report ────────────────────────────────────────────────────────────────────
from collections import Counter, defaultdict

bucket_counts = Counter(r['bucket'] for r in results)
bucket_entries = defaultdict(list)
for r in results:
    bucket_entries[r['bucket']].append(r)

print(f"\n=== REPORT ===")
print(f"Total entries classified: {len(results)}")
print()
for b in ['A','B','C','D','E']:
    count = bucket_counts[b]
    pct   = 100 * count / len(results)
    print(f"  Bucket {b}: {count:6,}  ({pct:.1f}%)")

for b in ['A','B','C','D','E']:
    bucket_name = {
        'A': 'Clean database IDs',
        'B': 'Species-prefixed symbols',
        'C': 'Bare gene symbols',
        'D': 'Descriptive names',
        'E': 'Noise / unmappable',
    }[b]
    print(f"\n--- Bucket {b}: {bucket_name} — 5 samples ---")
    for r in bucket_entries[b][:5]:
        print(f"  [{r['sub_type']:30}]  {r['entity'][:80]}")

# Full Bucket E listing
print(f"\n--- Bucket E: ALL entries ({bucket_counts['E']} total) ---")
e_entries = [r['entity'] for r in bucket_entries['E']]
for i, e in enumerate(e_entries[:50]):
    print(f"  {i+1:3}. {e}")
if len(e_entries) > 50:
    print(f"  ... ({len(e_entries) - 50} more)")

# Sub-type breakdown for Bucket A
print(f"\n--- Bucket A sub-types ---")
a_subtypes = Counter(r['sub_type'] for r in bucket_entries['A'])
for sub, cnt in a_subtypes.most_common():
    print(f"  {sub:30}  {cnt:6,}")

# Species prefix breakdown for Bucket B (top 20)
print(f"\n--- Bucket B top species prefixes ---")
b_prefixes = Counter(r['extracted_species_prefix'] for r in bucket_entries['B'])
for pfx, cnt in b_prefixes.most_common(20):
    print(f"  {pfx:5}  {cnt:6,}")
