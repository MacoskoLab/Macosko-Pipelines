import os
import gzip
import shutil
import argparse
import numpy as np
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser(description='Subset a diffusion matrix to the beads listed in a Puck csv')
    parser.add_argument("-i", "--in_dir", help="folder with matrix.csv.gz, sb1.txt.gz, sb2.txt.gz", type=str, default=".")
    parser.add_argument("-o", "--out_dir", help="output folder", type=str, default=None)
    parser.add_argument("-p", "--puck", help="Puck csv whose barcodes to keep (headerless sb,x,y)", type=str, default=None)
    parser.add_argument("-b", "--bead", help="which side the puck barcodes index", type=int, default=2)
    parser.add_argument("-c", "--chunksize", help="matrix rows per chunk", type=int, default=25_000_000)
    parser.add_argument("-l", "--compresslevel", help="gzip level for the outputs", type=int, default=6)

    args, unknown = parser.parse_known_args()
    [print(f"WARNING: unknown command-line argument {u}") for u in unknown]
    return args

# Load arguments
args = get_args()
in_dir = args.in_dir                 ; print(f"input directory = {in_dir}")
out_dir = args.out_dir               ; print(f"output directory = {out_dir}")
puck = args.puck                     ; print(f"puck = {puck}")
bead = args.bead                     ; print(f"bead = {bead}")
chunksize = args.chunksize           ; print(f"chunksize = {chunksize}")
compresslevel = args.compresslevel   ; print(f"compresslevel = {compresslevel}")

assert bead in [1, 2], "bead must be 1 or 2"
assert out_dir is not None, "--out_dir is required"
assert puck is not None, "--puck is required"
assert os.path.isfile(puck), f"{puck} not found"
for f in ['matrix.csv.gz', 'sb1.txt.gz', 'sb2.txt.gz']:
    assert os.path.isfile(os.path.join(in_dir, f)), f"{f} not found in {in_dir}"
assert os.path.abspath(in_dir) != os.path.abspath(out_dir), "out_dir must differ from in_dir"
os.makedirs(out_dir, exist_ok=True)
assert os.path.exists(out_dir)

sel = bead        # the side the puck barcodes index
oth = 3 - bead    # the other side, kept only where links survive
matrix_path = os.path.join(in_dir, 'matrix.csv.gz')

# Read the barcode lists - line i (1-based) is the barcode named by index i
def read_sb(bead):
    with gzip.open(os.path.join(in_dir, f"sb{bead}.txt.gz"), 'rt') as f:
        return f.read().split()

sb = {b: read_sb(b) for b in (1, 2)}
n = {b: len(sb[b]) for b in (1, 2)}
for b in (1, 2):
    assert len(set(sb[b])) == n[b], f"sb{b}.txt.gz has duplicate barcodes"
    print(f"sb{b}.txt.gz: {n[b]:,} barcodes")

# Read the barcodes to keep out of the puck
with open(puck) as f:
    puck_sb = [l.split(',', 1)[0] for l in f if l.strip()]
assert len(set(puck_sb)) == len(puck_sb), f"{puck} has duplicate barcodes"
missing = set(puck_sb) - set(sb[sel])
assert not missing, f"{len(missing)} puck barcodes are absent from sb{sel}.txt.gz (wrong puck or wrong --bead?)"
print(f"{puck}: {len(puck_sb):,} barcodes to keep")

want = set(puck_sb)
keep = {sel: np.fromiter((s in want for s in sb[sel]), dtype=bool, count=n[sel])}
assert keep[sel].sum() == len(puck_sb)

col = {b: f"sb{b}_index" for b in (1, 2)}

# Pass 1: find which rows survive, and which barcodes on the other side keep a link
print("Pass 1/2: scanning the matrix...")
seen_oth = np.zeros(n[oth] + 1, dtype=bool)
total_rows = kept_rows = 0
for chunk in pd.read_csv(matrix_path, chunksize=chunksize, dtype=np.int32):
    total_rows += len(chunk)
    i_sel = chunk[col[sel]].to_numpy()
    i_oth = chunk[col[oth]].to_numpy()
    assert i_sel.min() >= 1 and i_sel.max() <= n[sel], f"{col[sel]} out of range 1..{n[sel]}"
    assert i_oth.min() >= 1 and i_oth.max() <= n[oth], f"{col[oth]} out of range 1..{n[oth]}"
    m = keep[sel][i_sel - 1]
    kept_rows += int(m.sum())
    seen_oth[i_oth[m]] = True
    print(f"  {total_rows:,} rows scanned, {kept_rows:,} kept")
keep[oth] = seen_oth[1:]

n_new = {b: int(keep[b].sum()) for b in (1, 2)}
print(f"connections: {kept_rows:,} / {total_rows:,} kept ({100*kept_rows/total_rows:.1f}%)")
for b in (1, 2):
    print(f"sb{b}: {n_new[b]:,} / {n[b]:,} kept ({n[b]-n_new[b]:,} dropped)")
assert kept_rows > 0, "the selection kept no connections"

# Renumber both axes to a gapless 1..N. knn.py infers the sparse matrix shape from
# max(index) and never reads the barcode files, so a gap would misalign silently.
def build_remap(keep):
    remap = np.zeros(len(keep) + 1, dtype=np.int32)
    remap[1:][keep] = np.arange(1, int(keep.sum()) + 1, dtype=np.int32)
    return remap

remap = {b: build_remap(keep[b]) for b in (1, 2)}

# Pass 2: write the subset matrix. The input is sorted by (sb1_index, sb2_index) and
# the remaps are monotonic, so streaming it keeps that ordering.
out_matrix = os.path.join(out_dir, 'matrix.csv.gz')
print(f"Pass 2/2: writing {out_matrix}...")
written = 0
with gzip.open(out_matrix, 'wt', compresslevel=compresslevel) as f:
    for chunk in pd.read_csv(matrix_path, chunksize=chunksize, dtype=np.int32):
        i_sel = chunk[col[sel]].to_numpy()
        m = keep[sel][i_sel - 1]
        new_sel = remap[sel][i_sel[m]]
        new_oth = remap[oth][chunk[col[oth]].to_numpy()[m]]
        assert new_sel.all() and new_oth.all(), "a kept row maps to a dropped barcode"
        out = pd.DataFrame({col[sel]: new_sel, col[oth]: new_oth,
                            'umi': chunk['umi'].to_numpy()[m]})
        out.to_csv(f, index=False, header=(written == 0),
                   columns=['sb1_index', 'sb2_index', 'umi'], lineterminator='\n')
        written += len(out)
        print(f"  {written:,} / {kept_rows:,} rows written")
assert written == kept_rows, f"wrote {written} rows, expected {kept_rows}"

# Write the barcode lists in the new index order (headerless, one per line)
for b in (1, 2):
    path = os.path.join(out_dir, f"sb{b}.txt.gz")
    with gzip.open(path, 'wt', compresslevel=compresslevel) as f:
        f.write('\n'.join(s for s, k in zip(sb[b], keep[b]) if k) + '\n')
    print(f"wrote {path}: {n_new[b]:,} barcodes")

# helpers.py estimate_diameter() bins on the total barcode count, so a large enough
# cut can silently change the puck scale recon.py applies
def diameter_bin(total):
    for cutoff, micron in [(17e6, 70_000), (8.43e6, 40_000), (4.22e6, 30_000),
                           (1.63e6, 20_000), (0.3e6, 12_000)]:
        if total > cutoff:
            return micron
    return None

d_old, d_new = diameter_bin(n[1] + n[2]), diameter_bin(n_new[1] + n_new[2])

# Copy metadata.csv verbatim (terra/reconstruction.wdl greps *_barcodes_manual out of
# it for its cache check) and append the subset counts as new headerless key,value rows
meta_in = os.path.join(in_dir, 'metadata.csv')
if os.path.isfile(meta_in):
    meta_out = os.path.join(out_dir, 'metadata.csv')
    shutil.copyfile(meta_in, meta_out)
    with open(meta_in, 'rb') as f:
        needs_newline = os.path.getsize(meta_in) > 0 and f.read()[-1:] != b'\n'
    with open(meta_out, 'a') as f:
        if needs_newline:
            f.write('\n')
        for k, v in [('subset_puck', os.path.abspath(puck)),
                     ('subset_bead', bead),
                     ('subset_sb1', n_new[1]),
                     ('subset_sb2', n_new[2]),
                     ('subset_sb1_dropped', n[1] - n_new[1]),
                     ('subset_sb2_dropped', n[2] - n_new[2]),
                     ('subset_connections', kept_rows),
                     ('subset_connections_dropped', total_rows - kept_rows),
                     # terra/reconstruction.wdl reads these to auto-pass recon.py -D
                     ('subset_diameter_prev', d_old),
                     ('subset_diameter_new', d_new)]:
            f.write(f"{k},{v}\n")
    print(f"wrote {meta_out}: copied + subset_* keys")
else:
    print(f"WARNING: no metadata.csv in {in_dir}, not written")

# Warn about anything in in_dir that is now stale rather than silently copying it
for f in ['knn1.npz', 'knn2.npz', 'readumi_per_sb1.csv.gz', 'readumi_per_sb2.csv.gz']:
    if os.path.isfile(os.path.join(in_dir, f)):
        print(f"NOTE: not copying {f} - it is indexed by the old barcode set")

if d_new is None:
    print(f"WARNING: only {n_new[1]+n_new[2]:,} barcodes remain, below the smallest bin "
          f"estimate_diameter() knows; recon.py will exit unless you pass -D")
elif d_old == d_new:
    print(f"estimated diameter unchanged: {d_new} micron")
else:
    print(f"WARNING: estimated diameter changes {d_old} -> {d_new} micron; "
          f"pass -D {d_old} to recon.py to keep the original scale")

print(f"\nDone. Next:\n"
      f"  python knn.py -i {out_dir} -o {out_dir} -b {bead} -k 2\n"
      f"  python recon.py -i {out_dir} -o {out_dir} -b {bead}"
      + ("" if d_old == d_new else f" -D {d_old}"))
