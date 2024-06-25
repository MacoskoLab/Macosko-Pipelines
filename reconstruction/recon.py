import os
import sys
import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
from umap import UMAP
import umap.plot
from umap.umap_ import nearest_neighbors
# os.chdir("/home/nsachdev/reconstruction/1.2")
assert all(os.path.isfile(file) for file in ['matrix.csv.gz', 'sb1.txt.gz', 'sb2.txt.gz'])

df = pd.read_csv('matrix.csv.gz', compression='gzip', header=None, names=['sb1', 'sb2', 'umi'])
df.sb1 -= 1 # convert from 1- to 0-indexed
df.sb2 -= 1 # convert from 1- to 0-indexed

with gzip.open('sb1.txt.gz', 'rt') as f:
    sb1 = [line.strip() for line in f.readlines()]
with gzip.open('sb2.txt.gz', 'rt') as f:
    sb2 = [line.strip() for line in f.readlines()]

assert sorted(list(set(df.sb1))) == list(range(len(sb1)))
assert sorted(list(set(df.sb2))) == list(range(len(sb2)))

# puck = pd.read_csv('V15A_recon_loc_30000_80_0.2_fbmax1000.csv')

print(f"{len(sb1)} R1 barcodes")
print(f"{len(sb2)} R2 barcodes")
#print(pd.Series([bc[0:14] for bc in sb1]).isin(puck['V15A']).sum())
#print(pd.Series([bc[0:14] for bc in sb2]).isin(puck['V15A']).sum())
#print(f"{puck.shape[0]} puck barcodes")

var_name = globals()["os"]
for var_name in globals():
    if var_name[0] != "_":
        print(f"Variable '{var_name}': {sys.getsizeof(globals()[var_name])} bytes")



fig, axs = plt.subplots(3, 1, figsize=(7, 8))

# Plot 1: Histogram of UMI (Log Scale)
axs[0].hist(np.log10(df['umi']), bins=100, edgecolor='black')
axs[0].set_title('Histogram of UMI (Log Scale)')
axs[0].set_xlabel('log10 UMI')
axs[0].set_ylabel('Count')
axs[0].grid(True)

# Plot 2: Histogram of Summarized UMI by SB1
gdf_sb1 = df.groupby('sb1')['umi'].sum().reset_index()
meanumi_sb1 = np.log10(gdf_sb1['umi'].mean())
axs[1].hist(np.log10(gdf_sb1['umi']), bins=30, edgecolor='black')
axs[1].axvline(meanumi_sb1, color='red', linestyle='dashed', linewidth=2)
axs[1].text(meanumi_sb1, axs[1].get_ylim()[1]*0.9, f'Mean: {10**meanumi_sb1:.1f}', color='red', fontsize=12, ha='center')
axs[1].set_title('Histogram of Summarized UMI by SB1')
axs[1].set_xlabel('log10 UMI')
axs[1].set_ylabel('Frequency')
axs[1].grid(True)

# Plot 3: Histogram of Summarized UMI by SB2
gdf_sb2 = df.groupby('sb2')['umi'].sum().reset_index()
meanumi_sb2 = np.log10(gdf_sb2['umi'].mean())
axs[2].hist(np.log10(gdf_sb2['umi']), bins=30, edgecolor='black')
axs[2].axvline(meanumi_sb2, color='red', linestyle='dashed', linewidth=2)
axs[2].text(meanumi_sb2, axs[2].get_ylim()[1]*0.9, f'Mean: {10**meanumi_sb2:.1f}', color='red', fontsize=12, ha='center')
axs[2].set_title('Histogram of Summarized UMI by SB2')
axs[2].set_xlabel('log10 UMI')
axs[2].set_ylabel('Frequency')
axs[2].grid(True)

# Adjust layout to prevent overlap
plt.tight_layout()

# Save the plot as a PDF
plt.savefig('umi_histograms.pdf')


mat = coo_matrix((df['umi'], (df['sb1'], df['sb2'])))


knn = nearest_neighbors(mat,
                        n_neighbors=25,
                        metric="cosine",
                        metric_kwds=None,
                        angular=False,
                        random_state=None,
                        low_memory=False,
                        verbose=True
                       )

print(knn)
