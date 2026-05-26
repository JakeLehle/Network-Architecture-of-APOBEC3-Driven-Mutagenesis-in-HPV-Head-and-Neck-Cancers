import scanpy as sc
import matplotlib
matplotlib.use('Agg')

adata = sc.read_h5ad("/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_4/00_input/adata_final.h5ad")

print("=== SHAPE ===")
print(adata.shape)

print("\n=== .obs columns ===")
for col in adata.obs.columns:
    dtype = adata.obs[col].dtype
    nuniq = adata.obs[col].nunique()
    example = adata.obs[col].dropna().iloc[0] if len(adata.obs[col].dropna()) > 0 else 'NA'
    print(f"  {col:40s}  dtype={str(dtype):15s}  nunique={nuniq:6d}  example={example}")

print("\n=== .var columns ===")
for col in adata.var.columns:
    print(f"  {col}")

print("\n=== .obsm keys ===")
for k in adata.obsm.keys():
    print(f"  {k}  shape={adata.obsm[k].shape}")

print("\n=== .layers keys ===")
for k in adata.layers.keys():
    print(f"  {k}")

print("\n=== .uns keys (top level) ===")
for k in sorted(adata.uns.keys()):
    print(f"  {k}")
