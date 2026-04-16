#!/usr/bin/env python3
"""
Step05_Revised_HNSC_A3_vs_SBS2.py
====================================
Figure 1: three panels + supplemental.

Panel 1a: A3A+A3B vs SBS2 scatter (necessary but not sufficient)
Panel 1b: A3A vs A3B colored by SBS2 (visual gradient, depth-sorted)
Panel 1c: Box-and-whisker of SBS2 across 4 A3A/A3B quadrants + 2x2 heatmap
Supplemental: Zoomed view of low-A3 / high-SBS2 boundary region

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os, sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from scipy.stats import spearmanr, mannwhitneyu

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
SHARED_VCF = "/master/jlehle/SHARED/TCGA/VCF"

EXPRESSION_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_1", "TCGA_master_FPKM_UQ.tsv")
METADATA_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_1", "TCGA_sample_metadata_final.tsv")
ORIGINAL_MUT_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_1", "Mutation_Table_Tumors_TCGA.tsv")
NEW_COUNTS_PATH = os.path.join(SHARED_VCF, "SigProfiler_output", "TCGA_SBS_signature_counts.tsv")
MANIFEST_PATH = os.path.join(SHARED_VCF, "manifests", "TCGA_MuTect2_master_manifest.tsv")

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "FIG_1")
PANEL_DIR = os.path.join(OUTPUT_DIR, "FIGURE_1_PANELS")
TROUBLE_DIR = os.path.join(OUTPUT_DIR, "TROUBLESHOOTING")
os.makedirs(PANEL_DIR, exist_ok=True); os.makedirs(TROUBLE_DIR, exist_ok=True)

COLOR_SBS2_HIGH = "#ed6a5a"; COLOR_CREAM = "#f4f1bb"; COLOR_NORMAL = "#9bc1bc"
COLOR_DARK_GRAY = "#4D4D4D"; COLOR_BLACK = "#000000"
FONT_SIZE = 28
A3_ALL = ['APOBEC3A','APOBEC3B','APOBEC3C','APOBEC3D','APOBEC3F','APOBEC3G','APOBEC3H']

# =============================================================================
# HELPERS
# =============================================================================
report_lines = []
def log(msg=""):
    print(msg, flush=True); report_lines.append(str(msg))
def banner(title, char="="):
    log(""); log(char*80); log(f"  {title}"); log(char*80)
def save_fig(fig, name):
    for ext in ['pdf','png']:
        fig.savefig(os.path.join(PANEL_DIR, f"{name}.{ext}"), dpi=300, bbox_inches='tight')
    plt.close(fig); log(f"  Saved: {name}.pdf/.png")

def find_axis_break(values, n_check=10, min_gap_ratio=1.5):
    sv = np.sort(values)[::-1]
    if len(sv) < 3: return 0, 0, False
    best_r, best_i = 1.0, -1
    for i in range(min(n_check, len(sv)-1)):
        r = sv[i]/sv[i+1] if sv[i+1]>0 else 1.0
        if r > best_r: best_r, best_i = r, i
    if best_r > min_gap_ratio and best_i >= 0:
        bs = sv[best_i+1]*1.10; be = sv[best_i]*0.97
        log(f"    Axis break: ratio={best_r:.2f}, {bs:.1f}-{be:.1f}")
        return bs, be, True
    return 0, 0, False

def add_inner_axes(ax, x_lo, x_hi, y_lo, y_hi, y_color=COLOR_BLACK, lw=3, zorder=1):
    if x_lo <= 0 <= x_hi:
        ax.plot([0,0],[max(0,y_lo),y_hi], color=y_color, lw=lw, zorder=zorder, solid_capstyle='butt')
    ax.plot([max(0,x_lo),x_hi],[0,0], color=COLOR_BLACK, lw=lw, zorder=zorder, solid_capstyle='butt')

def add_break_markers(axL, axR, d=0.015):
    for kw in [dict(transform=axL.transAxes)]:
        axL.plot((1-d,1+d),(-d,+d), color=COLOR_BLACK, clip_on=False, lw=1.5, **kw)
        axL.plot((1-d,1+d),(1-d,1+d), color=COLOR_BLACK, clip_on=False, lw=1.5, **kw)
    for kw in [dict(transform=axR.transAxes)]:
        axR.plot((-d,+d),(-d,+d), color=COLOR_BLACK, clip_on=False, lw=1.5, **kw)
        axR.plot((-d,+d),(1-d,1+d), color=COLOR_BLACK, clip_on=False, lw=1.5, **kw)
    axR.spines['left'].set_visible(False); axR.tick_params(left=False)
    axL.spines['right'].set_visible(False)

def clean_spines(ax):
    for s in ax.spines.values(): s.set_visible(False)

# =============================================================================
# STEPS 1-5: DATA LOADING AND MATCHING
# =============================================================================
banner("STEP 1: Load expression matrix")
raw = pd.read_csv(EXPRESSION_PATH, sep='\t', header=None, low_memory=False)
log(f"  Raw shape: {raw.shape}")
gene_symbols = raw.iloc[1].values; data = raw.iloc[3:].copy()
data.columns = range(len(data.columns)); entity_ids = data[4].astype(str).values
a3_col_indices = {}
for gene in A3_ALL:
    m = np.where(gene_symbols[5:] == gene)[0]
    if len(m) > 0: a3_col_indices[gene] = m[0]+5; log(f"  Found {gene} at col {m[0]+5}")
a3_expr = pd.DataFrame({'Entity_ID': entity_ids})
for g, c in a3_col_indices.items():
    a3_expr[g] = pd.to_numeric(data[c].values, errors='coerce')

banner("STEP 2: Filter to HNSC tumors")
metadata = pd.read_csv(METADATA_PATH, sep='\t')
a3_expr['Project_ID'] = a3_expr['Entity_ID'].map(metadata.set_index('Entity_ID')['Project_ID'].to_dict())
a3_hnsc = a3_expr[a3_expr['Project_ID']=='TCGA-HNSC'].copy()
a3_hnsc = a3_hnsc[a3_hnsc['Entity_ID'].str[13:15].isin([f'{i:02d}' for i in range(1,10)])].copy()
a3_hnsc['rna_sample_type'] = a3_hnsc['Entity_ID'].str[13:15]
log(f"  HNSC tumors: {len(a3_hnsc)}")

banner("STEP 3: Load SigProfiler counts")
raw_ct = pd.read_csv(NEW_COUNTS_PATH, sep='\t'); fc = raw_ct.columns[0]
if any(raw_ct[fc].astype(str).str.startswith('SBS')):
    raw_ct = raw_ct.set_index(fc); new_ct = raw_ct.T.copy()
    new_ct.index.name='WES_Barcode'; new_ct = new_ct.reset_index()
    new_ct = new_ct[~new_ct['WES_Barcode'].str.contains('Cancer_Type|^$',na=True,regex=True)].copy()
    for c in [x for x in new_ct.columns if x.startswith('SBS')]:
        new_ct[c] = pd.to_numeric(new_ct[c], errors='coerce')
else: new_ct = raw_ct.rename(columns={fc:'WES_Barcode'})
manifest = pd.read_csv(MANIFEST_PATH, sep='\t')
hwes = set(manifest[manifest['Cancer_Type']=='HNSC']['Entity_ID'].astype(str))
nch = new_ct[new_ct['WES_Barcode'].isin(hwes)].copy()
nch = nch[nch['WES_Barcode'].str[13:15].isin([f'{i:02d}' for i in range(1,10)])].copy()
nch['wes_sample_type']=nch['WES_Barcode'].str[13:15]; nch['Case_ID']=nch['WES_Barcode'].str[:12]
log(f"  HNSC WES tumors: {len(nch)}")

banner("STEP 4: Barcode matching")
om = pd.read_csv(ORIGINAL_MUT_PATH, sep='\t', usecols=[
    'TCGA_Gene_Expression_Entity_ID','Mutation_Signature__File_Orginal_Entity_ID'])
cw = dict(zip(om['TCGA_Gene_Expression_Entity_ID'],om['Mutation_Signature__File_Orginal_Entity_ID']))
rids = set(a3_hnsc['Entity_ID'].values); wset = set(nch['WES_Barcode'].values)
dm = [{'RNA_Barcode':r,'WES_Barcode':cw[r],'match_source':'DIRECT'}
      for r in rids if r in cw and cw[r] in wset]
ddf = pd.DataFrame(dm); dp = set(d['RNA_Barcode'][:12] for d in dm)
log(f"  Direct: {len(ddf)}")
am = set(ddf['RNA_Barcode'].values) if len(ddf)>0 else set()
um = a3_hnsc[~a3_hnsc['Entity_ID'].isin(am)].copy()
cm = []
for _,r in um.iterrows():
    rb,rc,rs = r['Entity_ID'],r['Entity_ID'][:12],r['rna_sample_type']
    if rc in dp: continue
    cd = nch[(nch['Case_ID']==rc)&(nch['wes_sample_type']==rs)]
    if len(cd)>=1: cm.append({'RNA_Barcode':rb,'WES_Barcode':cd.iloc[0]['WES_Barcode'],'match_source':'CASE_ID'})
cdf = pd.DataFrame(cm); log(f"  Case_ID: {len(cdf)}")
amdf = pd.concat([ddf,cdf],ignore_index=True)
amdf['patient']=amdf['RNA_Barcode'].str[:12]
amdf['pri']=amdf['match_source'].map({'DIRECT':0,'CASE_ID':1}).fillna(99)
amdf = amdf.sort_values('pri').drop_duplicates(subset='patient',keep='first').drop(columns=['pri','patient'])
log(f"  Total: {len(amdf)}")

banner("STEP 5: Build final table")
amdf['SBS2'] = amdf['WES_Barcode'].map(nch.set_index('WES_Barcode')['SBS2'].to_dict())
final = a3_hnsc.merge(amdf, left_on='Entity_ID', right_on='RNA_Barcode', how='inner')
final['A3A']=final['APOBEC3A']; final['A3B']=final['APOBEC3B']
final['A3A_plus_A3B']=final['A3A']+final['A3B']
median_sbs2 = final['SBS2'].median()
final['outlier_low_A3_high_SBS2'] = (final['SBS2']>median_sbs2)&(final['A3A_plus_A3B']<1.0)

med_a3a = final['A3A'].median(); med_a3b = final['A3B'].median()
final['A3A_status'] = np.where(final['A3A']>=med_a3a,'HIGH','LOW')
final['A3B_status'] = np.where(final['A3B']>=med_a3b,'HIGH','LOW')
final['quadrant'] = final['A3B_status']+'_A3B_'+final['A3A_status']+'_A3A'
log(f"  Final: {len(final)} tumors, median SBS2={median_sbs2:.0f}")

x_all = final['A3A_plus_A3B'].values; y_all = final['SBS2'].values
lm = (x_all<=10)&(x_all>0)
safe_slope = np.max(y_all[lm]/x_all[lm])*1.05 if lm.sum()>0 else np.max(y_all)/10
X_THR = median_sbs2/safe_slope if safe_slope>0 else 10; y_thr = median_sbs2
regions = []
for xi,yi in zip(x_all,y_all):
    if xi<=X_THR and yi>safe_slope*max(xi,0.001): regions.append('teal')
    elif yi>y_thr: regions.append('coral')
    else: regions.append('cream')
final['region'] = regions
n_teal=regions.count('teal'); n_coral=regions.count('coral'); n_cream=regions.count('cream')
log(f"  Regions: teal={n_teal}, coral={n_coral}, cream={n_cream}")
final.to_csv(os.path.join(OUTPUT_DIR, "HNSC_A3_SBS2_matched_v2.tsv"), sep='\t', index=False)

# =============================================================================
# FIGURE SETUP
# =============================================================================
plt.rcParams.update({'font.size':FONT_SIZE,'axes.titlesize':FONT_SIZE,'axes.labelsize':FONT_SIZE,
                     'xtick.labelsize':FONT_SIZE-6,'ytick.labelsize':FONT_SIZE-6,'legend.fontsize':FONT_SIZE-8})
rcol = {'teal':COLOR_NORMAL,'coral':COLOR_SBS2_HIGH,'cream':COLOR_CREAM}
pcols = [rcol[r] for r in regions]
x_max=x_all.max()*1.05; y_max=y_all.max()*1.05

# =============================================================================
# PANEL 1a
# =============================================================================
banner("Panel 1a")
bs,be,hb = find_axis_break(x_all)
def draw_1a(ax, xl, xh, is_main=True):
    ax.fill([0,0,X_THR,X_THR],[0,y_max,y_max,y_thr],color=COLOR_NORMAL,alpha=0.5,zorder=0)
    ax.fill([X_THR,X_THR,x_max*2,x_max*2],[y_thr,y_max,y_max,y_thr],color=COLOR_SBS2_HIGH,alpha=0.5,zorder=0)
    ax.fill([0,X_THR,X_THR,x_max*2,x_max*2,0],[0,y_thr,y_thr,y_thr,0,0],color=COLOR_CREAM,alpha=0.5,zorder=0)
    add_inner_axes(ax,xl,xh,0,y_max,y_color=COLOR_NORMAL if is_main else COLOR_BLACK,lw=3,zorder=1)
    ax.plot([0,X_THR],[0,y_thr],'--',color=COLOR_DARK_GRAY,lw=1,alpha=0.7,zorder=2)
    ax.plot([X_THR,x_max*2],[y_thr,y_thr],'--',color=COLOR_DARK_GRAY,lw=1,alpha=0.7,zorder=2)
    ax.scatter(x_all,y_all,c=pcols,s=50,alpha=0.7,edgecolors=COLOR_BLACK,linewidths=0.4,rasterized=True,zorder=3)
    ax.set_xlim(xl,xh); ax.set_ylim(-y_max*0.02,y_max)

if hb:
    fig,(axM,axB) = plt.subplots(1,2,sharey=True,figsize=(16,10),gridspec_kw={'width_ratios':[3,1],'wspace':0.04})
    draw_1a(axM,-x_max*0.02,bs,True); draw_1a(axB,be,x_max,False)
    clean_spines(axM); clean_spines(axB); add_break_markers(axM,axB)
    axM.set_ylabel('SBS2 Weight (mutation count)',fontsize=FONT_SIZE)
    fig.text(0.5,0.01,'A3A + A3B Expression (FPKM-UQ)',ha='center',fontsize=FONT_SIZE)
    leg = [Patch(facecolor=COLOR_SBS2_HIGH,edgecolor=COLOR_BLACK,alpha=0.7,label=f'A3 + SBS2 HIGH (n={n_coral})'),
           Patch(facecolor=COLOR_CREAM,edgecolor=COLOR_BLACK,alpha=0.7,label=f'A3 + SBS2 LOW (n={n_cream})'),
           Patch(facecolor=COLOR_NORMAL,edgecolor=COLOR_BLACK,alpha=0.7,label=f'No A3 + High SBS2 (n={n_teal})')]
    axM.legend(handles=leg,loc='upper right',framealpha=0.9)
    plt.tight_layout(); save_fig(fig,"Panel_1a_A3sum_vs_SBS2")
else:
    fig,ax = plt.subplots(figsize=(14,10)); draw_1a(ax,-x_max*0.02,x_max,True); clean_spines(ax)
    ax.set_xlabel('A3A + A3B Expression (FPKM-UQ)',fontsize=FONT_SIZE)
    ax.set_ylabel('SBS2 Weight (mutation count)',fontsize=FONT_SIZE)
    leg = [Patch(facecolor=COLOR_SBS2_HIGH,edgecolor=COLOR_BLACK,alpha=0.7,label=f'A3 + SBS2 HIGH (n={n_coral})'),
           Patch(facecolor=COLOR_CREAM,edgecolor=COLOR_BLACK,alpha=0.7,label=f'A3 + SBS2 LOW (n={n_cream})'),
           Patch(facecolor=COLOR_NORMAL,edgecolor=COLOR_BLACK,alpha=0.7,label=f'No A3 + High SBS2 (n={n_teal})')]
    ax.legend(handles=leg,loc='upper right',framealpha=0.9)
    plt.tight_layout(); save_fig(fig,"Panel_1a_A3sum_vs_SBS2")

# =============================================================================
# PANEL 1b
# =============================================================================
banner("Panel 1b")
pdf = final.sort_values('SBS2',ascending=True).copy()
x2=pdf['A3A'].values; y2=pdf['A3B'].values; c2=pdf['SBS2'].values; c2l=np.log1p(c2)
y2m=y2.max()*1.05; bs2,be2,hb2 = find_axis_break(x2); x2m=x2.max()*1.05
def draw_1b(ax,xl,xh,is_main=True):
    add_inner_axes(ax,xl,xh,0,y2m,y_color=COLOR_BLACK,lw=3,zorder=0)
    sc = ax.scatter(x2,y2,c=c2l,cmap='magma',s=50,alpha=0.7,edgecolors=COLOR_BLACK,linewidths=0.3,rasterized=True,zorder=2)
    ax.set_xlim(xl,xh); ax.set_ylim(-y2m*0.02,y2m); return sc
if hb2:
    fig,(axM,axB) = plt.subplots(1,2,sharey=True,figsize=(16,10),gridspec_kw={'width_ratios':[3,1],'wspace':0.04})
    draw_1b(axM,-x2m*0.02,bs2,True); sc=draw_1b(axB,be2,x2m,False)
    clean_spines(axM); clean_spines(axB); add_break_markers(axM,axB)
    cb = plt.colorbar(sc,ax=axB,shrink=0.8,pad=0.02)
    ct = [t for t in [0,1,5,20,50,100,200,500,700] if t<=c2.max()]
    cb.set_ticks([np.log1p(t) for t in ct]); cb.set_ticklabels([str(t) for t in ct])
    cb.set_label('SBS2 Weight',fontsize=FONT_SIZE-6); cb.ax.tick_params(labelsize=FONT_SIZE-8)
    axM.set_ylabel('APOBEC3B Expression (FPKM-UQ)',fontsize=FONT_SIZE)
    fig.text(0.5,0.01,'APOBEC3A Expression (FPKM-UQ)',ha='center',fontsize=FONT_SIZE)
    plt.tight_layout(); save_fig(fig,"Panel_1b_A3A_vs_A3B_SBS2")
else:
    fig,ax = plt.subplots(figsize=(14,10)); sc=draw_1b(ax,-x2m*0.02,x2m,True); clean_spines(ax)
    cb = plt.colorbar(sc,ax=ax,shrink=0.8,pad=0.02)
    ct = [t for t in [0,1,5,20,50,100,200,500,700] if t<=c2.max()]
    cb.set_ticks([np.log1p(t) for t in ct]); cb.set_ticklabels([str(t) for t in ct])
    cb.set_label('SBS2 Weight',fontsize=FONT_SIZE-6)
    ax.set_xlabel('APOBEC3A Expression (FPKM-UQ)',fontsize=FONT_SIZE)
    ax.set_ylabel('APOBEC3B Expression (FPKM-UQ)',fontsize=FONT_SIZE)
    plt.tight_layout(); save_fig(fig,"Panel_1b_A3A_vs_A3B_SBS2")

# =============================================================================
# PANEL 1c: Box-and-whisker + 2x2 heatmap (UPDATED)
# =============================================================================
banner("Panel 1c: Boxplot + Heatmap")

quad_order = ['LOW_A3B_LOW_A3A','HIGH_A3B_LOW_A3A','LOW_A3B_HIGH_A3A','HIGH_A3B_HIGH_A3A']
quad_labels = ['A3B$^{-}$ A3A$^{-}$','A3B$^{+}$ A3A$^{-}$','A3B$^{-}$ A3A$^{+}$','A3B$^{+}$ A3A$^{+}$']
quad_colors = [COLOR_CREAM, '#e8c87a', '#e8a87a', COLOR_SBS2_HIGH]

fig = plt.figure(figsize=(20, 10))
gs = gridspec.GridSpec(1, 2, width_ratios=[2.5, 1.5], wspace=0.35)
ax_box = fig.add_subplot(gs[0])
ax_heat = fig.add_subplot(gs[1])

# --- Boxplot ---
box_data = [final[final['quadrant']==q]['SBS2'].values for q in quad_order]
medians = [np.median(d) for d in box_data]
ns = [len(d) for d in box_data]

bp = ax_box.boxplot(box_data, patch_artist=True, widths=0.6, showfliers=True,
                     flierprops=dict(marker='o', markersize=4, alpha=0.3,
                                     markerfacecolor=COLOR_DARK_GRAY, markeredgecolor=COLOR_DARK_GRAY))
for patch, color in zip(bp['boxes'], quad_colors):
    patch.set_facecolor(color); patch.set_alpha(0.75)
    patch.set_edgecolor(COLOR_BLACK); patch.set_linewidth(1.2)
for element in ['whiskers','caps']:
    for line in bp[element]: line.set_color(COLOR_DARK_GRAY); line.set_linewidth(1.2)
for line in bp['medians']:
    line.set_color(COLOR_BLACK); line.set_linewidth(2)

# Trend line connecting medians (diamonds)
ax_box.plot(range(1,5), medians, 'D-', color=COLOR_DARK_GRAY, markersize=10,
            markerfacecolor='white', markeredgecolor=COLOR_DARK_GRAY, markeredgewidth=2,
            linewidth=2, zorder=5)

# MEDIAN LABELS — positioned at ~50% of the y-axis range, above the box area
y_mid = y_all.max() * 0.50
for i, med in enumerate(medians):
    ax_box.text(i+1, y_mid, f'median = {med:.0f}',
                ha='center', fontsize=FONT_SIZE-8, fontweight='bold', color=COLOR_DARK_GRAY)

ax_box.set_xticklabels(quad_labels, fontsize=FONT_SIZE-6)
ax_box.set_ylabel('SBS2 Weight (mutation count)', fontsize=FONT_SIZE)
ax_box.tick_params(labelsize=FONT_SIZE-6)

# Significance bracket
p_both_vs_neither = mannwhitneyu(box_data[0], box_data[3], alternative='two-sided').pvalue
y_sig = max(np.percentile(box_data[3],99), np.percentile(box_data[0],99)) * 1.1
ax_box.plot([1, 1, 4, 4], [y_sig, y_sig*1.03, y_sig*1.03, y_sig],
            color=COLOR_DARK_GRAY, linewidth=1.5)
sig_text = f'p = {p_both_vs_neither:.2e}' if p_both_vs_neither < 0.001 else f'p = {p_both_vs_neither:.4f}'
ax_box.text(2.5, y_sig*1.05, sig_text, ha='center', fontsize=FONT_SIZE-8, color=COLOR_DARK_GRAY)

ax_box.spines['top'].set_visible(False); ax_box.spines['right'].set_visible(False)

# --- 2x2 Heatmap (wider, all black text) ---
# Rows = A3B (HIGH top, LOW bottom), Cols = A3A (LOW left, HIGH right)
heat_matrix = np.array([
    [medians[1], medians[3]],   # HIGH A3B row
    [medians[0], medians[2]],   # LOW A3B row
])
heat_n = np.array([
    [ns[1], ns[3]],
    [ns[0], ns[2]],
])

im = ax_heat.imshow(heat_matrix, cmap='YlOrRd', aspect='equal',
                     vmin=0, vmax=max(medians)*1.1)

# All text in BLACK
for i in range(2):
    for j in range(2):
        val = heat_matrix[i,j]
        n = heat_n[i,j]
        ax_heat.text(j, i, f'{val:.0f}\n(n={n})', ha='center', va='center',
                     fontsize=FONT_SIZE-4, fontweight='bold', color=COLOR_BLACK)

ax_heat.set_xticks([0,1]); ax_heat.set_xticklabels(['LOW','HIGH'], fontsize=FONT_SIZE-4)
ax_heat.set_yticks([0,1]); ax_heat.set_yticklabels(['HIGH','LOW'], fontsize=FONT_SIZE-4)
ax_heat.set_xlabel('A3A Expression', fontsize=FONT_SIZE-2)
ax_heat.set_ylabel('A3B Expression', fontsize=FONT_SIZE-2)

cb2 = plt.colorbar(im, ax=ax_heat, shrink=0.6, pad=0.08)
cb2.set_label('Median SBS2', fontsize=FONT_SIZE-6)
cb2.ax.tick_params(labelsize=FONT_SIZE-8)

plt.tight_layout()
save_fig(fig, "Panel_1c_Boxplot_Heatmap")

# =============================================================================
# SUPPLEMENTAL
# =============================================================================
banner("Supplemental: Low-A3 zoom")
hs = final[final['SBS2']>median_sbs2]; lo4 = hs.nsmallest(4,'A3A_plus_A3B')
if len(lo4)>0:
    fig,ax = plt.subplots(figsize=(10,8))
    zxm=5; zym=lo4['SBS2'].max()*1.3
    dx=np.linspace(0,min(X_THR,zxm),100); dy=safe_slope*dx
    ax.fill_between(dx,dy,zym,alpha=0.3,color=COLOR_NORMAL,zorder=0)
    ax.fill_between(dx,0,dy,alpha=0.3,color=COLOR_CREAM,zorder=0)
    if X_THR<zxm:
        ax.fill([X_THR,X_THR,zxm,zxm],[y_thr,zym,zym,y_thr],color=COLOR_SBS2_HIGH,alpha=0.3,zorder=0)
        ax.fill([X_THR,zxm,zxm,X_THR],[0,0,y_thr,y_thr],color=COLOR_CREAM,alpha=0.3,zorder=0)
    ax.plot(dx,dy,'--',color=COLOR_DARK_GRAY,lw=1.5,alpha=0.8,zorder=1)
    ax.axhline(y_thr,color=COLOR_DARK_GRAY,ls='--',lw=1.5,alpha=0.8,zorder=1)
    ax.plot([0,0],[0,zym],color=COLOR_NORMAL,lw=3,zorder=1,solid_capstyle='butt')
    ax.plot([0,zxm],[0,0],color=COLOR_BLACK,lw=3,zorder=1,solid_capstyle='butt')
    iz = final[(final['A3A_plus_A3B']<=zxm)&(final['SBS2']<=zym)]
    ax.scatter(iz['A3A_plus_A3B'],iz['SBS2'],c=[rcol[r] for r in iz['region']],
               s=80,alpha=0.8,edgecolors=COLOR_BLACK,linewidths=0.5,zorder=2)
    for _,r in lo4.iterrows():
        ax.annotate(r['Entity_ID'][:12],xy=(r['A3A_plus_A3B'],r['SBS2']),
                    xytext=(15,5),textcoords='offset points',fontsize=14,fontweight='bold',
                    color='#CC0000',arrowprops=dict(arrowstyle='->',color='#CC0000',lw=1.5))
    ax.set_xlim(-0.1,zxm); ax.set_ylim(-zym*0.05,zym)
    ax.set_xlabel('A3A + A3B Expression (FPKM-UQ)',fontsize=22)
    ax.set_ylabel('SBS2 Weight (mutation count)',fontsize=22)
    ax.set_title(f'Zoomed View: Low A3 Expression Region\n'
                 f'({len(lo4)} tumors with lowest A3 in SBS2 > median group labeled)',fontsize=18)
    ax.tick_params(labelsize=18); clean_spines(ax)
    plt.tight_layout(); save_fig(fig,"Supplemental_Low_A3_High_SBS2_Zoom")

# =============================================================================
# DIAGNOSTIC
# =============================================================================
banner("DIAGNOSTIC NUMBERS FOR TEXT")
log(f"  n tumors = {len(final)}")
log(f"  SBS2 range: 0-{final['SBS2'].max():.0f}, median={median_sbs2:.0f}")
log(f"  SBS2>0: {(final['SBS2']>0).sum()}/{len(final)} ({100*(final['SBS2']>0).mean():.1f}%)")
log(f"  Regions: teal={n_teal}, coral={n_coral}, cream={n_cream}")
log(f"  Slope={safe_slope:.4f}, X_THR={X_THR:.1f}")
log(f"  Median A3A={med_a3a:.2f}, A3B={med_a3b:.2f}")
log(f"\n  Quadrant SBS2 medians:")
for q,lab in zip(quad_order,['LOW/LOW','HIGH_A3B/LOW_A3A','LOW_A3B/HIGH_A3A','HIGH/HIGH']):
    sub=final[final['quadrant']==q]
    log(f"    {lab}: n={len(sub)}, median SBS2={sub['SBS2'].median():.1f}")
log(f"\n  Wilcoxon p-values:")
pairs = [('LOW_A3B_LOW_A3A','HIGH_A3B_LOW_A3A','A3B effect alone'),
         ('LOW_A3B_LOW_A3A','LOW_A3B_HIGH_A3A','A3A effect alone'),
         ('HIGH_A3B_LOW_A3A','HIGH_A3B_HIGH_A3A','Adding A3A to A3B'),
         ('LOW_A3B_LOW_A3A','HIGH_A3B_HIGH_A3A','Both vs neither')]
for q1,q2,lab in pairs:
    g1=final[final['quadrant']==q1]['SBS2']; g2=final[final['quadrant']==q2]['SBS2']
    if len(g1)>5 and len(g2)>5:
        _,p = mannwhitneyu(g1,g2,alternative='two-sided')
        log(f"    {lab}: {q1} (med={g1.median():.1f}) vs {q2} (med={g2.median():.1f}), p={p:.2e}")
log(f"\n  Low-A3 outliers (SBS2>median, A3A+A3B<1.0): {final['outlier_low_A3_high_SBS2'].sum()}")
for _,r in final[final['outlier_low_A3_high_SBS2']].iterrows():
    log(f"    {r['Entity_ID'][:12]}: A3A={r['A3A']:.2f}, A3B={r['A3B']:.2f}, SBS2={r['SBS2']:.0f}")
log(f"\n  Spearman correlations:")
for g in ['A3A','A3B','A3A_plus_A3B']:
    rho,p = spearmanr(final[g],final['SBS2'])
    log(f"    {g} vs SBS2: rho={rho:.4f}, p={p:.2e}")
rp = os.path.join(TROUBLE_DIR,"figure1_final_pipeline_report.txt")
with open(rp,'w') as f: f.write('\n'.join(report_lines))
log(f"\n  Report: {rp}")
banner("FIGURE 1 PIPELINE COMPLETE")
