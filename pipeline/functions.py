# ============================================================
# Helpers clean_tso_pipeline.ipynb
# ============================================================
import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, roc_auc_score, precision_recall_curve, average_precision_score
import pandas as pd
import pyranges as pr
from pptx.util import Inches, Pt
from pptx import Presentation
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

def assert_columns(df, required, name="df"):
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise KeyError(
            f"[{name}] Missing columns: {missing}\n"
            f"[{name}] Available columns: {df.columns.tolist()}\n"
            f"[{name}] Head:\n{df.head(3)}"
        )

def save_fig(fig, path, dpi=200):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

def safe_auc(y, s):
    # Handle degenerate cases where only one class present
    y = np.asarray(y)
    s = np.asarray(s)
    if len(np.unique(y)) < 2:
        return np.nan
    return roc_auc_score(y, s)

# ============================================================
# Score engineering (NO imbalance_score)
# ============================================================
def add_baf_scores(df):
    req = ["ratio_0-0.25", "ratio_0.25-0.5", "ratio_0.5-0.75", "ratio_0.75-1.0"]
    assert_columns(df, req, "add_baf_scores input")
    df = df.copy()

    df["outer_mass"] = df["ratio_0-0.25"] + df["ratio_0.75-1.0"]
    df["inner_mass"] = df["ratio_0.25-0.5"] + df["ratio_0.5-0.75"]
    df["skew"] = (df["ratio_0-0.25"] - df["ratio_0.75-1.0"]).abs()

    return df

# ============================================================
# Single metric plots (return figs)
# ============================================================
def fig_roc(df, label_col, score_col):
    d = df.dropna(subset=[label_col, score_col]).copy()
    y = d[label_col].to_numpy()
    s = d[score_col].to_numpy()

    fig = plt.figure(figsize=(4.4, 3.6))
    if len(np.unique(y)) < 2:
        plt.text(0.5, 0.5, "Only one class\nAUC=NA", ha="center", va="center")
        plt.xlim(0,1); plt.ylim(0,1)
        plt.title(f"ROC ({score_col})")
        plt.xlabel("FPR"); plt.ylabel("TPR")
        return fig

    fpr, tpr, _ = roc_curve(y, s, drop_intermediate=False)
    roc_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, lw=2, label=f"AUC={roc_auc:.3f}")
    plt.plot([0,1], [0,1], "--", lw=1)
    plt.xlabel("FPR"); plt.ylabel("TPR")
    plt.title(f"ROC ({score_col})")
    plt.legend(loc="lower right")
    return fig

def fig_pr(df, label_col, score_col):
    d = df.dropna(subset=[label_col, score_col]).copy()
    y = d[label_col].to_numpy()
    s = d[score_col].to_numpy()

    fig = plt.figure(figsize=(4.4, 3.6))
    if len(np.unique(y)) < 2:
        plt.text(0.5, 0.5, "Only one class\nAP=NA", ha="center", va="center")
        plt.xlim(0,1); plt.ylim(0,1)
        plt.title(f"PR ({score_col})")
        plt.xlabel("Recall"); plt.ylabel("Precision")
        return fig

    prec, rec, _ = precision_recall_curve(y, s)
    ap = average_precision_score(y, s)
    plt.plot(rec, prec, lw=2, label=f"AP={ap:.3f}")
    plt.xlabel("Recall"); plt.ylabel("Precision")
    plt.title(f"PR ({score_col})")
    plt.legend(loc="lower left")
    return fig

def fig_hist(df, label_col, score_col, neg_name="0", pos_name="1"):
    d = df.dropna(subset=[label_col, score_col]).copy()
    fig = plt.figure(figsize=(4.4, 3.6))
    plt.hist(d.loc[d[label_col] == 0, score_col], bins=30, alpha=0.6, label=neg_name)
    plt.hist(d.loc[d[label_col] == 1, score_col], bins=30, alpha=0.6, label=pos_name)
    plt.xlabel(score_col); plt.ylabel("count")
    plt.title(f"Score dist ({score_col})")
    plt.legend()
    return fig


def prepare_baf_origin(dataframe,baf_col='BAF',het_chr_col= "#CHR",baf_missing="./.",depth=0):
    """Extract and clean BAF dataframe from input."""
    test=dataframe.copy()
    baf_missing="./."
    test[baf_col] = test[baf_col].replace(baf_missing, np.nan)
    test["BAF_n"] = pd.to_numeric(test[baf_col], errors="coerce")
    test["POS"]   = pd.to_numeric(test["POS"], errors="coerce")
    test["CHROM"] = test[het_chr_col].astype(str)
    test2=test[test['DP']>0] #### previously was DP > 0
    baf_df = test2.dropna(subset=["CHROM","POS","BAF_n"])
    return baf_df

def baf_overlap_distribution(
    segments_df,
    baf_df,
    seg_chr="chr",
    seg_start="startpos",
    seg_end="endpos",
    baf_chr="CHROM",
    baf_pos="POS",
    baf_val="BAF_n",
):
    """
    Join SNP-level BAF to segments/bins and compute quartile distribution per segment.

    Returns
    -------
    hits : DataFrame
        SNP-in-segment table.
    seg_dist : DataFrame
        One row per segment with:
          - Chromosome, segment_id, startpos, endpos
          - n_baf_points
          - counts: 0-0.25, 0.25-0.5, 0.5-0.75, 0.75-1.0
          - ratios: ratio_0-0.25, ratio_0.25-0.5, ratio_0.5-0.75, ratio_0.75-1.0
        (Plus any extra columns originally in segments_df, carried through.)
    """

    seg = segments_df.copy()
    baf = baf_df.copy()

    seg.columns = seg.columns.str.strip()
    baf.columns = baf.columns.str.strip()

    # --- standardize chromosome
    seg["Chromosome"] = seg[seg_chr].astype(str).str.replace("^chr", "", regex=True)
    baf["Chromosome"] = baf[baf_chr].astype(str).str.replace("^chr", "", regex=True)

    # --- numeric coords + baf
    seg[seg_start] = pd.to_numeric(seg[seg_start], errors="coerce")
    seg[seg_end]   = pd.to_numeric(seg[seg_end], errors="coerce")
    baf[baf_pos]   = pd.to_numeric(baf[baf_pos], errors="coerce")
    baf[baf_val]   = pd.to_numeric(baf[baf_val], errors="coerce")

    seg = seg.dropna(subset=["Chromosome", seg_start, seg_end]).copy()
    baf = baf.dropna(subset=["Chromosome", baf_pos, baf_val]).copy()

    # keep BAF in [0,1]
    baf = baf[(baf[baf_val] >= 0) & (baf[baf_val] <= 1)].copy()

    # --- segment id
    seg = seg.sort_values(["Chromosome", seg_start, seg_end]).reset_index(drop=True)
    seg["segment_id"] = np.arange(len(seg))

    # --- 0-based half-open coordinates for PyRanges
    seg["Start"] = seg[seg_start].astype(int) - 1
    seg["End"]   = seg[seg_end].astype(int)

    baf["Start"] = baf[baf_pos].astype(int) - 1
    baf["End"]   = baf[baf_pos].astype(int)

    seg_pr = pr.PyRanges(seg[["Chromosome","Start","End","segment_id", seg_start, seg_end]])
    baf_pr = pr.PyRanges(baf[["Chromosome","Start","End", baf_pos, baf_val]])

    hits = baf_pr.join(seg_pr).df

    # Ensure columns exist even if empty join
    if hits.empty:
        # return empty seg_dist with expected ratio columns
        seg_base = seg.drop(columns=["Start","End"]).rename(columns={seg_start:"startpos", seg_end:"endpos"})
        for lab in ["0-0.25","0.25-0.5","0.5-0.75","0.75-1.0"]:
            seg_base[lab] = 0
            seg_base[f"ratio_{lab}"] = np.nan
        seg_base["n_baf_points"] = 0
        return hits, seg_base

    # --- assign quartile bin
    bins_edges = [0.0, 0.25, 0.5, 0.75, 1.0]
    labels = ["0-0.25","0.25-0.5","0.5-0.75","0.75-1.0"]
    hits["baf_bin"] = pd.cut(
        hits[baf_val],
        bins=bins_edges,
        labels=labels,
        include_lowest=True,
        right=True
    )

    # --- counts per segment x bin
    counts_long = (
        hits.groupby(["Chromosome","segment_id","baf_bin"], observed=True)
            .size().rename("n_points").reset_index()
    )

    counts_wide = (
        counts_long.pivot_table(
            index=["Chromosome","segment_id"],
            columns="baf_bin",
            values="n_points",
            fill_value=0,
            observed=True
        ).reset_index()
    )

    totals = (
        hits.groupby(["Chromosome","segment_id"], observed=True)
            .size().rename("n_baf_points").reset_index()
    )

    # --- base segment table (carry-through ALL original segment columns)
    seg_base = seg.drop(columns=["Start","End"]).copy()
    seg_base = seg_base.rename(columns={seg_start:"startpos", seg_end:"endpos"})
    seg_base["segment_len_bp"] = seg_base["endpos"].astype(int) - seg_base["startpos"].astype(int) + 1

    seg_dist = (
        seg_base.merge(totals, on=["Chromosome","segment_id"], how="left")
                .merge(counts_wide, on=["Chromosome","segment_id"], how="left")
    )

    seg_dist["n_baf_points"] = seg_dist["n_baf_points"].fillna(0).astype(int)

    # make sure all 4 bins exist + ratios
    for lab in labels:
        if lab not in seg_dist.columns:
            seg_dist[lab] = 0
        seg_dist[lab] = seg_dist[lab].fillna(0).astype(int)
        seg_dist[f"ratio_{lab}"] = np.where(
            seg_dist["n_baf_points"] > 0,
            seg_dist[lab] / seg_dist["n_baf_points"],
            np.nan
        )

    return hits, seg_dist


#######################
### functions to plot 
#######################

def add_title(slide, title, M=0.35, title_h=0.45):
    tx = slide.shapes.add_textbox(Inches(M), Inches(M), Inches(11.69 - 2*M), Inches(title_h))
    p = tx.text_frame.paragraphs[0]
    p.text = title
    p.font.size = Pt(18)

def add_picture(slide, path, x, y, w, h):
    slide.shapes.add_picture(path, Inches(x), Inches(y), width=Inches(w), height=Inches(h))

def add_table(slide, df, x, y, w, h, font_size=12):
    """
    df should be small. We'll render as a pptx table.
    """
    rows = df.shape[0] + 1
    cols = df.shape[1]
    table_shape = slide.shapes.add_table(rows, cols, Inches(x), Inches(y), Inches(w), Inches(h))
    table = table_shape.table

    # header
    for j, col in enumerate(df.columns):
        cell = table.cell(0, j)
        cell.text = str(col)
        for p in cell.text_frame.paragraphs:
            for r in p.runs:
                r.font.size = Pt(font_size)
                r.font.bold = True

    # body
    for i in range(df.shape[0]):
        for j in range(cols):
            cell = table.cell(i+1, j)
            val = df.iat[i, j]
            cell.text = "" if pd.isna(val) else str(val)
            for p in cell.text_frame.paragraphs:
                for r in p.runs:
                    r.font.size = Pt(font_size)

    return table_shape


def build_metrics_panel_slide(
    prs,
    title,
    track_path,
    metric_image_paths,  # dict score -> {"roc":..., "pr":..., "hist":...}
):
    """
    1 slide:
      Title
      Top: track
      Bottom: 3 rows (outer_mass, inner_mass, skew) x 3 cols (ROC/PR/HIST)
    """
    blank = prs.slide_layouts[6]
    slide = prs.slides.add_slide(blank)

    M = 0.35
    title_h = 0.45
    Gx = 0.15
    Gy = 0.12

    add_title(slide, title, M=M, title_h=title_h)

    # Track
    track_y = M + title_h
    track_h = 3.6
    add_picture(slide, track_path, x=M, y=track_y, w=11.69 - 2*M, h=track_h)

    # Grid for metrics
    # remaining height
    y0 = track_y + track_h + Gy
    remaining = 8.27 - y0 - M
    # 3 rows
    row_h = (remaining - 2*Gy) / 3.0
    # 3 cols
    col_w = (11.69 - 2*M - 2*Gx) / 3.0

    scores = ["outer_mass", "inner_mass", "skew"]
    for i, sc in enumerate(scores):
        imgs = metric_image_paths[sc]
        yy = y0 + i*(row_h + Gy)
        add_picture(slide, imgs["roc"],  x=M,              y=yy, w=col_w, h=row_h)
        add_picture(slide, imgs["pr"],   x=M + col_w + Gx, y=yy, w=col_w, h=row_h)
        add_picture(slide, imgs["hist"], x=M + 2*(col_w+Gx), y=yy, w=col_w, h=row_h)

    return slide

A4_W, A4_H = Inches(11.69), Inches(8.27)

def build_track_only_slide(prs, title, track_path):
    blank = prs.slide_layouts[6]
    slide = prs.slides.add_slide(blank)

    M = 0.35
    title_h = 0.45
    add_title(slide, title, M=M, title_h=title_h)

    add_picture(slide, track_path, x=M, y=M + title_h, w=11.69 - 2*M, h=8.27 - (M + title_h) - M)
    return slide

def build_table_slide(prs, title, table_df, model_auc):
    blank = prs.slide_layouts[6]
    slide = prs.slides.add_slide(blank)

    M = 0.35
    title_h = 0.45
    add_title(slide, title, M=M, title_h=title_h)

    # Add model AUC text
    box = slide.shapes.add_textbox(Inches(M), Inches(M + title_h), Inches(11.69 - 2*M), Inches(0.4))
    p = box.text_frame.paragraphs[0]
    p.text = f"Multivariate Logistic Regression AUC (in-sample): {model_auc:.3f}" if np.isfinite(model_auc) else "Multivariate Logistic Regression AUC: NA"
    p.font.size = Pt(14)

    # Format table for PPT
    tdf = table_df.copy()
    tdf["univariate_auc"] = tdf["univariate_auc"].map(lambda x: f"{x:.3f}" if pd.notna(x) else "NA")
    tdf["logreg_coef"] = tdf["logreg_coef"].map(lambda x: f"{x:+.3f}" if pd.notna(x) else "NA")

    # Table placement
    # Make it big and readable
    add_table(
        slide,
        tdf[["feature", "univariate_auc", "logreg_coef"]],
        x=M,
        y=M + title_h + 0.55,
        w=11.69 - 2*M,
        h=8.27 - (M + title_h + 0.75) - M,
        font_size=14
    )
    return slide

def std_chr(x):
    s = pd.Series(x).astype(str).str.strip()
    s = s.str.replace("^chr", "", regex=True)
    s = s.str.replace(r"\.0$", "", regex=True)  # "1.0" -> "1"
    s = s.str.upper().replace({"MT": "M"})
    return s

def make_fixed_bins(chr_sizes_df, bin_size_bp=10_000_000, chr_col="chr", len_col="chr_len_bp"):
    """
    Build fixed-width bins across each chromosome using chr_sizes_df.
    Output: Chromosome, bin_start (1-based), bin_end (1-based), bin_id
    """
    cs = chr_sizes_df.copy()
    cs[chr_col] = std_chr(cs[chr_col])
    cs[len_col] = pd.to_numeric(cs[len_col], errors="coerce")
    cs = cs.dropna(subset=[chr_col, len_col])

    rows = []
    for chrom, L in cs[[chr_col, len_col]].itertuples(index=False):
        L = int(L)
        # bins in 1-based inclusive coords
        for start in range(1, L + 1, bin_size_bp):
            end = min(start + bin_size_bp - 1, L)
            rows.append((chrom, start, end))

    bins = pd.DataFrame(rows, columns=["Chromosome", "bin_start", "bin_end"])
    bins["bin_id"] = np.arange(len(bins))
    return bins



def aggregate_bins_by_total_cn(bins_df, aggregation_method="mean"):
    """
    Aggregate bins within the same segment that have the same total CN.
    
    For each segment, if consecutive bins have the same total_cn, merge them
    and aggregate their BAF metrics using mean or median.
    
    Parameters
    ----------
    bins_df : DataFrame
        Bins with segment_id, total_cn, and BAF metrics
    aggregation_method : str
        'mean' or 'median' for aggregating BAF ratios
        
    Returns
    -------
    aggregated_df : DataFrame
        Merged bins with aggregated BAF metrics
    """
    df = bins_df.copy()
    
    # Sort by segment and position
    df = df.sort_values(["segment_id", "startpos"]).reset_index(drop=True)
    
    # Create grouping key: consecutive bins with same segment_id AND total_cn
    df["cn_group"] = (
        (df["segment_id"] != df["segment_id"].shift()) |
        (df["total_cn"] != df["total_cn"].shift())
    ).cumsum()
    
    # Aggregation function
    agg_func = np.mean if aggregation_method == "mean" else np.median
    
    # Aggregate
    agg_dict = {
        "Chromosome": "first",
        "startpos": "min",  # Take earliest start
        "endpos": "max",    # Take latest end
        "segment_id": "first",
        "original_seg_start": "first",
        "original_seg_end": "first",
        "total_cn": "first",
        "nMajor": "first",  # Assuming same within group
        "nMinor": "first",
        
        # BAF counts - sum them
        "0-0.25": "sum",
        "0.25-0.5": "sum",
        "0.5-0.75": "sum",
        "0.75-1.0": "sum",
        "n_baf_points": "sum",
        
        # Derived metrics - will recalculate
        "bin_actual_size": "sum",
        "bin_id": "first",  # Keep first bin_id as reference
    }
    
    aggregated = df.groupby("cn_group").agg(agg_dict).reset_index(drop=True)
    
    # Recalculate ratios based on aggregated counts
    labels = ["0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1.0"]
    for lab in labels:
        aggregated[f"ratio_{lab}"] = np.where(
            aggregated["n_baf_points"] > 0,
            aggregated[lab] / aggregated["n_baf_points"],
            np.nan
        )
    
    # Recalculate BAF features
    aggregated["outer_mass"] = aggregated["ratio_0-0.25"] + aggregated["ratio_0.75-1.0"]
    aggregated["inner_mass"] = aggregated["ratio_0.25-0.5"] + aggregated["ratio_0.5-0.75"]
    aggregated["skew"] = (aggregated["ratio_0-0.25"] - aggregated["ratio_0.75-1.0"]).abs()
    
    # Add metadata
    aggregated["bin_type"] = "aggregated_total_cn"
    aggregated["is_full_bin"] = False
    aggregated["segment_len_bp"] = aggregated["endpos"] - aggregated["startpos"] + 1
    
    print(f"Aggregated {len(df)} bins → {len(aggregated)} total-CN groups")
    
    return aggregated



def aggregate_bins_by_allelic_cn(bins_df, aggregation_method="mean"):
    """
    Aggregate bins within the same segment that have the same allelic CN (nMajor + nMinor).
    
    For each segment, if consecutive bins have the same nMajor AND nMinor, merge them
    and aggregate their BAF metrics.
    
    Parameters
    ----------
    bins_df : DataFrame
        Bins with segment_id, nMajor, nMinor, and BAF metrics
    aggregation_method : str
        'mean' or 'median' for aggregating BAF ratios
        
    Returns
    -------
    aggregated_df : DataFrame
        Merged bins with aggregated BAF metrics
    """
    df = bins_df.copy()
    
    # Sort by segment and position
    df = df.sort_values(["segment_id", "startpos"]).reset_index(drop=True)
    
    # Create grouping key: consecutive bins with same segment_id AND allelic CN
    df["allelic_group"] = (
        (df["segment_id"] != df["segment_id"].shift()) |
        (df["nMajor"] != df["nMajor"].shift()) |
        (df["nMinor"] != df["nMinor"].shift())
    ).cumsum()
    
    # Aggregation function
    agg_func = np.mean if aggregation_method == "mean" else np.median
    
    # Aggregate
    agg_dict = {
        "Chromosome": "first",
        "startpos": "min",
        "endpos": "max",
        "segment_id": "first",
        "original_seg_start": "first",
        "original_seg_end": "first",
        "nMajor": "first",
        "nMinor": "first",
        "total_cn": "first",
        
        # BAF counts - sum them
        "0-0.25": "sum",
        "0.25-0.5": "sum",
        "0.5-0.75": "sum",
        "0.75-1.0": "sum",
        "n_baf_points": "sum",
        
        # Derived metrics
        "bin_actual_size": "sum",
        "bin_id": "first",
    }
    
    aggregated = df.groupby("allelic_group").agg(agg_dict).reset_index(drop=True)
    
    # Recalculate ratios
    labels = ["0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1.0"]
    for lab in labels:
        aggregated[f"ratio_{lab}"] = np.where(
            aggregated["n_baf_points"] > 0,
            aggregated[lab] / aggregated["n_baf_points"],
            np.nan
        )
    
    # Recalculate BAF features
    aggregated["outer_mass"] = aggregated["ratio_0-0.25"] + aggregated["ratio_0.75-1.0"]
    aggregated["inner_mass"] = aggregated["ratio_0.25-0.5"] + aggregated["ratio_0.5-0.75"]
    aggregated["skew"] = (aggregated["ratio_0-0.25"] - aggregated["ratio_0.75-1.0"]).abs()
    
    # Add metadata
    aggregated["bin_type"] = "aggregated_allelic_cn"
    aggregated["is_full_bin"] = False
    aggregated["segment_len_bp"] = aggregated["endpos"] - aggregated["startpos"] + 1
    
    print(f"Aggregated {len(df)} bins → {len(aggregated)} allelic-CN groups")
    
    return aggregated


def old_make_bins_from_segments_linear_plus_right(
    segments_df,
    bin_size=10_000_000,
    seg_chr="chr",
    seg_start="startpos",
    seg_end="endpos",
    cn_cols=("nMajor", "nMinor", "total_cn"),
    min_segment_size=None
):
    """
    Create bins from segments using linear tiling PLUS right-anchored bin.
    
    Strategy:
    1. Create full bins from left
    2. If remainder < bin_size:
       - Create partial bin for remainder
       - ALSO create full bin from right
       - These two will overlap
    
    Example: 24 Mb segment, 10 Mb bins
      Bin 1: 0-10 Mb (full, left)
      Bin 2: 10-20 Mb (full, left)
      Bin 3: 20-24 Mb (partial, remainder)
      Bin 4: 14-24 Mb (full, right, overlaps with Bin 3)
    """
    seg = segments_df.copy()
    seg.columns = seg.columns.str.strip()
    
    seg["Chromosome"] = seg[seg_chr].astype(str).str.replace("^chr", "", regex=True)
    seg[seg_start] = pd.to_numeric(seg[seg_start], errors="coerce")
    seg[seg_end] = pd.to_numeric(seg[seg_end], errors="coerce")
    
    for col in cn_cols:
        if col in seg.columns:
            seg[col] = pd.to_numeric(seg[col], errors="coerce")
    
    seg = seg.dropna(subset=["Chromosome", seg_start, seg_end])
    seg = seg.sort_values(["Chromosome", seg_start]).reset_index(drop=True)
    seg["segment_id"] = np.arange(len(seg))
    
    bins_list = []
    bin_counter = 0
    
    for _, row in seg.iterrows():
        chrom = row["Chromosome"]
        s_start = int(row[seg_start])
        s_end = int(row[seg_end])
        seg_id = row["segment_id"]
        seg_len = s_end - s_start + 1
        
        cn_vals = {col: row[col] for col in cn_cols if col in seg.columns}
        
        # Step 1: Create full bins from left
        left_bins = []
        pos = s_start
        
        while pos + bin_size - 1 <= s_end:
            bin_end = pos + bin_size - 1
            
            left_bins.append({
                "Chromosome": chrom,
                "startpos": pos,
                "endpos": bin_end,
                "bin_id": bin_counter,
                "segment_id": seg_id,
                "original_seg_start": s_start,
                "original_seg_end": s_end,
                "bin_type": "full_left",
                "bin_position": "left",
                "is_full_bin": True,
                "bin_actual_size": bin_size,
                "overlap_type": "none",
                **cn_vals
            })
            
            bin_counter += 1
            pos = bin_end + 1
        
        # Step 2: Check for remainder
        remainder_start = pos
        remainder_len = s_end - remainder_start + 1
        
        if remainder_len > 0:
            # Create partial bin for remainder
            bins_list.extend(left_bins)
            
            bins_list.append({
                "Chromosome": chrom,
                "startpos": remainder_start,
                "endpos": s_end,
                "bin_id": bin_counter,
                "segment_id": seg_id,
                "original_seg_start": s_start,
                "original_seg_end": s_end,
                "bin_type": "partial_remainder",
                "bin_position": "remainder",
                "is_full_bin": False,
                "bin_actual_size": remainder_len,
                "overlap_type": "left_of_overlap" if remainder_len < bin_size else "none",
                **cn_vals
            })
            
            partial_bin_id = bin_counter
            bin_counter += 1
            
            # If remainder < bin_size, also create full bin from right
            if remainder_len < bin_size:
                right_start = s_end - bin_size + 1
                
                # Calculate overlap with partial bin
                overlap_start = max(right_start, remainder_start)
                overlap_end = s_end
                overlap_len = overlap_end - overlap_start + 1
                
                bins_list.append({
                    "Chromosome": chrom,
                    "startpos": right_start,
                    "endpos": s_end,
                    "bin_id": bin_counter,
                    "segment_id": seg_id,
                    "original_seg_start": s_start,
                    "original_seg_end": s_end,
                    "bin_type": "full_right",
                    "bin_position": "right",
                    "is_full_bin": True,
                    "bin_actual_size": bin_size,
                    "overlap_type": "right_of_overlap",
                    "overlap_with_bin_id": partial_bin_id,
                    "overlap_start": overlap_start,
                    "overlap_end": overlap_end,
                    "overlap_len": overlap_len,
                    **cn_vals
                })
                
                # Update partial bin with overlap info
                bins_list[-2]["overlap_with_bin_id"] = bin_counter
                bins_list[-2]["overlap_start"] = overlap_start
                bins_list[-2]["overlap_end"] = overlap_end
                bins_list[-2]["overlap_len"] = overlap_len
                
                bin_counter += 1
        else:
            # No remainder - segment is exactly divisible
            bins_list.extend(left_bins)
    
    bins_df = pd.DataFrame(bins_list)
    
    # Ensure overlap columns exist
    for col in ["overlap_with_bin_id", "overlap_start", "overlap_end", "overlap_len"]:
        if col not in bins_df.columns:
            bins_df[col] = np.nan
    
    return bins_df



def old_baf_distribution_on_segment_bins_with_aggregation(
    bins_df,
    baf_df,
    bin_chr="Chromosome",
    bin_start="startpos",
    bin_end="endpos",
    baf_chr="CHROM",
    baf_pos="POS",
    baf_val="BAF_n",
    aggregate_overlapping=True,
    remove_right_after_aggregation=True
):
    """
    Calculate BAF quartile distributions for segment bins.
    
    For overlapping partial_remainder + full_right pairs:
    - Computes weighted mean of BAF metrics
    - Assigns aggregated metrics to partial_remainder bin
    - Removes full_right bin
    """
    import pyranges as pr
    
    bins = bins_df.copy()
    baf = baf_df.copy()
    
    bins.columns = bins.columns.str.strip()
    baf.columns = baf.columns.str.strip()
    
    bins["Chromosome"] = bins[bin_chr].astype(str).str.replace("^chr", "", regex=True)
    baf["Chromosome"] = baf[baf_chr].astype(str).str.replace("^chr", "", regex=True)
    
    bins = bins.sort_values(["Chromosome", bin_start]).reset_index(drop=True)
    if "bin_id" not in bins.columns:
        bins["bin_id"] = np.arange(len(bins))
    
    baf = baf[[baf_chr, baf_pos, baf_val]].copy()
    baf["Chromosome"] = baf[baf_chr].astype(str).str.replace("^chr", "", regex=True)
    baf[baf_pos] = pd.to_numeric(baf[baf_pos], errors="coerce")
    baf[baf_val] = pd.to_numeric(baf[baf_val], errors="coerce")
    baf = baf.dropna(subset=[baf_pos, baf_val])
    baf = baf[(baf[baf_val] >= 0) & (baf[baf_val] <= 1)]
    
    bins["Start"] = pd.to_numeric(bins[bin_start], errors="coerce").astype(int) - 1
    bins["End"] = pd.to_numeric(bins[bin_end], errors="coerce").astype(int)
    
    baf["Start"] = baf[baf_pos].astype(int) - 1
    baf["End"] = baf[baf_pos].astype(int)
    
    bins_pr = pr.PyRanges(bins[["Chromosome", "Start", "End", "bin_id"]])
    baf_pr = pr.PyRanges(baf[["Chromosome", "Start", "End", baf_pos, baf_val]])
    
    hits = baf_pr.join(bins_pr).df
    
    bins_edges = [0.0, 0.25, 0.5, 0.75, 1.0]
    labels = ["0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1.0"]
    
    if len(hits):
        hits["baf_bin"] = pd.cut(
            hits[baf_val], 
            bins=bins_edges, 
            labels=labels,
            include_lowest=True, 
            right=True
        )
    else:
        hits["baf_bin"] = pd.Series(dtype="category")
    
    # ========== PHASE A: Calculate BAF metrics per bin ==========
    counts_long = (
        hits.groupby(["bin_id", "baf_bin"], observed=True)
            .size().rename("n_points").reset_index()
    )
    
    counts_wide = (
        counts_long.pivot_table(
            index="bin_id", 
            columns="baf_bin",
            values="n_points", 
            fill_value=0, 
            observed=True
        ).reset_index()
    )
    
    totals = (
        hits.groupby("bin_id", observed=True)
            .size().rename("n_baf_points").reset_index()
    )
    
    bins_base = bins[[bin_chr, bin_start, bin_end, "bin_id"]].copy()
    bins_base = bins_base.rename(columns={
        bin_chr: "Chromosome",
        bin_start: "startpos",
        bin_end: "endpos"
    })
    
    bins_dist = (
        bins_base.merge(totals, on="bin_id", how="left")
                 .merge(counts_wide, on="bin_id", how="left")
    )
    
    bins_dist["n_baf_points"] = bins_dist["n_baf_points"].fillna(0).astype(int)
    
    for lab in labels:
        if lab not in bins_dist.columns:
            bins_dist[lab] = 0
        bins_dist[lab] = bins_dist[lab].fillna(0).astype(int)
        bins_dist[f"ratio_{lab}"] = np.where(
            bins_dist["n_baf_points"] > 0,
            bins_dist[lab] / bins_dist["n_baf_points"],
            np.nan
        )
    
    bins_with_dist = bins.merge(
        bins_dist[["bin_id", "n_baf_points"] + labels + [f"ratio_{lab}" for lab in labels]],
        on="bin_id",
        how="left"
    )
    
    # ========== PHASE B: AGGREGATE OVERLAPPING BINS ==========
    # THIS IS THE WEIGHTED AVERAGING PHASE!
    
    if aggregate_overlapping and 'bin_type' in bins_with_dist.columns:
        print("\n" + "="*70)
        print("AGGREGATING OVERLAPPING BINS (weighted by bin length)")
        print("Keeping partial_remainder with aggregated metrics, removing full_right")
        print("="*70)
        
        # Identify overlapping pairs
        partial_bins = bins_with_dist[bins_with_dist['bin_type'] == 'partial_remainder'].copy()
        
        bins_to_remove = []
        
        for _, partial in partial_bins.iterrows():
            # Find corresponding full_right bin
            if 'overlap_with_bin_id' in partial and pd.notna(partial['overlap_with_bin_id']):
                right_bin_id = int(partial['overlap_with_bin_id'])
                
                # Get the full_right bin
                right_bin = bins_with_dist[bins_with_dist['bin_id'] == right_bin_id]
                
                if len(right_bin) == 0:
                    continue
                    
                right_bin = right_bin.iloc[0]
                partial_bin_id = partial['bin_id']
                
                print(f"\nAggregating pair:")
                print(f"  Partial (bin_id={partial_bin_id}): {partial['startpos']}-{partial['endpos']} ({partial['bin_actual_size']/1e6:.2f} Mb)")
                print(f"  Right (bin_id={right_bin_id}): {right_bin['startpos']}-{right_bin['endpos']} ({right_bin['bin_actual_size']/1e6:.2f} Mb)")
                
                # Weights based on bin size (length)
                weight_partial = partial['bin_actual_size']
                weight_right = right_bin['bin_actual_size']
                total_weight = weight_partial + weight_right
                
                print(f"  Weights: partial={weight_partial/1e6:.2f} Mb ({weight_partial/total_weight*100:.1f}%), right={weight_right/1e6:.2f} Mb ({weight_right/total_weight*100:.1f}%)")
                
                # ===== WEIGHTED MEAN CALCULATION =====
                aggregated_metrics = {}
                
                # For each quartile ratio, compute weighted mean
                for lab in labels:
                    ratio_col = f'ratio_{lab}'
                    
                    ratio_partial = partial[ratio_col] if pd.notna(partial[ratio_col]) else 0.0
                    ratio_right = right_bin[ratio_col] if pd.notna(right_bin[ratio_col]) else 0.0
                    
                    # Weighted mean by length
                    weighted_ratio = (ratio_partial * weight_partial + ratio_right * weight_right) / total_weight
                    aggregated_metrics[ratio_col] = weighted_ratio
                    
                    print(f"    {ratio_col}: ({ratio_partial:.3f} × {weight_partial/1e6:.2f} + {ratio_right:.3f} × {weight_right/1e6:.2f}) / {total_weight/1e6:.2f} = {weighted_ratio:.3f}")
                
                # Weighted n_baf_points
                n_baf_partial = partial['n_baf_points'] if pd.notna(partial['n_baf_points']) else 0
                n_baf_right = right_bin['n_baf_points'] if pd.notna(right_bin['n_baf_points']) else 0
                
                aggregated_n_baf = int((n_baf_partial * weight_partial + n_baf_right * weight_right) / total_weight)
                aggregated_metrics['n_baf_points'] = aggregated_n_baf
                
                # Back-calculate counts from weighted ratios
                for lab in labels:
                    ratio_col = f'ratio_{lab}'
                    count = int(aggregated_metrics[ratio_col] * aggregated_n_baf)
                    aggregated_metrics[lab] = count
                
                # Calculate derived metrics
                outer_mass = aggregated_metrics['ratio_0-0.25'] + aggregated_metrics['ratio_0.75-1.0']
                inner_mass = aggregated_metrics['ratio_0.25-0.5'] + aggregated_metrics['ratio_0.5-0.75']
                
                print(f"  Aggregated n_baf_points (weighted): {aggregated_n_baf}")
                print(f"  Aggregated outer_mass: {aggregated_metrics['ratio_0-0.25']:.3f} + {aggregated_metrics['ratio_0.75-1.0']:.3f} = {outer_mass:.3f}")
                print(f"  Aggregated inner_mass: {aggregated_metrics['ratio_0.25-0.5']:.3f} + {aggregated_metrics['ratio_0.5-0.75']:.3f} = {inner_mass:.3f}")
                
                # UPDATE PARTIAL_REMAINDER BIN with aggregated metrics
                idx = bins_with_dist[bins_with_dist['bin_id'] == partial_bin_id].index[0]
                
                for lab in labels:
                    bins_with_dist.at[idx, lab] = aggregated_metrics[lab]
                    bins_with_dist.at[idx, f'ratio_{lab}'] = aggregated_metrics[f'ratio_{lab}']
                
                bins_with_dist.at[idx, 'n_baf_points'] = aggregated_metrics['n_baf_points']
                
                # Mark as aggregated
                bins_with_dist.at[idx, 'is_aggregated'] = True
                bins_with_dist.at[idx, 'aggregated_from_bin_id'] = right_bin_id
                bins_with_dist.at[idx, 'aggregation_weight_partial'] = weight_partial
                bins_with_dist.at[idx, 'aggregation_weight_right'] = weight_right
                
                # Update bin_type
                bins_with_dist.at[idx, 'bin_type'] = 'partial_remainder_aggregated'
                
                # Mark FULL_RIGHT bin for removal
                if remove_right_after_aggregation:
                    bins_to_remove.append(right_bin_id)
                    print(f"  → Keeping partial bin {partial_bin_id} with aggregated metrics")
                    print(f"  → Removing full_right bin {right_bin_id}")
        
        # Remove full_right bins
        if remove_right_after_aggregation and bins_to_remove:
            print(f"\nRemoving {len(bins_to_remove)} full_right bins after aggregation")
            bins_with_dist = bins_with_dist[~bins_with_dist['bin_id'].isin(bins_to_remove)].copy()
            bins_with_dist = bins_with_dist.reset_index(drop=True)
        
        print("\n" + "="*70)
    
    # Add aggregation flags if not present
    if 'is_aggregated' not in bins_with_dist.columns:
        bins_with_dist['is_aggregated'] = False
    if 'aggregated_from_bin_id' not in bins_with_dist.columns:
        bins_with_dist['aggregated_from_bin_id'] = np.nan
    
    return hits, bins_with_dist,baf_df


################### new versions

def make_bins_from_segments_linear_plus_right(
    segments_df,
    bin_size=10_000_000,
    seg_chr="chr",
    seg_start="startpos",
    seg_end="endpos",
    cn_cols=("nMajor", "nMinor", "total_cn"),
    min_segment_size=None
):
    """
    Create bins from segments using linear tiling PLUS right-anchored bin.
    
    Strategy:
    1. Create full bins from left
    2. If remainder < bin_size AND segment_length > bin_size:
       - Create partial bin for remainder
       - ALSO create full bin from right
       - These two will overlap
    3. If segment_length <= bin_size:
       - Create only one partial bin (NO right-anchored bin)
    
    Example 1: 6 Mb segment, 10 Mb bins (NO right-anchored)
      Bin 1: 0-6 Mb (partial, no overlap)
    
    Example 2: 15 Mb segment, 10 Mb bins (YES right-anchored)
      Bin 1: 0-10 Mb (full, left)
      Bin 2: 10-15 Mb (partial, remainder)
      Bin 3: 5-15 Mb (full, right, overlaps with Bin 2)
    
    Example 3: 24 Mb segment, 10 Mb bins (YES right-anchored)
      Bin 1: 0-10 Mb (full, left)
      Bin 2: 10-20 Mb (full, left)
      Bin 3: 20-24 Mb (partial, remainder)
      Bin 4: 14-24 Mb (full, right, overlaps with Bin 3)
    """
    seg = segments_df.copy()
    seg.columns = seg.columns.str.strip()
    
    seg["Chromosome"] = seg[seg_chr].astype(str).str.replace("^chr", "", regex=True)
    seg[seg_start] = pd.to_numeric(seg[seg_start], errors="coerce")
    seg[seg_end] = pd.to_numeric(seg[seg_end], errors="coerce")
    
    for col in cn_cols:
        if col in seg.columns:
            seg[col] = pd.to_numeric(seg[col], errors="coerce")
    
    seg = seg.dropna(subset=["Chromosome", seg_start, seg_end])
    seg = seg.sort_values(["Chromosome", seg_start]).reset_index(drop=True)
    seg["segment_id"] = np.arange(len(seg))
    
    bins_list = []
    bin_counter = 0
    
    for _, row in seg.iterrows():
        chrom = row["Chromosome"]
        s_start = int(row[seg_start])
        s_end = int(row[seg_end])
        seg_id = row["segment_id"]
        seg_len = s_end - s_start + 1
        
        cn_vals = {col: row[col] for col in cn_cols if col in seg.columns}
        
        # Step 1: Create full bins from left
        left_bins = []
        pos = s_start
        
        while pos + bin_size - 1 <= s_end:
            bin_end = pos + bin_size - 1
            
            left_bins.append({
                "Chromosome": chrom,
                "startpos": pos,
                "endpos": bin_end,
                "bin_id": bin_counter,
                "segment_id": seg_id,
                "original_seg_start": s_start,
                "original_seg_end": s_end,
                "original_seg_length": seg_len,
                "bin_type": "full_left",
                "bin_position": "left",
                "is_full_bin": True,
                "bin_actual_size": bin_size,
                "overlap_type": "none",
                **cn_vals
            })
            
            bin_counter += 1
            pos = bin_end + 1
        
        # Step 2: Check for remainder
        remainder_start = pos
        remainder_len = s_end - remainder_start + 1
        
        if remainder_len > 0:
            # Create partial bin for remainder
            bins_list.extend(left_bins)
            
            bins_list.append({
                "Chromosome": chrom,
                "startpos": remainder_start,
                "endpos": s_end,
                "bin_id": bin_counter,
                "segment_id": seg_id,
                "original_seg_start": s_start,
                "original_seg_end": s_end,
                "original_seg_length": seg_len,
                "bin_type": "partial_remainder",
                "bin_position": "remainder",
                "is_full_bin": False,
                "bin_actual_size": remainder_len,
                "overlap_type": "left_of_overlap" if (remainder_len < bin_size and seg_len > bin_size) else "none",
                **cn_vals
            })
            
            partial_bin_id = bin_counter
            bin_counter += 1
            
            # CRITICAL DECISION: Only create right-anchored bin if segment is LARGE
            # If remainder < bin_size AND segment_length > bin_size: create full_right
            # If segment_length <= bin_size: do NOT create full_right
            if remainder_len < bin_size and seg_len > bin_size:
                right_start = s_end - bin_size + 1
                
                # Calculate overlap with partial bin
                overlap_start = max(right_start, remainder_start)
                overlap_end = s_end
                overlap_len = overlap_end - overlap_start + 1
                
                bins_list.append({
                    "Chromosome": chrom,
                    "startpos": right_start,
                    "endpos": s_end,
                    "bin_id": bin_counter,
                    "segment_id": seg_id,
                    "original_seg_start": s_start,
                    "original_seg_end": s_end,
                    "original_seg_length": seg_len,
                    "bin_type": "full_right",
                    "bin_position": "right",
                    "is_full_bin": True,
                    "bin_actual_size": bin_size,
                    "overlap_type": "right_of_overlap",
                    "overlap_with_bin_id": partial_bin_id,
                    "overlap_start": overlap_start,
                    "overlap_end": overlap_end,
                    "overlap_len": overlap_len,
                    **cn_vals
                })
                
                # Update partial bin with overlap info
                bins_list[-2]["overlap_with_bin_id"] = bin_counter
                bins_list[-2]["overlap_start"] = overlap_start
                bins_list[-2]["overlap_end"] = overlap_end
                bins_list[-2]["overlap_len"] = overlap_len
                
                bin_counter += 1
                
                print(f"  Segment {seg_id} ({seg_len/1e6:.2f} Mb): Created full_right bin for {remainder_len/1e6:.2f} Mb remainder")
            else:
                if seg_len <= bin_size:
                    print(f"  Segment {seg_id} ({seg_len/1e6:.2f} Mb): Small segment, no full_right bin created")
        else:
            # No remainder - segment is exactly divisible
            bins_list.extend(left_bins)
    
    bins_df = pd.DataFrame(bins_list)
    
    # Ensure overlap columns exist
    for col in ["overlap_with_bin_id", "overlap_start", "overlap_end", "overlap_len"]:
        if col not in bins_df.columns:
            bins_df[col] = np.nan
    
    # Ensure original_seg_length exists
    if "original_seg_length" not in bins_df.columns:
        bins_df["original_seg_length"] = np.nan
    
    return bins_df


def baf_distribution_on_segment_bins_with_aggregation(
    bins_df,
    baf_df,
    bin_chr="Chromosome",
    bin_start="startpos",
    bin_end="endpos",
    baf_chr="CHROM",
    baf_pos="POS",
    baf_val="BAF_n",
    aggregate_overlapping=True,
    remove_right_after_aggregation=True
):
    """
    Calculate BAF quartile distributions for segment bins.
    
    For overlapping partial_remainder + full_right pairs:
    - ONLY aggregates if original segment length > bin_size
    - Computes weighted mean of BAF metrics
    - Assigns aggregated metrics to partial_remainder bin
    - Removes full_right bin
    
    For small segments (length <= bin_size):
    - NO aggregation (there's no full_right bin anyway)
    """
    import pyranges as pr
    
    bins = bins_df.copy()
    baf = baf_df.copy()
    
    bins.columns = bins.columns.str.strip()
    baf.columns = baf.columns.str.strip()
    
    bins["Chromosome"] = bins[bin_chr].astype(str).str.replace("^chr", "", regex=True)
    baf["Chromosome"] = baf[baf_chr].astype(str).str.replace("^chr", "", regex=True)
    
    bins = bins.sort_values(["Chromosome", bin_start]).reset_index(drop=True)
    if "bin_id" not in bins.columns:
        bins["bin_id"] = np.arange(len(bins))
    
    baf = baf[[baf_chr, baf_pos, baf_val]].copy()
    baf["Chromosome"] = baf[baf_chr].astype(str).str.replace("^chr", "", regex=True)
    baf[baf_pos] = pd.to_numeric(baf[baf_pos], errors="coerce")
    baf[baf_val] = pd.to_numeric(baf[baf_val], errors="coerce")
    baf = baf.dropna(subset=[baf_pos, baf_val])
    baf = baf[(baf[baf_val] >= 0) & (baf[baf_val] <= 1)]
    
    bins["Start"] = pd.to_numeric(bins[bin_start], errors="coerce").astype(int) - 1
    bins["End"] = pd.to_numeric(bins[bin_end], errors="coerce").astype(int)
    
    baf["Start"] = baf[baf_pos].astype(int) - 1
    baf["End"] = baf[baf_pos].astype(int)
    
    bins_pr = pr.PyRanges(bins[["Chromosome", "Start", "End", "bin_id"]])
    baf_pr = pr.PyRanges(baf[["Chromosome", "Start", "End", baf_pos, baf_val]])
    
    hits = baf_pr.join(bins_pr).df
    
    bins_edges = [0.0, 0.25, 0.5, 0.75, 1.0]
    labels = ["0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1.0"]
    
    if len(hits):
        hits["baf_bin"] = pd.cut(
            hits[baf_val], 
            bins=bins_edges, 
            labels=labels,
            include_lowest=True, 
            right=True
        )
    else:
        hits["baf_bin"] = pd.Series(dtype="category")
    
    # ========== PHASE A: Calculate BAF metrics per bin ==========
    counts_long = (
        hits.groupby(["bin_id", "baf_bin"], observed=True)
            .size().rename("n_points").reset_index()
    )
    
    counts_wide = (
        counts_long.pivot_table(
            index="bin_id", 
            columns="baf_bin",
            values="n_points", 
            fill_value=0, 
            observed=True
        ).reset_index()
    )
    
    totals = (
        hits.groupby("bin_id", observed=True)
            .size().rename("n_baf_points").reset_index()
    )
    
    bins_base = bins[[bin_chr, bin_start, bin_end, "bin_id"]].copy()
    bins_base = bins_base.rename(columns={
        bin_chr: "Chromosome",
        bin_start: "startpos",
        bin_end: "endpos"
    })
    
    bins_dist = (
        bins_base.merge(totals, on="bin_id", how="left")
                 .merge(counts_wide, on="bin_id", how="left")
    )
    
    bins_dist["n_baf_points"] = bins_dist["n_baf_points"].fillna(0).astype(int)
    
    for lab in labels:
        if lab not in bins_dist.columns:
            bins_dist[lab] = 0
        bins_dist[lab] = bins_dist[lab].fillna(0).astype(int)
        bins_dist[f"ratio_{lab}"] = np.where(
            bins_dist["n_baf_points"] > 0,
            bins_dist[lab] / bins_dist["n_baf_points"],
            np.nan
        )
    
    bins_with_dist = bins.merge(
        bins_dist[["bin_id", "n_baf_points"] + labels + [f"ratio_{lab}" for lab in labels]],
        on="bin_id",
        how="left"
    )
    
    # ========== PHASE B: AGGREGATE OVERLAPPING BINS ==========
    # ONLY if original segment length > bin_size
    
    if aggregate_overlapping and 'bin_type' in bins_with_dist.columns:
        print("\n" + "="*70)
        print("AGGREGATING OVERLAPPING BINS (weighted by bin length)")
        print("ONLY for segments with original_seg_length > bin_size")
        print("Keeping partial_remainder with aggregated metrics, removing full_right")
        print("="*70)
        
        # Identify overlapping pairs
        partial_bins = bins_with_dist[bins_with_dist['bin_type'] == 'partial_remainder'].copy()
        
        bins_to_remove = []
        aggregated_count = 0
        skipped_count = 0
        
        for _, partial in partial_bins.iterrows():
            # Check if this partial bin should be aggregated
            # It should have overlap_with_bin_id AND original segment should be large
            if 'overlap_with_bin_id' not in partial or pd.isna(partial['overlap_with_bin_id']):
                # No overlap partner - this is from a small segment
                skipped_count += 1
                continue
            
            # Check original segment length
            if 'original_seg_length' in partial:
                orig_seg_len = partial['original_seg_length']
                bin_size = partial['bin_actual_size'] if 'bin_actual_size' in partial else 10_000_000
                
                # CRITICAL CHECK: Only aggregate if original segment > bin_size
                if orig_seg_len <= bin_size:
                    print(f"\n  Skipping aggregation for segment {partial['segment_id']}")
                    print(f"    Reason: Original segment length ({orig_seg_len/1e6:.2f} Mb) <= bin_size ({bin_size/1e6:.2f} Mb)")
                    skipped_count += 1
                    continue
            
            # Find corresponding full_right bin
            right_bin_id = int(partial['overlap_with_bin_id'])
            
            # Get the full_right bin
            right_bin = bins_with_dist[bins_with_dist['bin_id'] == right_bin_id]
            
            if len(right_bin) == 0:
                continue
                
            right_bin = right_bin.iloc[0]
            partial_bin_id = partial['bin_id']
            
            print(f"\n  Aggregating pair (segment {partial['segment_id']}, orig_len={partial.get('original_seg_length', 0)/1e6:.2f} Mb):")
            print(f"    Partial (bin_id={partial_bin_id}): {partial['startpos']}-{partial['endpos']} ({partial['bin_actual_size']/1e6:.2f} Mb)")
            print(f"    Right (bin_id={right_bin_id}): {right_bin['startpos']}-{right_bin['endpos']} ({right_bin['bin_actual_size']/1e6:.2f} Mb)")
            
            # Weights based on bin size (length)
            weight_partial = partial['bin_actual_size']
            weight_right = right_bin['bin_actual_size']
            total_weight = weight_partial + weight_right
            
            print(f"    Weights: partial={weight_partial/1e6:.2f} Mb ({weight_partial/total_weight*100:.1f}%), right={weight_right/1e6:.2f} Mb ({weight_right/total_weight*100:.1f}%)")
            
            # ===== WEIGHTED MEAN CALCULATION =====
            aggregated_metrics = {}
            
            # For each quartile ratio, compute weighted mean
            for lab in labels:
                ratio_col = f'ratio_{lab}'
                
                ratio_partial = partial[ratio_col] if pd.notna(partial[ratio_col]) else 0.0
                ratio_right = right_bin[ratio_col] if pd.notna(right_bin[ratio_col]) else 0.0
                
                # Weighted mean by length
                weighted_ratio = (ratio_partial * weight_partial + ratio_right * weight_right) / total_weight
                aggregated_metrics[ratio_col] = weighted_ratio
                
                print(f"      {ratio_col}: ({ratio_partial:.3f} × {weight_partial/1e6:.2f} + {ratio_right:.3f} × {weight_right/1e6:.2f}) / {total_weight/1e6:.2f} = {weighted_ratio:.3f}")
            
            # Weighted n_baf_points
            n_baf_partial = partial['n_baf_points'] if pd.notna(partial['n_baf_points']) else 0
            n_baf_right = right_bin['n_baf_points'] if pd.notna(right_bin['n_baf_points']) else 0
            
            aggregated_n_baf = int((n_baf_partial * weight_partial + n_baf_right * weight_right) / total_weight)
            aggregated_metrics['n_baf_points'] = aggregated_n_baf
            
            # Back-calculate counts from weighted ratios
            for lab in labels:
                ratio_col = f'ratio_{lab}'
                count = int(aggregated_metrics[ratio_col] * aggregated_n_baf)
                aggregated_metrics[lab] = count
            
            # Calculate derived metrics
            outer_mass = aggregated_metrics['ratio_0-0.25'] + aggregated_metrics['ratio_0.75-1.0']
            inner_mass = aggregated_metrics['ratio_0.25-0.5'] + aggregated_metrics['ratio_0.5-0.75']
            
            print(f"    Aggregated n_baf_points (weighted): {aggregated_n_baf}")
            print(f"    Aggregated outer_mass: {aggregated_metrics['ratio_0-0.25']:.3f} + {aggregated_metrics['ratio_0.75-1.0']:.3f} = {outer_mass:.3f}")
            print(f"    Aggregated inner_mass: {aggregated_metrics['ratio_0.25-0.5']:.3f} + {aggregated_metrics['ratio_0.5-0.75']:.3f} = {inner_mass:.3f}")
            
            # UPDATE PARTIAL_REMAINDER BIN with aggregated metrics
            idx = bins_with_dist[bins_with_dist['bin_id'] == partial_bin_id].index[0]
            
            for lab in labels:
                bins_with_dist.at[idx, lab] = aggregated_metrics[lab]
                bins_with_dist.at[idx, f'ratio_{lab}'] = aggregated_metrics[f'ratio_{lab}']
            
            bins_with_dist.at[idx, 'n_baf_points'] = aggregated_metrics['n_baf_points']
            
            # Mark as aggregated
            bins_with_dist.at[idx, 'is_aggregated'] = True
            bins_with_dist.at[idx, 'aggregated_from_bin_id'] = right_bin_id
            bins_with_dist.at[idx, 'aggregation_weight_partial'] = weight_partial
            bins_with_dist.at[idx, 'aggregation_weight_right'] = weight_right
            
            # Update bin_type
            bins_with_dist.at[idx, 'bin_type'] = 'partial_remainder_aggregated'
            
            # Mark FULL_RIGHT bin for removal
            if remove_right_after_aggregation:
                bins_to_remove.append(right_bin_id)
                print(f"    → Keeping partial bin {partial_bin_id} with aggregated metrics")
                print(f"    → Removing full_right bin {right_bin_id}")
            
            aggregated_count += 1
        
        # Remove full_right bins
        if remove_right_after_aggregation and bins_to_remove:
            print(f"\n  Summary:")
            print(f"    - Aggregated: {aggregated_count} pairs")
            print(f"    - Skipped (small segments): {skipped_count}")
            print(f"    - Removing {len(bins_to_remove)} full_right bins")
            bins_with_dist = bins_with_dist[~bins_with_dist['bin_id'].isin(bins_to_remove)].copy()
            bins_with_dist = bins_with_dist.reset_index(drop=True)
        
        print("\n" + "="*70)
    
    # Add aggregation flags if not present
    if 'is_aggregated' not in bins_with_dist.columns:
        bins_with_dist['is_aggregated'] = False
    if 'aggregated_from_bin_id' not in bins_with_dist.columns:
        bins_with_dist['aggregated_from_bin_id'] = np.nan
    
    return hits, bins_with_dist,baf_df



################### temporary functions for tessting 

def plot_allelic_cn_track_with_segment_loh(
    bins_df, 
    dp_filter, 
    n_positions, 
    gap=1_000_000
):
    """
    Plot ALLELIC CN (nMajor/nMinor) with segment-level LoH shown as colored X-axis regions.
    """
    df = bins_df.copy()
    assert_columns(df, ["Chromosome", "startpos", "endpos", "nMajor", "nMinor", "segment_id"], 
                   "plot_allelic_cn_track input")
    
    # Compute segment-level LoH
    segment_loh = (
        df.groupby('segment_id')
        .agg({
            'total_cn': 'first',
            'nMajor': 'first',
            'nMinor': 'first'
        })
    )
    segment_loh['is_segment_loh'] = (
        (segment_loh['total_cn'] > 0) &
        ((segment_loh['nMajor'] == 0) ^ (segment_loh['nMinor'] == 0))
    )
    
    # Map back to bins
    df = df.merge(
        segment_loh[['is_segment_loh']], 
        left_on='segment_id', 
        right_index=True, 
        how='left'
    )
    
    # Chromosome formatting
    df["chr"] = df["Chromosome"].astype(str).str.replace("^chr", "", regex=True)
    
    def chr_key(c):
        if c == "X": return 23
        if c == "Y": return 24
        if c in ("M", "MT"): return 25
        try: return int(c)
        except: return 10**9
    
    chroms = sorted(df["chr"].unique(), key=chr_key)
    chr_len = df.groupby("chr")["endpos"].max().reindex(chroms)
    offsets = (chr_len + gap).cumsum() - (chr_len + gap)
    offset_map = offsets.to_dict()
    
    df["gstart"] = df["startpos"] + df["chr"].map(offset_map)
    df["gend"] = df["endpos"] + df["chr"].map(offset_map)
    df = df.sort_values(["chr", "startpos"]).reset_index(drop=True)
    
    # Create figure with 2 panels
    fig, (ax_cn, ax_baf) = plt.subplots(
        2, 1, figsize=(14, 7), sharex=True,
        gridspec_kw={"height_ratios": [2.5, 1.5], "hspace": 0.05}
    )
    
    # ========== PANEL 1: CN TRACK ==========
    
    # Plot allelic CN
    dodge = 0.08
    maj = df["nMajor"].to_numpy(dtype=float)
    minr = df["nMinor"].to_numpy(dtype=float)
    
    y_major = np.where(maj == 0, 0.0, maj + dodge)
    y_minor = np.where(minr == 0, 0.0, minr - dodge)
    
    seg_major = np.stack([np.c_[df["gstart"], y_major], np.c_[df["gend"], y_major]], axis=1)
    seg_minor = np.stack([np.c_[df["gstart"], y_minor], np.c_[df["gend"], y_minor]], axis=1)
    
    ax_cn.add_collection(LineCollection(seg_major, colors="C0", linewidths=4, alpha=0.9, zorder=2))
    ax_cn.add_collection(LineCollection(seg_minor, colors="C1", linewidths=4, alpha=0.9, zorder=2))
    ax_cn.autoscale()
    ax_cn.set_ylabel("Allelic CN", fontsize=12)
    ax_cn.set_ylim(-0.5, 4)
    ax_cn.set_yticks([0, 1, 2, 3, 4])
    ax_cn.set_title(f"Allelic CN (nMajor/nMinor) | DP ≥ {dp_filter} - {n_positions} BAF positions", 
                    fontsize=13)
    
    # Legend
    ax_cn.legend(handles=[
        plt.Line2D([0], [0], color="C0", lw=4, label="nMajor"),
        plt.Line2D([0], [0], color="C1", lw=4, label="nMinor"),
    ], loc="upper right", fontsize=10)
    
    # ========== PANEL 2: BAF QUARTILES ==========
    
    centers = offsets + chr_len / 2
    ax_baf.set_xticks(centers.values)
    ax_baf.set_xticklabels(centers.index)
    
    bins = [
        ("0-0.25", "ratio_0-0.25"),
        ("0.25-0.5", "ratio_0.25-0.5"),
        ("0.5-0.75", "ratio_0.5-0.75"),
        ("0.75-1.0", "ratio_0.75-1.0"),
    ]
    ax_baf.set_ylim(0, 4)
    ax_baf.set_yticks(np.arange(0.5, 4.0, 1.0))
    ax_baf.set_yticklabels([b[0] for b in bins])
    ax_baf.set_ylabel("BAF bins", fontsize=12)
    ax_baf.set_xlabel("Genome coordinate (bp)", fontsize=12)
    
    for yline in [1, 2, 3]:
        ax_baf.axhline(yline, linewidth=0.8, alpha=0.3, zorder=1)
    
    # Plot BAF quartiles
    for _, r in df.iterrows():
        x0 = r["gstart"]
        w = r["gend"] - r["gstart"]
        for i, (_, ratio_col) in enumerate(bins):
            ratio = r.get(ratio_col, np.nan)
            if not np.isfinite(ratio) or ratio <= 0:
                continue
            ax_baf.add_patch(Rectangle(
                (x0, i), w, 1.0,
                facecolor="C2", edgecolor="none",
                alpha=float(ratio),
                zorder=1
            ))
    
    # ========== COLOR X-AXIS FOR LOH REGIONS ==========
    
    # Get axis limits
    xlim = ax_baf.get_xlim()
    
    # Color the X-axis spine for LoH regions
    for _, r in df.iterrows():
        if r["is_segment_loh"]:
            x0 = r["gstart"]
            w = r["gend"] - r["gstart"]
            
            # Add colored bar below X-axis
            ax_baf.add_patch(Rectangle(
                (x0, -0.3), w, 0.3,
                facecolor="lightgreen", 
                edgecolor="darkgreen",
                linewidth=1,
                clip_on=False,
                transform=ax_baf.transData,
                zorder=10
            ))
    
    # Add legend for LoH coloring
    loh_patch = Rectangle((0, 0), 1, 1, fc="lightgreen", ec="darkgreen", linewidth=1)
    ax_baf.legend(handles=[loh_patch], labels=['LoH segments'], 
                  loc='upper right', fontsize=9)
    
    plt.tight_layout()
    return fig


def plot_total_cn_track_with_segment_loh(
    bins_df, 
    dp_filter, 
    n_positions, 
    gap=1_000_000
):
    """
    Plot TOTAL CN with segment-level LoH shown as colored X-axis regions.
    """
    df = bins_df.copy()
    assert_columns(df, ["Chromosome", "startpos", "endpos", "total_cn", "segment_id"], 
                   "plot_total_cn_track input")
    
    # Compute segment-level LoH (total CN mode)
    segment_loh = (
        df.groupby('segment_id')
        .agg({'total_cn': 'first'})
    )
    segment_loh['is_segment_loh'] = (segment_loh['total_cn'] == 1)
    
    # Map back to bins
    df = df.merge(
        segment_loh[['is_segment_loh']], 
        left_on='segment_id', 
        right_index=True, 
        how='left'
    )
    
    # Chromosome formatting
    df["chr"] = df["Chromosome"].astype(str).str.replace("^chr", "", regex=True)
    
    def chr_key(c):
        if c == "X": return 23
        if c == "Y": return 24
        if c in ("M", "MT"): return 25
        try: return int(c)
        except: return 10**9
    
    chroms = sorted(df["chr"].unique(), key=chr_key)
    chr_len = df.groupby("chr")["endpos"].max().reindex(chroms)
    offsets = (chr_len + gap).cumsum() - (chr_len + gap)
    offset_map = offsets.to_dict()
    
    df["gstart"] = df["startpos"] + df["chr"].map(offset_map)
    df["gend"] = df["endpos"] + df["chr"].map(offset_map)
    df = df.sort_values(["chr", "startpos"]).reset_index(drop=True)
    
    # Create figure with 2 panels
    fig, (ax_cn, ax_baf) = plt.subplots(
        2, 1, figsize=(14, 7), sharex=True,
        gridspec_kw={"height_ratios": [2.5, 1.5], "hspace": 0.05}
    )
    
    # ========== PANEL 1: TOTAL CN TRACK ==========
    
    # Plot TOTAL CN
    y = df["total_cn"].to_numpy(dtype=float)
    segs = np.stack([np.c_[df["gstart"], y], np.c_[df["gend"], y]], axis=1)
    
    ax_cn.add_collection(LineCollection(segs, colors="purple", linewidths=4, alpha=0.9, zorder=2))
    ax_cn.autoscale()
    ax_cn.set_ylabel("Total CN", fontsize=12)
    
    ymax = max(4, int(np.nanmax(y)) if len(y) else 4)
    ax_cn.set_ylim(-0.5, ymax)
    ax_cn.set_yticks(list(range(0, ymax + 1)))
    ax_cn.set_title(f"Total CN (nMajor + nMinor) | DP ≥ {dp_filter} - {n_positions} BAF positions", 
                    fontsize=13)
    
    # Legend
    ax_cn.legend(handles=[plt.Line2D([0], [0], color="purple", lw=4, label="total_cn")],
                 loc="upper right", fontsize=10)
    
    # ========== PANEL 2: BAF QUARTILES ==========
    
    centers = offsets + chr_len / 2
    ax_baf.set_xticks(centers.values)
    ax_baf.set_xticklabels(centers.index)
    
    bins = [
        ("0-0.25", "ratio_0-0.25"),
        ("0.25-0.5", "ratio_0.25-0.5"),
        ("0.5-0.75", "ratio_0.5-0.75"),
        ("0.75-1.0", "ratio_0.75-1.0"),
    ]
    
    ax_baf.set_ylim(0, 4)
    ax_baf.set_yticks(np.arange(0.5, 4.0, 1.0))
    ax_baf.set_yticklabels([b[0] for b in bins])
    ax_baf.set_ylabel("BAF bins", fontsize=12)
    ax_baf.set_xlabel("Genome coordinate (bp)", fontsize=12)
    
    for yline in [1, 2, 3]:
        ax_baf.axhline(yline, linewidth=0.8, alpha=0.3, zorder=1)
    
    # Plot BAF quartiles
    for _, r in df.iterrows():
        x0 = r["gstart"]
        w = r["gend"] - r["gstart"]
        for i, (_, ratio_col) in enumerate(bins):
            ratio = r.get(ratio_col, np.nan)
            if not np.isfinite(ratio) or ratio <= 0:
                continue
            ax_baf.add_patch(Rectangle(
                (x0, i), w, 1.0,
                facecolor="C2", edgecolor="none",
                alpha=float(ratio),
                zorder=1
            ))
    
    # ========== COLOR X-AXIS FOR LOH REGIONS ==========
    
    # Get axis limits
    xlim = ax_baf.get_xlim()
    
    # Color the X-axis spine for LoH regions (different color for total CN mode)
    for _, r in df.iterrows():
        if r["is_segment_loh"]:
            x0 = r["gstart"]
            w = r["gend"] - r["gstart"]
            
            # Add colored bar below X-axis
            ax_baf.add_patch(Rectangle(
                (x0, -0.3), w, 0.3,
                facecolor="lightsalmon", 
                edgecolor="darkorange",
                linewidth=1,
                clip_on=False,
                transform=ax_baf.transData,
                zorder=10
            ))
    
    # Add legend for LoH coloring
    loh_patch = Rectangle((0, 0), 1, 1, fc="lightsalmon", ec="darkorange", linewidth=1)
    ax_baf.legend(handles=[loh_patch], labels=['LoH segments (CN=1)'], 
                  loc='upper right', fontsize=9)
    
    plt.tight_layout()
    return fig



def create_segment_bin_slides(
    prs,
    bins_with_metrics,
    dp_filter,
    n_positions,
    sample_name,
    file_type,
    bin_size_mb,
    outdir="ppt_imgs_segment_bins"
):
    """Create PowerPoint slides - SIMPLIFIED, NO AGGREGATIONS, NO CDF."""
    os.makedirs(outdir, exist_ok=True)
    
    # Validate segment_id exists
    if 'segment_id' not in bins_with_metrics.columns:
        print("  WARNING: segment_id not found, creating from groupby...")
        bins_with_metrics['segment_id'] = (
            bins_with_metrics.groupby(['Chromosome', 'nMajor', 'nMinor'])
            .ngroup()
        )
    
    # Define LOH labels (allelic mode)
    bins_with_metrics["loh_signal_allelic"] = (
        (bins_with_metrics["total_cn"] > 0) &
        ((bins_with_metrics["nMajor"] == 0) ^ (bins_with_metrics["nMinor"] == 0))
    ).astype(int)
    
    # Define LOH labels (total CN mode)
    bins_with_metrics["loh_signal_total"] = (
        bins_with_metrics["total_cn"] == 1
    ).astype(int)
    
    # Use allelic as default for backwards compatibility
    bins_with_metrics["loh_signal"] = bins_with_metrics["loh_signal_allelic"]
    
    # Add segment_len_bp if missing
    if "segment_len_bp" not in bins_with_metrics.columns:
        bins_with_metrics["segment_len_bp"] = (
            bins_with_metrics["endpos"] - bins_with_metrics["startpos"] + 1
        )
    
    prefix = f"{sample_name}_{file_type}_dp{dp_filter}_bin{bin_size_mb}Mb"
    
    # ========== DATA EXPLORATION SLIDES ==========
    print("  Creating Data Exploration Slides...")
    
    # Slide 1: BAF points vs segment length (scatter)
    baf_seg_path = os.path.join(outdir, f"{prefix}_baf_vs_length.png")
    save_fig(
        plot_baf_points_vs_segment_length(bins_with_metrics, sample_name, file_type, dp_filter),
        baf_seg_path
    )
    build_track_only_slide(
        prs,
        title=f"{sample_name} | {file_type} | BAF Points vs Segment Length | DP≥{dp_filter}",
        track_path=baf_seg_path
    )
    
    # Slide 2: Histogram of BAF points per bin - ALL BINS
    baf_hist_path = os.path.join(outdir, f"{prefix}_baf_per_bin_hist_all.png")
    save_fig(
        plot_baf_points_per_bin_histogram(bins_with_metrics, sample_name, file_type, dp_filter),
        baf_hist_path
    )
    build_track_only_slide(
        prs,
        title=f"{sample_name} | {file_type} | BAF Points per Bin (All) | DP≥{dp_filter}",
        track_path=baf_hist_path
    )
    
    # Slide 3: Histogram of BAF points per bin - LoH vs Non-LoH SEPARATE
    baf_loh_sep_path = os.path.join(outdir, f"{prefix}_baf_per_bin_loh_vs_nonloh.png")
    save_fig(
        plot_baf_points_per_bin_loh_vs_nonloh(bins_with_metrics, sample_name, file_type, dp_filter),
        baf_loh_sep_path
    )
    build_track_only_slide(
        prs,
        title=f"{sample_name} | {file_type} | BAF Points: LoH vs Non-LoH Bins | DP≥{dp_filter}",
        track_path=baf_loh_sep_path
    )
    
    # Slide 4: BAF points by bin type
    baf_bintype_path = os.path.join(outdir, f"{prefix}_baf_by_bintype.png")
    save_fig(
        plot_baf_points_by_bin_type(bins_with_metrics, sample_name, file_type, dp_filter),
        baf_bintype_path
    )
    build_track_only_slide(
        prs,
        title=f"{sample_name} | {file_type} | BAF Points by Bin Type | DP≥{dp_filter}",
        track_path=baf_bintype_path
    )
    
    # Slide 5: BAF points by LoH status (allelic) - OVERLAPPED
    baf_loh_allelic_path = os.path.join(outdir, f"{prefix}_baf_by_loh_allelic.png")
    save_fig(
        plot_baf_points_by_loh_status_allelic(bins_with_metrics, sample_name, file_type, dp_filter),
        baf_loh_allelic_path
    )
    build_track_only_slide(
        prs,
        title=f"{sample_name} | {file_type} | BAF Points by LoH Status (Allelic) | DP≥{dp_filter}",
        track_path=baf_loh_allelic_path
    )
    
    # Slide 6: BAF points by LoH status (total) - OVERLAPPED
    baf_loh_total_path = os.path.join(outdir, f"{prefix}_baf_by_loh_total.png")
    save_fig(
        plot_baf_points_by_loh_status_total(bins_with_metrics, sample_name, file_type, dp_filter),
        baf_loh_total_path
    )
    build_track_only_slide(
        prs,
        title=f"{sample_name} | {file_type} | BAF Points by LoH Status (Total CN) | DP≥{dp_filter}",
        track_path=baf_loh_total_path
    )
    
    # ========== SLIDE 7: Allelic Track with COLORED X-AXIS ==========
    print("  Creating Slide: Allelic track with colored X-axis...")
    track_allelic_path = os.path.join(outdir, f"{prefix}_bins_allelic_segment_loh_track.png")
    save_fig(
        plot_allelic_cn_track_with_segment_loh(bins_with_metrics, dp_filter, n_positions),
        track_allelic_path
    )
    
    build_track_only_slide(
        prs,
        title=f"{sample_name} | {file_type} | Allelic CN + Segment LoH | DP≥{dp_filter} | {bin_size_mb}Mb",
        track_path=track_allelic_path
    )
    
    # ========== SLIDE 8: Total CN Track with COLORED X-AXIS ==========
    print("  Creating Slide: Total CN track with colored X-axis...")
    track_total_path = os.path.join(outdir, f"{prefix}_bins_total_segment_loh_track.png")
    save_fig(
        plot_total_cn_track_with_segment_loh(bins_with_metrics, dp_filter, n_positions),
        track_total_path
    )
    
    build_track_only_slide(
        prs,
        title=f"{sample_name} | {file_type} | Total CN + Segment LoH (CN=1) | DP≥{dp_filter} | {bin_size_mb}Mb",
        track_path=track_total_path
    )
    
    # ========== SLIDES 9-10: ROC + PR with optimal threshold ==========
    print("  Creating Slide: ROC + PR for outer_mass...")
    roc_outer_path = os.path.join(outdir, f"{prefix}_outer_mass_roc.png")
    pr_outer_path = os.path.join(outdir, f"{prefix}_outer_mass_pr_optimal.png")
    hist_outer_path = os.path.join(outdir, f"{prefix}_outer_mass_hist.png")
    
    save_fig(plot_roc_curve(bins_with_metrics, "loh_signal", "outer_mass"), roc_outer_path)
    save_fig(plot_pr_curve_with_optimal_threshold(bins_with_metrics, "loh_signal", "outer_mass"), pr_outer_path)
    save_fig(plot_histogram(bins_with_metrics, "loh_signal", "outer_mass"), hist_outer_path)
    
    build_roc_pr_slide(
        prs,
        title=f"{sample_name} | {file_type} | outer_mass Performance | DP≥{dp_filter}",
        roc_path=roc_outer_path,
        pr_path=pr_outer_path,
        hist_path=hist_outer_path
    )
    
    print("  Creating Slide: ROC + PR for inner_mass...")
    roc_inner_path = os.path.join(outdir, f"{prefix}_inner_mass_roc.png")
    pr_inner_path = os.path.join(outdir, f"{prefix}_inner_mass_pr_optimal.png")
    hist_inner_path = os.path.join(outdir, f"{prefix}_inner_mass_hist.png")
    
    save_fig(plot_roc_curve(bins_with_metrics, "loh_signal", "inner_mass"), roc_inner_path)
    save_fig(plot_pr_curve_with_optimal_threshold(bins_with_metrics, "loh_signal", "inner_mass"), pr_inner_path)
    save_fig(plot_histogram(bins_with_metrics, "loh_signal", "inner_mass"), hist_inner_path)
    
    build_roc_pr_slide(
        prs,
        title=f"{sample_name} | {file_type} | inner_mass Performance | DP≥{dp_filter}",
        roc_path=roc_inner_path,
        pr_path=pr_inner_path,
        hist_path=hist_inner_path
    )
    
    # ========== SLIDE 11: Feature Table ==========
    print("  Creating Slide: Feature importance table...")
    feat_table_bins, model_auc_bins = compute_feature_tables(
        bins_with_metrics, 
        label_col="loh_signal"
    )
    
    build_table_slide(
        prs,
        title=f"{sample_name} | {file_type} | Feature Importance (No skew) | DP≥{dp_filter}",
        table_df=feat_table_bins,
        model_auc=model_auc_bins
    )
    
    print("  All slides created successfully!")
    
    return {"bins": bins_with_metrics}

def build_slide_from_figure(prs, title, fig, dpi=75):  # Changed from 150 to 75
    """
    Add matplotlib figure directly to PowerPoint without saving to disk.
    
    Parameters
    ----------
    prs : Presentation
        PowerPoint presentation object
    title : str
        Slide title
    fig : matplotlib.figure.Figure
        Matplotlib figure to add
    dpi : int
        Resolution (default: 75 for smaller files, 150 for higher quality)
    """
    from io import BytesIO
    
    # Create slide
    blank = prs.slide_layouts[6]
    slide = prs.slides.add_slide(blank)
    
    # Add title
    M = 0.35
    title_h = 0.45
    add_title(slide, title, M=M, title_h=title_h)
    
    # Save figure to memory buffer with lower DPI
    buffer = BytesIO()
    fig.savefig(buffer, format='png', dpi=dpi, bbox_inches='tight')
    buffer.seek(0)
    
    # Add image from buffer
    left = Inches(M)
    top = Inches(M + title_h)
    width = Inches(11.69 - 2*M)
    height = Inches(8.27 - (M + title_h) - M)
    
    slide.shapes.add_picture(buffer, left, top, width=width, height=height)
    
    # Close buffer and figure
    buffer.close()
    plt.close(fig)
    
    return slide

def create_roc_pr_hist_figure(df, label_col="loh_signal", score_col="outer_mass"):
    """
    Create 3-panel figure with ROC, PR, and Histogram.
    Returns matplotlib figure (not saved to disk).
    """
    fig, (ax_roc, ax_pr, ax_hist) = plt.subplots(1, 3, figsize=(18, 5))
    
    data = df[[label_col, score_col]].dropna()
    
    if len(data) == 0 or len(data[label_col].unique()) < 2:
        for ax in [ax_roc, ax_pr, ax_hist]:
            ax.text(0.5, 0.5, "Insufficient data", ha='center', va='center', fontsize=14)
        plt.tight_layout()
        return fig
    
    y_true = data[label_col].values
    y_score = data[score_col].values
    
    # Panel 1: ROC
    fpr, tpr, _ = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)
    
    ax_roc.plot(fpr, tpr, color='darkorange', lw=2, label=f'AUC = {roc_auc:.3f}')
    ax_roc.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--', label='Random')
    ax_roc.set_xlim([0.0, 1.0])
    ax_roc.set_ylim([0.0, 1.05])
    ax_roc.set_xlabel('False Positive Rate', fontsize=12)
    ax_roc.set_ylabel('True Positive Rate', fontsize=12)
    ax_roc.set_title(f'ROC Curve - {score_col}', fontsize=13)
    ax_roc.legend(loc="lower right")
    ax_roc.grid(alpha=0.3)
    
    # Panel 2: Precision-Recall
    precision, recall, thresholds = precision_recall_curve(y_true, y_score)
    pr_auc = auc(recall, precision)
    
    f1_scores = 2 * (precision * recall) / (precision + recall + 1e-10)
    optimal_idx = np.argmax(f1_scores)
    optimal_threshold = thresholds[optimal_idx] if optimal_idx < len(thresholds) else thresholds[-1]
    optimal_precision = precision[optimal_idx]
    optimal_recall = recall[optimal_idx]
    optimal_f1 = f1_scores[optimal_idx]
    
    ax_pr.plot(recall, precision, color='blue', lw=2, label=f'AUC = {pr_auc:.3f}')
    ax_pr.plot(optimal_recall, optimal_precision, 'ro', markersize=10, 
              label=f'F1-optimal (t={optimal_threshold:.3f})')
    ax_pr.set_xlim([0.0, 1.0])
    ax_pr.set_ylim([0.0, 1.05])
    ax_pr.set_xlabel('Recall', fontsize=12)
    ax_pr.set_ylabel('Precision', fontsize=12)
    ax_pr.set_title(f'Precision-Recall - {score_col}', fontsize=13)
    ax_pr.legend(loc="lower left")
    ax_pr.grid(alpha=0.3)
    
    # Panel 3: Histogram
    loh = data[data[label_col] == 1][score_col]
    non_loh = data[data[label_col] == 0][score_col]
    
    bins = np.linspace(data[score_col].min(), data[score_col].max(), 30)
    ax_hist.hist(non_loh, bins=bins, alpha=0.5, label='Non-LoH', color='blue', edgecolor='black')
    ax_hist.hist(loh, bins=bins, alpha=0.5, label='LoH', color='red', edgecolor='black')
    ax_hist.set_xlabel(score_col, fontsize=12)
    ax_hist.set_ylabel('Count', fontsize=12)
    ax_hist.set_title(f'Distribution - {score_col}', fontsize=13)
    ax_hist.legend()
    ax_hist.grid(alpha=0.3)
    
    plt.tight_layout()
    return fig


def create_input_data_scarHRD_justOdd(segments_df, sample):
    """
    Infer allele-specific copy numbers (A_cn, B_cn) from total CN and BAF metrics.
    
    Uses only outer_mass and mid_mass:
    - outer_mass: BAF near 0 or 1 (LOH signal)
    - mid_mass: BAF near 0.5 (balanced heterozygosity)
    
    Parameters
    ----------
    segments_df : DataFrame
        Must contain: Chromosome, event_start0, event_end0, totalCN, outer_mass, mid_mass
    sample : str
        Sample identifier
        
    Returns
    -------
    DataFrame with columns: SampleID, Chromosome, Start_position, End_position, 
                           total_cn, A_cn, B_cn, ploidy
    """
    
    df = segments_df.copy()
    
    # ADD THIS CHECK AT THE START
    if len(df) == 0:
        print(f"⚠️  WARNING: No segments for sample {sample}")
        return pd.DataFrame(columns=[
            "SampleID", "Chromosome", "Start_position", "End_position",
            "total_cn", "A_cn", "B_cn", "ploidy"
        ])
    
    # Ensure totalCN is numeric
    df["totalCN"] = pd.to_numeric(df["totalCN"], errors="coerce")
    
    # Create comparison columns (CORRECTED)
    df['Outer>mid_mass'] = df['outer_mass'] > df['mid_mass']
    df['Outer<mid_mass'] = df['outer_mass'] < df['mid_mass']
    df['Outer=mid_mass'] = df['outer_mass'] == df['mid_mass']  
    
    # Find which metric is maximum
    cols = ['outer_mass', 'mid_mass']
    df["max_metric"] = df[cols].replace(0, np.nan).idxmax(axis=1)
    
    # Create masks for different CN types
    is_even = (df["totalCN"] % 2 == 0) & df["totalCN"].notna() #### par
    is_odd = (df["totalCN"] % 2 == 1) & df["totalCN"].notna()  #### impar 
    is_loh1 = (df["totalCN"] == 1)
    is_cn0 = (df["totalCN"] == 0)
    is_odd_gt1 = is_odd & (df["totalCN"] > 1)
    
    # Initialize output columns
    df["A_cn"] = np.nan
    df["B_cn"] = np.nan
    
    # =========================
    # EVEN totalCN
    # =========================
    print(f"Processing EVEN CN segments: {is_even.sum()}")
    
    # outer_mass dominant → Complete LOH (n:0) ##### NEW ---> IF IS EVEN EVERYTHING IS /2
    mask = is_even # & df["max_metric"].eq("outer_mass") & (df["outer_mass"] != 0)
    df.loc[mask, "A_cn"] = df["totalCN"] / 2
    df.loc[mask, "B_cn"] = df["totalCN"] / 2
    #print(f"  EVEN + outer_mass → n:0 : {mask.sum()} segments")
    print(f"  EVEN → n/2 FOR EACH ALLELE : {mask.sum()} segments")
    # mid_mass dominant → Balanced (n/2:n/2)
    #mask = is_even & df["max_metric"].eq("mid_mass") & (df["mid_mass"] != 0)
    #df.loc[mask, "A_cn"] = df.loc[mask, "totalCN"] / 2
    #df.loc[mask, "B_cn"] = df.loc[mask, "totalCN"] / 2
    #print(f"  EVEN + mid_mass → n/2:n/2 : {mask.sum()} segments")
    
    # =========================
    # ODD totalCN == 1 (true LOH)
    # =========================
    print(f"\nProcessing CN=1 segments: {is_loh1.sum()}")
    df.loc[is_loh1, "A_cn"] = 1
    df.loc[is_loh1, "B_cn"] = 0
    
    # =========================
    # ODD totalCN == 0 (deletion)
    # =========================
    print(f"Processing CN=0 segments: {is_cn0.sum()}")
    df.loc[is_cn0, "A_cn"] = 0
    df.loc[is_cn0, "B_cn"] = 0
    
    # =========================
    # ODD totalCN > 1
    # =========================
    print(f"\nProcessing ODD CN > 1 segments: {is_odd_gt1.sum()}")
    
    # outer_mass dominant → LOH with amplification (n:0)
    mask = is_odd_gt1 & df["max_metric"].eq("outer_mass") & (df["outer_mass"] != 0)
    df.loc[mask, "A_cn"] = df.loc[mask, "totalCN"]
    df.loc[mask, "B_cn"] = 0
    print(f"  ODD>1 + outer_mass → n:0 : {mask.sum()} segments")
    
    # mid_mass dominant → As balanced as possible (floor(n/2):ceil(n/2))
    mask = is_odd_gt1 & df["max_metric"].eq("mid_mass") & (df["mid_mass"] != 0)
    df.loc[mask, "A_cn"] = np.floor(df.loc[mask, "totalCN"] / 2)
    df.loc[mask, "B_cn"] = np.ceil(df.loc[mask, "totalCN"] / 2)
    print(f"  ODD>1 + mid_mass → floor(n/2):ceil(n/2) : {mask.sum()} segments")
    
    # ===================================================
    # When max_metric == NaN and totalCN is Even
    # ===================================================
    mask = is_even & df["max_metric"].isna()
    df.loc[mask, "A_cn"] = df.loc[mask, "totalCN"] / 2
    df.loc[mask, "B_cn"] = df.loc[mask, "totalCN"] / 2
    print(f"  Even + No metric → n/2:n/2 : {mask.sum()} segments")

    # ===================================================
    # When max_metric == NaN and totalCN is Odd
    # ===================================================
    mask = is_odd & df["max_metric"].isna() 
    df.loc[mask, "A_cn"] = df.loc[mask, "totalCN"] - 1
    df.loc[mask, "B_cn"] = 1
    print(f"  Odd + no metric → A=n-1; B=1 : {mask.sum()} segments")

    # =========================
    # Final formatting
    # =========================
    df["A_cn"] = df["A_cn"].round().astype("Int64")
    df["B_cn"] = df["B_cn"].round().astype("Int64")
    df["ploidy"] = df["totalCN"].mean()

    # Prepare output columns
    output_df = pd.DataFrame({
        "SampleID": sample,
        "Chromosome": df["Chromosome"],
        "Start_position": df["start"],
        "End_position": df["end"],
        "total_cn": df["totalCN"],
        "A_cn": df["A_cn"],
        "B_cn": df["B_cn"],
        "ploidy": df["ploidy"]
    })
    
    print(f"\nTotal segments processed: {len(output_df)}")
    
    # FIX THE ERROR HERE - Check if dataframe is empty before accessing
    if len(df) > 0:
        print(f"Ploidy: {df['ploidy'].iloc[0]:.2f}")
    else:
        print(f"Ploidy: N/A (no segments)")
    
    return output_df

####### improving feature for embedding analysis
"""
Enhanced feature engineering functions for BAF analysis.
"""

import numpy as np
import pandas as pd
from scipy.stats import skew, kurtosis


def calculate_enhanced_features_per_segment(baf_group):
    """
    Calculate enhanced features for a single segment.
    
    Parameters
    ----------
    baf_group : DataFrame
        BAF positions for one segment (one segment_id, one dp_filter)
        Must contain: BAF_corrected, BAF_raw, DP, REFcounts, ALTcounts
    
    Returns
    -------
    dict : Enhanced features for this segment
    """
    
    if len(baf_group) == 0:
        return {
            # Position count
            'n_baf_points': 0,
            
            # Variance features
            'baf_std': 0,
            'baf_cv': 0,
            'baf_iqr': 0,
            'baf_mad': 0,
            'baf_stability_score': 0,
            
            # Coverage features
            'mean_dp': 0,
            'median_dp': 0,
            'min_dp': 0,
            'max_dp': 0,
            'dp_std': 0,
            'dp_cv': 0,
            'dp_range': 0,
            'dp_iqr': 0,
            'fraction_low_dp': 0,
            'fraction_high_dp': 0,
            'fraction_medium_dp': 0,
            
            # Allelic ratio features
            'mean_allelic_ratio': 0.5,
            'allelic_ratio_std': 0,
            'fraction_extreme_ratio': 0,
            'fraction_balanced_ratio': 0,
            'max_allelic_imbalance': 0,
            'mean_ref_counts': 0,
            'mean_alt_counts': 0,
            'total_ref_counts': 0,
            'total_alt_counts': 0,
            
            # Spatial pattern features
            'baf_autocorrelation': 0,
            'baf_trend': 0,
            'n_transitions': 0,
            'transition_rate': 0,
            'max_consecutive_loh': 0,
            'max_consecutive_het': 0,
            
            # Distribution shape features
            'baf_skewness': 0,
            'baf_kurtosis': 0,
            'bimodality_coefficient': 0,
            'n_peaks': 0,
            
            # Purity correction features
            'mean_correction_magnitude': 0,
            'max_correction_magnitude': 0,
            'correction_consistency': 0,
        }
    
    # Extract data
    baf_corrected = baf_group['BAF_corrected'].values
    baf_raw = baf_group['BAF_raw'].values
    dp_vals = baf_group['DP'].values
    ref_counts = baf_group['REFcounts'].values
    alt_counts = baf_group['ALTcounts'].values
    
    n_positions = len(baf_corrected)
    
    features = {}
    
    # ========================================
    # 1. POSITION COUNT
    # ========================================
    
    features['n_baf_points'] = n_positions
    
    # ========================================
    # 2. BAF VARIANCE FEATURES
    # ========================================
    
    baf_std = np.std(baf_corrected)
    baf_mean = np.mean(baf_corrected)
    baf_median = np.median(baf_corrected)
    baf_cv = baf_std / (baf_mean + 0.01)
    baf_iqr = np.percentile(baf_corrected, 75) - np.percentile(baf_corrected, 25)
    baf_mad = np.median(np.abs(baf_corrected - baf_median))
    baf_stability = 1 / (1 + baf_std)
    
    features['baf_std'] = baf_std
    features['baf_cv'] = baf_cv
    features['baf_iqr'] = baf_iqr
    features['baf_mad'] = baf_mad
    features['baf_stability_score'] = baf_stability
    
    # ========================================
    # 3. COVERAGE FEATURES
    # ========================================
    
    mean_dp = np.mean(dp_vals)
    median_dp = np.median(dp_vals)
    min_dp = np.min(dp_vals)
    max_dp = np.max(dp_vals)
    dp_std = np.std(dp_vals)
    dp_cv = dp_std / (mean_dp + 0.01)
    dp_range = max_dp - min_dp
    dp_iqr = np.percentile(dp_vals, 75) - np.percentile(dp_vals, 25)
    
    fraction_low = (dp_vals < 5).sum() / n_positions
    fraction_high = (dp_vals > 50).sum() / n_positions
    fraction_medium = ((dp_vals >= 10) & (dp_vals <= 30)).sum() / n_positions
    
    features['mean_dp'] = mean_dp
    features['median_dp'] = median_dp
    features['min_dp'] = min_dp
    features['max_dp'] = max_dp
    features['dp_std'] = dp_std
    features['dp_cv'] = dp_cv
    features['dp_range'] = dp_range
    features['dp_iqr'] = dp_iqr
    features['fraction_low_dp'] = fraction_low
    features['fraction_high_dp'] = fraction_high
    features['fraction_medium_dp'] = fraction_medium
    
    # ========================================
    # 4. ALLELIC RATIO FEATURES
    # ========================================
    
    total = ref_counts + alt_counts
    allelic_ratio = np.where(total > 0, alt_counts / total, 0.5)
    
    mean_ratio = np.mean(allelic_ratio)
    ratio_std = np.std(allelic_ratio)
    fraction_extreme = ((allelic_ratio < 0.1) | (allelic_ratio > 0.9)).sum() / n_positions
    fraction_balanced = ((allelic_ratio >= 0.4) & (allelic_ratio <= 0.6)).sum() / n_positions
    max_imbalance = np.max(np.abs(allelic_ratio - 0.5))
    
    mean_ref = np.mean(ref_counts)
    mean_alt = np.mean(alt_counts)
    total_ref = np.sum(ref_counts)
    total_alt = np.sum(alt_counts)
    
    features['mean_allelic_ratio'] = mean_ratio
    features['allelic_ratio_std'] = ratio_std
    features['fraction_extreme_ratio'] = fraction_extreme
    features['fraction_balanced_ratio'] = fraction_balanced
    features['max_allelic_imbalance'] = max_imbalance
    features['mean_ref_counts'] = mean_ref
    features['mean_alt_counts'] = mean_alt
    features['total_ref_counts'] = total_ref
    features['total_alt_counts'] = total_alt
    
    # ========================================
    # 5. SPATIAL PATTERN FEATURES
    # ========================================
    
    if n_positions >= 3:
        # Autocorrelation
        if n_positions > 1:
            autocorr = np.corrcoef(baf_corrected[:-1], baf_corrected[1:])[0, 1]
            autocorr = 0 if np.isnan(autocorr) else autocorr
        else:
            autocorr = 0
        
        # Trend
        positions = np.arange(n_positions)
        if n_positions > 2:
            trend = np.polyfit(positions, baf_corrected, 1)[0]
        else:
            trend = 0
        
        # Transitions
        is_loh = (baf_corrected < 0.2) | (baf_corrected > 0.8)
        is_het = (baf_corrected >= 0.4) & (baf_corrected <= 0.6)
        
        transitions = 0
        for i in range(n_positions - 1):
            if is_loh[i] != is_loh[i+1]:
                transitions += 1
        
        transition_rate = transitions / n_positions
        
        # Max consecutive runs
        max_loh_run = 0
        current_run = 0
        for val in is_loh:
            if val:
                current_run += 1
                max_loh_run = max(max_loh_run, current_run)
            else:
                current_run = 0
        
        max_het_run = 0
        current_run = 0
        for val in is_het:
            if val:
                current_run += 1
                max_het_run = max(max_het_run, current_run)
            else:
                current_run = 0
        
        features['baf_autocorrelation'] = autocorr
        features['baf_trend'] = trend
        features['n_transitions'] = transitions
        features['transition_rate'] = transition_rate
        features['max_consecutive_loh'] = max_loh_run
        features['max_consecutive_het'] = max_het_run
    else:
        features['baf_autocorrelation'] = 0
        features['baf_trend'] = 0
        features['n_transitions'] = 0
        features['transition_rate'] = 0
        features['max_consecutive_loh'] = 0
        features['max_consecutive_het'] = 0
    
    # ========================================
    # 6. DISTRIBUTION SHAPE FEATURES
    # ========================================
    
    if n_positions >= 5:
        baf_skew = skew(baf_corrected)
        baf_kurt = kurtosis(baf_corrected)
        
        # Bimodality coefficient
        n = n_positions
        bimodality = (baf_skew**2 + 1) / (baf_kurt + 3 * ((n-1)**2 / ((n-2)*(n-3))))
        
        # Number of peaks
        hist, _ = np.histogram(baf_corrected, bins=20, range=(0, 1))
        peaks = 0
        for i in range(1, len(hist)-1):
            if hist[i] > hist[i-1] and hist[i] > hist[i+1]:
                if hist[i] > 0.1 * hist.max():
                    peaks += 1
        
        features['baf_skewness'] = baf_skew
        features['baf_kurtosis'] = baf_kurt
        features['bimodality_coefficient'] = bimodality
        features['n_peaks'] = peaks
    else:
        features['baf_skewness'] = 0
        features['baf_kurtosis'] = 0
        features['bimodality_coefficient'] = 0
        features['n_peaks'] = 0
    
    # ========================================
    # 7. PURITY CORRECTION FEATURES
    # ========================================
    
    correction_mag = np.abs(baf_corrected - baf_raw)
    mean_correction = np.mean(correction_mag)
    max_correction = np.max(correction_mag)
    
    correction_direction = baf_corrected - baf_raw
    consistency = np.abs(np.mean(np.sign(correction_direction)))
    
    features['mean_correction_magnitude'] = mean_correction
    features['max_correction_magnitude'] = max_correction
    features['correction_consistency'] = consistency
    
    return features


def add_enhanced_features_to_merged(merged_df, baf_temp):
    """
    Add enhanced features to merged dataframe.
    
    Parameters
    ----------
    merged_df : DataFrame
        Merged segments dataframe
    baf_temp : DataFrame
        Position-level BAF data
    
    Returns
    -------
    DataFrame with enhanced features
    """
    
    print("  ⭐ Adding enhanced features...")
    
    enhanced_features_list = []
    
    # Group by segment_id and dp_filter
    grouped = baf_temp.groupby(['segment_id', 'dp_filter'])
    total_groups = len(grouped)
    
    for idx, ((seg_id, dp_filter), group) in enumerate(grouped, 1):
        if idx % 1000 == 0:
            print(f"    Progress: {idx}/{total_groups} segments...")
        
        # Calculate features (only 1 argument!)
        features = calculate_enhanced_features_per_segment(group)
        
        # Add identifiers
        features['segment_id'] = seg_id
        features['dp_filter'] = dp_filter
        
        enhanced_features_list.append(features)
    
    # Convert to DataFrame
    enhanced_df = pd.DataFrame(enhanced_features_list)
    
    # Merge with main dataframe
    merged_enhanced = merged_df.merge(
        enhanced_df,
        on=['segment_id', 'dp_filter'],
        how='left'
    )
    
    n_new_features = len(enhanced_df.columns) - 2
    print(f"  ✓ Added {n_new_features} enhanced features")
    
    return merged_enhanced

def add_length_normalized_features(df):
    """
    Add length-normalized features.
    
    Parameters
    ----------
    df : DataFrame
        Segments dataframe with n_baf_points (or n_positions) and segment_len_bp
    
    Returns
    -------
    DataFrame with normalized features
    """
    
    print("  ⭐ Adding length-normalized features...")
    
    df = df.copy()
    
    # ========================================
    # HANDLE COLUMN NAME VARIATIONS
    # ========================================
    
    # Check which column exists for position count
    if 'n_baf_points' in df.columns:
        n_positions_col = 'n_baf_points'
    elif 'n_positions' in df.columns:
        n_positions_col = 'n_positions'
    else:
        print("  ✗ ERROR: Neither 'n_baf_points' nor 'n_positions' column found!")
        print(f"  Available columns: {df.columns.tolist()}")
        raise KeyError("Missing position count column")
    
    print(f"  Using position count column: '{n_positions_col}'")
    
    # ========================================
    # SEGMENT LENGTH CONVERSIONS
    # ========================================
    
    df['segment_len_kb'] = df['segment_len_bp'] / 1000
    df['segment_len_mb'] = df['segment_len_bp'] / 1_000_000
    
    # ========================================
    # POSITION DENSITY ⭐ KEY FEATURE
    # ========================================
    
    df['positions_per_kb'] = df[n_positions_col] / (df['segment_len_kb'] + 0.001)
    df['positions_per_mb'] = df[n_positions_col] / (df['segment_len_mb'] + 0.001)
    df['positions_per_bp'] = df[n_positions_col] / (df['segment_len_bp'] + 1)
    
    # ========================================
    # READ DENSITY
    # ========================================
    
    if 'mean_dp' in df.columns:
        df['estimated_total_reads'] = df['mean_dp'] * df[n_positions_col]
        df['reads_per_kb'] = df['estimated_total_reads'] / (df['segment_len_kb'] + 0.001)
        df['reads_per_mb'] = df['estimated_total_reads'] / (df['segment_len_mb'] + 0.001)
    
    # ========================================
    # LOG TRANSFORMS
    # ========================================
    
    df['log_segment_len'] = np.log10(df['segment_len_bp'] + 1)
    df['log_n_positions'] = np.log10(df[n_positions_col] + 1)
    
    if 'mean_dp' in df.columns:
        df['log_mean_dp'] = np.log10(df['mean_dp'] + 1)
    
    # ========================================
    # COVERAGE UNIFORMITY
    # ========================================
    
    EXPECTED_DENSITY = 0.1
    df['expected_positions'] = df['segment_len_bp'] * EXPECTED_DENSITY
    df['coverage_ratio'] = df[n_positions_col] / (df['expected_positions'] + 1)
    
    # ========================================
    # CONFIDENCE SCORE
    # ========================================
    
    density_score = np.clip(df['positions_per_kb'] / 5.0, 0, 1)
    absolute_score = np.clip(df[n_positions_col] / 10.0, 0, 1)
    
    if 'mean_dp' in df.columns:
        dp_score = np.clip(df['mean_dp'] / 20.0, 0, 1)
    else:
        dp_score = 1.0
    
    length_kb = df['segment_len_kb']
    length_score = np.where(
        length_kb < 1,
        length_kb,
        np.where(length_kb > 1000, 1000 / length_kb, 1.0)
    )
    
    df['confidence_score'] = (density_score * absolute_score * dp_score * length_score) ** 0.25
    
    # ========================================
    # SEGMENT SIZE CATEGORIES
    # ========================================
    
    df['segment_size_category'] = pd.cut(
        df['segment_len_kb'],
        bins=[0, 1, 10, 100, 1000, np.inf],
        labels=['tiny', 'small', 'medium', 'large', 'huge']
    )
    
    print(f"  ✓ Added length-normalized features")
    
    return df