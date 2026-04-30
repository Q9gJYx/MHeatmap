# %% [markdown]
# # Using mheatmap for visualizing unsupervised cell type clustering on Tabula Sapiens datasets
# Generates visualizations for all available datasets with three heatmap types: raw, mosaic, and spectral

# %%
# Additional dependencies: uv pip install scanpy

# %%
import scanpy as sc

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.sparse import issparse
import seaborn as sns
from mheatmap import mosaic_heatmap, spectral_permute

# %%
# Configuration
SCRIPT_DIR = Path(__file__).parent
DATA_DIR = Path("D:/Projects/Tabula_Sapien/data")

# All available Tabula Sapiens datasets
# The datasets can be downloaded from: https://figshare.com/articles/dataset/Tabula_Sapiens_v2/27921984
DATASETS = {
    "bone_marrow": DATA_DIR / "TS_Bone_Marrow_v1.h5ad",
    "vasculature": DATA_DIR / "TS_Vasculature.h5ad",
    "stomach": DATA_DIR / "Stomach_TSP1_30_version2d_10X_smartseq_scvi_Nov122024.h5ad",
    "spleen": DATA_DIR / "Spleen_TSP1_30_version2d_10X_smartseq_scvi_Nov122024.h5ad",
}


# %%
# Utility functions for automatic param finding (from bio.py)
def compute_modularity_components(A, labels):
    """Compute f1 (intra-cluster edge density) and f2 (degree concentration)"""
    if issparse(A):
        A = A.tocsr()
    m = A.sum() / 2.0
    if m == 0:
        return 0.0, 1.0
    k = np.array(A.sum(axis=1)).flatten()
    unique_clusters = np.unique(labels)
    e_cc_values, a_c2_values = [], []
    for cid in unique_clusters:
        idx = np.where(labels == cid)[0]
        if len(idx) == 0:
            continue
        A_c = A[idx, :][:, idx]
        e_cc = A_c.sum() / (2.0 * m)
        e_cc_values.append(e_cc)
        a_c = k[idx].sum() / (2.0 * m)
        a_c2_values.append(a_c**2)
    return float(np.sum(e_cc_values)), float(np.sum(a_c2_values))


def evaluate_resolution(adata, A, resolution, key_prefix="leiden_res", random_state=42):
    """Run Leiden and compute f1/f2 for a given resolution"""
    key_name = f"{key_prefix}_{resolution:.6f}"
    sc.tl.leiden(
        adata,
        adjacency=A,
        resolution=resolution,
        key_added=key_name,
        random_state=random_state,
        flavor="igraph",
    )
    labels = adata.obs[key_name].astype(int).values
    n_clusters = np.unique(labels).size
    f1, f2 = compute_modularity_components(A, labels)
    return {
        "resolution": resolution,
        "f1": f1,
        "f2": f2,
        "n_clusters": int(n_clusters),
        "labels": labels,
        "key_name": key_name,
    }


def evaluate_optim(coords):
    """Evaluate distance from diagonal for knee point detection"""
    start = np.array([0.0, 0.0])
    end = np.array([1.0, 1.0])
    line_vec = end - start
    line_len = np.linalg.norm(line_vec)
    line_unitvec = line_vec / line_len
    vec_to_point = coords - start
    proj_length = np.dot(vec_to_point, line_unitvec)
    proj_point = start + proj_length * line_unitvec
    dist = np.linalg.norm(coords - proj_point)
    return dist


def evaluate_improvement(point_D, point_A, point_B, tolerance=1e-4):
    """Evaluate if the new point (D) expands the Pareto front"""
    f1_A, f2_A = point_A["f1"], point_A["f2"]
    f1_B, f2_B = point_B["f1"], point_B["f2"]
    f1_D, f2_D = point_D["f1"], point_D["f2"]
    delta_f1 = f1_B - f1_A
    delta_f2 = f2_B - f2_A
    if abs(delta_f1) < 1e-10 and abs(delta_f2) < 1e-10:
        return False, 0.0
    a = delta_f2
    b = -delta_f1
    c = delta_f1 * f2_A - delta_f2 * f1_A
    distance = abs(a * f1_D + b * f2_D + c) / np.sqrt(a**2 + b**2)
    is_between = min(f1_A, f1_B) <= f1_D <= max(f1_A, f1_B)
    is_better = (distance > tolerance) and is_between
    return is_better


def search_pareto_front(
    adata,
    A,
    delta_min,
    delta_max,
    max_depth=5,
    tolerance=1e-4,
    random_state=42,
    verbose=False,
):
    """Complete Pareto front search algorithm"""
    all_points = []
    point_min = evaluate_resolution(adata, A, delta_min, random_state)
    point_max = evaluate_resolution(adata, A, delta_max, random_state)
    all_points.extend([point_min, point_max])

    def recursive_search(point_A, point_B, depth=0):
        if depth >= max_depth:
            return []
        if abs(point_B["resolution"] - point_A["resolution"]) < tolerance:
            return []
        delta_new = (point_B["resolution"] + point_A["resolution"]) / 2
        try:
            point_new = evaluate_resolution(adata, A, delta_new, random_state)
        except Exception:
            return []
        is_better = evaluate_improvement(point_new, point_A, point_B, tolerance)
        new_points = []
        if is_better:
            new_points.append(point_new)
            new_points.extend(recursive_search(point_A, point_new, depth + 1))
            new_points.extend(recursive_search(point_new, point_B, depth + 1))
        return new_points

    new_points = recursive_search(point_min, point_max, depth=0)
    all_points.extend(new_points)
    return sorted(all_points, key=lambda p: p["resolution"])


def find_knee_point(pareto_points):
    """Find the knee point on the Pareto front using distance from diagonal"""
    f1_arr = np.array([p["f1"] for p in pareto_points])
    f2_arr = np.array([p["f2"] for p in pareto_points])
    coords = np.column_stack([f1_arr, f2_arr])
    start = np.array([0.0, 0.0])
    end = np.array([1.0, 1.0])
    line_vec = end - start
    line_len = np.linalg.norm(line_vec)
    line_unitvec = line_vec / line_len
    distances = [evaluate_optim(point) for point in coords]
    knee_idx = np.argmax(distances)
    return pareto_points[knee_idx]


def plot_raw_heatmap(cmat, output_path):
    """Plot a standard matplotlib heatmap without reordering or mosaic layout"""
    fig, ax = plt.subplots(figsize=(cmat.shape[1] * 0.5, cmat.shape[0] * 0.4))
    sns.heatmap(
        cmat,
        annot=False,
        cmap="Blues",
        ax=ax,
        cbar=False,
        mask=(cmat == 0),
        xticklabels=False,
        yticklabels=False,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor("black")
        spine.set_linewidth(2)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {output_path.name}")


def plot_mosaic_heatmap(cmat, output_path):
    """Plot using mheatmap's mosaic_heatmap"""
    ax = mosaic_heatmap(
        cmat,
        annot=False,
        cmap="Blues",
        cbar=False,
        mask=(cmat == 0),
        xticklabels=False,
        yticklabels=False,
    )
    plt.gcf().set_size_inches(cmat.shape[1] * 0.8, cmat.shape[0] * 0.5)
    ax.set_xlabel("")
    ax.set_ylabel("")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor("black")
        spine.set_linewidth(2)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {output_path.name}")


def plot_spectral_heatmap(cmat, output_path):
    """Plot spectrally reordered confusion matrix using mheatmap"""
    row_labels = np.asarray(cmat.index)
    col_labels = np.asarray(cmat.columns)
    result = spectral_permute(
        cmat.to_numpy(), row_labels, mode="tw", col_labels=col_labels
    )
    reordered_cmat = pd.DataFrame(
        result["reordered_matrix"],
        index=result["reordered_row_labels"],
        columns=result["reordered_col_labels"],
    )
    ax = mosaic_heatmap(
        reordered_cmat,
        annot=False,
        cmap="Blues",
        cbar=False,
        mask=(reordered_cmat == 0),
        xticklabels=False,
        yticklabels=False,
    )
    plt.gcf().set_size_inches(
        reordered_cmat.shape[1] * 0.8, reordered_cmat.shape[0] * 0.5
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor("black")
        spine.set_linewidth(2)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {output_path.name}")


def plot_spectral_heatmap_plain(cmat, output_path):
    """Plot spectrally reordered confusion matrix using standard seaborn heatmap (no mosaic layout)"""
    row_labels = np.asarray(cmat.index)
    col_labels = np.asarray(cmat.columns)
    result = spectral_permute(
        cmat.to_numpy(), row_labels, mode="tw", col_labels=col_labels
    )
    reordered_cmat = pd.DataFrame(
        result["reordered_matrix"],
        index=result["reordered_row_labels"],
        columns=result["reordered_col_labels"],
    )
    fig, ax = plt.subplots(
        figsize=(reordered_cmat.shape[1] * 0.5, reordered_cmat.shape[0] * 0.4)
    )
    sns.heatmap(
        reordered_cmat,
        annot=False,
        cmap="Blues",
        ax=ax,
        cbar=False,
        mask=(reordered_cmat == 0),
        xticklabels=False,
        yticklabels=False,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor("black")
        spine.set_linewidth(2)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {output_path.name}")


def plot_combined_heatmap(cmat, output_path, dataset_name):
    """Stack all 4 heatmap visualizations: Original, TW, Mosaic, MHeatmap (TW+Mosaic)"""
    fig, axes = plt.subplots(4, 1, figsize=(cmat.shape[1] * 0.8, cmat.shape[0] * 2.5))
    fig.suptitle(dataset_name, fontsize=24, fontweight="bold", y=1.02)

    # 1. Raw heatmap (Original)
    sns.heatmap(
        cmat,
        annot=False,
        cmap="Blues",
        ax=axes[0],
        cbar=False,
        mask=(cmat == 0),
        xticklabels=False,
        yticklabels=False,
    )
    axes[0].set_title("Original", fontsize=36)
    axes[0].set_xlabel("")
    axes[0].set_ylabel("")
    for spine in axes[0].spines.values():
        spine.set_visible(True)
        spine.set_edgecolor("black")
        spine.set_linewidth(2)

    # 2. Spectral plain heatmap (TW)
    row_labels = np.asarray(cmat.index)
    col_labels = np.asarray(cmat.columns)
    result_tw = spectral_permute(
        cmat.to_numpy(), row_labels, mode="tw", col_labels=col_labels
    )
    tw_cmat = pd.DataFrame(
        result_tw["reordered_matrix"],
        index=result_tw["reordered_row_labels"],
        columns=result_tw["reordered_col_labels"],
    )
    sns.heatmap(
        tw_cmat,
        annot=False,
        cmap="Blues",
        ax=axes[1],
        cbar=False,
        mask=(tw_cmat == 0),
        xticklabels=False,
        yticklabels=False,
    )
    axes[1].set_title("TW", fontsize=36)
    axes[1].set_xlabel("")
    axes[1].set_ylabel("")
    for spine in axes[1].spines.values():
        spine.set_visible(True)
        spine.set_edgecolor("black")
        spine.set_linewidth(2)

    # 3. Mosaic heatmap (no reordering)
    mosaic_heatmap(
        cmat,
        annot=False,
        cmap="Blues",
        ax=axes[2],
        cbar=False,
        mask=(cmat == 0),
        xticklabels=False,
        yticklabels=False,
    )
    axes[2].set_title("Mosaic", fontsize=36)
    axes[2].set_xlabel("")
    axes[2].set_ylabel("")
    for spine in axes[2].spines.values():
        spine.set_visible(True)
        spine.set_edgecolor("black")
        spine.set_linewidth(2)

    # 4. MHeatmap (TW + Mosaic)
    ax = mosaic_heatmap(
        tw_cmat,
        annot=False,
        cmap="Blues",
        ax=axes[3],
        cbar=False,
        mask=(tw_cmat == 0),
        xticklabels=False,
        yticklabels=False,
    )
    axes[3].set_title("MHeatmap (TW+Mosaic)", fontsize=36)
    axes[3].set_xlabel("")
    axes[3].set_ylabel("")
    for spine in axes[3].spines.values():
        spine.set_visible(True)
        spine.set_edgecolor("black")
        spine.set_linewidth(2)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {output_path.name}")


def plot_pareto_front(candidates, optim, output_path, dataset_name):
    """Plot the Pareto front with knee point"""
    fig, ax = plt.subplots(figsize=(6, 6))
    f1_arr = [p["f1"] for p in candidates]
    f2_arr = [p["f2"] for p in candidates]
    ax.plot(f1_arr, f2_arr, "b-o", markersize=4)
    ax.plot(
        optim["f1"],
        optim["f2"],
        "r*",
        markersize=15,
        label=f"Knee point (res={optim['resolution']:.3f})",
    )
    ax.set_xlabel("f1 (intra-cluster edge density)")
    ax.set_ylabel("f2 (degree concentration)")
    ax.set_title(f"Pareto Front: {dataset_name}")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {output_path.name}")


def plot_umap_comparison(adata, gt_col, output_path, dataset_name):
    """Plot UMAP colored by Leiden clusters and ground truth"""
    if "X_umap" not in adata.obsm:
        print("No UMAP embedding available")
        return
    print("Plotting UMAP comparison...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    umap_df = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"])
    umap_df["leiden"] = adata.obs["leiden_clusters"].values
    for cluster in sorted(umap_df["leiden"].unique()):
        mask = umap_df["leiden"] == cluster
        axes[0].scatter(
            umap_df.loc[mask, "UMAP1"],
            umap_df.loc[mask, "UMAP2"],
            label=cluster,
            alpha=0.6,
            s=2,
        )
    axes[0].set_title("UMAP colored by Leiden clusters")
    axes[0].legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=6, ncol=2)
    axes[0].set_xlabel("UMAP1")
    axes[0].set_ylabel("UMAP2")
    umap_df["cell_type"] = adata.obs[gt_col].values
    for cell_type in sorted(umap_df["cell_type"].unique()):
        mask = umap_df["cell_type"] == cell_type
        axes[1].scatter(
            umap_df.loc[mask, "UMAP1"],
            umap_df.loc[mask, "UMAP2"],
            label=cell_type,
            alpha=0.6,
            s=2,
        )
    axes[1].set_title(f"UMAP colored by ground truth ({gt_col})")
    axes[1].legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=6, ncol=2)
    axes[1].set_xlabel("UMAP1")
    axes[1].set_ylabel("UMAP2")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {output_path.name}")


# %%
# Process each dataset
for dataset_name, data_path in DATASETS.items():
    print(f"\n{'='*60}")
    print(f"Processing: {dataset_name}")
    print(f"{'='*60}")

    output_dir = SCRIPT_DIR / dataset_name
    output_dir.mkdir(exist_ok=True)

    # Load data
    print(f"Loading data from {data_path}...")
    try:
        adata = sc.read_h5ad(data_path)
    except Exception as e:
        print(f"Error loading {data_path}: {e}")
        continue
    print(f"Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

    # Ensure neighbors graph exists
    if "neighbors" not in adata.uns or "connectivities" not in adata.obsp:
        print("Recomputing neighbors graph...")
        if "X_pca" not in adata.obsm:
            print("Computing PCA...")
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            sc.pp.pca(adata, n_comps=50)
        sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15)
    else:
        print("Using precomputed neighbors graph")

    A = adata.obsp["connectivities"]
    print(f"Adjacency matrix shape: {A.shape}")

    # Run Leiden with automatic param finding via Pareto front
    print("Running Pareto front search for optimal resolution...")
    candidates = search_pareto_front(adata, A, 0.5, 1.5, 10)
    print(f"Evaluated {len(candidates)} resolution candidates")
    optim = find_knee_point(candidates)
    print(
        f"Optimal resolution: {optim['resolution']:.4f} -> {optim['n_clusters']} clusters"
    )

    sc.tl.leiden(
        adata,
        adjacency=A,
        key_added="leiden_clusters",
        random_state=42,
        flavor="igraph",
        resolution=optim["resolution"],
    )
    print(
        f"Leiden clustering complete: {adata.obs['leiden_clusters'].nunique()} clusters"
    )

    # Ground truth comparison
    gt_col = "cell_ontology_class"
    print(f"Ground truth classes: {adata.obs[gt_col].nunique()}")
    print(f"Ground truth distribution:\n{adata.obs[gt_col].value_counts()}")

    # Build confusion matrix
    cmat = pd.crosstab(adata.obs[gt_col], adata.obs["leiden_clusters"])
    print(f"Confusion matrix shape: {cmat.shape}")

    # Plot all heatmap types
    print("\nGenerating visualizations...")

    # 1. Raw heatmap (standard matplotlib, no reordering)
    plot_raw_heatmap(
        cmat,
        output_dir / "confusion_matrix_raw.png",
    )

    # 2. Mosaic heatmap (mheatmap's proportional layout)
    plot_mosaic_heatmap(
        cmat,
        output_dir / "confusion_matrix_mosaic.png",
    )

    # 3. Spectral reordered heatmap (mheatmap's spectral_permute)
    plot_spectral_heatmap(
        cmat,
        output_dir / "confusion_matrix_spectral.png",
    )

    # 3b. Spectral reordered heatmap (plain seaborn, no mosaic layout)
    plot_spectral_heatmap_plain(
        cmat,
        output_dir / "confusion_matrix_spectral_plain.png",
    )

    # 4b. Combined 4-panel heatmap visualization
    plot_combined_heatmap(
        cmat, output_dir / "confusion_matrix_combined.png", dataset_name=dataset_name
    )

    # 5. Pareto front plot
    plot_pareto_front(candidates, optim, output_dir / "pareto_front.png", dataset_name)

    # 6. UMAP comparison
    plot_umap_comparison(
        adata, gt_col, output_dir / "umap_comparison.png", dataset_name
    )

    print(f"\nCompleted {dataset_name}. Output saved to: {output_dir}")

print("\n" + "=" * 60)
print("=== All Tabula Sapiens Datasets Processing Complete ===")
print("=" * 60)
