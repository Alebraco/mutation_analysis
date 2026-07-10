"""
Plotting utilities for mutation analysis pipeline.

Includes:
- Bubble plot (mutation summary)
- Mutation spectrum plot
- Time trajectory plot
- Allele Freq Plot by Group
- Parallel Mutation Plot
- Zoomed genome plot (Kosterlitz-style)
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import to_hex, to_rgb

import os
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch, Rectangle, Polygon


def _validate_columns(df, required_cols):
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns: {missing}. "
            f"Available columns are: {list(df.columns)}"
        )


from .utils import parse_line_label as _parse_line_label
from .utils import get_strain_columns


def _sort_day_labels(day_values):
    def key_func(x):
        m = re.search(r"(\d+)", str(x))
        return int(m.group(1)) if m else 999999

    return sorted(day_values, key=key_func)


def _adjust_color_lightness(color, factor):
    rgb = np.array(to_rgb(color))
    if factor >= 1:
        adjusted = rgb + (1 - rgb) * (factor - 1)
    else:
        adjusted = rgb * factor
    adjusted = np.clip(adjusted, 0, 1)
    return to_hex(adjusted)


def _build_color_map(groups, days, base_colors=None):
    if base_colors is None:
        default_palette = ["#70AD47", "#8E63CE", "#ED7D31", "#5B9BD5", "#C0504D", "#4BACC6"]
        base_colors = {g: default_palette[i % len(default_palette)] for i, g in enumerate(groups)}

    n_days = len(days)
    if n_days == 1:
        factors = [1.0]
    elif n_days == 2:
        factors = [1.25, 0.8]
    else:
        factors = np.linspace(1.35, 0.65, n_days)

    color_map = {}
    for g in groups:
        for d, f in zip(days, factors):
            color_map[(g, d)] = _adjust_color_lightness(base_colors[g], float(f))
    return color_map


def _nice_round(x):
    x = float(x)
    if x <= 0:
        return 0
    if x < 10:
        return int(round(x))
    if x < 50:
        return int(round(x / 5) * 5)
    if x < 200:
        return int(round(x / 10) * 10)
    if x < 1000:
        return int(round(x / 50) * 50)
    return int(round(x / 100) * 100)


def _get_size_legend_values(values, n_levels=3):
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]

    if arr.size == 0:
        return [1, 2, 3]

    quantiles = np.linspace(0.2, 0.8, n_levels)
    vals = np.quantile(arr, quantiles)
    vals = sorted(set(_nice_round(v) for v in vals if _nice_round(v) > 0))

    if len(vals) < 2:
        vals = sorted(set([
            _nice_round(np.min(arr)),
            _nice_round(np.median(arr)),
            _nice_round(np.max(arr))
        ]))
        vals = [v for v in vals if v > 0]

    return vals if vals else [1, 2, 3]


def _compute_size_scale(values, target_max_area=700.0):
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr) & (arr > 0)]

    if arr.size == 0:
        return 1.0

    vmax = np.max(arr)
    if vmax == 0:
        return 1.0

    return target_max_area / vmax


def plot_mutation_bubble(
    input_file,
    output_file=None,
    x_col="Nonsynonymous",
    y_col="Average Frequency",
    size_col="Total Mutations",
    line_col="Line",
    title="Mutation Frequency and Nonsynonymous Proportion",
    panel_by_group=False,
    base_colors=None,
    add_reference_line=True,
    reference_x=2 / 3,
    figsize=(12, 7),
    dpi=200,
    show=True,
):
    """
    Plot mutation summary bubble plot from a CSV/TSV/Excel file.

    Parameters
    ----------
    input_file : str
        Path to mutation summary file.
    output_file : str, optional
        If provided, save plot to this path.
    x_col, y_col, size_col, line_col : str
        Column names in the mutation summary file.
    panel_by_group : bool
        If True, create one panel per group.
    """

    input_path = Path(input_file)
    suffix = input_path.suffix.lower()

    if suffix == ".csv":
        df = pd.read_csv(input_path)
    elif suffix == ".tsv":
        df = pd.read_csv(input_path, sep="\t")
    elif suffix in {".xlsx", ".xls"}:
        df = pd.read_excel(input_path)
    else:
        raise ValueError(f"Unsupported input file format: {suffix}")

    _validate_columns(df, [line_col, x_col, y_col, size_col])

    parsed = df[line_col].apply(lambda x: pd.Series(_parse_line_label(x)))
    parsed.columns = ["Day", "Group", "Replicate"]
    df[["Day", "Group", "Replicate"]] = parsed

    df = df.dropna(subset=[x_col, y_col, size_col]).copy()

    if df["Group"].notna().any():
        groups = list(dict.fromkeys(df["Group"].dropna().astype(str)))
    else:
        groups = ["All"]
        df["Group"] = "All"

    if df["Day"].notna().any():
        days = _sort_day_labels(list(dict.fromkeys(df["Day"].dropna().astype(str))))
    else:
        days = ["All"]
        df["Day"] = "All"

    color_map = _build_color_map(groups, days, base_colors=base_colors)
    df["point_color"] = df.apply(
        lambda row: color_map.get((str(row["Group"]), str(row["Day"])), "#A5A5A5"),
        axis=1,
    )

    size_scale = _compute_size_scale(df[size_col])
    df["plot_size"] = df[size_col].astype(float) * size_scale
    size_values = _get_size_legend_values(df[size_col])

    if panel_by_group:
        n_panels = len(groups)
        fig, axes = plt.subplots(
            1, n_panels,
            figsize=figsize,
            dpi=dpi,
            sharey=True if n_panels > 1 else False
        )
        if n_panels == 1:
            axes = [axes]

        for ax, group in zip(axes, groups):
            sub = df[df["Group"].astype(str) == str(group)]

            ax.scatter(
                sub[x_col],
                sub[y_col],
                s=sub["plot_size"],
                c=sub["point_color"],
                alpha=0.75,
                edgecolors="white",
                linewidths=1.0,
            )

            if add_reference_line:
                ax.axvline(x=reference_x, color="gray", linestyle="--", linewidth=1)

            ax.set_title(str(group), fontsize=13)
            ax.set_xlabel("Proportion of Nonsynonymous Mutations")
            ax.grid(True, linestyle="--", alpha=0.25)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

        axes[0].set_ylabel("Average Mutation Frequency")

        color_handles = []
        for g in groups:
            for d in days:
                key = (g, d)
                if key in color_map:
                    color_handles.append(
                        Line2D(
                            [0], [0],
                            marker="o",
                            color="w",
                            label=f"{g} {d}",
                            markerfacecolor=color_map[key],
                            markersize=8,
                        )
                    )

        fig.legend(
            handles=color_handles,
            title="Population / Day",
            loc="center left",
            bbox_to_anchor=(0.92, 0.7),
            frameon=False,
        )

        size_handles = [
            plt.scatter([], [], s=v * size_scale, color="gray", alpha=0.5, label=f"{int(v)}")
            for v in size_values
        ]
        fig.legend(
            handles=size_handles,
            title=size_col,
            loc="center left",
            bbox_to_anchor=(0.92, 0.3),
            frameon=False,
        )

        fig.suptitle(title, y=1.02, fontsize=14)
        plt.tight_layout(rect=[0, 0, 0.88, 1])

    else:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

        ax.scatter(
            df[x_col],
            df[y_col],
            s=df["plot_size"],
            c=df["point_color"],
            alpha=0.75,
            edgecolors="white",
            linewidths=1.0,
        )

        if add_reference_line:
            ax.axvline(x=reference_x, color="gray", linestyle="--", linewidth=1)

        ax.set_xlabel("Proportion of Nonsynonymous Mutations", fontsize=12)
        ax.set_ylabel("Average Mutation Frequency", fontsize=12)
        ax.set_title(title, fontsize=14, pad=12)

        ax.grid(True, linestyle="--", alpha=0.25)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        color_handles = []
        for g in groups:
            for d in days:
                key = (g, d)
                if key in color_map:
                    color_handles.append(
                        Line2D(
                            [0], [0],
                            marker="o",
                            color="w",
                            label=f"{g} {d}",
                            markerfacecolor=color_map[key],
                            markersize=8,
                        )
                    )

        leg1 = ax.legend(
            handles=color_handles,
            title="Population / Day",
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
            frameon=False,
        )
        ax.add_artist(leg1)

        size_handles = [
            plt.scatter([], [], s=v * size_scale, color="gray", alpha=0.5, label=f"{int(v)}")
            for v in size_values
        ]
        ax.legend(
            handles=size_handles,
            title=size_col,
            loc="lower left",
            bbox_to_anchor=(1.02, 0.0),
            frameon=False,
        )

        plt.tight_layout()

    if output_file is not None:
        fig.savefig(output_file, dpi=300, bbox_inches="tight")

    if show:
        plt.show()

    return fig


def plot_mutation_spectrum(
    input_file,
    output_file=None,
    line_col="Line",
    figsize=(10, 6),
    dpi=200,
    show=True,
):
    """
    Plot mutation spectrum (stacked bar plot).
    Aggregates replicates by Group + Day.
    """

    from pathlib import Path

    input_path = Path(input_file)
    suffix = input_path.suffix.lower()

    if suffix == ".csv":
        df = pd.read_csv(input_path)
    elif suffix == ".tsv":
        df = pd.read_csv(input_path, sep="\t")
    elif suffix in {".xlsx", ".xls"}:
        df = pd.read_excel(input_path)
    else:
        raise ValueError(f"Unsupported input file format: {suffix}")

    # 必要列
    spectrum_cols = [
        "Nonsynonymous",
        "Synonymous",
        "Intergenic",
        "Nonsense",
        "Noncoding"
    ]

    _validate_columns(df, [line_col] + spectrum_cols)

    # 解析 Line
    parsed = df[line_col].apply(lambda x: pd.Series(_parse_line_label(x)))
    parsed.columns = ["Day", "Group", "Replicate"]
    df[["Day", "Group", "Replicate"]] = parsed

    df = df.dropna(subset=spectrum_cols)

    # 如果没有解析成功，fallback
    if df["Group"].isna().all():
        df["Group"] = "All"
    if df["Day"].isna().all():
        df["Day"] = "All"

    #  核心：aggregate replicates
    agg = (
        df.groupby(["Group", "Day"])[spectrum_cols]
        .mean()
        .reset_index()
    )

    # 排序 Day（D60 → D120 → D180）
    if agg["Day"].nunique() > 1:
        day_order = _sort_day_labels(list(agg["Day"].unique()))
        agg["Day"] = pd.Categorical(agg["Day"], categories=day_order, ordered=True)

    agg = agg.sort_values(["Group", "Day"])

    # x轴标签
    agg["Label"] = agg["Group"].astype(str) + "-" + agg["Day"].astype(str)

    
    color_map = {
        "Nonsynonymous": "#4CAF50",  # green
        "Synonymous": "#2196F3",     # blue
        "Intergenic": "#FF9800",     # orange
        "Nonsense": "#F44336",       # red
        "Noncoding": "#9E9E9E"       # gray
    }

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    bottom = np.zeros(len(agg))

    # stacked bar
    for col in spectrum_cols:
        ax.bar(
            agg["Label"],
            agg[col],
            bottom=bottom,
            label=col,
            color=color_map.get(col, None),
            edgecolor="white"
        )
        bottom += agg[col].values

    # 样式
    ax.set_ylabel("Proportion", fontsize=12)
    ax.set_title("Mutation Spectrum", fontsize=14, pad=10)

    ax.set_ylim(0, 1)
    ax.grid(axis="y", linestyle="--", alpha=0.3)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.xticks(rotation=45, ha="right")

    ax.legend(title="Mutation Type", frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")

    plt.tight_layout()

    if output_file is not None:
        fig.savefig(output_file, dpi=300, bbox_inches="tight")

    if show:
        plt.show()

    return fig


def plot_time_trajectory(
    input_file,
    y_col="Total Mutations",
    output_file=None,
    show=True,
):


    input_path = Path(input_file)
    suffix = input_path.suffix.lower()
    if suffix == ".csv":
        df = pd.read_csv(input_path)
    elif suffix == ".tsv":
        df = pd.read_csv(input_path, sep="\t")
    elif suffix in {".xlsx", ".xls"}:
        df = pd.read_excel(input_path)
    else:
        raise ValueError(f"Unsupported input file format: {suffix}")

    parsed = df["Line"].apply(lambda x: pd.Series(_parse_line_label(x)))
    parsed.columns = ["Day", "Group", "Replicate"]
    df[["Day", "Group", "Replicate"]] = parsed

    group_parsed = not df["Group"].isna().all()
    if not group_parsed:
        df["Group"] = "All"
    else:
        unlabeled = df.loc[df["Group"].isna(), "Line"].tolist()
        if unlabeled:
            print(f"Time trajectory plot: {len(unlabeled)} sample(s) had no parsable "
                  f"day/treatment and are grouped as 'unassigned': {unlabeled}")
            df["Group"] = df["Group"].fillna("unassigned")

    if df["Day"].isna().all():
        sample_values = df["Line"].dropna().head(3).tolist()
        raise ValueError(
            "Could not parse day information from the 'Line' column. "
            "Each label must contain a timepoint in the form 'd<number>' (e.g. 'd60')"
            "followed by optional condition and replicate (e.g. 'd120-me1', 'sm-d120-me1-p')."
            f"The labels found in 'Line' were: {sample_values}"
        )

    df["Day_num"] = df["Day"].str.extract(r"(\d+)").astype(float)

    fig, ax = plt.subplots(figsize=(8, 6), dpi=200)

    palette = plt.cm.tab10.colors
    unique_groups = sorted(df["Group"].dropna().unique())
    colors = {g: palette[i % len(palette)] for i, g in enumerate(unique_groups)}

    for group in df["Group"].unique():
        sub = df[df["Group"] == group]

        ax.scatter(
            sub["Day_num"],
            sub[y_col],
            color=colors[group],
            alpha=0.5,
            s=40,
            label=group if group_parsed else None,
        )

        mean_df = sub.groupby("Day_num")[y_col].mean().reset_index()

        ax.plot(
            mean_df["Day_num"],
            mean_df[y_col],
            color=colors[group],
            linewidth=2.5,
            marker="o",
            label=None,
        )

    ax.set_xlabel("Time (Day)")
    ax.set_ylabel(y_col)
    ax.set_title(f"Time Trajectory of {y_col}")

    ax.grid(True, linestyle="--", alpha=0.3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if group_parsed:
        handles, labels = ax.get_legend_handles_labels()
        handles.append(Line2D([0], [0], color="#666666", linewidth=2.5, marker="o", label="mean"))
        labels.append("mean")
        ax.legend(handles, labels, frameon=False)

    plt.tight_layout()

    if output_file:
        fig.savefig(output_file, dpi=300, bbox_inches="tight")

    if show:
        plt.show()

    return fig

def plot_allele_distribution(
    input_file,
    output_file=None,
    figsize=(12, 6),
    dpi=200,
    show=True,
):
    """
    Plot allele frequency distribution from cleaned_data.
    Creates two panels: one for ME and one for PE, with overlaid histograms by Day.
    """

    input_path = Path(input_file)
    suffix = input_path.suffix.lower()

    if suffix == ".csv":
        df = pd.read_csv(input_path)
    elif suffix == ".tsv":
        df = pd.read_csv(input_path, sep="\t")
    elif suffix in {".xlsx", ".xls"}:
        df = pd.read_excel(input_path)
    else:
        raise ValueError(f"Unsupported input file format: {suffix}")

    required_cols = ["seq_id", "position", "mutation"]
    _validate_columns(df, required_cols)

    # All remaining columns are sample columns containing allele frequencies
    value_cols = [col for col in df.columns if col not in required_cols]

    # Wide -> long
    long_df = df.melt(
        id_vars=required_cols,
        value_vars=value_cols,
        var_name="Sample",
        value_name="Frequency"
    )

    # Keep valid numeric frequencies only
    long_df["Frequency"] = pd.to_numeric(long_df["Frequency"], errors="coerce")
    long_df = long_df.dropna(subset=["Frequency"]).copy()

    # Parse sample labels into Day / Group / Replicate
    parsed = long_df["Sample"].apply(lambda x: pd.Series(_parse_line_label(x)))
    parsed.columns = ["Day", "Group", "Replicate"]
    long_df[["Day", "Group", "Replicate"]] = parsed

    # Normalize for consistency
    long_df["Group"] = long_df["Group"].astype(str).str.lower()
    long_df["Day"] = long_df["Day"].astype(str).str.upper()

    # Keep only rows where Group was successfully parsed
    parsed_df = long_df[long_df["Group"] != "none"].copy()

    if parsed_df.empty:
        # Fall back: treat all samples as one group, no day breakdown
        parsed_df = long_df.copy()
        parsed_df["Group"] = "all"
        parsed_df["Day"] = "all"

    long_df = parsed_df

    days = _sort_day_labels(list(dict.fromkeys(long_df["Day"].dropna().astype(str))))
    groups = sorted(long_df["Group"].dropna().unique())

    palette = plt.cm.tab10.colors
    day_colors = {day: palette[i % len(palette)] for i, day in enumerate(days)}

    n = len(groups)
    fig, axes = plt.subplots(1, n, figsize=figsize, dpi=dpi, sharey=True)
    axes = [axes] if n == 1 else list(axes)

    for ax, group in zip(axes, groups):
        group_df = long_df[long_df["Group"] == group]

        for day in days:
            sub = group_df[group_df["Day"] == day]
            if sub.empty:
                continue

            ax.hist(
                sub["Frequency"],
                bins=30,
                alpha=0.5,
                label=day,
                color=day_colors.get(day, None),
                edgecolor="white",
            )

        ax.set_title(group.upper(), fontsize=13)
        ax.set_xlabel("Allele Frequency", fontsize=12)
        ax.grid(True, linestyle="--", alpha=0.25)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[0].set_ylabel("Count", fontsize=12)
    for ax in axes:
        ax.legend(
        title="Day",
        loc="upper right",
        fontsize=8,
        title_fontsize=9,
        frameon=False
        )

    fig.suptitle("Allele Frequency Distribution by Group", fontsize=14, y=1.02)
    plt.tight_layout()

    if output_file is not None:
        fig.savefig(output_file, dpi=300, bbox_inches="tight")

    if show:
        plt.show()

    return fig


def plot_parallel_mutation_heatmap(
    parallel_csv,
    top_n=15,
    day_order=("d60", "d120", "d180"),
    condition_order=("me", "pe"),
    output_file=None,
    show=True,
):
    """
    Create a gene-level parallel mutation heatmap.

    Parameters
    ----------
    parallel_csv : str
        Path to gene_parallel_mutations.csv
    top_n : int, default=15
        Number of top genes to display, ranked by strain_count
    day_order : tuple
        Desired order of day blocks
    condition_order : tuple
        Desired order of conditions within each day block

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object for saving outside the function
    """

    # -----------------------------
    # 1. Read data
    # -----------------------------
    df = pd.read_csv(parallel_csv)

    required_cols = {"gene", "strain_count"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Metadata columns that are not strain columns
    meta_cols = {"gene", "description", "shared_strains", "strain_count"}
    strain_cols = [c for c in df.columns if c not in meta_cols]

    if not strain_cols:
        raise ValueError("No strain columns detected in gene_parallel_mutations.csv")

    # -----------------------------
    # 2. Parse strain information
    # -----------------------------
    # Accepts e.g. SM-D120-ME1-P, kz19-d60-me1, or d180-me1 (no prefix, no suffix).
    pattern = re.compile(r"^(?:.*?-)?(d\d+)-([A-Za-z]+)(\d+)(?:-.*)?$", re.IGNORECASE)

    strain_info = []
    for col in strain_cols:
        match = pattern.match(col)
        if match:
            day = match.group(1).lower()
            condition = match.group(2).lower()
            replicate = int(match.group(3))
        else:
            day = "unknown"
            condition = "unknown"
            replicate = 999

        strain_info.append({
            "strain": col,
            "day": day,
            "condition": condition,
            "replicate": replicate
        })
     
    strain_info_df = pd.DataFrame(strain_info)
    strain_info_df["day"] = strain_info_df["day"].astype(str).str.lower()
    strain_info_df["condition"] = strain_info_df["condition"].astype(str).str.lower()

    has_day = bool((strain_info_df["day"] != "unknown").any())
    has_condition = bool((strain_info_df["condition"] != "unknown").any())

    # Auto-discover conditions present in the data but not listed in condition_order
    discovered_conditions = [
        c for c in dict.fromkeys(strain_info_df["condition"])
        if c != "unknown" and c.lower() not in {x.lower() for x in condition_order}
    ]
    condition_order = tuple(condition_order) + tuple(discovered_conditions)

    # Ranking for sorting
    day_rank = {d.lower(): i for i, d in enumerate(day_order)}
    cond_rank = {c.lower(): i for i, c in enumerate(condition_order)}

    strain_info_df["day_rank"] = strain_info_df["day"].map(day_rank).fillna(999)
    strain_info_df["cond_rank"] = strain_info_df["condition"].map(cond_rank).fillna(999)

    strain_info_df = strain_info_df.sort_values(
        by=["day_rank", "cond_rank", "replicate", "strain"]
    ).reset_index(drop=True)

    ordered_strains = strain_info_df["strain"].tolist()

    # -----------------------------
    # 3. Select top genes
    # -----------------------------
    df = df.sort_values(by=["strain_count", "gene"], ascending=[False, True]).copy()
    df_top = df.head(top_n).copy()

    heatmap_df = df_top.set_index("gene")[ordered_strains].copy()
    heatmap_df = heatmap_df.apply(pd.to_numeric, errors="coerce").fillna(0)
    heatmap_df = (heatmap_df > 0).astype(int)

    matrix = heatmap_df.values
    n_rows, n_cols = matrix.shape

    # -----------------------------
    # 4. Create figure
    # -----------------------------
    fig_width = max(12, n_cols * 0.35)
    fig_height = max(6, n_rows * 0.45 + 1.5)

    fig = plt.figure(figsize=(fig_width, fig_height))

    row_specs = []
    if has_day:
        row_specs.append(("day", 0.3))
    if has_condition:
        row_specs.append(("cond", 0.3))
    row_specs.append(("main", 4))

    gs = fig.add_gridspec(
        len(row_specs), 1,
        height_ratios=[h for _, h in row_specs],
        hspace=0.05,
    )
    axes = {name: fig.add_subplot(gs[i]) for i, (name, _) in enumerate(row_specs)}
    ax_day = axes.get("day")
    ax_cond = axes.get("cond")
    ax = axes["main"]
    top_ax = axes[row_specs[0][0]]

    # -----------------------------
    # 5. Annotation bars
    # -----------------------------
    day_palette = {d.lower(): i for i, d in enumerate(day_order)}
    day_palette["unknown"] = len(day_order)
    cond_palette = {c.lower(): i for i, c in enumerate(condition_order)}
    cond_palette["unknown"] = len(condition_order)

    _day_base = plt.cm.Blues(np.linspace(0.3, 0.8, max(len(day_order), 1)))
    day_colors = [to_hex(c) for c in _day_base] + ["#d9d9d9"]
    # Pinned colors for known conditions; any new condition falls back to the
    # remaining Set2 colors so its color stays stable regardless of column order.
    fixed_cond_colors = {
        "me": "#66c2a5",
        "pe": "#fc8d62",
        "mpm": "#8da0cb",
        "pmp": "#e78ac3",
    }
    _pinned = set(fixed_cond_colors.values())
    _auto_base = [to_hex(c) for c in plt.cm.Set2.colors if to_hex(c) not in _pinned]
    cond_colors = []
    _auto_i = 0
    for c in condition_order:
        key = c.lower()
        if key in fixed_cond_colors:
            cond_colors.append(fixed_cond_colors[key])
        else:
            cond_colors.append(_auto_base[_auto_i % len(_auto_base)])
            _auto_i += 1
    cond_colors.append("#d9d9d9")  # unknown
    binary_colors = ["#f2f2f2", "#08306b"]

    day_cmap = ListedColormap(day_colors)
    cond_cmap = ListedColormap(cond_colors)
    binary_cmap = ListedColormap(binary_colors)

    annotation_axes = []
    if ax_day is not None:
        day_vals = np.array([
            day_palette.get(d, day_palette["unknown"])
            for d in strain_info_df["day"]
        ]).reshape(1, -1)
        ax_day.imshow(day_vals, aspect="auto", cmap=day_cmap, interpolation="none",
                      vmin=0, vmax=len(day_colors) - 1)
        annotation_axes.append((ax_day, "Day"))

    if ax_cond is not None:
        cond_vals = np.array([
            cond_palette.get(c, cond_palette["unknown"])
            for c in strain_info_df["condition"]
        ]).reshape(1, -1)
        ax_cond.imshow(cond_vals, aspect="auto", cmap=cond_cmap, interpolation="none",
                       vmin=0, vmax=len(cond_colors) - 1)
        annotation_axes.append((ax_cond, "Condition"))

    for a, label in annotation_axes:
        a.set_xticks([])
        a.set_yticks([0])
        a.set_yticklabels([label], fontsize=10)
        for spine in a.spines.values():
            spine.set_visible(False)

    # -----------------------------
    # 6. Main heatmap
    # -----------------------------
    ax.imshow(matrix, aspect="auto", cmap=binary_cmap, interpolation="none",vmin=0, vmax=1)

    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(ordered_strains, rotation=90, fontsize=8)

    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(heatmap_df.index, fontsize=9)

    ax.set_xlabel("Strains", fontsize=11)
    ax.set_ylabel("Genes", fontsize=11)
    top_ax.set_title(
        f"Parallel Mutation Heatmap (Top {len(heatmap_df.index)} Genes)",
        fontsize=13, pad=10,
    )

    # Cell gridlines
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="#c7c7c7", linestyle="-", linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    # -----------------------------
    # 7. Boundary lines
    # -----------------------------
    boundary_axes = [ax] + [a for a, _ in annotation_axes]
    # Condition boundaries must not split the Day bar (a single day stays solid).
    cond_boundary_axes = [ax] + [a for a, label in annotation_axes if label != "Day"]

    # Dashed red boundaries between adjacent conditions within each day
    if has_day and has_condition:
        for day in dict.fromkeys(strain_info_df["day"]):
            day_subset = strain_info_df[strain_info_df["day"] == day]
            positions = day_subset.index.to_list()
            conds = day_subset["condition"].to_list()
            for i in range(1, len(positions)):
                if conds[i] != conds[i - 1]:
                    boundary_x = positions[i] - 0.5
                    for axis in cond_boundary_axes:
                        axis.axvline(x=boundary_x, color="red", linestyle="--", linewidth=1)

    # Solid black boundaries between days
    if has_day:
        present_days = list(dict.fromkeys(strain_info_df["day"]))
        for day in present_days[:-1]:
            day_subset = strain_info_df[strain_info_df["day"] == day]
            if day_subset.empty:
                continue
            boundary_x = day_subset.index.max() + 0.5
            for axis in boundary_axes:
                axis.axvline(x=boundary_x, color="black", linestyle="-", linewidth=2)
    # -----------------------------
    # 8. Legends
    # -----------------------------
    legend_specs = []

    if has_day:
        day_handles = [
            Patch(facecolor=day_colors[day_palette[d]], edgecolor="none", label=d)
            for d in day_order if d in strain_info_df["day"].values
        ]
        if day_handles:
            legend_specs.append(("Day", day_handles))

    if has_condition:
        cond_handles = [
            Patch(facecolor=cond_colors[cond_palette[c]], edgecolor="none", label=c.upper())
            for c in condition_order if c in strain_info_df["condition"].values
        ]
        if cond_handles:
            legend_specs.append(("Condition", cond_handles))

    mutation_handles = [
        Patch(facecolor=binary_colors[0], edgecolor="black", label="0"),
        Patch(facecolor=binary_colors[1], edgecolor="black", label="1"),
    ]
    legend_specs.append(("Mutation", mutation_handles))

    y_anchor = 1.00
    for i, (title, handles) in enumerate(legend_specs):
        leg = ax.legend(
            handles=handles,
            title=title,
            loc="upper left",
            bbox_to_anchor=(1.02, y_anchor),
            frameon=False,
        )
        if i < len(legend_specs) - 1:
            ax.add_artist(leg)
        y_anchor -= 0.05 * (len(handles) + 1) + 0.05

    if output_file is not None:
        fig.savefig(output_file, dpi=300, bbox_inches="tight")

    if show:
        plt.show()

    return fig


_STRUCT_COLORS = {
    "snp": "#d62728",
    "ins": "#000000",
    "del": "#ffffff",
    "other": "#7f7f7f",
}

# Match breseq patterns for different structural mutation types
_SNP_RE = re.compile(r"^[ACGTN]\s*(?:→|->|>)\s*[ACGTN]$", re.IGNORECASE)
_DEL_RE = re.compile(r"^Δ\s*([\d,]+)\s*bp$", re.IGNORECASE)
_INS_RE = re.compile(r"^\+\s*([\d,]+)\s*bp$")
_INS_SEQ_RE = re.compile(r"^\+\s*([ACGTN]+)$", re.IGNORECASE)
_REPEAT_RE = re.compile(
    r"^\(\s*([ACGTN]+)\s*\)\s*(\d+)\s*(?:→|->|>)\s*(\d+)$", re.IGNORECASE
)
_POSITION_RE = re.compile(r"^\s*([\d,]+)")


def _parse_position(value):
    """
    Parse a breseq position cell into an int, handling thousands separators
    and the ":N" sub-position suffix used for insertions.
    """
    if pd.isna(value):
        return None
    m = _POSITION_RE.match(str(value))
    if not m:
        return None
    return int(m.group(1).replace(",", ""))


def _classify_structural_mutation(mutation):
    """
    Classify a `mutation` string into a structural type (snp/ins/del/other)
    """
    s = str(mutation).strip()
    if _SNP_RE.match(s):
        return "snp", 0
    m = _DEL_RE.match(s)
    if m:
        return "del", int(m.group(1).replace(",", ""))
    m = _INS_RE.match(s)
    if m:
        return "ins", int(m.group(1).replace(",", ""))
    m = _INS_SEQ_RE.match(s)
    if m:
        return "ins", len(m.group(1))
    m = _REPEAT_RE.match(s)
    if m:
        unit, before, after = len(m.group(1)), int(m.group(2)), int(m.group(3))
        delta = (after - before) * unit
        return ("ins", delta) if delta > 0 else ("del", -delta)
    return "other", 0


def plot_zoomed_genome(
    clean_csv,
    reference_path,
    companion_path=None,
    ancestor=None,
    output_file=None,
    min_gap="gene",
    gene_split_gap=10000,
    section_buffer=1,
    shared_only=True,
    min_strains=0.25,
    figsize=(15, 8),
    dpi=200,
    show=True,
):
    """
    Genome plot adapted from Kosterlitz (https://github.com/livkosterlitz/Breseq_genome_plots)

    Whole genome is a bar on top, the regions with mutations are
    magnified into a strain table below (marks are colored by 
    mutation type, rows grouped and colored by day+condition).
    A variant histogram with counts is at the bottom.
    Parallel mutations are represented by aligned marks
    and wide histogram bars.

    Parameters
    ----------
    clean_csv : str
        Path to cleaned_data.csv
    reference_path : str
        Reference genome (.gbk/.fasta/.gff). 
        Used only for contig lengths and coordinates 
    companion_path : str, optional
        Companion FASTA for a bare GFF reference.
    ancestor : str, optional
        Ancestor column name, excluded from the strain rows.
    min_gap : int or "gene", default "gene"
        How to group mutations into sections.
        "gene" = one section per gene
        int = new region when mutations are more than {min_gap} bp apart
    gene_split_gap : int, default 10000
        Only used when min_gap="gene". A gene name that reappears more than
        this many bases away is treated as a separate locus.
    section_buffer : int, default 1
        Bases added to each end of a region.
    shared_only : bool, default True
        Keep only sections hit by more than one strain (parallel regions).
        Set False to display every single region hit by mutations.
    min_strains : float or int, default 0.25
        Minimum number of strains a region must be mutated in to be kept,
        used to focus the plot on parallel and drop isolated mutations.
        Overrides shared_only.
        - a fraction in (0, 1) is a proportion of strains
          (0.25 keeps regions hit by more than 1/4 of strains)
          the default adapts to datasets of any size
        - an integer is an absolute strain count
        - None uses the default shared_only (more or equal to 2 strains).

    Output
    -------
    fig : matplotlib.figure.Figure or None
        None if there is nothing to plot.
    """
    # Read data and reference

    df = pd.read_csv(clean_csv)
    _validate_columns(df, ["seq_id", "position", "mutation", "gene"])

    from .simulation.reference_loader import load_sequences
    from .simulation.mutation_model import build_contig_pointers

    seq = load_sequences(reference_path, companion_path)
    pointer, n_bases = build_contig_pointers(seq)
    if not n_bases:
        print("Zoomed genome plot skipped: reference has zero length.")
        return None

    strain_cols = get_strain_columns(df, ancestor)
    if not strain_cols:
        print("Zoomed genome plot skipped: no strain columns found.")
        return None

    # Convert data from wide to long format (one row per strain per mutation)
    records = []
    unmapped = set()
    missing_contig = 0
    unparsed_pos = 0
    for _, row in df.iterrows():
        if pd.isna(row["seq_id"]):
            missing_contig += 1
            continue
        contig = str(row["seq_id"])
        if contig not in pointer:
            unmapped.add(contig)
            continue
        pos = _parse_position(row["position"])
        if pos is None:
            unparsed_pos += 1
            continue
        gpos = pointer[contig][0] + pos
        mtype, _size = _classify_structural_mutation(row["mutation"])
        gene = row["gene"]
        gene = str(gene).strip() if pd.notna(gene) and str(gene).strip() else None
        for strain in strain_cols:
            cell = row[strain]
            if pd.isna(cell) or str(cell).strip() in ("", "?"):
                continue
            day, group, replicate = _parse_line_label(strain)
            records.append({
                "strain": strain,
                "day": day,
                "group": group,
                "replicate": replicate,
                "gpos": gpos,
                "mut": str(row["mutation"]).strip(),
                "gene": gene if gene is not None else f"intergenic@{gpos}",
                "mtype": mtype,
                "size": _size,
            })

    if unmapped:
        print(f"Zoomed genome plot: {len(unmapped)} contig(s) in the table were "
              f"not found in the reference and were skipped: {sorted(unmapped)}")
    if missing_contig:
        print(f"Zoomed genome plot: {missing_contig} row(s) with no seq_id were skipped.")
    if unparsed_pos:
        print(f"Zoomed genome plot: {unparsed_pos} row(s) with an unparsable "
              f"position were skipped.")
    if not records:
        print("Zoomed genome plot skipped: no mutations to plot.")
        return None

    long_df = pd.DataFrame(records)

    # Group mutations into regions
    gene_mode = isinstance(min_gap, str) and min_gap.lower() == "gene"
    long_df["gend"] = long_df["gpos"] + long_df["size"].where(
        long_df["mtype"] == "del", 0)

    site_gene = {}
    for _, rec in long_df.iterrows():
        site_gene.setdefault(rec["gpos"], rec["gene"])
        site_gene.setdefault(rec["gend"], rec["gene"])

    pos2cid, cid = {}, 0
    prev_pos, prev_gene = None, None
    for pos in sorted(site_gene):
        gene = site_gene[pos]
        if prev_pos is not None:
            if gene_mode:
                split = gene != prev_gene or (pos - prev_pos) > gene_split_gap
            else:
                split = (pos - prev_pos) > int(min_gap)
            if split:
                cid += 1
        pos2cid[pos] = cid
        prev_pos, prev_gene = pos, gene

    long_df["section_key"] = long_df["gpos"].map(pos2cid)
    long_df["section_end"] = long_df["gend"].map(pos2cid)

    sites_by_cid, strains_by_cid, gene_by_cid = {}, {}, {}
    for pos, c in pos2cid.items():
        sites_by_cid.setdefault(c, []).append(pos)
        gene_by_cid.setdefault(c, site_gene[pos])
    for _, rec in long_df.iterrows():
        for c in range(rec["section_key"], rec["section_end"] + 1):
            strains_by_cid.setdefault(c, set()).add(rec["strain"])

    regions = []
    for c, sites in sites_by_cid.items():
        regions.append({
            "key": c,
            "init": min(sites) - section_buffer,
            "fin": max(sites) + section_buffer,
            "n_strains": len(strains_by_cid.get(c, ())),
            "label": gene_by_cid[c] if gene_mode else "",
        })

    n_experiment_strains = len(strain_cols)
    if min_strains is None:
        threshold = 2 if shared_only else 1
    elif 0 < min_strains < 1:
        threshold = max(2, int(round(min_strains * n_experiment_strains)))
    else:
        threshold = int(min_strains)
    if threshold > 1:
        regions = [r for r in regions if r["n_strains"] >= threshold]
        if not regions:
            print(f"Zoomed genome plot skipped: no regions hit by at least "
                  f"{threshold} strains (lower min_strains or set shared_only=False).")
            return None

    regions.sort(key=lambda r: r["init"])
    keep_keys = {r["key"] for r in regions}
    long_df = long_df[long_df["section_key"].isin(keep_keys)].copy()

    # Each variant site is represented as a column in the plot.
    kept_sites = sorted(pos for pos, c in pos2cid.items() if c in keep_keys)
    n_sites = len(kept_sites)
    site_idx = {pos: i for i, pos in enumerate(kept_sites)}

    # Column width: equal per site, but adjusted to plots with only a few sites
    # to use narrow cells (marks read as vertical bars) instead of wide blocks

    MAX_COL_W = 0.025
    unit_w = 1.0 / n_sites if n_sites else 1.0
    capped = unit_w > MAX_COL_W
    col_w = min(unit_w, MAX_COL_W)
    table_w = n_sites * col_w
    table_x0 = 0.0

    def col_x(i):
        return table_x0 + i * col_w

    for s in regions:
        idxs = [site_idx[p] for p in sites_by_cid[s["key"]] if p in site_idx]
        s["new_init"] = col_x(min(idxs))
        s["new_fin"] = col_x(max(idxs) + 1)
        s["gb_init"] = table_x0 + (s["init"] / n_bases) * table_w
        s["gb_fin"] = table_x0 + (s["fin"] / n_bases) * table_w
    sec_by_key = {s["key"]: s for s in regions}

    # Order treatments, then order strains within a treatment by variant count
    strain_counts = long_df.groupby("strain").size().to_dict()
    strain_first = long_df.groupby("strain")["gpos"].min().to_dict()

    def strain_sort_key(strain):
        day, group, rep = _parse_line_label(strain)
        day_num = int(re.search(r"\d+", day).group()) if day else 10 ** 9
        return (day_num, group or "~", strain_counts[strain],
                strain_first[strain], strain)

    ordered_strains = sorted(set(long_df["strain"]), key=strain_sort_key)

    parsed = {s: _parse_line_label(s) for s in ordered_strains}
    groups = [g for g in dict.fromkeys(p[1] for p in parsed.values()) if g]
    days = _sort_day_labels([d for d in dict.fromkeys(p[0] for p in parsed.values()) if d])
    color_map = _build_color_map(groups, days) if groups and days else {}

    unlabeled = [s for s in ordered_strains if not (parsed[s][0] and parsed[s][1])]
    if unlabeled:
        print(f"Zoomed genome plot: {len(unlabeled)} strain(s) had no parsable "
              f"day/treatment and are grouped as 'unassigned': {unlabeled}")

    def treatment_label(strain):
        day, group, _ = parsed[strain]
        if day and group:
            return f"{day}-{group}"
        return "unassigned"

    def treatment_color(strain):
        day, group, _ = parsed[strain]
        return color_map.get((group, day), "#c7c7c7")

    # Plot figure
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    section_colors = ["#e9e9e9", "#d0d0d0"]
    grid_color = "#b0b0b0"
    n_str = len(ordered_strains)

    # Vertical bands
    gene_width = 0.03
    mag_spacer = 0.05
    tab_spacer = 0.025
    row_h = (1.0 - gene_width - mag_spacer - tab_spacer) / (n_str + 1)
    y_gene_bot = 1.0 - gene_width
    y_sec_top = y_gene_bot - mag_spacer
    y_sec_bot = y_sec_top - row_h
    y_table_top = y_sec_bot - tab_spacer
    y_table_bot = y_table_top - n_str * row_h

    # Horizontal bar panel on the right
    hist_spacer, hist_len = 0.02, 0.5
    bx0 = table_x0 + table_w + hist_spacer
    bx1 = bx0 + hist_len
    y_axis = y_table_bot - 0.02

    # Genome bar
    ax.add_patch(Rectangle((table_x0, y_gene_bot), table_w, gene_width,
                           facecolor="ivory", edgecolor="black", lw=0.6))
    for i, s in enumerate(regions):
        gb_w = max(s["gb_fin"] - s["gb_init"], 0.0005)
        edge = "black" if gb_w > 0.002 else "none"
        ax.add_patch(Rectangle((s["gb_init"], y_gene_bot), gb_w, gene_width,
                               facecolor="black" if gb_w <= 0.002 else section_colors[i % 2],
                               edgecolor=edge, lw=0.4))
    ax.text(table_x0, 1.005, "0 Mb", ha="left", va="bottom", fontsize=10)
    ax.text(table_x0 + table_w, 1.005, f"{n_bases / 1e6:.1f} Mb",
            ha="right", va="bottom", fontsize=10)

    # Magnification shapes
    poly_colors = ["#9a9a9a", "#7a7a7a"]
    for i, s in enumerate(regions):
        poly = Polygon(
            [(s["gb_init"], y_gene_bot), (s["gb_fin"], y_gene_bot),
             (s["new_fin"], y_sec_top), (s["new_init"], y_sec_top)],
            closed=True, facecolor=poly_colors[i % 2], alpha=0.7,
            edgecolor="none",
        )
        ax.add_patch(poly)

    # Region header bar and optional gene labels
    label_sections = len(regions) <= 60
    for i, s in enumerate(regions):
        w = s["new_fin"] - s["new_init"]
        ax.add_patch(Rectangle((s["new_init"], y_sec_bot), w, row_h,
                               facecolor=section_colors[i % 2],
                               edgecolor="none" if w <= 0.002 else grid_color, lw=0.5))
        if label_sections and w > 0.01 and s["label"]:
            ax.text((s["new_init"] + s["new_fin"]) / 2.0, y_sec_top + 0.005,
                    s["label"], rotation=90, ha="center", va="bottom",
                    fontsize=6, color="#333333")

    # Mutation table: one row per strain, colored by treatment.
    records_by_strain = {s: g for s, g in long_df.groupby("strain")}
    point_frac = 0.65 if capped else 1.0
    mark_w = col_w * point_frac
    for r, strain in enumerate(ordered_strains):
        row_top = y_table_top - r * row_h
        row_bot = row_top - row_h
        ax.add_patch(Rectangle((table_x0, row_bot), table_w, row_h,
                               facecolor=treatment_color(strain), alpha=0.35,
                               edgecolor="none"))
        grp = records_by_strain.get(strain)
        if grp is None:
            continue
        for _, rec in grp.iterrows():
            i0 = site_idx.get(rec["gpos"])
            if i0 is None:
                continue
            color = _STRUCT_COLORS.get(rec["mtype"], _STRUCT_COLORS["other"])
            if rec["mtype"] == "del" and rec["size"] > 0:
                i1 = max(site_idx.get(rec["gend"], i0), i0)
                x0 = col_x(i0)
                mw = col_x(i1 + 1) - x0
            else:
                x0 = col_x(i0) + (col_w - mark_w) / 2.0
                mw = mark_w
            ax.add_patch(Rectangle((x0, row_bot), mw, row_h, facecolor=color,
                                   edgecolor="none", zorder=3))

    px_per_site = figsize[0] * dpi / (1.72 * max(n_sites, 1))
    if px_per_site >= 4:
        for i in range(1, n_sites):
            x = col_x(i)
            ax.plot([x, x], [y_table_bot, y_sec_top],
                    color="#dddddd", lw=0.25, zorder=4)
    elif px_per_site < 2:
        print(f"Zoomed genome plot: {n_sites} site columns are too narrow to "
              f"separate; internal grid omitted (raise figsize/dpi or min_strains).")
    for s in regions:
        ax.plot([s["new_init"], s["new_init"]], [y_table_bot, y_sec_top],
                color=grid_color, lw=0.5, zorder=5)
    ax.plot([table_x0 + table_w, table_x0 + table_w], [y_table_bot, y_sec_top],
            color=grid_color, lw=0.5, zorder=5)
    for r in range(n_str + 1):
        y = y_table_top - r * row_h
        ax.plot([table_x0, table_x0 + table_w], [y, y],
                color=grid_color, lw=0.4, zorder=5)

    # Treatment labels
    block_start = 0
    for r in range(1, n_str + 1):
        boundary = (r == n_str) or (treatment_label(ordered_strains[r]) != treatment_label(ordered_strains[block_start]))
        if boundary:
            top_y = y_table_top - block_start * row_h
            bot_y = y_table_top - r * row_h
            ax.text(table_x0 - 0.01, (top_y + bot_y) / 2.0,
                    treatment_label(ordered_strains[block_start]),
                    ha="right", va="center", fontsize=8)
            block_start = r

    # One bar per strain, number of variants it carries
    max_count = max(strain_counts[s] for s in ordered_strains) or 1
    unit = hist_len / max_count
    for r, strain in enumerate(ordered_strains):
        row_top = y_table_top - r * row_h
        ax.add_patch(Rectangle((bx0, row_top - row_h), strain_counts[strain] * unit,
                               row_h, facecolor=treatment_color(strain), alpha=0.35,
                               edgecolor="black", lw=0.6))

    step = 10 ** max(0, int(np.ceil(np.log10(max_count / 13.0))))
    axis_end = bx0 + hist_len + 0.02
    ax.annotate("", xy=(axis_end, y_axis), xytext=(bx0, y_axis),
                arrowprops=dict(arrowstyle="->", color="black", lw=0.8))
    for tick in range(step, max_count + 1, step):
        tx = bx0 + tick * unit
        ax.plot([tx, tx], [y_axis, y_axis - 0.02], color="black", lw=0.8)
        ax.text(tx, y_axis - 0.03, str(tick), ha="center", va="top", fontsize=9)
    ax.text(bx0 + hist_len * 0.5, y_axis - 0.08,
            "number of variants", ha="center", va="top", fontsize=11)

    # Legend for mutation types
    legend_handles = [
        Patch(facecolor=_STRUCT_COLORS["snp"], edgecolor="black", label="SNP"),
        Patch(facecolor=_STRUCT_COLORS["ins"], edgecolor="black", label="insertion"),
        Patch(facecolor=_STRUCT_COLORS["del"], edgecolor="black", label="deletion"),
        Patch(facecolor=_STRUCT_COLORS["other"], edgecolor="black", label="other"),
    ]
    # Legend below the "number of variants" label, centered on the figure
    ax.legend(handles=legend_handles, loc="upper center", ncol=4,
              bbox_to_anchor=(0.5, y_axis - 0.12),
              bbox_transform=ax.get_yaxis_transform(),
              frameon=False, handlelength=1.4, columnspacing=1.6)

    ax.set_title(f"Zoomed genome plot: {n_str} strains",
                 fontsize=11, fontweight="bold")
    ax.set_xlim(-0.15, axis_end + 0.05)
    ax.set_ylim(y_axis - 0.20, 1.06)
    ax.axis("off")

    if output_file is not None:
        fig.savefig(output_file, dpi=300, bbox_inches="tight")

    if show:
        plt.show()

    return fig
