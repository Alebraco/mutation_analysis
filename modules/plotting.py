"""
Plotting utilities for mutation analysis pipeline.

Includes:
- Bubble plot (mutation summary)
- Mutation spectrum plot
- Time trajectory plot
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Optional, Dict, Tuple, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import to_hex, to_rgb


def _validate_columns(df: pd.DataFrame, required_cols: List[str]) -> None:
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns: {missing}. "
            f"Available columns are: {list(df.columns)}"
        )


def _parse_line_label(line: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Parse labels like:
    - sm-d120-me1-p
    - d120-me1
    Returns: (day, group, replicate)
    """
    if pd.isna(line):
        return None, None, None

    match = re.search(r"([dD]\d+)-([A-Za-z]+)(\d+)", str(line))
    if not match:
        return None, None, None

    day = match.group(1)
    group = match.group(2)
    replicate = match.group(3)
    return day, group, replicate


def _sort_day_labels(day_values: List[str]) -> List[str]:
    def key_func(x: str):
        m = re.search(r"(\d+)", str(x))
        return int(m.group(1)) if m else 999999

    return sorted(day_values, key=key_func)


def _adjust_color_lightness(color: str, factor: float) -> str:
    rgb = np.array(to_rgb(color))
    if factor >= 1:
        adjusted = rgb + (1 - rgb) * (factor - 1)
    else:
        adjusted = rgb * factor
    adjusted = np.clip(adjusted, 0, 1)
    return to_hex(adjusted)


def _build_color_map(groups: List[str], days: List[str], base_colors: Optional[Dict[str, str]] = None):
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


def _nice_round(x: float) -> int:
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


def _get_size_legend_values(values, n_levels: int = 3):
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


def _compute_size_scale(values, target_max_area: float = 700.0) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr) & (arr > 0)]

    if arr.size == 0:
        return 1.0

    vmax = np.max(arr)
    if vmax == 0:
        return 1.0

    return target_max_area / vmax


def plot_mutation_bubble(
    input_file: str,
    output_file: Optional[str] = None,
    x_col: str = "Nonsynonymous",
    y_col: str = "Average Frequency",
    size_col: str = "Total Mutations",
    line_col: str = "Line",
    title: str = "Mutation Frequency and Nonsynonymous Proportion",
    panel_by_group: bool = False,
    base_colors: Optional[Dict[str, str]] = None,
    add_reference_line: bool = True,
    reference_x: float = 2 / 3,
    figsize: tuple = (12, 7),
    dpi: int = 200,
    show: bool = True,
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
    input_file: str,
    output_file: Optional[str] = None,
    line_col: str = "Line",
    figsize: tuple = (10, 6),
    dpi: int = 200,
    show: bool = True,
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
    input_file: str,
    y_col: str = "Total Mutations",
    output_file: Optional[str] = None,
    show: bool = True,
):


    df = pd.read_excel(input_file)

    # parse
    parsed = df["Line"].apply(lambda x: pd.Series(_parse_line_label(x)))
    parsed.columns = ["Day", "Group", "Replicate"]
    df[["Day", "Group", "Replicate"]] = parsed

    df["Day_num"] = df["Day"].str.extract(r"(\d+)").astype(float)

    fig, ax = plt.subplots(figsize=(8, 6), dpi=200)

    colors = {"ME": "#4CAF50", "PE": "#8E63CE"}

    for group in df["Group"].unique():
        sub = df[df["Group"] == group]

        # scatter（所有replicate）
        ax.scatter(
            sub["Day_num"],
            sub[y_col],
            color=colors[group],
            alpha=0.5,
            s=40,
            label=f"{group} (replicates)"
        )

        #  mean line
        mean_df = sub.groupby("Day_num")[y_col].mean().reset_index()

        ax.plot(
            mean_df["Day_num"],
            mean_df[y_col],
            color=colors[group],
            linewidth=2.5,
            marker="o",
            label=f"{group} (mean)"
        )

    ax.set_xlabel("Time (Day)")
    ax.set_ylabel(y_col)
    ax.set_title(f"Time Trajectory of {y_col}")

    ax.grid(True, linestyle="--", alpha=0.3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.legend(frameon=False)

    plt.tight_layout()

    if output_file:
        fig.savefig(output_file, dpi=300, bbox_inches="tight")

    if show:
        plt.show()

    return fig

def plot_allele_distribution(
    input_file: str,
    output_file: Optional[str] = None,
    figsize: tuple = (12, 6),
    dpi: int = 200,
    show: bool = True,
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

    # Keep only rows with parsed Group/Day
    long_df = long_df[
        long_df["Group"].isin(["me", "pe"]) &
        long_df["Day"].str.match(r"D\d+", na=False)
    ].copy()

    if long_df.empty:
        raise ValueError("No valid allele frequency data found after parsing sample labels.")

    days = _sort_day_labels(list(dict.fromkeys(long_df["Day"].dropna().astype(str))))
    groups = ["me", "pe"]

    day_colors = {
        "D60": "#5B9BD5",
        "D120": "#ED7D31",
        "D180": "#70AD47",
    }

    fig, axes = plt.subplots(1, 2, figsize=figsize, dpi=dpi, sharey=True)

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
















