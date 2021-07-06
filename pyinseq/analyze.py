#!/usr/bin/env python
"""
Analyzes resulting output
"""
import numpy as np
import pandas as pd
import seaborn as sn
import matplotlib.pyplot as plt


def read_sites_file(sample: str, settings: "runner.Settings") -> pd.DataFrame:
    return pd.read_csv(f"{settings.path}{sample}_sites.txt", sep="\t")


def read_summary_table(settings: "runner.Settings") -> pd.DataFrame:
    return pd.read_csv(settings.summary_table, sep="\t")


def t_fifty(sample: str, settings: "runner.Settings") -> pd.DataFrame:
    """T50: The minimum number of transposon insertion sites in the sample that
    account for at least 50% of the samples's reads.
    """
    # calculate sum of counts in df
    df = read_sites_file(sample, settings)
    return np.sum(df.total_counts.sort_values().cumsum() > (df.total_counts.sum() / 2))


def spearman_correlation(samples: dict, settings: "runner.Settings") -> pd.DataFrame:
    """Calculates spearman correlation among samples and plots heatmap"""
    df = read_summary_table(settings)
    corr_df = df[list(samples.keys())].corr("spearman")
    corr_df.to_csv(f"{settings.analysis_path}correlation.csv", sep="\t")
    plot_heatmap(corr_df)
    return


def plot_insertions(samples: dict, settings: "runner.Settings") -> pd.DataFrame:
    """Plot insertions across each contig/chromosome."""
    for sample in samples.keys():

        sites_df = read_sites_file(sample, settings)

        for i, group in sites_df.groupby("contig"):
            fig, ax = plt.subplots()
            group.plot(
                kind="scatter",
                x="nucleotide",
                y="cpm",
                color="blue",
                alpha=0.5,
                title=str(i),
                ax=ax,
            )
            fig.savefig(
                f"{settings.figures_path}insertions_scatter_{sample}_{str(i)}.pdf",
                format="pdf",
            )
            return


def plot_heatmap(df):
    pass


def plot_MA(sample, settings):
    pass


def plot_pairwise(samples, settings):
    pass


def df_raw_counts(samples: dict, settings: "runner.Settings") -> pd.DataFrame:
    """Returns total raw counts from samples in a DataFrame."""
    # Results DataFrame with desired columns
    result_df = pd.DataFrame(
        columns=["contig", "nucleotide", "left_counts", "right_counts", "total_counts"]
    )
    col = result_df.columns[2:5]

    # Loop over samples and add cell values to result_df
    for sample in samples.keys():
        df = pd.read_csv(f"{settings.path}{sample}_sites.txt", sep="\t").drop(
            "cpm", axis=1
        )
        for i, row in df.iterrows():
            if row["nucleotide"] in result_df.nucleotide.tolist():
                for c in col:
                    previous = result_df.loc[
                        row["nucleotide"] == result_df["nucleotide"], c
                    ].item()
                    result_df.loc[row["nucleotide"] == result_df["nucleotide"], c] = (
                        previous + row[c]
                    )
            else:
                result_df = result_df.append(row)
    result_df.to_csv(
        f"{settings.analysis_path}{settings.experiment}_raw_counts.csv", sep="\t"
    )
    return


def df_cpm_counts(samples: dict, settings: "runner.Settings") -> pd.DataFrame:
    """Returns total cpm counts from samples in a DataFrame."""
    # Results DataFrame with desired columns
    result_df = pd.DataFrame(columns=["contig", "nucleotide", "cpm"])

    # Loop over samples and add cell values to result_df
    for sample in samples.keys():
        df = pd.read_csv(
            f"{settings.path}{sample}_sites.txt", sep="\t", usecols=[0, 1, 5]
        )
        for i, row in df.iterrows():
            if row["nucleotide"] in result_df.nucleotide.tolist():
                previous = result_df.loc[
                    row["nucleotide"] == result_df["nucleotide"], "cpm"
                ].item()
                result_df.loc[row["nucleotide"] == result_df["nucleotide"], "cpm"] = (
                    previous + row["cpm"]
                )
            else:
                result_df = result_df.append(row)
    result_df.to_csv(
        f"{settings.analysis_path}{settings.experiment}_cpm_counts.csv", sep="\t"
    )
    return


def df_gene_counts(samples: dict, settings: "runner.Settings") -> pd.DataFrame:
    """Return total count of hits per gene."""
    # Results DataFrame with desired columns
    results_df = pd.DataFrame(
        columns=[
            "contig",
            "left_counts",
            "right_counts",
            "total_counts",
            "cpm",
            "locus_tag",
        ]
    )
    col = results_df.columns[1:5]

    # Loop over samples and add cell values to result_df
    for sample in samples.keys():
        # TODO: caclulate average of three_primeness of gene across samples
        # Add gene index to results_df
        df = pd.read_csv(f"{settings.path}{sample}_genes.txt", sep="\t").drop(
            ["nucleotide", "three_primeness"], axis=1
        )
        for i, row in df.iterrows():
            if row["locus_tag"] in results_df.locus_tag.tolist():
                for c in col:
                    previous = results_df.loc[
                        results_df["locus_tag"] == row["locus_tag"], c
                    ].item()
                    results_df.loc[results_df["locus_tag"] == row["locus_tag"], c] = (
                        previous + row[c]
                    )
            else:
                results_df = results_df.append(row)
    results_df.to_csv(
        f"{settings.analysis_path}{settings.experiment}_gene_counts.csv", sep="\t"
    )

    return


def main():
    """Start here."""
    samples = f"results/{'example02'}/samples.yml"
    with open(samples, "r") as f:
        organize_samples(f)


if __name__ == "__main__":
    main()
