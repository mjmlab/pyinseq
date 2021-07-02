#!/usr/bin/env python3

"""

Count and process the bowtie hits at each position in each sample

"""

import csv

# Module Imports
from pyinseq.logger import pyinseq_logger

logger = pyinseq_logger.logger


def map_sites(sample, settings):
    """Map insertions to nucleotide sites."""
    # Placeholder for dictionary of mapped reads in format:
    # {(contig, position) : [left_counts, right_counts]}
    map_dict = {}
    # overall_total = denominator for cpm calculation
    overall_total = 0
    bowtie_file = settings.path.joinpath(f"{sample}_bowtie.txt")
    sites_file = settings.path.joinpath(f"{sample}_sites.txt")
    logger.info(
        f"Sample {sample}: Tally site data from {bowtie_file.name} into {sites_file.name} file."
    )
    with open(bowtie_file, "r") as fi:
        for line in fi:
            bowtie_data = line.rstrip().split("\t")
            # Calculate transposon insertion point = insertion_NT
            contig, insertion_NT, read_length = (
                str(bowtie_data[2]),
                int(bowtie_data[3]),
                len(bowtie_data[4]),
            )
            if bowtie_data[1] == "+":  # positive strand read
                insertion_NT = insertion_NT + read_length - 1
                map_dict.setdefault((contig, insertion_NT), [0, 0])[0] += 1  # left
            else:  # negative strand read
                insertion_NT = insertion_NT + 1
                map_dict.setdefault((contig, insertion_NT), [0, 0])[1] += 1  # right
            overall_total += 1
    # write tab-delimited of contig/nucleotide/left_counts/right_counts/total_counts/cpm
    # use the total_counts as the denominator for cpm calculation
    with open(sites_file, "a") as fo:
        writer = csv.writer(fo, delimiter="\t", dialect="excel", lineterminator="\n")
        header_entry = (
            "contig",
            "nucleotide",
            "left_counts",
            "right_counts",
            "total_counts",
            "cpm",
        )
        writer.writerow(header_entry)
        sample_dict = {sample: {"site hits": 0}}
        for insertion in sorted(map_dict):
            left_counts = map_dict[insertion][0]
            right_counts = map_dict[insertion][1]
            total_counts = map_dict[insertion][0] + map_dict[insertion][1]
            cpm = float(1e6) * total_counts / overall_total
            row_entry = (
                insertion[0],
                insertion[1],
                left_counts,
                right_counts,
                total_counts,
                cpm,
            )
            writer.writerow(row_entry)
            sample_dict[sample]["site hits"] += total_counts
    # Summarize site data into io
    pyinseq_logger.logger_io.write(
        f"- {sample}: {overall_total} aligned reads mapped to sites\n"
    )
    # Return a dict where {sample: {sites: 100} }
    return sample_dict


def map_genes(sample, settings):
    """Maps insertions to genes

    1. Writes a csv file listing the gene for each insertion with the tabs:
    contig, nucleotide, Lcounts, Rcounts, totalCounts, cpm, threePrimeness, locus_tag
    2. Returns a dictionary of aggregate counts per gene (for any gene with at
    least one hit)

    Insertions that map to multiple genes result in multiple lines in the CSV file
    and are counted for both genes in the returned dictionary.

    ThreePrimeness = insertion location in gene (5' end = 0.0, 3' end = 1.0)
    All insertions are written to the file but only ones <= disruption threshold
    are counted in the dictionary.
    """
    # List of tuples of genome features
    genome = ftt_lookup(settings.organism, settings.experiment)
    # list of tuples of each mapped insertion to be immediately written per insertion
    mapped_hit_list = []
    # Dictionary with running total of cpm per gene; keys are genes, values are aggregate cpm
    # only hits in the first part of the gene are added to the count, as defined
    # by the disruption threshold.
    # if disruption = 1.0 then every hit in the gene is included
    gene_dict = {}
    sites_file = settings.path.joinpath(f"{sample}_sites.txt")
    genes_file = settings.path.joinpath(f"{sample}_genes.txt")
    second_message = f"Transposon hits that disrupt the 5-prime-most {round(settings.disruption * 100)}% of each gene are tallied"
    if settings.disruption == 1:
        second_message = (
            "Transposon hits that disrupt any position in a gene are tallied"
        )
    logger.info(
        f"Sample {sample}: Tally site data from {sites_file.name} to gene-level data in {genes_file.name}. {second_message}"
    )
    with open(sites_file, "r", newline="") as csvfileR:
        sites_reader = csv.reader(csvfileR, delimiter="\t")
        next(sites_reader, None)  # skip the headers
        sites_in_genes = 0
        for line in sites_reader:
            contig, nucleotide, left_counts, right_counts, total_counts, cpm = line[0:6]
            nucleotide = int(nucleotide)
            cpm = float(cpm)
            # Used to save previous feature for intergenic calling
            previous_feature = ""
            for gene in genome:
                (
                    locus,
                    start,
                    end,
                    strand,
                    length,
                    pid,
                    gene,
                    locus_tag,
                    code,
                    cog,
                    product,
                ) = (
                    gene[0],
                    int(gene[1]),
                    int(gene[2]),
                    gene[3],
                    int(gene[4]),
                    *gene[5:11],
                )
                # contig from insertion; locus from lookup table
                if contig == locus:
                    if start <= nucleotide <= end:
                        sites_in_genes += int(total_counts)
                        # 0.0 = 5'end ; 1.0 = 3'end
                        if strand == "+":
                            three_primeness = (1 + nucleotide - start) / length
                        else:
                            three_primeness = (end - nucleotide) / length
                        mappedHit = (
                            contig,
                            nucleotide,
                            left_counts,
                            right_counts,
                            total_counts,
                            cpm,
                            three_primeness,
                            locus_tag,
                        )
                        # Checks minimum count and max ratio
                        if int(total_counts) >= settings.min_counts and (
                            min(int(left_counts), int(right_counts))
                            * settings.max_ratio
                        ) >= max(int(left_counts), int(right_counts)):
                            mapped_hit_list.append(mappedHit)
                            # Filter based on location in the gene
                            if three_primeness <= settings.disruption:
                                # Add to the total for that gene --
                                # Single-element list (rather than integer) so
                                # that it is subscriptable to add cpm counts
                                gene_dict.setdefault(locus_tag, [0])[0] += cpm
                previous_feature = locus_tag
    # Write individual insertions to *_genes.txt
    sample_dict = {sample: {"gene hits": 0}}
    with open(genes_file, "w", newline="") as csvfileW:
        headers = (
            "contig",
            "nucleotide",
            "left_counts",
            "right_counts",
            "total_counts",
            "cpm",
            "three_primeness",
            "locus_tag",
        )
        mapped_gene_writer = csv.writer(csvfileW, delimiter="\t", lineterminator="\n")
        mapped_gene_writer.writerow(headers)
        for hit in mapped_hit_list:
            mapped_gene_writer.writerow(hit)
            sample_dict[sample]["gene hits"] += int(hit[4])
    # Summarize site data into io
    pyinseq_logger.logger_io.write(
        f"- {sample}: {sites_in_genes} of mapped reads fall in genes\n"
    )
    # Return aggregated insertions by gene (filtered on 5'-3' threshold)
    return sample_dict


def build_gene_table(organism, sample_dict, gene_mappings, experiment=""):
    """
    For each entry in a feature table (.ftt) list the summary of hits
    for each sample in the experiment
    """
    logger.info("Aggregate gene mapping from all samples into the summary_data_table")

    gene_table = ftt_lookup(organism, experiment)

    # Header will be extended in the future to
    # list the experiment and barcode of each sample of interest
    header = [
        "Contig",
        "Start",
        "End",
        "Strand",
        "Length",
        "PID",
        "Gene",
        "Synonym",
        "Code",
        "COG",
        "Product",
    ]

    # Add header row to ftt file
    gene_table.insert(0, header)

    # Get individual mapped hits
    # mappedHitList = mapToGene(organism, experiment)

    # current column, sample being matched
    current_column = len(gene_table[0]) - 1

    for sample in sample_dict:
        # Add the new sample name as a new column the table
        current_column += 1
        gene_table[0].append(sample)
        # Fill the rest of the new column with 0 as the count for each gene
        for row in gene_table[1:]:
            row.append(0)
        # add the sample's results to the building gene_table
        mapped_genes = gene_mappings[sample]
        # for each row in the table, *try* from mapped_genes. Add if found.
        for gene in mapped_genes:
            for i, f in enumerate(gene_table):
                hit_locus_tag = gene
                ftt_locus_tag = f[7]
                # matches based on locus_tag.
                if hit_locus_tag == ftt_locus_tag:
                    gene_table[i][current_column] += mapped_genes[gene][0]
        with open(f"results/{experiment}/summary_gene_table.txt", "w") as fo:
            writer = csv.writer(
                fo, delimiter="\t", dialect="excel", lineterminator="\n"
            )
            writer.writerows(gene_table)

    # Summarize build gene table step
    pyinseq_logger.logger_io.write(
        f"- Gene table contains {len(sample_dict)} samples (columns) "
        f"for {len(gene_table) - 1} genes (rows) \n"
    )
    return


def ftt_lookup(organism, experiment=""):
    """Import the ftt file and process as a dictionary of lookup values
    indexed on Synonym (i.e., Locus Tag)
    {'VF_0001': {'locus': 'CP000020', 'start': ...},
        'VF_0002': {'locus': 'CP000020', 'start': ...}}
    """
    ftt_list = []
    with open(
        f"results/{experiment}/genome_lookup/{organism}.ftt", newline=""
    ) as csv_file:
        ftt_reader = csv.reader(csv_file, delimiter="\t")
        for line in ftt_reader:
            # ignore header row
            if line[0] != ("Locus"):
                # Locus, Location_Start, Location_End, Strand, Length, PID,
                # Gene, Synonym, Code, COG, Product
                feature_data = line[0:11]
                ftt_list.append(feature_data)
    return ftt_list


def main():
    """Start here."""
    pass


if __name__ == "__main__":
    main()
