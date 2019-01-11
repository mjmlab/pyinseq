#!/usr/bin/env python3
"""
Counts the bowtie hits at each position in each sample

"""

import csv


def map_sites(sample, samples_dict, settings):
    """Map insertions to nucleotide sites."""
    # Placeholder for dictionary of mapped reads in format:
    # {(contig, position) : [left_counts, right_counts]}
    map_dict = {}
    # overall_total = denominator for cpm calculation
    overall_total = 0
    cpm = 0
    bowtie_file = settings.path + sample + "_bowtie.txt"
    sites_file = settings.path + sample + "_sites.txt"
    with open(bowtie_file, "r") as fi:
        for line in fi:
            bowtie_data = line.rstrip().split("\t")
            # Calculate transposon insertion point = transposonNT
            contig, insertion_NT, read_length = (
                str(bowtie_data[2]),
                int(bowtie_data[3]),
                len(bowtie_data[4]),
            )
            if bowtie_data[1] == "+":  # positive strand read
                insertion_NT = insertion_NT + read_length - 1
                map_dict.setdefault((contig, insertion_NT), [0, 0])[0] += 1  # Lcount
            else:  # negative strand read
                insertion_NT = insertion_NT + 1
                map_dict.setdefault((contig, insertion_NT), [0, 0])[1] += 1  # Rcount
            overall_total += 1
    # write tab-delimited of contig/nucleotide/left_counts/right_counts/total_counts/cpm
    # use the index totalCounts as the denominator for cpm calculation
    with open(sites_file, "a") as fo:
        writer = csv.writer(fo, delimiter="\t", dialect="excel")
        header_entry = (
            "contig",
            "nucleotide",
            "left_counts",
            "right_counts",
            "total_counts",
            "cpm",
        )
        writer.writerow(header_entry)
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
    return map_dict


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
    sites_file = settings.path + sample + "_sites.txt"
    genes_file = settings.path + sample + "_genes.txt"
    # TODO(minimum counts and maximum ratio)\
    # min_counts = settings.min_counts
    # max_ratio = settings.max_ratio
    with open(sites_file, "r", newline="") as csvfileR:
        sites_reader = csv.reader(csvfileR, delimiter="\t")
        next(sites_reader, None)  # skip the headers
        for line in sites_reader:
            contig, nucleotide, left_counts, right_counts, total_counts, cpm = line[0:6]
            nucleotide = int(nucleotide)
            cpm = float(cpm)
            # Used to save previous feature for intergenic calling
            previous_feature = ""
            for gene in genome:
                locus, start, end, strand, length, pid, gene, locus_tag, code, cog, product = (
                    gene[0],
                    int(gene[1]),
                    int(gene[2]),
                    *gene[3:11],
                )
                # contig from insertion; locus from lookup table
                if contig == locus:
                    # TODO: simplify nucleotide location match with gene
                    if nucleotide >= start:
                        if nucleotide <= end:
                            # 0.0 = 5'end ; 1.0 = 3'end
                            # TODO: Should featureEnd have +1 added?
                            if strand == "+":
                                three_primeness = (nucleotide - start) / (end - start)
                            if strand == "-":
                                three_primeness = (end - nucleotide) / (end - start)
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
        mapped_gene_writer = csv.writer(csvfileW, delimiter="\t")
        mapped_gene_writer.writerow(headers)
        for hit in mapped_hit_list:
            mapped_gene_writer.writerow(hit)
    # Return aggregated insertions by gene (filtered on 5'-3' threshold)
    return gene_dict


def build_gene_table(organism, sample_dict, gene_mappings, experiment=""):
    """
       For each entry in a feature table (.ftt) list the summary of hits
       for each sample in the experiment
    """

    # TODO(Bring back in the header row in the future. Use it here; ignore it for previous steps)
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
        # mapped_genes = gene_mappings.get(sampleName)
        # TODO(simplify this:)
        # for each row in the table, *try* from mapped_genes. Add if found.
        for gene in mapped_genes:
            for i, f in enumerate(gene_table):
                hit_locus_tag = gene
                ftt_locus_tag = f[7]
                # matches based on locusTag.
                # In future should I instead create an index field in the .ftt?
                if hit_locus_tag == ftt_locus_tag:
                    gene_table[i][current_column] += mapped_genes[gene][0]
        with open(f"results/{experiment}/summary_gene_table.txt", "w") as fo:
            writer = csv.writer(fo, delimiter="\t", dialect="excel")
            writer.writerows(gene_table)


def ftt_lookup(organism, experiment=""):
    """Import the ftt file and process as a dictionary of lookup values
       indexed on Synonym (i.e., Locus Tag)
       {'VF_0001': {'locus': 'CP000020', 'start': ...},
           'VF_0002': {'locus': 'CP000020', 'start': ...}}
    """

    # TODO: Error checking when generating the ftt file that locus tags are \
    # unique and complete.
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
