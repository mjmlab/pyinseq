#!/usr/bin/env python3
"""

GenBank conversion utilities for the pyinseq pipeline.

Files are written to genome_lookup/ directory for the EXPERIMENT:
    EXPERIMENT/genome_lookup/genome.fna

# gbk2fna()
Convert GenBank sequence to fasta sequence
Multilocus GenBank converts to one multifasta GenBank file
Locus headers are the fasta headers
Maintains original newlines (typically leaving up to 60 nucleotides per line)

# gbk2table()
Convert GenBank to feature table
Format similar to .ptt and .rnt files except:
- full tabular (locus as a field)
- start and end positions as separate fields
- includes the following features:
    CDS
    rRNA
    tRNA
    misc_RNA
- Multilocus GenBank converts to multi-.ftt file
'Unlike .ptt files that show the number of amino acids as 'length'

Optional convert to GFF3 format also.

"""

import csv
import logging
import re
import sys

logger = logging.getLogger("pyinseq")


def write_to_file(outfile, row_data):
    with open(outfile, "w") as fo:
        writer = csv.writer(fo, delimiter="\t", lineterminator="\n")
        for row in row_data:
            writer.writerow(row)


def gbk2fna(infile, organism, output_directory=""):
    """Convert genbank format to fna format."""
    with open(infile, "r") as fi:
        fna_rows = []
        n_nt = {}
        dna_seq = False  # in the DNA sequence of the file
        for i, line in enumerate(fi):

            # Don't parse blank lines
            if line.strip():
                parts = line.split()

                # Locus (replicon) as header
                if parts[0] == "LOCUS":
                    locus = parts[1]
                    n_nt[locus] = 0
                    fna_rows.append((f">{locus}",))

                # DNA Sequence
                if parts[0] == "//":
                    dna_seq = False
                if dna_seq:
                    sequence = "".join(n for n in line.strip() if n.isalpha())
                    fna_rows.append((f"{sequence}",))
                    n_nt[locus] += len(sequence.strip())
                if parts[0] == "ORIGIN":
                    dna_seq = True

        outfile = f"{output_directory}{organism}.fna"
        write_to_file(outfile, fna_rows)

        return {"fasta": fna_rows, "nucleotides": n_nt}


def gbk2table(infile, organism, output_directory="", gff=False):
    """Convert genbank format to feature table format (.ftt, and optionally .gff)."""
    with open(infile, "r") as fi:
        ftt_rows = [
            (
                "Locus",
                "Location_Start",
                "Location_End",
                "Strand",
                "Length",
                "PID",
                "Gene",
                "Synonym",
                "Code",
                "COG",
                "Product",
            )
        ]

        # Collect/write fasta sequence and collect contig lengths
        fasta = gbk2fna(infile, organism, output_directory)
        fna_rows = fasta["fasta"]
        fna_nucleotides = fasta["nucleotides"]

        if gff:
            gff_rows = [("##gff-version 3",)]
            for contig in fna_nucleotides:
                gff_rows.append(
                    (f"##sequence-region {contig} 1 {fna_nucleotides[contig]}",)
                )

        # Initialize variables
        features = False  # in the FEATURES section of the GenBank file
        new_feature = False  # collecting data for a new feature
        parse_types = ["CDS", "tRNA", "rRNA", "misc_RNA"]
        type = "CDS"
        strand = "+"
        length = 0
        protein_id = "-"
        gene = "-"
        locus_tag = "-"
        code = "-"
        cog = "-"
        product = "-"
        product_append = False  # append the current line to product

        for i, line in enumerate(fi):

            # Don't parse blank lines
            if line.strip():
                parts = line.split()

                # PRINT HEADER FOR THE LOCUS
                # 2 LINES:
                # LOCUS <tab> locus name
                # Location <tab> Strand etc...
                if parts[0] == "LOCUS":
                    locus = parts[1]

                if features:
                    # Print line before go on to next feature
                    # (gene, COG, protein id not required)
                    # Reset flags/defaults
                    if line[5:21].rstrip():
                        if new_feature:
                            # if locus_tag:
                            if not product_append:
                                ftt_rows.append(
                                    (
                                        locus,
                                        first,
                                        last,
                                        strand,
                                        str(length),
                                        protein_id,
                                        gene,
                                        locus_tag,
                                        code,
                                        cog,
                                        product,
                                    )
                                )
                                if gff:
                                    gff_id = (
                                        f"ID:{locus_tag}"
                                        if gene == "-"
                                        else f"ID:{locus_tag};Name:{gene}"
                                    )
                                    gff_rows.append(
                                        (
                                            locus,
                                            ".",
                                            type,  # need to extract the type of feature
                                            first,
                                            last,
                                            ".",
                                            strand,
                                            ".",
                                            gff_id,
                                        )
                                    )

                                new_feature = False

                    if line[5:21].rstrip() in parse_types:
                        new_feature = True  # Feature that should be written
                        type = line[5:21].rstrip()
                        protein_id = "-"
                        gene = "-"
                        locus_tag = "-"
                        code = "-"
                        cog = "-"
                        product = "-"

                        # NOTES ABOUT FEATURES
                        # 1. At ends of contigs greater than/less than signs
                        #    (> / <) are removed.
                        # 2. Complicated features use only the outer bounds
                        #    join(481257..481331,481333..482355) uses 481257..482355
                        location = re.search(r"(\d+)\.+.*\.(\d+)", parts[1])
                        first = location.group(1)
                        last = location.group(2)
                        try:
                            length = int(last) - int(first) + 1
                        except AttributeError:
                            error_complex_feature = (
                                f"PyINSeq Error: Complex feature coordinates at or near {locus_tag} "
                                "in GenBank file. Additional attention required."
                            )
                            print(error_complex_feature)
                            exit(0)
                        strand = "-" if parts[1].startswith("complement") else "+"

                    if "/protein_id=" in parts[0]:
                        protein_id = parts[0][13:-1]

                    if "/gene=" in parts[0]:
                        gene = parts[0][7:-1]

                    if "/locus_tag=" in parts[0]:
                        locus_tag = parts[0][12:-1]

                    # Multi-line product description
                    if product_append:
                        product = product + " " + line.strip()
                        if product.count('"') != 2:
                            product_append = True
                        # TODO(Error handling if not exactly 2 parentheses)
                        if product.count('"') == 2:
                            product = product.strip('"')
                            product_append = False

                    if "/product=" in parts[0]:
                        product = line.strip()[9:]
                        if product.count('"') != 2:
                            product_append = True
                        if product.count('"') == 2:
                            product = product.strip('"')

                if parts[0] == "ORIGIN":
                    features = False  # Not in FEATURES any more
                if parts[0] == "FEATURES":
                    features = True

    outfile = f"{output_directory}{organism}.ftt"
    write_to_file(outfile, ftt_rows)
    if gff:
        gff_rows.append(("##FASTA",))
        gff_rows = gff_rows + fna_rows
        outfile = f"{output_directory}{organism}.gff"
        write_to_file(outfile, gff_rows)


def main():
    """Start here."""
    input_file = sys.argv[1]
    organism = sys.argv[2]
    # gbk2fna(input_file, organism)
    gbk2table(input_file, organism, gff=True)


if __name__ == "__main__":
    main()
