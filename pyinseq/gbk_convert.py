#!/usr/bin/env python3
"""

GenBank conversion utilities for the pyinseq pipeline.

# gbk2fna()
Convert GenBank sequence to fasta sequence
Multilocus GenBank converts to one multifasta GenBank file
Locus headers are the fasta headers
Maintains original newlines (typically leaving up to 60 nucleotides per line)
File is written to temp/ directory for the EXPERIMENT:
    EXPERIMENT/temp/genome.fna

# gbk2ftt()
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

"""

import csv
import logging
import os
import re
import sys

logger = logging.getLogger("pyinseq")


def gbk2fna(infile, organism, output_directory=""):
    """Convert genbank format to fna format."""
    with open(infile, "r") as fi:
        out_file = f"{output_directory}{organism}.fna"
        print(f"  Nucleotide file output file: {out_file}")
        if not os.path.exists(f"{output_directory}"):
            print(f"Error: {output_directory} directory was not created.")
            exit(1)
        with open(out_file, "w") as fo:
            dna_seq = False  # in the DNA sequence of the file
            for i, line in enumerate(fi):

                # Don't parse blank lines
                if line.strip():
                    parts = line.split()

                    # Locus (replicon) as header
                    if parts[0] == "LOCUS":
                        locus = parts[1]
                        fo.write(f">{locus}\n")

                    # DNA Sequence
                    if parts[0] == "//":
                        dna_seq = False
                    if dna_seq:
                        sequence = "".join(
                            n for n in line.strip() if n.isalpha())
                        fo.write(f"{sequence}\n")
                    if parts[0] == "ORIGIN":
                        dna_seq = True


def gbk2ftt(infile, organism, output_directory=""):
    """Convert genbank format to ptt-like ftt format."""
    with open(infile, "r") as fi:
        out_file = f"{output_directory}{organism}.ftt"
        print(f"  Feature table output file: {out_file}")
        with open(out_file, "w") as fo:
            writer = csv.writer(fo, delimiter="\t", dialect="excel")
            header = (
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
            writer.writerow(header)

            # Initialize variables
            features = False  # in the FEATURES section of the GenBank file
            new_feature = False  # collecting data for a new feature
            parse_types = ["CDS", "tRNA", "rRNA", "misc_RNA"]
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
                                    output = (
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
                                    writer.writerow(output)
                                    new_feature = False

                        if line[5:21].rstrip() in parse_types:
                            new_feature = True  # Feature that should be written
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
                            location = re.search(
                                r"(\d+)\.+.*\.(\d+)", parts[1])
                            first = location.group(1)
                            last = location.group(2)
                            try:
                                location_raw = f"{first}..{last}"
                                length = int(last) - int(first) + 1
                            except AttributeError:
                                error_complex_feature = (
                                    f"PyINSeq Error: Complex feature coordinates at or near {locus_tag} "
                                    "in GenBank file. Additional attention required."
                                )
                                print(error_complex_feature)
                                exit(0)
                            strand = "-" if parts[1].startswith(
                                "complement") else "+"

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


def main():
    """Start here."""
    input_file = sys.argv[1]
    organism = sys.argv[2]
    # gbk2fna(inputfile, organism)
    gbk2ftt(input_file, organism)


if __name__ == "__main__":
    main()
