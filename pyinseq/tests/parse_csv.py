#!/usr/bin/env python
import csv

def read_data(data):
    with open(data, 'r') as f:
        data = [row for row in csv.reader(f.read().splitlines(), delimiter='\t')]
    return data

data = 'results/example01/summary_gene_table.txt'
read_data(data)
