#!/usr/bin/env python
'''
Analyzes resulting output
'''

from collections import OrderedDict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import yaml


def read_sites_file(sample, settings):
    return pd.read_csv(settings.path + sample + '_sites.txt', sep='\t')


def nfifty(sample, settings):
    '''N50: number of insertions that account for >=50% of counts'''
    # calculate sum of counts in d
    df = read_sites_file(sample, settings)
    return np.sum(df.total_counts.sort_values().cumsum() > (df.total_counts.sum() / 2))


def plot_insertions(sample, settings):
    '''Plot insertions across each contig/chromosome.'''
    df = read_sites_file(sample, settings)
    for i, group in df.groupby('contig'):
        plt.figure()
        group.plot(kind='scatter', x='nucleotide', y='cpm', color='blue', alpha=0.5, title=str(i))
        plt.savefig(settings.path + 'insertions_scatter_{0}_{1}.pdf'.format(sample, str(i)), format='pdf')



def df_raw_counts(samples):
    pass

def df_cpm_counts(samples):
    pass

def df_gene_counts(samples):
    pass

def pearson_correlation():
    pass

def main():
    '''Start here.'''
    samples = 'results/{}/samples.yml'.format('example02')
    with open(samples, 'r') as f:
        organize_samples(f)

if __name__ == "__main__":
    main()
