"""
Code not in use in the pipeline

Not needed or deprecated.
"""

def multifasta2dict(fna):
    """
    Returns a dictionary of fasta headers and sequences

    Input file in multifasta format
    >header1
    agctag
    >header2
    gattta

    Function returns the fasta data as a dictionary
    {'header1': 'agctag', 'header2': 'gattta'}
    """
    d = {}
    with open(fna, 'r') as fi:
        header = ''
        firstLine = False # next line is first line of sequence in contig
        for line in fi:
            line = line.rstrip()
            if line.startswith('>'): # header line; therefore new contig
                header = line[1:] # new header
                firstLine = True
            else: # sequence line
                if firstLine:
                    d[header] = line # create dictionary entry
                    firstLine = False
                # would this work in one step?
                else:
                    appendedSequence = d[header] + line
                    d[header] = appendedSequence
        return d



def taSites(fna):
    """ Enumerate TA dinucleotides in a fasta nucleotide file

    Returns a dictionary of 5'-TA positions for each contig in a multifasta file.
    Assumes circular genome; i.e. will check if the last [-1] nucleotide is a T
    and the first [0] is an A

    Function calls multifasta2dict() and returns the nucleotide positions
    of TA sites in each contig as dictionary values:
    {'header1': [20, 22, 24, ... 1401104], 'header2': [5, 16, 24, ... 39910]}
    """
    di = multifasta2dict(fna)
    do = {}
    for header in di:
        sequence = di[header]
        do[header] = []
        i = 0   # nucleotide index
        # Checks the linear molecule (i.e. all except the last base)
        while i < len(sequence)-1:
            if sequence[i:i+2] == 'ta':
                do[header].append(i+1)
            i += 1
        # Checks last base circular around to first base
        if i == len(sequence)-1:
            if (sequence[-1] + sequence[0]) == 'ta':
                do[header].append(i+1)
    return do


def createSamplesDirectory():
    """
    Create the samples directory if needed
    """

    # Create /pyinseq/samples/ path if does not exist already
    try:
        os.makedirs('samples/')
    except:
        pass

##  Not filtering now, showing all results
##  Will re-implement filtering later, as a formula in another function
def filterSortCounts(experiment, sample):
    """ Filter for min 1 L read, 1 R read and maximum 10-fold L/R differential

    Also sort by the first four fields

    Data that do not pass this filter are treated as artefacts and not
    used in subsequent steps, such as data normalization.

    """
    with open('{0}/{1}_output_bowtie_mapped.txt'.format(experiment, sample), 'r') as fi:
        li = []
        for line in fi:
            # contig/nucleotide/Lcount/Rcount/TotalCount
            tupi = tuple(line.rstrip().split('\t'))
            li.append(tupi)
        #print(li)

        # FILTER COUNT DATA
        lo = []
        for l in li:
            contig, nucleotide, Lcounts, Rcounts, totalCounts = l[0:5]
            lo.append(l)
            # TODO: RETURN TO FILTERING BELOW
            """# minimum 1 read in each direction.
            if int(Lcounts) >= 1 and int(Rcounts) >= 1:
                # maximum 10-fold L/R differential
                # L=1 R=10 ok, but not L=1 R=11
                if not (11 * min(Lcounts, Rcounts)) < (totalCounts):
                    lo.append(l)"""

        # SORT COUNT DATA
        loSorted = sorted(lo, key=itemgetter(0,1,2,3))

        # Write sorted/filtered data to tab-delimited file
        # experiment/sample/contig/nucleotide/Lcount/Rcount/totalCount
        with open('{0}/{1}_output_bowtie_mapped_filtered.txt'.format(experiment, sample), 'w') as fo:
            for e in loSorted:
                for x in e:
                    fo.write('{0}\t'.format(x.strip()))
                fo.write('\n')


def normalizeCpm(experiment, sample):
    """ Normalize every sample to 1E6 CPM
    Returns a list of tuples:
    [
    ('contig1', '999401', 80.00268809031984)
    ]
    contig, nucleotide, normalized_counts
    """
    # Total count of filtered, mapped reads for each sample
    countsDict = {}
    # list of tuples of each mapped insertion, counted
    with open('{0}/genome_lookup/{1}_bowtie_mapped.txt'.format(experiment, sample), newline='') as csvfile:
        readsReader = csv.reader(csvfile, delimiter='\t')
        for line in readsReader:
            contig, nucleotide, totalCounts = line[0], line[1], line[4]


    for tupi in counts:
        #Data for each insertion
        experiment = tupi[0]
        barcode = tupi[1]
        readCount = int(tupi[6])
        # Add to the total counts for that sample (barcode)
        ddTotalCountBySample[(experiment, barcode)] += readCount

    # TODO: Write counts for logging
    """
    for entry in ddTotalCountBySample:
        experiment, barcode = entry
        print('{exp}\t{bc}\t{counts}'.
            format(exp=experiment, bc=barcode, counts=str(ddTotalCountBySample[entry])))
    """

    # Transform the original counts data (total only) by dividing each count by
    # the total for the sample and multiplying by 1E6.
    # resulting dataset normCountsAll. For each tuple:
    #[0] = experiment
    #[1] = barcode
    #[2] = contig
    #[3] = TAnucleotide
    #[4] = countTotal (from tupi[6] totalCounts)
    #  Note that Left and Right counts are not carried through.
    normCountsAll = []
    for tupi in counts:
        # note: can add back rawCounts if it would be valuable
        # to have them in the output
        rawCounts = int(tupi[6])
        sample = (tupi[0], tupi[1])
        totalSampleCounts = ddTotalCountBySample[sample]
        normCounts = float(1E6) * rawCounts / totalSampleCounts
        newTup = (tupi[0], tupi[1], tupi[2], tupi[3], normCounts)
        normCountsAll.append(newTup)
    return normCountsAll

# This nicely converts the genome FTT into a dictionary indexed on locus_tag
# However, it is not ordered...
def fttLookup(organism, experiment=''):
    """
    Import the ftt file and process as a dictionary of lookup values
    indexed on Synonym (i.e., Locus Tag)
    {'VF_0001': {'locus': 'CP000020', 'start': ...},
        'VF_0002': {'locus': 'CP000020', 'start': ...}}

    """
    # TODO: Error checking when generating the ftt file that locus tags are \
    # unique and complete.
    fttDict = {}
    with open('{0}/genome_lookup/{1}.ftt'.format(experiment, organism), newline='') as csvfile:
        fttreader = csv.reader(csvfile, delimiter='\t')
        for line in fttreader:
            # ignore header row
            if line[0] != ('Locus'):
                Locus, Location_Start, Location_End, Strand, Length, PID, \
                    Gene, Synonym, Code, COG, Product = \
                line[0], line[1], line[2], line[3], line[4], line[5], \
                    line[6], line[7], line[8], line[9], line[10]
                fttDict[Synonym] = {
                    'locus': Locus,
                    'start': Location_Start,
                    'end': Location_End,
                    'strand': Strand,
                    'length': Length,
                    'pid': PID,
                    'gene': Gene,
                    'locus_tag': Synonym,
                    'code': Code,
                    'cog': COG,
                    'product': Product
                    }
    return fttDict
