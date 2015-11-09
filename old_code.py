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
