# Once the pipeline is operational, all of these items will get transitioned to either Issues or the Wiki. These are my semi-organized ongoing notes and I haven't had a chance to sort through all of them yet.

# Highest priority / bugs
- [ ] Add command line flags (move away from argv)
- [ ] Remove need to name the organism; instead call the files genome.ftt etc. Parse the organism/contig names from the GenBank file and FNA file and add it to the log information.

# Documentation
- [ ] Add documentation on the `README.md` on how to install and run the software
- [ ] Name all of the filters/steps and use consistently in the software: e.g., sequence quality filter; minimum reads per insertion filter; left/right bias filter (need a better name here!)
- [ ] See run_info.json example here: http://pachterlab.github.io/kallisto/starting.html

# Sofware design/overall planning
- [ ] Need to do read/writes more efficiently overall. Chunks? Buffers?
- [ ] Would it be better to attach gene information (locus tag, annotation, % into gene) when processing the bowtie results? Then perform subsequent operations on a table that contains all of this information. **Maybe append the gene-level INSeq results in columnar form to an ftt-based table that also has the contig listed on each row?**
- [ ] How to keep `experiment` as a global variable that doesn't need to be passed to every module? Or is it a better idea to pass it each time so that the modules can be better adapted to other uses (both now and in the future)?
- [ ] Consider `pandas` dataframes to organize the data. Could process each individual data set on their own and then use `pandas.merge` to combine data when appropriate (e.g., combine different samples).
- [ ] List sample name instead of (in addition to?) barcode sequence in id line?
- [ ] Seems like key will be to not go through data so many times. Need to extract/filter/array data all at once instead in so many pieces.
- [ ] Use Python read/write csv instead of manual.
- [ ] Consider parsing on fixed width rather than on space-delimited for LOCUS, etc
- [ ] Add code to do left vs. right transposon ends (and plan new version with that in mind)

# `config.py`
- [ ] generate `bowtie-build` path from bowtie path
- [ ] Write a script to provide directions to the user, to get input back to populate this file, and to print the contents of this file?

# `processMapping.py`
- [ ] `filterSortCounts()` - allow the user to specify cutoff parameters
- [ ] `filterSortCounts()` - write a file of insertions that were discarded by the cutoffs

# `mapToGene()`
- [ ] Truncated product names. Is that from here or from `gbk2ftt()` ?
- [ ] Would be better to call the `gbk2ftt()` function directly and return the list from there. Don't repeat that work. That currently only writes to file though.

# `assign.py`
- [ ] Clean up code to be more readable; i.e., nreads = ((i+1)/4)

# `demultiplex.py`
- [ ] Use StringIO and buffering to avoid doing so many `write`s. Check speed of an optional demultiplexing step.

# map_reads.py

# gbkconvert.py
- [ ] Extract COG from /note field in Refseq GBK files?

# Low priority/next steps
- [ ] Multiprocessing/multithreading
