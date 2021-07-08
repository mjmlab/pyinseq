# Changelog

## [0.3.0] - 2021-07-08
### Fixed
- Refactor demultiplex.py for readability and faster execution
- Logger is now compatible with snakemake logging module
- Fixed small parsing bugs in gbk_convert.py

### Changed
- Remove third party folder which contained 'bowtie' executables
- file paths are now handled by pathlib
- New settings class that is compatible with snakemake
- The map_reads script is now a stand-alone script
- Settings class which is used by snakemake and pyinseq
- Pyinseq and snakemake now use `pathlib` for handling filepaths
- Removed readthedocs website
- Removed docs/
- Removed third-party/

### Added
- Snakemake can now execute workflows using Snakefile in pyinseq
- Additional Snakefiles for genomprep and demultiplex of pyinseq
- Conda virtual environment files that allow conda to install `bowtie` software during runtime
- Parser.py script for parsing command line arguments
- Include '--test' option for running pytest on pyinseq
- New pytests for testing workflows in pyinseq
- Pyinseq will now summarize sample information and snakemake output
- New user guide for pyinseq


## [0.2.1] - 2021-06-02
### Fixed
- `pyinseq` alone brings up the help documentation
- Small fix to the `three_primeness` calculation. 
  A minimum of 3 reads is now required per site, and a Left:Right max read ratio of 10-fold to be tallied.
  
### Changed
- Only Python 3.6 and 3.7 are supported.
- `screed` module is used for opening/writing fastq files.

### Added
- `pyinseq genomeprep` subcommand will prepare genome files for pyinseq run. Also checks GenBank files before running.
- Added T50 calculation for sites files.
- Added progress bar for `demultiplex` function and for `writing` reads to sample files.
- `test_script.py` now compares directories and files from `pyinseq` runs to the expected output.
- Parameter `--min_counts`: minimum number of reads at a site required to be tallied. Default is 3
- Parameter `--max_ratio`: max ratio allowed between left and right reads around a TA insertion site. Default is 10-fold.
- Parameter `--transposon_seq`: define transposon sequence that is found at end of reads to help in demultiplexing. Default is ACAGGTTG
- Parameter `--barcode_length`: length of barcode index used to demultiplex reads into samples, allows for 4 - 16 nt. Default is 4.
- Parameter `--gff3`: enables `pyinseq` to write gff3 version of genome files.


## [0.2.0] - 2017-07-16
### Added
- Documentation hosting `http://pyinseq.readthedocs.io/`.
- `pyinseq demultiplex` command (including --notrim option).
- Improved logging of messages during the run.
- `log.txt` file to record log output to file.


## [0.1.0] - 2016-12-13
### Added
- Initial release.
- `pyinseq` command demultiplexes and maps samples.
- Automated testing using pytest, on TravisCI
