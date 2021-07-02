#!/usr/bin/env python3

"""

Tests for running pyinseq pipelines

"""

from .test_utils import runscript, datadir, compare_directories


def test_pyinseq_script_no_args(datadir, tmpdir):
    args = []
    status, out, err = runscript("pyinseq", args, directory=str(tmpdir))
    assert out.split("\n")[1][0:5] == "usage"


def test_pyinseq_script(datadir, tmpdir):
    input_fn = datadir("input/example01.fastq")
    sample_fn = datadir("input/example01.txt")
    gb_fn = datadir("input/ES114v2.gb")
    output_name = "example_pyinseq"
    expected_output = datadir("output_pyinseq")
    output_dir = tmpdir.join("results/example_pyinseq")

    args = [
        "-i",
        input_fn,
        "-s",
        sample_fn,
        "-g",
        gb_fn,
        "-e",
        output_name,
        "-d",
        "0.9",
        "--additional_params",
        "--use-conda",
    ]
    status, out, err = runscript("pyinseq", args, directory=str(tmpdir))
    assert not status

    # Compare directory outputs
    compare_directories(
        expected_output,
        output_dir,
        ignore=[
            "E001_01_bowtie.txt",
            "E001_02_bowtie.txt",
            ".DS_Store",
            "log.txt",
            "summary_log.txt",
        ],
    )


def test_pyinseq_script_unique_transp_barcode_length(datadir, tmpdir):
    input_fn = datadir("input/example03.fastq")
    sample_fn = datadir("input/example03.txt")
    gb_fn = datadir("input/ES114v2.gb")
    output_name = "example_pyinseq_unique_transposon_barcode"
    expected_output = datadir("output_pyinseq_unique_transposon_barcode_length")
    output_dir = tmpdir.join("results/example_pyinseq_unique_transposon_barcode")

    args = [
        "-i",
        input_fn,
        "-s",
        sample_fn,
        "-g",
        gb_fn,
        "-e",
        output_name,
        "-d",
        "0.9",
        "--barcode_length",
        "9",
        "--transposon_seq",
        "TCGCACGG",
        "--additional_params",
        "--use-conda",
    ]
    status, out, err = runscript("pyinseq", args, directory=str(tmpdir))
    assert not status

    # Compare directory outputs
    compare_directories(
        expected_output,
        output_dir,
        ignore=[
            "E001_01_bowtie.txt",
            "E001_02_bowtie.txt",
            ".DS_Store",
            "log.txt",
            "summary_log.txt",
        ],
    )


def test_pyinseq_script_unique_transp_barcode_length(datadir, tmpdir):
    input_fn = datadir("input/example03.fastq")
    sample_fn = datadir("input/example03.txt")
    gb_fn = datadir("input/ES114v2.gb")
    output_name = "example_pyinseq_unique_transposon_barcode"
    expected_output = datadir("output_pyinseq_unique_transposon_barcode_length")
    output_dir = tmpdir.join("results/example_pyinseq_unique_transposon_barcode")

    args = [
        "-i",
        input_fn,
        "-s",
        sample_fn,
        "-g",
        gb_fn,
        "-e",
        output_name,
        "-d",
        "0.9",
        "--barcode_length",
        "9",
        "--transposon_seq",
        "TCGCACGG",
    ]
    status, out, err = runscript("pyinseq", args, directory=str(tmpdir))

    assert not status
    dcmp = filecmp.dircmp(
        expected_output,
        output_dir,
        ignore=["E001_01_bowtie.txt", "E001_02_bowtie.txt", ".DS_Store", "log.txt"],
    )

    # check that log file is created from pyinseq
    assert "log.txt" in os.listdir(output_dir)

    # checks that files are same in both directories
    assert not dcmp.left_only and not dcmp.right_only

    # check files to see if content differs
    assert not dcmp.diff_files
    assert not dcmp.funny_files  # Check for files that cannot be compared

    # because subdirs is a dict keyed by subdir name
    # with dircmp objects as values
    for subdcmp in dcmp.subdirs.values():
        assert not subdcmp.diff_files
        assert not subdcmp.funny_files
        assert not subdcmp.left_only and not subdcmp.right_only


def test_pyinseq_demultiplex_script(datadir, tmpdir):
    input_fn = datadir("input/example01.fastq")
    sample_fn = datadir("input/example01.txt")
    output_name = "test_demultiplex"
    expected_output = datadir("output_demultiplex")
    output_dir = tmpdir.join(f"results/{output_name}")

    args = ["demultiplex", "-i", input_fn, "-s", sample_fn, "-e", output_name]
    status, out, err = runscript("pyinseq", args, directory=str(tmpdir))
    assert not status

    # Compare directory outputs
    compare_directories(
        expected_output, output_dir, ignore=[".DS_Store", "log.txt", "summary_log.txt"]
    )


def test_pyinseq_demultiplex_notrim_script(datadir, tmpdir):
    input_fn = datadir("input/example01.fastq")
    sample_fn = datadir("input/example01.txt")
    output_name = "test_demultiplex_notrim"
    expected_output = datadir("output_demultiplex_notrim")
    output_dir = tmpdir.join(f"results/{output_name}")

    args = [
        "demultiplex",
        "-i",
        input_fn,
        "-s",
        sample_fn,
        "-e",
        output_name,
        "--notrim",
    ]
    status, out, err = runscript("pyinseq", args, directory=str(tmpdir))
    assert not status

    # Compare directory outputs
    compare_directories(
        expected_output, output_dir, ignore=[".DS_Store", "log.txt", "summary_log.txt"]
    )


def test_pyinseq_genomeprep_script(datadir, tmpdir):
    output_name = "test_genomeprep"
    expected_output = datadir("output_genomeprep")
    output_dir = tmpdir.join(f"results/{output_name}")
    gb_fn = datadir("input/ES114v2.gb")

    args = [
        "genomeprep",
        "-e",
        output_name,
        "-g",
        gb_fn,
        "--additional_params",
        "--use-conda",
    ]
    status, out, err = runscript("pyinseq", args, directory=str(tmpdir))
    assert not status

    # Compare directory outputs
    compare_directories(expected_output, output_dir, ignore=[".DS_Store", "log.txt"])


def test_pyinseq_genomeprep_gff_script(datadir, tmpdir):
    output_name = "test_genomeprep_gff"
    expected_output = datadir("output_genomeprep_gff")
    output_dir = tmpdir.join(f"results/{output_name}")
    gb_fn = datadir("input/ES114v2.gb")

    args = [
        "genomeprep",
        "-e",
        output_name,
        "-g",
        gb_fn,
        "--gff",
        "--additional_params",
        "--use-conda",
    ]
    status, out, err = runscript("pyinseq", args, directory=str(tmpdir))
    assert not status

    # Compare directory outputs
    compare_directories(expected_output, output_dir, ignore=[".DS_Store", "log.txt"])
