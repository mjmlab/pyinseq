## Makefile to intall `pyinseq` in the path and run the test script

install:
	python setup.py install

test:
	py.test --cov=./pyinseq/tests

clean_dump:
	rm -rf pyinseq/tests/dump/*
	rm -rf pyinseq/tests/dump/.snakemake

get_dag:
	snakemake --dag -s ../../workflows/PyinseqWorkflow/Snakefile --cores 2 --configfile test_pyinseq-config.yaml | dot -Tpng > pyinseq-snakemake_dag.png
