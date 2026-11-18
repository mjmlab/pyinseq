## Makefile to install `pyinseq` and run the test suite

install:
	python -m pip install .

install-dev:
	python -m pip install -e ".[dev]"

test:
	python -m pytest --cov=pyinseq

lint:
	black --check pyinseq/

clean_dump:
	rm -rf pyinseq/tests/dump/*
	rm -rf pyinseq/tests/dump/.snakemake

get_dag:
	snakemake --dag -s ../../workflows/PyinseqWorkflow/Snakefile --cores 2 --configfile test_pyinseq-config.yaml | dot -Tpng > pyinseq-snakemake_dag.png
