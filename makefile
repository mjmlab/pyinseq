## Makefile to intall `pyinseq` in the path and run the test script

install:
	python setup.py install

test:
	py.test --cov=./pyinseq/tests

clean_dump:
	rm -rf pyinseq/tests/dump/*

get_dag:
	snakemake --dag -s ../../workflows/PyinseqSnakemake/Snakefile --cores 2 --configfile test_pyinseq-snakemake_args_pyinseq-snake_config.yaml | dot -Tpng > pyinseq-snakemake_dag.png
