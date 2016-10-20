## Makefile to intall `pyinseq` in the path and run the test script

install:
	python setup.py install

test:
	py.test --cov=./ pyinseq/tests
