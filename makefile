## Makefile to run the example script and the test script

install:
	python setup.py install

test:
	py.test pyinseq/tests/

