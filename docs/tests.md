# Automated tests

For each pull request [Travis CI](https://travis-ci.org/mandel01/pyinseq) runs some basic tests to ensure that the `pyinseq` produces expected results from the `example01` dataset.

# Manual tests

You can run the tests manually by first running the `example01` data:

```
python run.py -i data/example/example01.fastq -s data/example/example01.txt -g data/example/ES114v2.gb -e example01   
```

Then run the tests:

```
python -m pytest pyinseq/tests/
```
