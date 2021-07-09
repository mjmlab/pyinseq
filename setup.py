#!/usr/bin/env python3

import os
import sys
from pathlib import Path

try:
    from setuptools import *
except ImportError:
    from distribute_setup import use_setuptools

    use_setuptools()
finally:
    from setuptools import *

from glob import glob

if sys.version_info < (3, 6):
    print >>sys.stderr, "ERROR: pyinseq requires python 3.6 or greater"
    sys.exit()

__version__ = open(os.path.join("pyinseq", "VERSION")).read().strip()

SCRIPTS = glob("scripts/*")

with open("README.md", "r") as fh:
    long_description = fh.read()


# Collect paths to Snakefile's and ENVS
snake_data = [
            "pyinseq/workflows/PyinseqWorkflow/Snakefile",
            "pyinseq/workflows/DemultiplexWorkflow/Snakefile",
            "pyinseq/workflows/GenomeprepWorkflow/Snakefile",
            "pyinseq/envs/bowtie.yaml",
]

def main():
    setup(
        name="pyinseq",
        version=__version__,
        description="Analysis of transposon insertion sequencing (INSeq) data in Python",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/mjmlab/pyinseq",
        author="Mark J. Mandel,Emanuel Burgos,Camille Scott,Loren Velasquez,Benjamin K. Johnson",
        author_email="mandel01@gmail.com",
        license="BSD",
        packages=find_packages(),
        include_package_data=True,
        package_dir={'pyinseq': 'pyinseq'},
        package_data={'pyinseq': snake_data},
        scripts=SCRIPTS,
        setup_requires=["pytest-runner"],
        tests_require=["pytest"],
        python_requires='!=3.8.*,>=3.6.*',
        classifiers=[
            "Development Status :: 4 - Beta",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: BSD License",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: POSIX",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
        zip_safe=False,
    )


if __name__ == "__main__":
    main()
