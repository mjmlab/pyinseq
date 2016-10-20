#!/usr/bin/env python3

import sys, platform, os

try:
    from setuptools import *
except ImportError:
    from distribute_setup import use_setuptools
    use_setuptools()
finally:
    from setuptools import *

from glob import glob

if sys.version_info < (3, 5):
    print >> sys.stderr, "ERROR: pyinseq requires python 3.5 or greater"
    sys.exit()

__version__ = open(os.path.join('pyinseq', 'VERSION')).read().strip()

SCRIPTS = glob('scripts/*')

def main():
    setup(  name = 'pyinseq',
            version = __version__,
            description = 'Analysis of transposon insertion sequencing (INSeq) data in Python',
            url = 'https://github.com/mandel01/pyinseq',
            author = 'Mark J. Mandel',
            author_email = 'mandel01@gmail.com',
            license = 'BSD',
            packages = find_packages(),
            scripts = SCRIPTS,
            setup_requires = ['pytest-runner'],
            tests_require = ['pytest'],
            install_requires = ['matplotlib>=1.5.0',
                                'seaborn>=0.6.0'
                                'numpy>=1.10.0',
                                'pandas>=0.18.1',
                                'pytest>=2.8.1',
                                'PyYAML>=3.11',
                                'regex>=2016.6.5',
                                'screed>=0.9'],
            zip_safe = False,
            include_package_data = True )


if __name__ == "__main__":
    main()
