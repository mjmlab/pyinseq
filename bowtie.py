#!/usr/bin/env python
"""
Build bowtie index for specified genome

Specify exact format of input .fna and .ptt files here

What about RNA genes?

"""

import subprocess

def bowtie_build():
    subprocess.check_call(['/Users/markmandel/bowtie-1.1.1/bowtie-build'])


# ===== Start here ===== #

def main():
    bowtie_build()

if __name__ == "__main__":
    main()
