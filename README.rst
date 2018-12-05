====
fnfi
====

fnfi is a collection of tools for detection of structural variants, originally designed for identifying telomere fusions.


Installation
------------
Install using::

    $ python setup.py install
    # Or
    $ pip install .

Requires python>=2.7, c++11 compatible compiler. For macosx, minimum version is 10.9. Python packages needed: click,
numpy, pandas, pysam, pybedtools, natsort, networkx, scikit-learn, ncls.

Recommended aligners: bwa mem, lastal

Usage
-----
Sub commands are::

    $ fnfi run
    $ fnfi find-reads
    $ fnfi align
    $ fnfi call-events

Basic usage command is::

    $ fufi run reference.fa your.bam > results.csv

For help use::

    $ fufi --help
    $ fufi find-reads|run|align|call-events --help

For more extensive documentation refer to the manual.