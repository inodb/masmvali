#!/usr/bin/env python
from setuptools import setup, find_packages
import sys, os

version = '0.1'

setup(name='masmvali',
      version=version,
      description="Metagenomic ASseMbly VALIdation",
      long_description="""MASMVALI (Metagenomic ASseMbly VALIdation) is an 
      assembly validator for metagenomics with known reference genomes.""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='Python Scilifelab, Metagenomic, Assembly, Validation',
      author='Ino de Bruijn, Johannes Alneberg', 
      author_email='ino.debruijn@scilifelab.se',
      url='https://github.com/inodb/masmvali',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=False,
      install_requires=['cython>=0.19.1',
                        'numpy>=1.7.1',
                        'biopython>=1.62b',
                        'pysam>=0.6',
                        'sh>=1.09',
                        'decorator>=3.4.0'],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )

