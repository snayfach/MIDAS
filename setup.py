try:
	from setuptools import setup
except:
	from distutils.core import setup
import os

setup(
	name = 'MIDAS',
	version = '1.3.0',
	description = 'An integrated pipeline and for estimating species and strain-level genomic variation from metagenomic data',
	license = 'GPL',
	author = 'Stephen Nayfach',
	author_email='snayfach@gmail.com',
	url='https://github.com/snayfach/MIDAS',
	install_requires = ['biopython >= 1.62', 'numpy >= 1.7.0', 'pysam >= 0.8.1', 'pandas >= 0.17.1']
)
