#!/usr/bin/python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

__version__ = '0.0.2'

import os
import sys
import subprocess
import platform

def download(url, outpath, progress=True):
	print("Downloading: %s" % url)
	command = "curl %s > %s" % (url, outpath)
	subprocess.call(command, shell=True)

def decompress(tar, file, remove=True):
	print("Decompressing: %s" % tar)
	c = 'tar -zxvf %s %s' % (tar, file)
	p = subprocess.Popen(c, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = p.communicate()
	os.remove(tar)

if __name__ == '__main__':
	
	url_base = "http://lighthouse.ucsf.edu/phylocnv"
	script_dir = os.path.dirname(os.path.abspath(__file__))
	main_dir = os.path.dirname(script_dir)
	os.chdir(main_dir)
	
	# examples
	file = "example.tar.gz"
	download('%s/%s' % (url_base, file), file, progress=True)
	decompress("example.tar.gz", "example")
	
	# reference database
	refdb_dir = '%s/ref_db' % main_dir
	if not os.path.isdir(refdb_dir): os.mkdir(refdb_dir)
	os.chdir(refdb_dir)
	files = ["README.txt", "annotations.txt", "membership.txt", "marker_genes.tar.gz", "genome_clusters.tar.gz", "ontologies.tar.gz"]
	for file in files:
		download('%s/%s' % (url_base, file), file, progress=True)
	decompress("marker_genes.tar.gz", "marker_genes")
	decompress("genome_clusters.tar.gz", "genome_clusters")
	decompress("ontologies.tar.gz", "ontologies")

