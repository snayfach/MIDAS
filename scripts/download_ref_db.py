#!/usr/bin/python

# PhyloCNV - estimation of single-nucleotide-variants and gene-copy-number from shotgun sequence data
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3

__version__ = '0.0.2'

import os
import sys
import urllib2
import subprocess

def download(url, outpath, progress=True):
	openurl = urllib2.urlopen(url)
	outfile = open(outpath, 'wb')
	file_size = int(openurl.info().getheaders("Content-Length")[0])
	print("Downloading: {0} Bytes: {1}".format(url, file_size))
	file_size_dl = 0
	block_sz = 8192
	last_percent = 0.0
	while True:
		# read
		buffer = openurl.read(block_sz)
		if not buffer:
			break
		# write
		outfile.write(buffer)
		# print progress
		file_size_dl += len(buffer)
		p = float(file_size_dl)/file_size
		if progress and (p - last_percent) > 0.0005:
			status = r"{0}  [{1:.2%}]".format(file_size_dl, p)
			status = status + chr(8)*(len(status)+1)
			sys.stdout.write(status)
			last_percent = p
	outfile.close()

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



