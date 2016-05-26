#!/usr/bin/env python

# MIDAS: Metagenomic Intra-species Diversity Analysis System
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import os, sys, subprocess, shutil

def download(url, progress=True):
	print("Downloading: %s" % url)
	command = "curl --retry 20 --speed-time 60 --speed-limit 10000 %s -O" % url
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
	download(os.path.join(url_base, file), progress=True)
	decompress(file, file.split('.tar.gz')[0])

	# reference database
	file = "midas_db_v1.0.tar.gz"
	download(os.path.join(url_base, file), progress=True)
	decompress(file, file.split('.tar.gz')[0])
	print("Renaming: %s to ref_db" % (file.split('.tar.gz')[0]))
	shutil.move(file.split('.tar.gz')[0], "ref_db")




