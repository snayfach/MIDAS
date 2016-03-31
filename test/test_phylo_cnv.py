#!/usr/bin/env python

import unittest
import tempfile
import shutil
import os
import subprocess

def check_exit_code(command):
	""" run shell command & return unix exit code """
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
	return(process.returncode)

class ImportDependencies(unittest.TestCase):
	""" test that all dependencies can be imported """
	def setUp(self):
		try:
			import numpy
			import pysam
			import microbe_census
			import Bio.SeqIO
			import phylo_cnv
			self.success = True
		except Exception:
			self.success = False
	def test_dependencies(self):
		self.assertTrue(self.success,
		msg="""\n\nOne or more dependencies failed to import.\nMake sure that dependencies have been properly installed""")

class HelpText(unittest.TestCase):
	""" check help text for all scripts """
	def setUp(self):
		self.scripts = [
			'../scripts/merge_genes.py -h',
			'../scripts/merge_species.py -h',
			'../scripts/merge_snps.py -h',
			'../scripts/pairwise_distances.py -h',
			'../scripts/run_phylo_cnv.py -h',
			'../scripts/run_phylo_cnv.py species -h',
			'../scripts/run_phylo_cnv.py genes -h',
			'../scripts/run_phylo_cnv.py snps -h']
		self.exit_codes = [check_exit_code(_) for _ in self.scripts]
	def test_help_text(self):
		self.assertTrue(all([_ == 0 for _ in self.exit_codes]))

class RunSpecies(unittest.TestCase):
	""" test run_phylo_cnv.py species """
	def setUp(self):
		self.command = '../scripts/run_phylo_cnv.py species ./test -1 ./test.fq.gz -n 100 --remove_temp'
	def test_help_text(self):
		self.assertTrue(check_exit_code(self.command)==0)
	def tearDown(self):
		shutil.rmtree('test')

class RunCNVs(unittest.TestCase):
	""" test run_phylo_cnv.py genes """
	def setUp(self):
		self.command = '../scripts/run_phylo_cnv.py genes ./test -1 ./test.fq.gz -n 100 --sp_id 57955 --remove_temp'
	def test_help_text(self):
		self.assertTrue(check_exit_code(self.command)==0)
	def tearDown(self):
		shutil.rmtree('test')

class RunSNPs(unittest.TestCase):
	""" test run_phylo_cnv.py snps """
	def setUp(self):
		self.command = '../scripts/run_phylo_cnv.py snps ./test -1 ./test.fq.gz -n 100 --sp_id 57955 --remove_temp'
	def test_help_text(self):
		self.assertTrue(check_exit_code(self.command)==0)
	def tearDown(self):
		shutil.rmtree('test')

if __name__ == '__main__':
	unittest.main()