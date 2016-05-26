#!/usr/bin/env python

import unittest
import shutil
import os
import subprocess

def run(command):
	""" run shell command & return unix exit code """
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
	return(process.returncode)

class ImportDependencies(unittest.TestCase):
	""" test that all dependencies can be imported """
	def setUp(self):
		self.failures = []
		try: import numpy
		except Exception: self.failures.append('numpy')
		try: import pandas
		except Exception: self.failures.append('pandas')
		try: import pysam
		except Exception: self.failures.append('pysam')
		try: import midas
		except Exception: self.failures.append('midas')
		try: import Bio.SeqIO
		except Exception: self.failures.append('Bio.SeqIO')
		
	def test_dependencies(self):
		self.assertTrue(len(self.failures)==0,
		msg="""\n\nThe following dependencies failed to import: %s.\nMake sure that dependencies have been properly installed""" % str(self.failures))

class HelpText(unittest.TestCase):
	""" check help text for all scripts """
	def setUp(self):
		self.scripts = [
			'../scripts/run_midas.py -h',
			'../scripts/run_midas.py species -h',
			'../scripts/run_midas.py genes -h',
			'../scripts/run_midas.py snps -h',
			'../scripts/merge_midas.py -h',
			'../scripts/merge_midas.py species -h',
			'../scripts/merge_midas.py genes -h',
			'../scripts/merge_midas.py snps -h']
		self.exit_codes = [run(_) for _ in self.scripts]
	
	def test_help_text(self):
		self.assertTrue(all([_ == 0 for _ in self.exit_codes]))

class RunSpecies(unittest.TestCase):
	""" test run_midas.py species """
	def setUp(self):
		self.command = '../scripts/run_midas.py species ./sample -1 ./test.fq.gz -n 100'
	def test_help_text(self):
		error = "Failed to execute the command: %s " % self.command
		self.assertTrue(run(self.command)==0, msg=error)

class RunGenes(unittest.TestCase):
	""" test run_midas.py genes """
	def setUp(self):
		self.command = '../scripts/run_midas.py genes ./sample -1 ./test.fq.gz -n 100 --species_id 57955'
	def test_help_text(self):
		error = "Failed to execute the command: %s " % self.command
		self.assertTrue(run(self.command)==0, msg=error)

class RunSNPs(unittest.TestCase):
	""" test run_midas.py snps """
	def setUp(self):
		self.command = '../scripts/run_midas.py snps ./sample -1 ./test.fq.gz -n 100 --species_id 57955'
	def test_help_text(self):
		error = "Failed to execute the command: %s " % self.command
		self.assertTrue(run(self.command)==0, msg=error)

class MergeSpecies(unittest.TestCase):
	""" test merge_midas.py species """
	def setUp(self):
		self.retcodes = []
		self.retcodes.append(run('../scripts/run_midas.py species ./sample -1 ./test.fq.gz -n 100'))
		self.retcodes.append(run('../scripts/merge_midas.py species ./species -i ./sample -t list'))
	def test_help_text(self):
		error = "Failed to execute the command: merge_midas.py species "
		self.assertTrue(sum(self.retcodes)==0, msg=error)

class MergeGenes(unittest.TestCase):
	""" test merge_midas.py species """
	def setUp(self):
		self.retcodes = []
		self.retcodes.append(run('../scripts/run_midas.py genes ./sample -1 ./test.fq.gz -n 100 --species_id 57955'))
		self.retcodes.append(run('../scripts/merge_midas.py genes ./genes -i ./sample -t list --species_id 57955 --sample_depth 0.0'))
	def test_help_text(self):
		error = "Failed to execute the command: merge_midas.py genes "
		self.assertTrue(sum(self.retcodes)==0, msg=error)

class MergeSNPs(unittest.TestCase):
	""" test merge_midas.py species """
	def setUp(self):
		self.retcodes = []
		self.retcodes.append(run('../scripts/run_midas.py snps ./sample -1 ./test.fq.gz -n 100 --species_id 57955'))
		self.retcodes.append(run('../scripts/merge_midas.py snps ./snps -i ./sample -t list --species_id 57955 --sample_depth 0.0 --max_sites 100 --fract_cov 0.0'))
	def test_help_text(self):
		error = "Failed to execute the command: merge_midas.py snps "
		self.assertTrue(sum(self.retcodes)==0, msg=error)

if __name__ == '__main__':
	unittest.main()
	shutil.rmtree('test')