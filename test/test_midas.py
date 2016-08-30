#!/usr/bin/env python

import unittest
import shutil
import os
import subprocess
from distutils.version import StrictVersion

def run(command):
	""" run shell command & return unix exit code """
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
	return(process.returncode)

class CheckEnv(unittest.TestCase):
	""" check environmental variables """
	def setUp(self):
		self.path_contents = []
		for _ in os.environ['PATH'].strip(':').split(':'):
			if os.path.isdir(_): self.path_contents += os.listdir(_)
		self.python_contents = []
		for _ in os.environ['PYTHONPATH'].strip(':').split(':'):
			if os.path.isdir(_): self.python_contents += os.listdir(_)
		
	def test_dependencies(self):
		self.assertTrue(
			'run_midas.py' in self.path_contents,
			msg="""\n\n'run_midas.py' not found in PATH environmental variable.\nMake sure '/path/to/MIDAS/scripts' has been added to your PATH"""
			)
		self.assertTrue(
			'midas' in self.python_contents,
			msg="""\n\n'midas' not found in PYTHONPATH environmental variable.\nMake sure '/path/to/MIDAS' has been added to your PYTHONPATH"""
			)

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

class CheckVersions(unittest.TestCase):
	""" check version numbers for dependencies """
	def setUp(self):
		self.modules = ['numpy', 'pandas', 'pysam', 'Bio.SeqIO']
		self.installeds = [module.__version__ for module in map(__import__, self.modules)]
		self.requireds = ['1.7.0', '0.17.1', '0.8.1', '1.6.2']
		
	def test_dependencies(self):
		for module, installed, required in zip(self.modules, self.installeds, self.requireds):
			if len(installed.split('.')) > 3:
				installed = '.'.join(installed.split('.')[0:3])
			self.assertTrue(
				StrictVersion(installed) >= StrictVersion(required),
			msg="""\n\nImported library '%s %s' is out of date. Required version is >= %s""" % (module, installed, required) )

class HelpText(unittest.TestCase):
	""" check help text for all scripts """
	def setUp(self):
		self.scripts = [
			'run_midas.py -h',
			'run_midas.py species -h',
			'run_midas.py genes -h',
			'run_midas.py snps -h',
			'merge_midas.py -h',
			'merge_midas.py species -h',
			'merge_midas.py genes -h',
			'merge_midas.py snps -h']
		self.exit_codes = [run(_) for _ in self.scripts]
	def test_help_text(self):
		error = "\n\nFailed to run MIDAS.\nMake sure you've appended /path/to/MIDAS/scripts to your PATH and /path/to/MIDAS to your PYTHONPATH variables"
		self.assertTrue(all([_ == 0 for _ in self.exit_codes]), msg=error)

class RunSpecies(unittest.TestCase):
	""" test run_midas.py species """
	def setUp(self):
		self.command = 'run_midas.py species ./sample -1 ./test.fq.gz -n 100'
	def test_help_text(self):
		error = "\n\nFailed to execute the command: %s " % self.command
		self.assertTrue(run(self.command)==0, msg=error)

class RunGenes(unittest.TestCase):
	""" test run_midas.py genes """
	def setUp(self):
		self.command = 'run_midas.py genes ./sample -1 ./test.fq.gz -n 100 --species_id 57955'
	def test_help_text(self):
		error = "\n\nFailed to execute the command: %s " % self.command
		self.assertTrue(run(self.command)==0, msg=error)

class RunSNPs(unittest.TestCase):
	""" test run_midas.py snps """
	def setUp(self):
		self.command = 'run_midas.py snps ./sample -1 ./test.fq.gz -n 100 --species_id 57955'
	def test_help_text(self):
		error = "\n\nFailed to execute the command: %s " % self.command
		self.assertTrue(run(self.command)==0, msg=error)

class MergeSpecies(unittest.TestCase):
	""" test merge_midas.py species """
	def setUp(self):
		self.retcodes = []
		self.retcodes.append(run('run_midas.py species ./sample -1 ./test.fq.gz -n 100'))
		self.retcodes.append(run('merge_midas.py species ./species -i ./sample -t list'))
	def test_help_text(self):
		error = "\n\nFailed to execute the command: merge_midas.py species "
		self.assertTrue(sum(self.retcodes)==0, msg=error)

class MergeGenes(unittest.TestCase):
	""" test merge_midas.py species """
	def setUp(self):
		self.retcodes = []
		self.retcodes.append(run('run_midas.py genes ./sample -1 ./test.fq.gz -n 100 --species_id 57955'))
		self.retcodes.append(run('merge_midas.py genes ./genes -i ./sample -t list --species_id 57955 --sample_depth 0.0'))
	def test_help_text(self):
		error = "\n\nFailed to execute the command: merge_midas.py genes "
		self.assertTrue(sum(self.retcodes)==0, msg=error)

class MergeSNPs(unittest.TestCase):
	""" test merge_midas.py species """
	def setUp(self):
		self.retcodes = []
		self.retcodes.append(run('run_midas.py snps ./sample -1 ./test.fq.gz -n 100 --species_id 57955'))
		self.retcodes.append(run('merge_midas.py snps ./snps -i ./sample -t list --species_id 57955 --sample_depth 0.0 --max_sites 100 --fract_cov 0.0'))
	def test_help_text(self):
		error = "\n\nFailed to execute the command: merge_midas.py snps "
		self.assertTrue(sum(self.retcodes)==0, msg=error)

if __name__ == '__main__':
	unittest.main()
	shutil.rmtree('test')