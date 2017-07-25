#!/usr/bin/env python

import unittest
import shutil
import os
import subprocess
import sys
from distutils.version import StrictVersion

def run(command):
	""" run shell command & return unix exit code """
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
	return(err, process.returncode)

class _01_CheckEnv(unittest.TestCase):
	def setUp(self):
		self.path_contents = []
		for _ in os.environ['PATH'].strip(':').split(':'):
			if os.path.isdir(_): self.path_contents += os.listdir(_)
		self.python_contents = []
		for _ in os.environ['PYTHONPATH'].strip(':').split(':'):
			if os.path.isdir(_): self.python_contents += os.listdir(_)
		
	def test_class(self):
		self.assertTrue(
			'run_midas.py' in self.path_contents,
			msg="""\n\n'run_midas.py' not found in PATH environmental variable.\nMake sure '/path/to/MIDAS/scripts' has been added to your PATH:\nexport PATH=$PATH:/path/to/MIDAS/scripts"""
			)
		self.assertTrue(
			'midas' in self.python_contents,
			msg="""\n\n'midas' not found in PYTHONPATH environmental variable.\nMake sure '/path/to/MIDAS' has been added to your PYTHONPATH:\nexport PYTHONPATH=$PYTHONPATH:/path/to/MIDAS"""
			)
		self.assertTrue(
			'MIDAS_DB' in os.environ,
			msg="""\n\n'MIDAS_DB' environmental variable not set.\nSet this variable and rerun the test:\nexport MIDAS_DB=/path/to/midas_db_v1.1"""
			)

class _02_ImportDependencies(unittest.TestCase):
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
		
	def test_class(self):
		self.assertTrue(len(self.failures)==0,
		msg="""\n\nThe following dependencies failed to import: %s.\nMake sure that dependencies have been properly installed""" % str(self.failures))

class _03_CheckVersions(unittest.TestCase):
	def setUp(self):
		self.modules = ['numpy', 'pandas', 'pysam', 'Bio.SeqIO']
		self.installeds = [module.__version__ for module in map(__import__, self.modules)]
		self.requireds = ['1.7.0', '0.17.1', '0.8.1', '1.6.2']
		
	def test_class(self):
		for module, installed, required in zip(self.modules, self.installeds, self.requireds):
			if len(installed.split('.')) > 3:
				installed = '.'.join(installed.split('.')[0:3])
			self.assertTrue(
				StrictVersion(installed) >= StrictVersion(required),
			msg="""\n\nImported library '%s %s' is out of date. Required version is >= %s""" % (module, installed, required) )

class _04_HelpText(unittest.TestCase):
	def test_class(self):
		commands = [
			'run_midas.py -h',
			'run_midas.py species -h',
			'run_midas.py genes -h',
			'run_midas.py snps -h',
			'merge_midas.py -h',
			'merge_midas.py species -h',
			'merge_midas.py genes -h',
			'merge_midas.py snps -h']
		for cmd in commands:
			err, code = run(cmd)
			self.assertTrue(code==0, msg=err)

class _05_RunSpecies(unittest.TestCase):
	def test_class(self):
		command = 'run_midas.py species ./sample -1 ./test.fq.gz -n 100'
		err, code = run(command)
		self.assertTrue(code==0, msg=err)

class _06_RunGenes(unittest.TestCase):
	def test_class(self):
		command = 'run_midas.py genes ./sample -1 ./test.fq.gz -n 100 --species_id Bacteroides_vulgatus_57955'
		err, code = run(command)
		self.assertTrue(code==0, msg=err)
		
class _07_RunSNPs(unittest.TestCase):
	def test_class(self):
		command = 'run_midas.py snps ./sample -1 ./test.fq.gz -n 100 --species_id Bacteroides_vulgatus_57955'
		err, code = run(command)
		self.assertTrue(code==0, msg=err)
		
class _08_MergeSpecies(unittest.TestCase):
	def test_class(self):
		command = 'merge_midas.py species ./species -i ./sample -t list'
		err, code = run(command)
		self.assertTrue(code==0, msg=err)

class _09_MergeGenes(unittest.TestCase):
	def test_class(self):
		command = 'merge_midas.py genes ./genes -i ./sample -t list --species_id Bacteroides_vulgatus_57955 --sample_depth 0.0'
		err, code = run(command)
		self.assertTrue(code==0, msg=err)

class _10_MergeSNPs(unittest.TestCase):
	def test_class(self):
		command = 'merge_midas.py snps ./snps -i ./sample -t list --species_id Bacteroides_vulgatus_57955 --all_samples --all_sites --max_sites 10000'
		err, code = run(command)
		self.assertTrue(code==0, msg=err)

class _11_SNPdiversity(unittest.TestCase):
	def test_class(self):
		command = 'snp_diversity.py snps/Bacteroides_vulgatus_57955'
		err, code = run(command)
		self.assertTrue(code==0, msg=err)

class _12_CallConsensus(unittest.TestCase):
	def test_class(self):
		command = 'call_consensus.py snps/Bacteroides_vulgatus_57955'
		err, code = run(command)
		self.assertTrue(code==0, msg=err)		
	
class _13_CompareGeneContent(unittest.TestCase):
	def test_class(self):
		command = 'compare_genes.py genes/Bacteroides_vulgatus_57955'
		err, code = run(command)
		self.assertTrue(code==0, msg=err)

class _14_QueryByCompound(unittest.TestCase):
	def test_class(self):
		command = 'query_by_compound.py -i sample -t list -c C00312'
		err, code = run(command)
		self.assertTrue(code==0, msg=err)
		
class _15_BuildDB(unittest.TestCase):
	def test_class(self):
		command = 'tar -zxvf genomes.tar.gz'
		err, code = run(command)
		command = 'build_midas_db.py genomes genomes.mapfile db --threads 10'
		err, code = run(command)
		self.assertTrue(code==0, msg=err)

if __name__ == '__main__':
	try:
		dir_name = os.path.dirname(os.path.abspath(__file__))
		os.chdir(dir_name)
		unittest.main(exit=False)
		for dir in ['sample', 'species', 'genes', 'snps', 'genomes', 'db']:
			shutil.rmtree(dir)
	except:
		print("")
		for dir in ['sample', 'species', 'genes', 'snps', 'genomes', 'db']:
			if os.path.exists(dir): shutil.rmtree(dir)
