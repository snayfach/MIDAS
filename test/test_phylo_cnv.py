#!/usr/bin/python

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
			'../scripts/annotate_genes.py -h',
			'../scripts/annotate_snps.py -h',
			'../scripts/merge_genes.py -h',
			'../scripts/merge_species.py -h',
			'../scripts/merge_snps.py -h',
			'../scripts/pairwise_distances.py -h',
			'../scripts/run_phylo_cnv.py -h',
			'../scripts/run_phylo_cnv.py species -h',
			'../scripts/run_phylo_cnv.py genes -h',
			'../scripts/run_phylo_cnv.py snvs -h']
		self.exit_codes = [check_exit_code(_) for _ in self.scripts]
	def test_help_text(self):
		self.assertTrue(all([_ == 0 for _ in self.exit_codes]))

class RunSpecies(unittest.TestCase):
	""" test run_phylo_cnv.py species """
	def setUp(self):
		self.command = '../scripts/run_phylo_cnv.py species -1 ./test.fq.gz -o ./species.txt -n 100'
	def test_help_text(self):
		self.assertTrue(check_exit_code(self.command)==0)
	def tearDown(self):
		os.remove('species.txt')

class RunCNVs(unittest.TestCase):
	""" test run_phylo_cnv.py genes """
	def setUp(self):
		self.command = '../scripts/run_phylo_cnv.py genes -1 ./test.fq.gz -o ./test -n 100 --gc_id 57955'
	def test_help_text(self):
		self.assertTrue(check_exit_code(self.command)==0)
	def tearDown(self):
		shutil.rmtree('test')

class RunSNVs(unittest.TestCase):
	""" test run_phylo_cnv.py snvs """
	def setUp(self):
		self.command = '../scripts/run_phylo_cnv.py snvs -1 ./test.fq.gz -o ./test -n 100 --gc_id 57955'
	def test_help_text(self):
		self.assertTrue(check_exit_code(self.command)==0)
	def tearDown(self):
		shutil.rmtree('test')

#class FileType(unittest.TestCase):
#	""" check whether filetype is correctly determined """
#
#	def setUp(self): # create test files in tmp directory
#		self.dir = tempfile.mkdtemp()
#		self.values = [['file.fq', 'fastq',
#		               ('@HWUSI-EAS574_102539073:1:100:10000:12882/1',
#						'AGCTCTTCCAGCGATACAATACCATCGTTCCTTCGGTAGCATC',
#						'+HWUSI-EAS574_102539073:1:100:10000:12882/1',
#						'GGGGFGFFGGAGDFG=EDEEDEBEEEEEEEEEEEAB?B?BEEE')],
#					   ['file.fa', 'fasta',
#		               ('>HWUSI-EAS574_102539073:1:100:10000:12882/1',
#						'AGCTCTTCCAGCGATACAATACCATCGTTCCTTCGGTAGCATC')],
#					   ['file.txt', None,
#		               ('some random text file')]]
#		for name, type, values in self.values:
#			inpath = os.path.join(self.dir, name)
#			infile = open(inpath, 'w')
#			for value in values: infile.write(value+'\n')
#
#	def test_detect_filetype(self):
#		# fastq
#		inpath = os.path.join(self.dir, self.values[0][0])
#		type = microbe_census.auto_detect_file_type(inpath)
#		self.assertEqual(type, self.values[0][1])
#		# fasta
#		inpath = os.path.join(self.dir, self.values[1][0])
#		type = microbe_census.auto_detect_file_type(inpath)
#		self.assertEqual(type, self.values[1][1])
#		# neither
#		inpath = os.path.join(self.dir, self.values[2][0])
#		with self.assertRaises(SystemExit):
#			microbe_census.auto_detect_file_type(inpath)
#
#	def tearDown(self): # clean up tmp directory
#		shutil.rmtree(self.dir)

if __name__ == '__main__':
	unittest.main()