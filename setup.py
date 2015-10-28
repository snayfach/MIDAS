from distutils.core import setup
from distutils.command import install_lib
from distutils import log
import os

class my_install_lib(install_lib.install_lib):
	""" Make all packaged binaries executable """
	def run(self):
		install_lib.install_lib.run(self)
		for fn in self.get_outputs():
			if fn.split('/')[-2] == 'bin':
				mode = ((os.stat(fn).st_mode) | 0555) & 07777
				log.info("changing mode of %s to %o", fn, mode)
				os.chmod(fn, mode)

setup(
	name = 'PhyloCNV',
	version = '0.0.2',
	description = 'An integrated pipeline and for estimating species and strain-level genomic variation from metagenomic data',
	packages = ['phylo_cnv'],
	package_data={'phylo_cnv': ['data/*', 'bin/Darwin/*', 'bin/Linux/*', 'example/*']},
	#scripts=['scripts/run_phylo_cnv.py', 'scripts/merge_snps.py', 'scripts/merge_genes.py'],
	license = 'GPL',
	author = 'Stephen Nayfach',
	author_email='snayfach@gmail.com',
	url='https://github.com/snayfach/PhyloCNV',
	install_requires = ['biopython', 'numpy', 'pysam', 'MicrobeCensus'],
	cmdclass={'install_lib':my_install_lib}
)
