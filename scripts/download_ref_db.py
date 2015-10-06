import os
import urllib
import subprocess
import time

script_dir = os.path.dirname(os.path.abspath(__file__))
main_dir = os.path.dirname(script_dir)
refdb_dir = '%s/ref_db' % main_dir
if not os.path.isdir(refdb_dir): os.mkdir(refdb_dir)
os.chdir(refdb_dir)

urllib.urlretrieve ("http://lighthouse.ucsf.edu/phylocnv/README.txt", "README.txt")
urllib.urlretrieve ("http://lighthouse.ucsf.edu/phylocnv/annotations.txt", "annotations.txt")
urllib.urlretrieve ("http://lighthouse.ucsf.edu/phylocnv/membership.txt", "membership.txt")

print("\nDownloading marker genes database")
start = time.time()
urllib.urlretrieve ("http://lighthouse.ucsf.edu/phylocnv/marker_genes.tar.gz", "marker_genes.tar.gz")
elapsed = time.time() - start
print( '%s seconds') % round(elapsed,2)

print("\nDecompressing marker genes database")
start = time.time()
p = subprocess.Popen('tar -zxvf marker_genes.tar.gz marker_genes', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = p.communicate()
os.remove('marker_genes.tar.gz')
elapsed = time.time() - start
print( '%s seconds') % round(elapsed,2)

print("\nDownloading genomes database")
start = time.time()
urllib.urlretrieve ("http://lighthouse.ucsf.edu/phylocnv/genome_clusters.tar.gz", "genome_clusters.tar.gz")
elapsed = time.time() - start
print( '%s seconds') % round(elapsed,2)

print("\nDecompressing genomes database")
start = time.time()
p = subprocess.Popen('tar -zxvf genome_clusters.tar.gz genome_clusters', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = p.communicate()
os.remove('genome_clusters.tar.gz')
elapsed = time.time() - start
print( '%s seconds') % round(elapsed,2)
