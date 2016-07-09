# Alternative installations

## CentOS 7
If the user has a brand new install of CentOS 7, the following instructions will work to install the prerequisites of MIDAS:

```
sudo yum install git gcc-c++ epel-release gcc-gfortran libgfortran -y
sudo yum install python-pip -y
sudo pip install --upgrade pip
sudo yum install python-devel python-nose numpy -y
sudo yum install blas-devel lapack-devel atlas-devel zlib-devel libcurl-devel -y
```
Then , assuming you have downloaded MIDAS already,
```
cd /path/to/MIDAS
sudo python setup.py install 
```

These instructions have been tested on a CentOS 7 VM in the Google Cloud Engine (GCE).

# #Ubuntu 16.04LTS
If the user has a brand new install of Ubuntu 16.04, the following instructions will work to install the pre-requisites:
```
sudo apt-get install gcc cpp python-pip python-nose python-setuptools python-numpy python-pandas python-pysam python-biopython
sudo apt-get install zlib1g-dev
pip install --upgrade pip # for some reason pip is not at the bleeding edge version ...
pip install --upgrade pysam # ubuntu does not have the right version of pysam 
```

These instructions have been tested on a Ubuntu 16.04LTS VM in the Google Cloud
