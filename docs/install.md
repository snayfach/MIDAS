### Download the software

MIDAS is a command-line tool that runs on Unix or Mac OSX

**Clone the latest version from GitHub (recommended)**   
`git clone https://github.com/snayfach/MIDAS`  
This makes it easy to update the software in the future using `git pull` as bugs are fixed and features are added

**Or, download the source code**: 
https://github.com/snayfach/MIDAS/archive/master.zip  
And unpack the project: `unzip MIDAS-master.zip`

### Install dependencies

The following Python modules should be pre-installed on your system
(tested version numbers are indicated in parenthesis):

* Numpy (>=1.7.0)
* BioPython (>=1.6.2)
* Pysam (>=0.8.1)
* Pandas (>=0.17.1)

The dependencies can be installed by running the setup.py script:  
`python MIDAS/setup.py install` or  
`sudo python MIDAS/setup.py install` to install as a superuser, or  
`python MIDAS/setup.py install --user` to install locally  

If you have pip installed, you can install these dependencies by running the following command:

`pip install -U numpy biopython pysam pandas` or  
`sudo pip install -U numpy biopython pysam pandas` to install as a superuser, or  
`pip install --user numpy biopython pysam pandas` to install locally  

### Dependencies: Alternate methods
[More detailed instructions for installing the dependencies on brand new CentOS and Ubuntu installs] (https://github.com/snayfach/MIDAS/blob/master/docs/install_.md)

### Update your environmental variables

This will enable you to run MIDAS scripts:  
`export PYTHONPATH=$PYTHONPATH:/path/to/MIDAS`  
`export PATH=$PATH:/path/to/MIDAS/scripts` 

Be sure to replace '/path/to' with the appropriate path on your system  
Add these commands to your .bash_profile to avoid entering the commands in the future

### Testing
You should be able to enter these commands into your terminal without getting an error:  
`run_midas.py -h`

For more complete testing, run:   
`cd /path/to/MIDAS`
`python test/test_midas.py -f`

Be sure to replace '/path/to' with the appropriate path on your system.
### Update MIDAS
Move to installation directory, pull the latest version, and install:  
`cd MIDAS`  
`git pull`  
`python setup.py install`


## Next step
[Download the reference database] (https://github.com/snayfach/MIDAS/blob/master/docs/ref_db.md)
