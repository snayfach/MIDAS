### Requirements

Requires Unix or Mac OSX operating system

Python dependencies:   

* Numpy (>=1.7.0)
* BioPython (>=1.6.2)
* Pysam (>=0.8.1)
* Pandas (>=0.17.1)

Dependencies are installed via setup.py  
Tested version numbers are indicated in parenthesis. Other versions may also work.

### Download the software


**Clone the latest version from GitHub (recommended)**   
`git clone https://github.com/snayfach/PhyloCNV`  
This makes it easy to update the software in the future using `git pull` as bugs are fixed and features are added

**Or, download the source code**: 
https://github.com/snayfach/PhyloCNV/archive/master.zip  
And unpack the project: `unzip PhyloCNV-master.zip`

### Install the software

**Run setup.py to install python dependencies**  
`python MIDAS/setup.py install` or  
`sudo python MIDAS/setup.py install` to install as a superuser or  
`python MIDAS/setup.py install --user` to install locally (~/.local)
[Read more] (https://github.com/snayfach/PhyloCNV/blob/master/docs/ref_db.md)

**Update your environmental variables**  
`export PYTHONPATH=$PYTHONPATH:/path/to/MIDAS`  
`export PATH=$PATH:/path/to/MIDAS/scripts` 

Be sure to replace '/path/to' with the appropriate path on your system  
Add these commands to your .bash_profile to avoid entering the commands in the future

### Testing
You should be able to enter the command into your terminal without getting an error:  
`run_midas.py -h`

For more complete testing, run:   
`python tests/test_midas.py`

### Update MIDAS
Move to installation directory, pull the latest version, and install:  
`cd MIDAS`  
`git pull`  
`python setup.py install`


## Next step
[Download the reference database] (https://github.com/snayfach/MIDAS/blob/master/docs/ref_db.md)
