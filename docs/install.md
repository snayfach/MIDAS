### Download the software

Options for downloading the software:

**Clone the latest version from GitHub (recommended)**   
`git clone https://github.com/snayfach/PhyloCNV`  
This makes it easy to update the software in the future using `git pull` as bugs are fixed and features are added

**Download the source code**: 
https://github.com/snayfach/PhyloCNV/archive/master.zip  
And unpack the project: `unzip PhyloCNV-master.zip`

### Install the software

Options for installing the software:

**Run setup.py (recommended)**  
`python PhyloCNV/setup.py install` or  
`sudo python PhyloCNV/setup.py install` to install as a superuser or  
`python PhyloCNV/setup.py install --user` to install locally (~/.local)  

**Manual install**  

Alternative to setup.py. First, check that required python libraries are installed. You should be able to enter the following command in the python interpreter without getting an error:  
`>>> import Bio.SeqIO`  
`>>> import numpy`  
`>>> import pysam`  
`>>> import microbe_census`

Update your PYTHONPATH environmental variable (replace '/path/to' with the appropriate path on your system):  
`export PYTHONPATH=$PYTHONPATH:/path/to/PhyloCNV` 

Add this to your .bash_profile to avoid entering the command in the future (you will need to open a new shell for this to take effect):  
`echo -e "\nexport PYTHONPATH=\$PYTHONPATH:/path/to/PhyloCNV" >> ~/.bash_profile` 

### Update your PATH (optional)

All of PhyloCNV's scripts are located in the directory PhyloCNV/scripts.  Optionally, add this directory to your PATH (replace '/path/to' with the appropriate path on your system):  
`export PATH=$PATH:/path/to/PhyloCNV/scripts` 

Optionally, add this to your .bash_profile to avoid entering the command in the future (you will need to open a new shell for this to take effect):  
`echo -e "\nexport PATH=\$PATH:/path/to/PhyloCNV/scripts" >> ~/.bash_profile`  

### Testing
You should be able to enter the command into your terminal without getting an error:  
`run_phylo_cnv.py -h`

### Update PhyloCNV
Move to installation directory, pull the latest version, and install:  
`cd PhyloCNV`  
`git pull`  
`python setup.py install`


## Next step
[Download the reference database] (https://github.com/snayfach/PhyloCNV/blob/master/docs/ref_db.md)
