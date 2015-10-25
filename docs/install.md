### Download the software
Clone the latest version from GitHub: `git clone https://github.com/snayfach/PhyloCNV`
This makes it easy to update the software as bugs are fixed and features are added

Alternatively, download the latest version of the software: https://github.com/snayfach/PhyloCNV/archive/master.zip
And unpack the project: `unzip master.zip`

### Install the software
Run setup.py to install:  
`python PhyloCNV/setup.py install` or  
`sudo python PhyloCNV/setup.py install` to install as a superuser or
`python PhyloCNV/setup.py install --user` to install locally (~/.local)  

Alternatively, you can manually install the software.
First, check that required python libraries are installed. You should be able to enter the following command in the python interpreter without getting an error:  
`>>> import Bio.SeqIO`  
`>>> import numpy`  
`>>> import pysam`  
`>>> import microbe_census`

Next, add the following to your PYTHONPATH environmental variable:  
`export PYTHONPATH=$PYTHONPATH:/path/to/PhyloCNV` or  
`echo -e "\nexport PYTHONPATH=\$PYTHONPATH:/path/to/PhyloCNV" >> ~/.bash_profile` to avoid entering the command in the future

Finally, add the scripts directory to your PATH environmental variable:  
`export PATH=$PATH:/path/to/PhyloCNV/scripts` or  
`echo -e "\nexport PATH=\$PATH:/path/to/PhyloCNV/scripts" >> ~/.bash_profile` to avoid entering the command in the future

### Testing
You should be able to enter the command into your terminal without getting an error:  
`run_phylo_cnv.py -h`

## Next step
[Download the reference database] (https://github.com/snayfach/PhyloCNV/blob/master/docs/ref_db.md)
