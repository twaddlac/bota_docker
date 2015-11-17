BOTA
===========

BOTA: Predicting Bacteria origined T-cell antigens

One line pitcher
===========
BOTA inspects all possible peptides in a genomic or metagenomic sequence and predicts its possibility of being presented on host T-cells.

Install
===========
use git:

    git clone https://bitbucket.org/luo-chengwei/BOTA

use hg:

    hg clone https://bitbucket.org/luo-chengwei/BOTA

You can also download the zip archive of it from bitbucket repository: 

https://bitbucket.org/luo-chengwei/BOTA

Dependencies
============

* Python-2.7 or above

+ Python libraries:

>BioPython

+ Third party pipelines: 

>standalone psortb v3.0+ (http://www.psort.org/)

>standalone HMMTOP (http://www.enzim.hu/hmmtop)

>standalone NetMHCII 2.0+ (http://www.cbs.dtu.dk/services/NetMHCII)

You don't have to install the package, just call BOTA.py from wherever you put the whole folder. 

Usage
===========

The basic BOTA analysis runs as below:

    python BOTA.py [options] -i <input_fasta_sequences> -o <output_directory>
    
Below is a detailed usage of BOTA.py:

  Options:

    --version             show program's version number and exit

    -h, --help            show this help message and exit

  
Required options:

  These options are required to run BOTA, and may be supplied in any order.

    -i FILE, --infile=FILE
                        Input microbial genome file(s) in fasta format (if
                        multiple, separate by coma).
    -o DIR, --outdir=DIR
                        The output directory of BOTA. If it doesn't exist,
                        BOTA will create it.

  
Optional parameters:

These options are optional, and may be supplied in any order.

    -t INT, --num_proc=INT
                        Number of processor for BOTA to use [default: 1; set 0
                        if you want to use all CPUs available].
    -m STRING, --mode=STRING
                        Mode of running BOTA, either "single" or "meta"
                        (default: single).
    --prodigal=DIR      The directory to prodigal binary, specify if not in
                        ENV (http://prodigal.ornl.gov/).
    --blat=DIR          The directory to blat binary, specify if not in ENV
                        (https://genome.ucsc.edu/FAQ/FAQblat.html).
    -p DIR, --psort=DIR
                        The directory to PSort, specify if not in ENV
                        (http://www.psort.org/).
    -n STRING, --netMHCIIpan=STRING
                        The directory to netMHCIIpan binary, specify if not in
                        ENV (http://www.cbs.dtu.dk/services/NetMHCII).
                        
Interpret output
===========

In the project/output directory you will find a folder called "results", in which you can find the below files and directories:

Tutorial
=====================

Below is an example to walk you through what a real BOTA analysis should look like.

