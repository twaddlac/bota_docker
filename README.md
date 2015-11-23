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
    --loci=STRING       Loci selection, either "human" or "mouse" (default:
                        human).
    --allele=STRING     Specify the allele you are interested in; if multiple,
                        separate them using coma.
                        For a full list of allele, type "python BOTA.py
                        --list_allele.
    --list_allele       Print full list of allele and leave.
    
    --min_pep_length=INT
                        The minimum length of epitode to consider (default:
                        12).
    --prodigal=DIR      The directory to prodigal binary, specify if not in
                        ENV (http://prodigal.ornl.gov/).
    --blat=DIR          The directory to blat binary, specify if not in ENV
                        (https://genome.ucsc.edu/FAQ/FAQblat.html).
    --psort=DIR         The directory to PSort, specify if not in ENV
                        (http://www.psort.org/).
    --hmmtop=DIR        The directory to HMMTOP, specify if not in ENV
                        (http://www.enzim.hu/hmmtop/).
    --netMHCII=STRING   The directory to netMHCII binary, specify if not in
                        ENV http://www.cbs.dtu.dk/services/NetMHCII).
                        
By running:

    python BOTA.py --list_allele
    
script will print out all possible alleles to select from as below:

    HLA-DRB10101
    HLA-DRB10301
    HLA-DRB10401
    HLA-DRB10404
    HLA-DRB10405
    HLA-DRB10701
    HLA-DRB10802
    HLA-DRB10901
    HLA-DRB11101
    HLA-DRB11302
    HLA-DRB11501
    HLA-DRB30101
    HLA-DRB40101
    HLA-DRB50101
    HLA-DPA10103-DPB10401
    HLA-DPA10103-DPB10201
    HLA-DPA10201-DPB10101
    HLA-DPA10201-DPB10501
    HLA-DPA10103
    HLA-DPB10301_DPB10401
    HLA-DPA10301-DPB10402
    HLA-DQA10101-DQB10501
    HLA-DQA10102-DQB10602
    HLA-DQA10301-DQB10302
    HLA-DQA10401-DQB10402
    HLA-DQA10501-DQB10201
    HLA-DQA10501-DQB10301
    H-2-IAb
    H-2-IAd
    
For instance, if you have a genome that is the EGD-e strain of Listeria monocytogenes in fastA format: "L.monocytogenes_EGD-e.fa", and you want to find out the epitodes that could be presented by mouse "H-2-IAb", you could run:

    python pipeline/BOTA.py -i test/L.monocytogenes_EGD-e.fa -o test/Listeria.BOTA/ -t 2 --allele H-2-IAb
    
Alternatively, if you are interested in all alleles from mouse MHC-II then you could run:

    python pipeline/BOTA.py -i test/L.monocytogenes_EGD-e.fa -o test/Listeria.BOTA/ -t 2 --loci mouse

The output of the eiptopes predicted would be in the file: test/Listeria.BOTA/L.monocytogenes_EGD-e.epitopes.out

Interpret output
===========

The output from the previous example looks like this:

    #epitope	core	gene	allele	start	end	affinity	cell_location
    AGAYSGAHLNPAVTI	YSGAHLNPA	ENA|AL591824|AL591824.1_1170|1194120-1194824|1	H-2-IAb	55	70	70.1000	CytoplasmicMembrane
    KVTFAGTVVQPIVEA	FAGTVVQPI	ENA|AL591824|AL591824.1_258|282755-284227|1	H-2-IAb	265	280	84.1000	Cellwall
    DVHYTWSIPNSTNVK	YTWSIPNST	ENA|AL591824|AL591824.1_156|157089-159470|1	H-2-IAb	65	80	22.3250	Cellwall
    DRVFTAGAPISSTNF	FTAGAPISS	ENA|AL591824|AL591824.1_156|157089-159470|1	H-2-IAb	243	258	49.1500	Cellwall
    TYFFTSPYARAVAVA	FTSPYARAV	ENA|AL591824|AL591824.1_2810|2879906-2880205|-1	H-2-IAb	40	55	21.1500	CytoplasmicMembrane
    MKPYNPTTALVVAED	YNPTTALVV	ENA|AL591824|AL591824.1_1671|1717193-1722328|-1	H-2-IAb	751	766	114.9250	Cellwall
    ......
    

Tutorial
=====================

Below is an example to walk you through what a real BOTA analysis should look like.

[TBD]

