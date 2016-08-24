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

After this step, you will need to download the latest Pfam-A database, just follow below url:

ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release

and download "Pfam-A.hmm.gz".

You will then gunzip it:

    gunzip Pfam-A.hmm.gz

Then one last thing, you will need to prepare it by running hmmpress:

    hmmpress Pfam-A.hmm
    
With the right dependencies installed and you are ready to go!


Dependencies
============

* Python-2.7 or above

+ Python libraries:

>BioPython

>NetworkX

+ Third party pipelines: 

>standalone psortb v3.0+ (http://www.psort.org/)

>standalone HMMTOP (http://www.enzim.hu/hmmtop)

>standalone HMMER v3.1b2 (http://www.hmmer.org)

You don't have to install the package, just call BOTA.py from wherever you put the whole folder. 


Input
===========

The required input files are the .fa file of the genome and the .gff file (GFF3 format) of the annotations.

You can run multiple genomes in one project, you will just need to specify each of them in the config file.

Usage
===========

The basic BOTA analysis runs as below:

    python BOTA.py [options] -i <input_fasta_sequences> -o <output_directory>
    
Below is a detailed usage of BOTA.py:

    Usage: BOTA.py -c/--config <config_file> -o/--outdir <output directory> [options]

    BOTA: Bacteria-Origin T-cell Antigen predictor


The configuration file format follows:

    # this a comment line
    hmmscan='hmmscan' # if alread in the PATH, then just leave it as 'hmmscan', if not, specify the path
    hmmtop='hmmtop' # the same as hmmscan
    psort='psort'	# the same as hmmscan
    [genome_name] # you need the squared bracket to contain whatever you want to call the genome
    fna	/path/to/genomics/fna/file  # this is a compulsory field
    gff	/path/to/gff/file  # optional. if not supplied, we will do prodigal protein-coding gene predictions
    hmmtop /path/to/hmmtop/file # optional. if not supplied, we will do hmmtop calculation.
    hmmscan /path/to/text-hmmscan/file # optional
    psort /path/to/psort/file # optional
    alleles list_of_alleles_separated_by_comma # Optional. you can also supply human or mouse to select all available alleles
              # if you don't specify, default to all alleles.
	gram # Optional. specify the organism is 'P', gram-positive, 'N', gram-negative, or 'A', achaea; if not specified, BOTA
			will try to determine it.

Add --help to see a full list of required and optional arguments to run BOTA

Additional information can also be found at:
 
    https://bitbucket.org/luo-chengwei/bota/wiki

Options:
 
   --version             show program's version number and exit

    -h, --help            show this help message and exit

Required options:

    These options are required to run BOTA, and may be supplied in any order.

    -c FILE, --config=FILE   The configuration file to define a project.
    -o DIR, --outdir=DIR  The output directory of BOTA. If it doesn't exist, BOTA will create it.

Optional parameters:

    These options are optional, and may be supplied in any order.

    -t INT, --nproc=INT    Number of processor for BOTA to use [default: 1; set 0 if you want to use all CPUs available].
  
Below is a list of HLA alleles that BOTA will support:

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

    BOTA.py -c Listeria_monocytogenes_EGD-e/config -o Listeria_monocytogenes_EGD-e/BOTA/ -t 2
    

The output of the eiptopes predicted would be in the file: test/Listeria.BOTA/L.monocytogenes_EGD-e.epitopes.out

Interpret output
===========

The output from the previous example looks like this:

    #peptide	gene_name	chr_acc	gene_start	gene_stop	strand	pep_start	pep_stop	score
    FSSATLNSA	TA05_25490	JXUN01000256.1	6793	9210	-	1111	1140	0.695987
    ELGALSLSA	TA05_25490	JXUN01000256.1	6793	9210	-	1195	1224	0.758599
    WPAGGLASA	TA05_19735	JXUN01000201.1	15072	17129	+	1399	1428	0.697193
    ISLALAAPSYAAEA	TA05_13960	JXUN01000151.1	17693	19804	+	43	87	0.758729
    FSVAAAMES	TA05_09645	JXUN01000116.1	12753	14144	-	889	918	0.826708
    FSVAYASQA	TA05_03005	JXUN01000059.1	24619	25944	-	412	441	0.957796
    AQGVVTAPAQNSTVAVA	nlpD	JXUN01000196.1	6459	7586	-	544	597	0.918355
    VTAPVTAPAVSTT	nlpD	JXUN01000196.1	6459	7586	-	682	723	0.961591
    ASKPTITYS	nlpD	JXUN01000196.1	6459	7586	-	592	621	0.646375
    TTEPTASST	nlpD	JXUN01000196.1	6459	7586	-	715	744	0.860312
    ......
    

Tutorial
=====================

Below is an example to walk you through what a real BOTA analysis should look like.

[TBD]

Citation
=====================

**Graham DB<sup>\+</sup>**, **Luo C<sup>\+</sup>**, Abelin JG, Matar CG, Conway KL, Lefkovith A, Jasso GJ, Causer K, Mani DR, Carr SA, and Xavier RJ, Antigen discovery
by MHCII peptidomics reveals biochemical features of immunodominance. [In Review]

