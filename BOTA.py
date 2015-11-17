#!/usr/bin/env python
# -*- coding: utf-8 -*- 

__author__ = 'Chengwei Luo (cluo@broadinstitute.org)'
__version__ = '0.1.0'
__date__ = 'Oct 2015'

"""
BOTA: Bacteria-Origin T-cell Antigen predictor

Copyright(c) 2015 Chengwei Luo (luo.chengwei@gatech.edu)

	This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>

https://bitbucket.org/luo-chengwei/BOTA

This script is part of the metaHGT package,

for help for BOTA, type:
python BOTA.py --help

"""

USAGE = \
"""Usage: %prog <required_parameters> [options]

BOTA: Bacteria-Origin T-cell Antigen predictor

Add --help to see a full list of required and optional
arguments to run BOTA

Additional information can also be found at:
https://bitbucket.org/luo-chengwei/bota/wiki

If you use BOTA in your work, please cite it as:
<BOTA citation TBD>

Copyright: Chengwei Luo, Broad Institute of MIT and Harvard, 2015
"""


import sys, os, re, glob, shutil, tempfile
from optparse import OptionParser, OptionGroup
from operator import itemgetter
from time import ctime, time
import multiprocessing as mp
from subprocess import call, PIPE, Popen

from Bio import SeqIO
from Bio.Seq import Seq

sys.path.append(os.path.basename(sys.argv[0])+'/algorithms/')

def which(program):
	def is_exe(fpath):
		return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

	fpath, fname = os.path.split(program)
	if fpath:
		if is_exe(program):
			return os.path.abspath(program)
	else:
		for path in os.environ["PATH"].split(os.pathsep):
			path = path.strip('"')
			exe_file = os.path.join(path, program)
			if is_exe(exe_file):
				return os.path.abspath(exe_file)
	return None

def run_prodigal(args):
	prodigal, sample_name, infile, genome_outdir, mode = args
	tmp_faa = '%s/%s.faa.tmp' % (genome_outdir, sample_name)
	tmp_ffn = '%s/%s.ffn.tmp' % (genome_outdir, sample_name)
	gbk = '%s/%s.gbk' % (genome_outdir, sample_name)
	cmds = [prodigal, '-a', tmp_faa, '-d', tmp_ffn, \
			'-f', 'gbk', '-i', infile, '-o', gbk, '-p', mode]
	call(cmds, stderr = PIPE, stdout = PIPE)
	# modify ffn and faa files
	faa = open('%s/%s.faa' % (genome_outdir, sample_name), 'w')
	ffn = open('%s/%s.ffn' % (genome_outdir, sample_name), 'w')
	for record in SeqIO.parse(tmp_faa, 'fasta'):
		x = record.description.split(' # ')
		tag = '%s|%s-%s|%s' % (x[0], x[1], x[2], x[3])
		faa.write('>%s\n%s\n' % (tag, record.seq))
	faa.close()
	for record in SeqIO.parse(tmp_ffn, 'fasta'):
		x = record.description.split(' # ')
		tag = '%s|%s-%s|%s' % (x[0], x[1], x[2], x[3])
		ffn.write('>%s\n%s\n' % (tag, record.seq))
	ffn.close()
	# remove the tmp files
	os.unlink(tmp_faa)
	os.unlink(tmp_ffn)
	sys.stdout.write('  [%s] Peptide prediction finished.\n' % sample_name)
	
	return 0

def	run_blat(args):
	blat, db, infile, outfile = args
	call([blat, '-prot', db, infile, '-out=blast8', outfile], stderr = PIPE, stdout = PIPE)
	return 0
	
def call_Gram(blat_file):
	thrs = {'k': 30, 'p': 40, 'c':50, 'o': 60}
	best_matches = {}
	for line in open(blat_file, 'rb'):
		cols = line.split('\t')
		gene, hit, identity, evalue, score = cols[0], cols[1], float(cols[2]), float(cols[-2]), float(cols[-1])
		if evalue > 1e-6: continue
		if gene not in best_matches: best_matches[gene] = [hit, identity, score]
		elif score > best_matches[gene][-1]: best_matches[gene] = [hit, identity, score]
	
	N = len(best_matches)
	
	X = 'N'
	origins = {}
	for x in best_matches:
		origin = best_matches[x][0].split('|')[0]
		if origin not in origins: origins[origin] = []
		origins[origin].append(best_matches[x][1])
	
	for x in origins:
		xn = len(origins[x])
		x_avg = sum(origins[x])/len(origins[x])
		tax_level = x.split('__')[0]
		if x_avg < thrs[tax_level]: continue
		if xn < 0.5*N: continue
		if x == 'k__Archaea': X = 'A'
		else: X = 'P'
		
	return X

def run_psort(args):
	psort, subdir, infile, outfile, gram = args
	try: os.mkdir(subdir)
	except: return 1
	# split infile into multiple files with 10 entries
	file_ind = 0
	ofh = None
	for i, record in enumerate(SeqIO.parse(infile, 'fasta')):
		if i % 10 == 0:
			if ofh != None: ofh.close()
			ofh = open(subdir+'/%i.faa' % (1+i/10), 'wb')
		ofh.write('>%s\n%s\n' % (record.name, record.seq))
	ofh.close()
	
	xfh = open(outfile, 'wb')
	for sub_infile in glob.glob(subdir+'/*.faa'):
		sub_outfile = sub_infile.replace('.faa', '.psort')
		if gram == 'P':
			cmd = [psort, '-p', '-o', 'terse', sub_infile]
		elif gram == 'A':
			cmd = [psort, '-a', '-o', 'terse', sub_infile]
		else:
			cmd = [psort, '-n', '-o', 'terse', sub_infile]
		with open(sub_outfile, 'wb') as ofh:
			call(cmd, stdout = ofh, stderr = ofh)
	for sub_outfile in glob.glob(subdir+'/*.psort'):
		for line_ind, line in enumerate(open(sub_outfile, 'rb')):
			if line_ind == 0: continue
			xfh.write(line)
	xfh.close()
		
	return 0
	
def main(argv = sys.argv[1:]):
	parser = OptionParser(usage = USAGE, version="Version: " + __version__)
	
	# Required arguments
	requiredOptions = OptionGroup(parser, "Required options",
		"These options are required to run BOTA, and may be supplied in any order.")
	
	requiredOptions.add_option("-i", "--infile", type = "string", metavar = "FILE",
		help = "Input microbial genome file(s) in fasta format (if multiple, separate by coma).")
		
	requiredOptions.add_option("-o", "--outdir", type = "string", metavar = "DIR",
		help = "The output directory of BOTA. If it doesn't exist, BOTA will create it.")
							
	parser.add_option_group(requiredOptions)

	# Optional arguments that need to be supplied if not the same as default
	optOptions = OptionGroup(parser, "Optional parameters",
		"These options are optional, and may be supplied in any order.")

	optOptions.add_option("-t", "--num_proc", type = "int", default = 1, metavar = 'INT',
		help = "Number of processor for BOTA to use [default: 1; set 0 if you want to use all CPUs available].")
	
	optOptions.add_option("-m", "--mode", type = "string", default = 'single', metavar = 'STRING',
		help = "Mode of running BOTA, either \"single\" or \"meta\" (default: single).")
	
	optOptions.add_option("--prodigal", type = "string", metavar = "DIR", default = 'prodigal',
		help = "The directory to prodigal binary, specify if not in ENV (http://prodigal.ornl.gov/).")
	
	optOptions.add_option("--blat", type = "string", metavar = "DIR", default = 'blat',
		help = "The directory to blat binary, specify if not in ENV (https://genome.ucsc.edu/FAQ/FAQblat.html).")
	
	optOptions.add_option("-p", "--psort", type = "string", metavar = "DIR", default = 'psort',
		help = "The directory to PSort, specify if not in ENV (http://www.psort.org/).")
	
	optOptions.add_option("-n", "--netMHCII", type = "string", metavar = "STRING", default = 'netMHCII',
		help = "The directory to netMHCII binary, specify if not in ENV http://www.cbs.dtu.dk/services/NetMHCII).")
	
							
	parser.add_option_group(optOptions)
	
	(options, args) = parser.parse_args(argv)
	
	if options.infile == None:
		parser.error('[FATAL]: No input fasta files specified.\n')
		exit(1)
	try: 
		infiles_st = options.infile.split(',')
	except: 
		sys.stderr.write('[FATAL]: Error in parsing input file string.\n')
		exit(1)
	
	infiles = []
	for infile_st in infiles_st:
		if not os.path.exists(infile_st):
			sys.stderr.write('[FATAL]: Error in locating %s.\n' % infile_st)
			exit(1)
		sample_name = '.'.join(os.path.basename(infile_st).split('.')[:-1])
		infiles.append([sample_name, os.path.abspath(infile_st)])
	
	if options.outdir == None:
		parser.error('[FATAL]: No output directory specified.\n')
		exit(1)
	if not os.path.exists(options.outdir):
		try: 
			os.mkdir(options.outdir)
		except:
			sys.stderr.write('[FATAL]: Error in creating output directory: %s.\n' % options.outdir)
	outdir = os.path.abspath(options.outdir)
	
	if options.mode not in ['single', 'meta']:
		sys.stderr.write('[FATAL]: the mode has to be either \"single\" or \"meta\", you supplied: %s\n' % options.mode)
		exit(1)
	
	nproc = 1
	if options.num_proc < 0:
		sys.stderr.write('[WARNING]: Cannot set number of CPUs to be %i as you supplied, will use 1 CPU instead.\n' % options.num_proc)
	elif options.num_proc > mp.cpu_count():
		sys.stderr.write('[WARNING]: Cannot set number of CPUs to be %i, which exceeds the max number of CPUs available (%i), \
			will use %i CPUs instead.\n' (options.num_proc, mp.cpu_count(), mp.cpu_count()))
		nproc = mp.cpu_count()
	else: nproc = options.num_proc
	
	# test 3rd party exe
	if not which(options.prodigal):
		sys.stderr.write('[FATAL]: Cannot locate Prodigal in ENV nor in the directory you supplied: %s.\n' % options.prodigal)
		exit(1)
	else: options.prodigal = which(options.prodigal)
	
	if not which(options.blat):
		sys.stderr.write('[FATAL]: Cannot locate blat in ENV nor in the directory you supplied: %s.\n' % options.blat)
		exit(1)
	else: options.blat = which(options.blat)
	
	if not which(options.psort):
		sys.stderr.write('[FATAL]: Cannot locate psort in ENV nor in the directory you supplied: %s.\n' % options.psort)
		exit(1)
	else: options.psort = which(options.psort)
	
	# test db files
	current_dir = os.path.dirname(os.path.realpath(__file__))
	db_dir = os.path.join(current_dir, 'db/')
	db_files = [os.path.join(db_dir, f) for f in ['Gram.faa',]]
	for db_file in db_files:
		if not os.path.exists(db_file):
			sys.stderr.write('[FATAL]: Cannot locate DB file %s\n' % os.path.basename(db_file))
			exit(1)
	
	# predict protein coding genes first
	sys.stdout.write('Now predicting peptides first in genomes...\n')
	gene_prediction_dir = os.path.join(outdir, 'genes')
	if not os.path.exists(gene_prediction_dir): os.mkdir(gene_prediction_dir)
	prodigal_cmds = []
	for sample_name, infile in infiles:
		prodigal_cmds.append([options.prodigal, sample_name, infile, gene_prediction_dir, options.mode])
	
	pool = mp.Pool(nproc)
	pool.map_async(run_prodigal, prodigal_cmds)
	pool.close()
	pool.join()
	
	sys.stdout.write('Peptide prediction done!\n\n')
	
	# profile genomes' taxonomy using AMPHORA, predict its taxonomy, Gram-positive/negative
	sys.stdout.write('Now predicting genomes Gram positive or negative...\n')
	gram_dir = os.path.join(outdir, 'gram/')
	if not os.path.exists(gram_dir): os.mkdir(gram_dir)
	blat_cmds = []
	for sample_name, infile in infiles:
		blat_file = os.path.join(gram_dir, '%s.blat' % sample_name)
		faa_file = os.path.join(gene_prediction_dir, '%s.faa' % sample_name)
		if not os.path.exists(blat_file):
			blat_cmds.append([options.blat, db_files[0], faa_file, blat_file])
	if len(blat_cmds) > 0:
		pool = mp.Pool(nproc)
		pool.map_async(run_blat, blat_cmds)
		pool.close()
		pool.join()
	
	Grams = {}
	sys.stdout.write('    # sample          Gram-positive/negative/Archaea [P/N/A]\n')
	for sample_name, infile in infiles:
		blat_file = os.path.join(gram_dir, '%s.blat' % sample_name)
		gram = call_Gram(blat_file)
		Grams[sample_name] = gram
		sys.stdout.write('    %s\t%s\n' % (sample_name, gram))
	sys.stdout.write('Done!\n\n')
	
	# filter the peptides down to candidate list
	sys.stdout.write('Now predicting locations of peptides...\n')
	select_candidate_peptide_cmds = []
	psort_cmds = []
	psort_dir = os.path.join(outdir, 'psort/')
	if not os.path.exists(psort_dir): os.mkdir(psort_dir)
	for sample_name, infile in infiles:
		genome_file = infile
		faa_file = gene_prediction_dir+'/'+sample_name + '.faa'
		psort_output = psort_dir+'/'+sample_name+'.psort'
		if not os.path.exists(psort_output) or os.stat(psort_output).st_size == 0:
			subdir = os.path.join(psort_dir, sample_name)
			psort_cmds.append([options.psort, subdir, faa_file, psort_output, Grams[sample_name]])
	
	if len(psort_cmds) > 0:
		pool = mp.Pool(nproc)
		pool.map_async(run_psort, psort_cmds)
		pool.close()
		pool.join()
	sys.stdout.write('Done!\n\n')	
	
		
	#	candidate_list_file = outdir+'/'+sample_name + '.candidate_peptides.txt'
	#	select_candidate_peptide([faa_file, infile, candidate_list_file])
		#select_candidate_peptide_cmds.append([faa_file, infile, candidate_list_file])

if __name__ == '__main__': main()
