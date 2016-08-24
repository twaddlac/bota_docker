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
"""Usage: %prog -c/--config <config_file> -o/--outdir <output directory> [options]

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
	
Add --help to see a full list of required and optional
arguments to run BOTA

Additional information can also be found at:
https://bitbucket.org/luo-chengwei/bota/wiki

If you use BOTA in your work, please cite it as:
<BOTA citation TBD>

Copyright: Chengwei Luo, Broad Institute of MIT and Harvard, 2015
"""


import sys, os, re, glob, shutil, cPickle, random
from optparse import OptionParser, OptionGroup
from operator import itemgetter
from time import ctime, time
import multiprocessing as mp
from subprocess import call, PIPE, Popen
from itertools import groupby

import networkx
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

# KERAS init
from keras.models import Sequential, model_from_json
#from keras.layers.core import Dense, Dropout, Activation
#from keras.optimizers import SGD, Adam, RMSprop
#from keras.utils import np_utils


sys.path.append(os.path.basename(sys.argv[0])+'/algorithms/')

class Peptide:
	def __init__(self):
		self.seq = None
		self.gene = None
		self.start = -1
		self.end = -1
		self.scores = []
	
	def printInfo(self):
		sys.stdout.write('%s\t%i-%i\n%s' % (self.gene, self.start, self.end, self.seq))

class GeneInfo:
	def __init__(self):
		self.name = None
		self.seq = None
		self.struct = None
		self.cellular_loc = None
		self.PW_score = []
		self.domains = []
	
	def get_PWscore(self, pos):
		if len(self.PW_score) == 0:
			return [None, None]
		s, n = self.PW_score[pos-1]
		return [s, n]
		
	def select_candidate(self):
		X = []
		for i in range(len(self.seq)): X.append(1)
		for ind, a in enumerate(list(self.hmmtop)):
			if a != 'O': X[ind] = 0
		domains = sorted(self.domains, key=lambda x: x[0])
		for s, e in domains:
			if e-s+1 < 30:
				for i in range(max(0, s-8), min(e+8, len(X)-1)): X[i] = 0
		if len(domains) >1:
			for (s1, e1), (s2, e2) in zip(domains[:-1], domains[1:]):
				if s2 - e1 < 20:
					s, e = sorted([s2, e1])
					for i in range(s, e):
						X[i-1] = 0
						
		
		candidates = []
		ind = 1
		for a in [list(j) for i, j in groupby(X)]:
			s = ind
			e = ind+len(a)-1
			scores = self.PW_score[s-1:e]
			ind+=len(a)
			if len(scores) == 0: continue
			sx, nx = sorted(scores, key=lambda x: x[0]+float(x[1])/1e13, reverse=True)[0]
			if a[0] == 1 and e-s > 8 and sx+nx/1e13 > 5e-11:
				scores = self.PW_score[s-1:e]
				P = Peptide()
				P.seq = self.seq[s-1:e]
				P.start = s
				P.end = e
				P.gene = self.name
				P.scores = self.PW_score[s-1:e]
				candidates.append(P)
		return candidates
			
		
class GenomeProj:
	"""The base class to hold a BOTA project's sub genome project information"""
	def __init__(self):
		self.name=None
		self.fna=None
		self.blat=None
		self.hmmscan = None
		self.hmmtop=None
		self.psort=None
		self.alleles=[]
		self.gram=None

class ProjInfo:
	"""The base class to hold a BOTA project information"""
	def __init__(self):
		""" base init method """
		self.prodigal = None
		self.psort = None
		self.hmmscan = None
		self.hmmtop = None
		self.hmmtop_arch = [None, None]
		self.pfam = None
		self.HLA_db = None
		self.gram_db = None
		self.pwm_db = None
		self.genomes = {}
		self.models={}
		
	def init_from_config(self, config):
		"""Initiate a BOTA project information for configuration file.
		
		Parameters
		----------
		config: the path to the configuration file with format detailed in usage.
		
		"""
		# init the db files
		script_path = os.path.dirname(os.path.realpath(__file__))
		db_path = os.path.join(script_path, 'db')
		if not os.path.exists(db_path):
			sys.stderr.write('[FATAL] Cannot locate the db directory.\n')
			exit(1)
		self.hmmtop_arch=[os.path.join(db_path, 'hmmtop.arch'), os.path.join(db_path, 'hmmtop.psv')]
		self.HLA_db=os.path.join(db_path, 'HLA.db')
		self.gram_db=os.path.join(db_path, 'Gram.faa')
		self.pwm_db = os.path.join(db_path, 'PWMatrix.pkl')
		
		if not os.path.exists(self.hmmtop_arch[0]):
			sys.stderr.write('[FATAL] Cannot locate HMMTOP arch file hmmtop.arch in db/.\n')
		if not os.path.exists(self.hmmtop_arch[1]):
			sys.stderr.write('[FATAL] Cannot locate HMMTOP parameters file hmmtop.psv in db/.\n')
		if not os.path.exists(self.HLA_db):
			sys.stderr.write('[FATAL] Cannot locate HLA allele DB in db/.\n')
		if not os.path.exists(self.gram_db):
			sys.stderr.write('[FATAL] Cannot locate the DB for Gram typing in db/.\n')
		if not os.path.exists(self.pwm_db):
			sys.stderr.write('[FATAL] Cannot locate the DB for PW matrix in db/.\n')
		
		# check pfam, if not there, then download it
		if self.pfam == None:
			for ext in ['h3f', 'h3i', 'h3m', 'h3p']:
				pfam_file = os.path.join(db_path, 'Pfam-A.hmm.%s' % ext)
				if not os.path.exists(pfam_file):
					sys.stderr.write('[FATAL] Cannot locate Pfam-A.hmm.%s\n' % ext)
					exit(1)
			self.pfam=os.path.join(db_path, 'Pfam-A.hmm')
		
		# get all alleles
		all_alleles = {'human':[], 'mouse':[]}
		for line in open(self.HLA_db, 'rb'):
			cols = line[:-1].split(',')
			if line[:3] == 'HLA': all_alleles['human'] += cols
			else: all_alleles['mouse'] += cols
		allele_set = set(all_alleles['human']+all_alleles['mouse'])
			
		genome_id = ''
		missing_execs = []
		for line in open(config, 'rb'):
			if line[0] == '#': continue # the comment line
			if re.search('blat|hmmscan|prodigal|hmmtop|psort|gram|alleles\=\'.+\'', line[:-1]) != None:
				try:
					k, v = re.search('(.+)\=\'(.+)\'', line.rstrip()).group(1, 2)
				except:
					sys.stderr.write('[FATAL] Error in parsing the config file genome file section.\n')
					exit(1)
				av = which(v)
				if av != None:
					if k == 'blat': self.blat = av
					elif k == 'psort': self.psort = av
					elif k == 'prodigal': self.prodigal = av
					elif k == 'hmmscan': self.hmmscan=av
					elif k == 'hmmtop': self.hmmtop=av
					else: raise ValueError
				else:
					sys.stderr.write('[FATAL] Error in parsing config. %s path %s does not seem to work.\n' % (k, v))
					#exit(1)	
			elif re.search('^\[.+\]', line.rstrip()) != None:
				genome_id = re.search('^\[(.+)\]', line.rstrip()).group(1)
				self.genomes[genome_id] = GenomeProj()
				self.genomes[genome_id].name=genome_id
			elif re.search('^gff|fna|psort|hmmtop|hmmscan|alleles|gram\t.+$', line.rstrip()) != None:
				k, v = line.rstrip().split('\t')
				if k != 'alleles' and not os.path.exists(v):
					sys.stderr.write('[FATAL] cannot locate %s file for genome=\'%s\', you supplied %s.\n' % (k, genome_id, v))
					exit(1)
				elif k == 'fna': self.genomes[genome_id].fna=v
				elif k == 'gff': self.genomes[genome_id].gff=v
				elif k == 'hmmscan': self.genomes[genome_id].hmmscan=v
				elif k == 'hmmtop': self.genomes[genome_id].hmmtop=v
				elif k == 'psort': self.genomes[genome_id].psort=v
				elif k == 'gram':
					if v in list('APV'): self.genomes[genome_id].gram=v
					else:
						sys.stderr.write('[FATAL] for genome %s Gram stain type you supplied: %s. \
							It has to be A, P, or N.\n' % (genome_id, v))
						exit(1)
				elif k == 'alleles':
					alleles = v.split(',')
					if len(alleles) == 1 and alleles[0] == 'human':
						self.genomes[genome_id].alleles=all_alleles['human']
					elif len(alleles) == 1 and alleles[0] == 'mouse':
						self.genomes[genome_id].alleles=all_alleles['mouse']
					elif len(alleles) == 1 and alleles[0] in allele_set:
						self.genomes[genome_id].alleles=alleles
					elif len(all_alleles) > 1:
						for allele in all_alleles:
							if allele not in allele_set:
								sys.stderr.write('[FATAL] You supplied %s allele, not in the allele list.\n' % allele)
								sys.stderr.write('This is a list of all alleles available:\n')
								sys.stderr.write('[human]\n')
								for allele in all_alleles['human']: sys.stderr.write(allele+'\n')
								sys.stderr.write('[mouse]\n')
								for allele in all_alleles['mouse']: sys.stderr.write(allele+'\n')
								sys.stderr.write('\n')
								exit(1)
							else: 
								self.genomes[genome_id].alleles.append(allele)
							
				else:
					sys.stderr.write('[FATAL] Illegal keyword %s found. \
					Specify files with \'gff\', \'fna\', \'hmmscan\', \'alleles\', \'gram\', or \'hmmtop\'.\n' % (k))
			elif line[:-1] == '': pass
			else:
				sys.stderr.write('[FATAL] Error in parsing config, no keyword found.\n')
				exit(1)
				
		# refresh alleles, add models
		model_path = os.path.join(script_path, 'models')
		if not os.path.exists(model_path):
			sys.stderr.write('[FATAL] Cannot locate the DNN model path.\n')
			exit(1)
		alleles_needed = set()
		for genome_id in self.genomes:
			if len(self.genomes[genome_id].alleles) == 0:
				self.genomes[genome_id].alleles = list(allele_set)
			for allele in self.genomes[genome_id].alleles:
				alleles_needed.add(allele)
		
		for allele in alleles_needed:
			model_arch = os.path.join(model_path, '%s.model_arch.json' % allele)
			model_weights = os.path.join(model_path, '%s.model_weights.h5' % allele)
			if not os.path.exists(model_arch):
				sys.stderr.write('[FATAL] Cannot locate the DL architecture for allele: %s\n' % (allele))
				exit()
			if not os.path.exists(model_weights):
				sys.stderr.write('[FATAL] Cannot locate the DL weigths for allele: %s\n' % (allele))
				exit()
			self.models[allele] = [model_arch, model_weights]
		
		return 0
	
	def print_projInfo(self, stream=sys.stdout):
		""" by default prints to sys.stdout; 
		can be redirected to other streams such as a file handle."""
		stream.write('################ Project Configurations ##############\n')
		stream.write('[3rd party programs]\n')
		stream.write('  prodigal: %s\n' % self.prodigal)
		stream.write('  hmmscan: %s\n' % self.hmmscan)
		stream.write('  hmmtop: %s\n' % self.hmmtop)
		stream.write('  psort: %s\n' % self.psort)
		stream.write('  \n')
		for genome_id in self.genomes:
			stream.write('[%s]\n' % genome_id)
			stream.write('  fna=%s\n' % self.genomes[genome_id].fna)
			stream.write('  gff=%s\n' % self.genomes[genome_id].gff)
			stream.write('  hmmtop=%s\n' % self.genomes[genome_id].hmmtop)
			stream.write('  hmmscan=%s\n' % self.genomes[genome_id].hmmscan)
			stream.write('  psort=%s\n' % self.genomes[genome_id].psort)
			stream.write('	Gram=%s\n' % self.genomes[genome_id].gram)
			stream.write('  [Alleles]\n')
			for allele in self.genomes[genome_id].alleles:
				stream.write('    %s\n' % allele)
			stream.write('\n')
		return 0

def which(program):
	"""Tests if a program is executable or exists
	
	Parameters
	----------
	program: the path to the executable
	
	Returns
	-------
	abspath: the absolute path of the program if it is executable and the system can
				locate it; otherwise None.
	"""
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

def parse_config(config_file):
	"""Parses the configuration file, returns the ProjInfo object
	
	Parameters
	----------
	config_file: the configuration file formatted as detailed in usage
	
	Returns
	-------
	proj_info: the initiated ProjInfo object.
	
	"""
	proj_info = ProjInfo()
	proj_info.init_from_config(config_file)
	return proj_info

def extract_gff_features(genome_proj, faa_file):
	gff = genome_proj.gff
	fna = genome_proj.fna
	seqs={}
	for rec in SeqIO.parse(fna, 'fasta'):
		seqs[rec.name]=str(rec.seq)
	
	gene_dict = {}
	cds_set = []
	for line in open(gff, 'rb'):
		if line[0]=='#': continue
		cols = line[:-1].split('\t')
		if cols[2] == 'CDS':
			try: parent_gene = re.search('Parent\=([^;]+)\;', cols[-1]).group(1)
			except: continue
			tag = '%s|%s-%s|%s' % (cols[0], cols[3], cols[4], cols[6])
			try:
				s, e = int(cols[3]), int(cols[4])
				nt=Seq(seqs[cols[0]][s-1:e])
				if cols[6] == '-': aa = nt.reverse_complement().translate(stop_symbol='')
				else: aa=nt.translate(stop_symbol='')
			except: contine
			cds_set.append([parent_gene, tag, aa])
		elif cols[2] == 'gene':
			try: gene_name = re.search('Name\=([^;]+)\;', cols[-1]).group(1)
			except: continue
			try: geneID = re.search('ID\=([^;]+)\;', cols[-1]).group(1)
			except: continue
			gene_dict[geneID] = gene_name
	
	faa = open(faa_file, 'wb')
	for parent_gene, tag, aa in cds_set:
		try: gene_name = gene_dict[parent_gene]
		except: continue
		xtag = '%s|gene=%s' % (tag, gene_name)
		faa.write('>%s\n%s\n' % (xtag, aa))
	faa.close()

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

def run_hmmscan(args):
	hmmscan, pfam_db, faa, hmmscan_out, nproc = args
	cmd = [hmmscan, '-o', hmmscan_out+'.temp', '--domtblout', hmmscan_out, '--domE', '0.0001', 
			'--cpu', str(nproc), pfam_db, faa]
	p1 = Popen(cmd, stdout=PIPE, stderr=PIPE)
	out, err = p1.communicate()
	os.unlink(hmmscan_out+'.temp')
	return 0

	
def call_Gram(args):
	blat, db, faa, outdir, genomeID = args
	blat_file = os.path.join(outdir, '%s.blat' % genomeID)
	
	blat_cmd = [blat, '-prot', db, faa, '-out=blast8', blat_file]
	p1 = Popen(blat_cmd, stdout=PIPE, stderr=PIPE)
	out, err = p1.communicate()
	
	thrs = {'k': 30, 'p': 40, 'c':50, 'o': 60}
	best_matches = {}
	for line in open(blat_file, 'rb'):
		cols = line.split('\t')
		gene, hit, identity, evalue, score = cols[0], cols[1], float(cols[2]), float(cols[-2]), float(cols[-1])
		if evalue > 1e-6: continue
		if gene not in best_matches: best_matches[gene] = [hit, identity, score]
		elif score > best_matches[gene][-1]: best_matches[gene] = [hit, identity, score]
	os.unlink(blat_file)
	
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

def run_indi_psort(args):
	psort, infile, outfile, gram = args
	if gram == 'P':
		cmd = [psort, '-p', '-o', 'terse', infile]
	elif gram == 'A':
		cmd = [psort, '-a', '-o', 'terse', infile]
	else:
		cmd = [psort, '-n', '-o', 'terse', infile]
	with open(outfile, 'wb') as ofh:
		p1=Popen(cmd, stdout=ofh, stderr=PIPE)
		p1.communicate()
	return 0
		

def run_psort(args):
	psort, genome_dir, faa_file, gram, outfile, nproc = args
	subdir = os.path.join(genome_dir, 'psort_sub')
	try: os.mkdir(subdir)
	except: return 1
	# split infile into multiple files with 10 entries
	file_ind = 0
	ofh = None
	for i, record in enumerate(SeqIO.parse(faa_file, 'fasta')):
		if i % 10 == 0:
			if ofh != None: ofh.close()
			ofh = open(subdir+'/%i.faa' % (1+i/10), 'wb')
		ofh.write('>%s\n%s\n' % (record.name, record.seq))
	ofh.close()
	
	cmds = []
	for sub_infile in glob.glob(subdir+'/*.faa'):
		sub_outfile = sub_infile.replace('.faa', '.psort')
		cmds.append([psort, sub_infile, sub_outfile, gram])
	
	
	pool=mp.Pool(nproc)
	pool.map_async(run_indi_psort, cmds)
	pool.close()
	pool.join()
	
	#for cmd in cmds: print cmd; run_indi_psort(cmd)
		
	xfh = open(outfile, 'wb')
	for sub_outfile in glob.glob(subdir+'/*.psort'):
		for line_ind, line in enumerate(open(sub_outfile, 'rb')):
			if line_ind == 0: continue
			xfh.write(line)
	xfh.close()
	shutil.rmtree(subdir)
	return 0

def run_indi_hmmtop(args):
	hmmtop, infile, outfile = args
	cmd = [hmmtop, '-if=%s' % infile, '-of=%s' % outfile, '-pl']
	p1=Popen(cmd, stdout = PIPE, stderr = PIPE)
	p1.communicate()
	return 0

def run_hmmtop(args):
	hmmtop, faa_file, outfile, genome_dir, nproc = args
	subdir = os.path.join(genome_dir, 'hmmtop_sub')
	try: os.mkdir(subdir)
	except: return 1
	
	# split infile into multiple files with 10 entries
	file_ind = 0
	ofh = None
	for i, record in enumerate(SeqIO.parse(faa_file, 'fasta')):
		if i % 10 == 0:
			if ofh != None: ofh.close()
			ofh = open(subdir+'/%i.faa' % (1+i/10), 'wb')
		ofh.write('>%s\n%s\n' % (record.name, record.seq))
	ofh.close()
	
	cmds = []
	for sub_infile in glob.glob(subdir+'/*.faa'):
		sub_outfile = sub_infile.replace('.faa', '.hmmtop')
		cmds.append([hmmtop, sub_infile, sub_outfile])
	
	pool=mp.Pool(nproc)
	pool.map_async(run_indi_hmmtop, cmds)
	pool.close()
	pool.join()
	
	#for cmd in cmds: run_indi_hmmtop(cmd)
	
	xfh = open(outfile, 'wb')
	for sub_outfile in glob.glob(subdir+'/*.hmmtop'):
		for line_ind, line in enumerate(open(sub_outfile, 'rb')): xfh.write(line)
	xfh.close()
	shutil.rmtree(subdir)
	return 0
	
def convert_hmmtop_output(infile):
	models = {}
	for line in open(infile, 'rb'):
		if line[0] == '>':
			gene, orientation, n_st = line[:-1].split()[2:5]
			models[gene] = {'seq': '', 'pred': ''}
			continue
		cols = [x for x in line[:-1].split() if x != '']
		if 'seq' not in cols and 'pred' not in cols: continue
		elif 'seq' in cols:
			for s in cols[1:-1]: models[gene]['seq']+=s
		elif 'pred' in cols:
			for s in cols[1:]: models[gene]['pred']+=s
	return models

def extract_structs(seq, pred):
	X = [[k,len(list(g))] for k, g in groupby(pred)]
	start_ind = 0
	end_ind = 0
	structs = {}
	for symbol, L in X:
		if symbol not in structs: structs[symbol] = []
		end_ind += L
		seq_seg = seq[start_ind: end_ind]
		structs[symbol].append([seq_seg, start_ind+1, end_ind])
		start_ind+=L
	return structs

def pw_scoring(pw_dict, score, pep):
	# cal absolute score first
	S = 1
	for x in score: S*=x
	# cal relative score
	pep_lst = list(pep)
	N = 0
	for ind in range(100):
		Sx = 1
		random.shuffle(pep_lst)
		for aa, pos_ind in zip(pep_lst, range(9)):
			try: Sx*=pw_dict[(aa, pos_ind)]
			except:
				aa = random.sample(list('ARNDCQEGHILKMFPSTWYV'), 1)[0]
				Sx*=pw_dict[(aa, pos_ind)]
		if Sx < S: N+=1
	return (S, N)
		

def run_pwmscore(args):
	faa_file, pw_dict, loc_dict, pwm_file, nproc = args
	# load all the sequences, and score them
	seq_scores = {}
	for rec_ind, rec in enumerate(SeqIO.parse(faa_file, 'fasta')):
		try: loc = loc_dict[rec.name]
		except: continue
		if loc not in ['Cellwall', 'Extracellular', 'OuterMembrane']: continue
		seq_scores[rec.name] = []
		seq = str(rec.seq).replace('*', '', 1000)
		for s_ind in range(len(seq)-8):
			peptide = seq[s_ind:s_ind+9]
			score = []
			for aa, pos_ind in zip(list(peptide), range(9)):
				try: score.append(pw_dict[(aa, pos_ind)])
				except:
					aa = random.sample(list('ARNDCQEGHILKMFPSTWYV'), 1)[0]
					score.append(pw_dict[(aa, pos_ind)])
			score_x = pw_scoring(pw_dict, score, peptide)
			seq_scores[rec.name].append(score_x)
	cPickle.dump(seq_scores, open(pwm_file, 'wb'))
	return 0
	 

def integrate_data(genomeID, genome_dir, genome_pkl):
	psort_file = os.path.join(genome_dir, '%s.psort' % genomeID)
	hmmtop_file = os.path.join(genome_dir, '%s.hmmtop' % genomeID)
	hmmscan_file = os.path.join(genome_dir, '%s.hmmscan' % genomeID)
	pwm_file = os.path.join(genome_dir, '%s.pwm' % genomeID)
	faa_file = os.path.join(genome_dir, '%s.faa' % genomeID)
	# make data
	gene_info_dict = {}
	# load ID and seq
	for rec in SeqIO.parse(faa_file, 'fasta'):
		g = GeneInfo()
		g.name = rec.name
		g.seq = str(rec.seq)
		gene_info_dict[rec.name] = g
	# load hmmtop
	models = convert_hmmtop_output(hmmtop_file)
	for geneID in models:
		pred_struct = models[geneID]['pred']
		gene_info_dict[geneID].hmmtop = pred_struct
	# load loc
	for line in open(psort_file, 'rb'):
		cols = line.rstrip().split('\t')
		geneID = cols[0].rstrip(' ')
		loc = cols[1]
		gene_info_dict[geneID].cellular_loc = loc
	# load domains
	for line in open(hmmscan_file, 'rb'):
		if line[0] == '#': continue
		cols = [x for x in line.rstrip().split(' ') if x != '']
		try:geneID = cols[3]
		except: print cols
		start, end = int(cols[19]), int(cols[20])
		gene_info_dict[geneID].domains.append([start, end])
	# load PW scores
	seq_scores = cPickle.load(open(pwm_file, 'rb'))
	for geneID in seq_scores:
		gene_info_dict[geneID].PW_score = seq_scores[geneID]
	
	cPickle.dump(gene_info_dict, open(genome_pkl, 'wb'))
	return 0
		
def transform_seq(seq, aa_dict):
	D = []
	for ind in range(0, len(seq)-9):
		subseq = seq[ind:ind+9]
		d = []
		for ind, aa in enumerate(list(subseq)):
			try:aa_ind = aa_dict[aa]
			except: aa_ind = random.randint(0, 19)
			d.append(aa_ind)
		D.append(d)
	D = np.array(D).astype('float32')
	return D
		
def main(argv = sys.argv[1:]):
	parser = OptionParser(usage = USAGE, version="Version: " + __version__)
	
	# Required arguments
	requiredOptions = OptionGroup(parser, "Required options",
		"These options are required to run BOTA, and may be supplied in any order.")
	
	requiredOptions.add_option("-c", "--config", type = "string", metavar = "FILE",
		help = "The configuration file to define a project.")
		
	requiredOptions.add_option("-o", "--outdir", type = "string", metavar = "DIR",
		help = "The output directory of BOTA. If it doesn't exist, BOTA will create it.")
							
	parser.add_option_group(requiredOptions)

	# Optional arguments that need to be supplied if not the same as default
	optOptions = OptionGroup(parser, "Optional parameters",
		"These options are optional, and may be supplied in any order.")

	optOptions.add_option("-t", "--nproc", type = "int", default = 1, metavar = "INT",
		help = "Number of processor for BOTA to use [default: 1; set 0 if you want to use all CPUs available].")
	
	parser.add_option_group(optOptions)
	
	(options, args) = parser.parse_args(argv)
	
	# load allele first
	if options.config == None:
		parser.error('[FATAL] No input configuration file specified.\n')
		exit(1)
	if options.outdir == None:
		parser.error('[FATAL] No output directory specified.\n')
		exit(1)
	if not os.path.exists(options.outdir):
		try: 
			os.mkdir(options.outdir)
		except:
			sys.stderr.write('[FATAL] Error in creating output directory:%s\n' % options.outdir)
			exit()
	
	logfile = os.path.join(options.outdir, 'project.log')
	logfh = open(logfile, 'wb')	
	project_info = parse_config(options.config)
	project_info.print_projInfo()
	project_info.print_projInfo(stream=logfh)
	logfh.close()
	
	sys.stdout.write('################# Data Preparation ################\n')
	# copy hmmtop db files
	current_dir = os.getcwd()
	for o in project_info.hmmtop_arch:
		t = os.path.join(current_dir, os.path.basename(o))
		if not os.path.exists(t): shutil.copyfile(o, t)
	
	# iterate through genomes
	for genomeID in project_info.genomes:
		genome_pkl = os.path.join(options.outdir, '%s.pkl' % genomeID)
		if os.path.exists(genome_pkl): continue
		genome_proj = project_info.genomes[genomeID]
		genome_dir = os.path.join(options.outdir, genomeID)
		if not os.path.exists(genome_dir): os.mkdir(genome_dir)
		
		# ID the AA sequences and genes
		sys.stdout.write('  [%s] Now extracting A.A. sequences.\n' % genomeID)
		faa_file = os.path.join(genome_dir, '%s.faa' % genomeID)
		if not os.path.exists(faa_file):
			if genome_proj.gff!=None:
				extract_gff_features(genome_proj, faa_file)
			else:
				prodigal = project_info.prodigal
				args = [prodigal, genomeID, genome_proj.fna, genome_dir, 'single']
				run_prodigal(args)
		
		# get Gram stain type
		sys.stdout.write('  [%s] Now deciding Gram stain type.\n' % genomeID)
		gram_file = os.path.join(genome_dir, '%s.gram' % genomeID)
		if not os.path.exists(gram_file):
			if genome_proj.gram == None:
				blat = project_info.blat
				db = project_info.gram_db
				genome_proj.gram = call_Gram([blat, db, faa_file, genome_dir, genomeID])
			with open(gram_file, 'wb') as gfh: gfh.write('%s\n' % genome_proj.gram)
		else:
			with open(gram_file, 'rb') as gfh: genome_proj.gram=gfh.readline().rstrip('\n')
		
		
		# hmmscan for domains
		sys.stdout.write('  [%s] HMMscan against Pfam-A for domain ID.\n' % genomeID)
		hmmscan_out = os.path.join(genome_dir, '%s.hmmscan' % genomeID)
		if not os.path.exists(hmmscan_out):
			pfam_db = project_info.pfam
			hmmscan = project_info.hmmscan
			args = [hmmscan, pfam_db, faa_file, hmmscan_out, options.nproc]
			run_hmmscan(args)
		
		# psort
		sys.stdout.write('  [%s] Determining protein subcellular location using PSORT.\n' % genomeID)
		psort_out = os.path.join(genome_dir, '%s.psort' % genomeID)
		if not os.path.exists(psort_out):
			psort = project_info.psort
			gram = genome_proj.gram
			psort_args = [psort, genome_dir, faa_file, gram, psort_out, options.nproc]
			run_psort(psort_args)
		
		# hmmtop
		sys.stdout.write('  [%s] Now predicting transmembrane structures using HMMTOP...\n' % genomeID)
		hmmtop_file = os.path.join(genome_dir, '%s.hmmtop' % genomeID)
		if not os.path.exists(hmmtop_file):
			hmmtop = project_info.hmmtop
			hmmtop_args = [project_info.hmmtop, faa_file, hmmtop_file, genome_dir, options.nproc]
			run_hmmtop(hmmtop_args)
		
		# PW Matrix
		sys.stdout.write('  [%s] Now PW Matrix scoring...\n' % genomeID)
		pwm_file = os.path.join(genome_dir, '%s.pwm' % genomeID)
		if not os.path.exists(pwm_file):
			loc_dict = {}
			# load the location information
			psort_out = os.path.join(genome_dir, '%s.psort' % genomeID)
			for line in open(psort_out, 'rb'):
				cols = line.rstrip().split('\t')
				geneID = cols[0].rstrip(' ')
				loc = cols[1]
				loc_dict[geneID] = loc
				
			aa_list, pwmatrix = cPickle.load(open(project_info.pwm_db, 'rb'))
			pw_dict = {}
			for aa_ind, aa in enumerate(aa_list):
				for pos_ind in range(9):
					pw_dict[(aa, pos_ind)] = pwmatrix[aa_ind, pos_ind]
			pwmscore_args = [faa_file, pw_dict, loc_dict, pwm_file, options.nproc]
			run_pwmscore(pwmscore_args)
			
		# make input data per GenomeID
		sys.stdout.write('  [%s] Integrating all data for DNN module...\n' % genomeID)
		# make data into genome_pkl
		integrate_data(genomeID, genome_dir, genome_pkl)
	
	# clean up		
	current_dir = os.getcwd()
	for o in project_info.hmmtop_arch:
		t = os.path.join(current_dir, os.path.basename(o))
		os.unlink(t)
	
	sys.stdout.write('################# DNN Model Prediction ################\n')
	# load models
	models = {}
	for genomeID in project_info.genomes:
		for allele in project_info.genomes[genomeID].alleles:
			model_arch, model_weights = project_info.models[allele]
			models[allele] = model_from_json(open(model_arch).read())
			models[allele].load_weights(model_weights)
	AAs = 'ARNDCQEGHILKMFPSTWYV'
	aa_dict = {}
	for ind, aa in enumerate(list(AAs)):
		aa_dict[aa] = ind
	
	for genomeID in project_info.genomes:
		genome_pkl = os.path.join(options.outdir, '%s.pkl' % genomeID)
		genome_data = cPickle.load(open(genome_pkl, 'rb'))
		
		candidates = []
		for geneID in genome_data:
			G = genome_data[geneID]
			if G.cellular_loc not in ['CellWall', 'Cellwall', 'Extracellular', 'OuterMembrane']: continue
			candidates += G.select_candidate()
		
		R = [] # store the results
		for pep in candidates:
			#if pep.gene != 'JXUN01000121.1|29853-30995|+|gene=TA05_10230': continue
			D = transform_seq(pep.seq, aa_dict)
			for allele in models:
				model = models[allele]
				try:
					Z = model.predict_classes(D, batch_size=128, verbose=0)
					S = map(itemgetter(1), model.predict(D, batch_size=128, verbose=0))
				except:
					continue
				
				ZX = []
				G = networkx.Graph()
				locs = []
				for ind, x in enumerate(Z):
					sx, nx = pep.scores[ind]
					#y = sx+nx/1e13
					if x==1 and sx>1e-10 and nx>=75:
						ZX.append(1)
						locs.append(ind)
					else: ZX.append(0)
				
				
				G.add_nodes_from(locs)
				for l1 in locs:
					for l2 in locs:
						if abs(l2-l1) <= 6: G.add_edge(l1, l2)
				for x in list(networkx.connected_components(G)):
					left, right = min(x), max(x)
					seq = pep.seq[left:right+9]
					z = sum(S[left:right+9])/(9+right-left)
					#z = max(S[left:right+9])
					if z < 0.6: continue
					s_coord = left+pep.start
					e_coord = right+9+pep.start
					chr = pep.gene.split('|')[0]
					R.append([allele, pep.gene, chr, s_coord, e_coord, seq, z]) 
				
		peptide_output = os.path.join(options.outdir, '%s.peptides.txt' % genomeID)
		pfh = open(peptide_output, 'wb')
		pfh.write('#peptide\tgene_name\tchr_acc\tgene_start\tgene_stop\tstrand\tpep_start\tpep_stop\tscore\n')
		for allele, gene, chr, s, e, seq, score in R:
			gene_start, gene_stop, strand, geneID = re.search('(\d+)\-(\d+)\|(.{1})\|gene\=(.+)$', gene).group(1,2,3,4)
			pep_start = s*3-2
			pep_stop = e*3
			pfh.write('%s\t%s\t%s\t%s\t%s\t%s\t%i\t%i\t%.6f\n' % \
				(seq, geneID, chr, gene_start, gene_stop, strand, pep_start, pep_stop, score))
		pfh.close()
		sys.stdout.write('  [%s] %i peptides predicted.\n' % (genomeID, len(R)))
	sys.stdout.write('Done.\n\n')
	return 0
	
	
if __name__ == '__main__': main()
