import os, sys, re, glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO

ofh = open('gram_positive.faa', 'w')
rec_ind = {}
for fasta in glob.glob('gram_positive/*/*.ffn'):
	x = os.path.dirname(fasta).split('/')[-1]
	if x not in rec_ind: rec_ind[x]=0
	for record in SeqIO.parse(fasta, 'fasta'):
		rec_ind[x] += 1
		protSeq = str(record.seq.translate()).replace('*','', 100)
		tag = '%s|%i' % (x, rec_ind[x])
		ofh.write('>%s\n%s\n' % (tag, protSeq))
ofh.close()