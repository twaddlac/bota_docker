import sys, os, re, glob
import random

species = {'p__Actinobacteria':[], 'p__Chloroflexi':[], 'o__Thermales':[],
			 'c__Bacilli':[], 'c__Clostridia':[], 'k__Archaea':[]}
for line in open('gsize.db', 'rb'):
	cols = line[:-1].split('\t')[0].split('|')
	sp = cols[-1]
	for x in cols:
		if x in species: species[x].append([sp, ''])

files = {}
for bz2 in glob.glob('ref_db/*.bz2'):
	sp = bz2.split('.')[-4]
	files[sp] = bz2

for x in species:
	for i, (sp, y) in enumerate(species[x]):
		if sp in files:
			y1 = files[sp]
			species[x][i][1] = y1
			
for x in species:
	if not os.path.exists('gram_positive/%s' % x): os.mkdir('gram_positive/%s' % x)
	files = []
	for sp, file in species[x]:
		if not os.path.exists(file): continue
		files.append(file)
	for file in random.sample(files, min(len(files), 20)):
		os.system('cp %s gram_positive/%s/' % (file, x))