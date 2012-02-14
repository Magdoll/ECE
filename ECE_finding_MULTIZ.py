import os,re,sys
import time, gzip
import numpy as np
from itertools import combinations,groupby
from collections import defaultdict
from cPickle import *

#import psyco
#psyco.full()

DIRNAME = './data/' #DIRNAME = '/projects/instr/multiz44way/'
OUTPUT_DIR = './results/'

species44 = ['hg18', 'panTro2', 'gorGor1', 'ponAbe2', 'rheMac2', 'calJac1', 'tarSyr1', 'micMur1', 'otoGar1', 'tupBel1', 'mm9', 'rn4', 'dipOrd1', 'cavPor3', 'speTri1', 'oryCun1', 'ochPri2', 'vicPac1', 'turTru1', 'bosTau4', 'equCab2', 'felCat3', 'canFam2', 'myoLuc1', 'pteVam1', 'eriEur1', 'sorAra1', 'loxAfr2', 'proCap1', 'echTel1', 'dasNov2', 'choHof1', 'monDom4', 'ornAna1', 'galGal3', 'taeGut1', 'anoCar1', 'xenTro2', 'tetNig1', 'fr2', 'gasAcu1', 'oryLat2', 'danRer5', 'petMar1']

species44_common_name = {\
		'hg18'   : 'Human       ',\
		'panTro2': 'Chimp       ',\
		'gorGor1': 'Gorilla     ',\
		'ponAbe2': 'Orangutan   ',\
		'rheMac2': 'Rhesus      ',\
		'calJac1': 'Marmoset    ',\
		'tarSyr1': 'Tarsier     ',\
		'micMur1': 'Mouse Lemur ',\
		'otoGar1': 'Bushbaby    ',\
		'tupBel1': 'TreeShrew   ',\
		'mm9'    : 'Mouse       ',\
		'rn4'    : 'Rat         ',\
		'dipOrd1': 'Kangaroo Rat',\
		'cavPor3': 'Guinea Pig  ',\
 		'speTri1': 'Squirrel    ',\
		'oryCun1': 'Rabbit      ',\
		'ochPri2': 'Pika        ',\
		'vicPac1': 'Alpaca      ',\
		'turTru1': 'Dolphin     ',\
		'bosTau4': 'Cow         ',\
		'equCab2': 'Horse       ',\
		'felCat3': 'Cat         ',\
		'canFam2': 'Dog         ',\
		'myoLuc1': 'Microbat    ',\
		'pteVam1': 'Megabat     ',\
		'eriEur1': 'Hedgehog    ',\
		'sorAra1': 'Shrew       ',\
		'loxAfr2': 'Elephant    ',\
		'proCap1': 'Rock hyrax  ',\
		'echTel1': 'Tenrec      ',\
		'dasNov2': 'Armadillo   ',\
		'choHof1': 'Sloth       ',\
		'monDom4': 'Opossum     ',\
		'ornAna1': 'Platypus    ',\
		'galGal3': 'Chicken     ',\
		'taeGut1': 'Zebra finch ',\
		'anoCar1': 'Lizard      ',\
		'xenTro2': 'X.tropicalis',\
		'tetNig1': 'Tetraodon   ',\
		'fr2'    : 'Fugu        ',\
		'gasAcu1': 'Stickleback ',\
		'oryLat2': 'Medaka      ',\
		'danRer5': 'Zebrafish   ',\
		'petMar1': 'Lamrepy     '}

# ######################################################################################
# SETTINGS
#
# m = 40   (at least 40 species)
# c = 0.8  (at least 80% conserved columns)
# s = 100  (at least 100 columns)
# 
# (1) the version of human genome is hg18 so the id for chromosome X would be hg18.chrX
# (2) gap characters are always -
# ######################################################################################

def read_maf_n_filter1(filename):
	"""
	Step 1: reads the .maf (or .maf.gz) file ( which is <filename> stored under <DIRNAME> )
	        discard blocks that are both:
			 (a) have < 40 species
			 AND
			 (b) < 50 bp 

			 and merge blocks that pass the above if they are separated by no more than 50 bp

	Returns: stats1 --- a list of discarded blocks ( # of rows, # of cols )
	         all_passed_blocks --- a list of passed-filter blocks, each block
			            is a dict of dict { species_id --> {'id','start','nogap_len','len','strand','seq'} }
	"""
	stats_filter1 = [] # contains, for each block discarded because it had less than 40 species
	                       # the number of columns in that block
	all_passed_blocks = []

	if filename.endswith('.gz'):
		f = gzip.open( os.path.join(DIRNAME, filename) )
	else:
		f = open(      os.path.join(DIRNAME, filename) )
	while 1:
		line = f.readline().strip()
		if line.startswith('##eof'):
			break
		if line.startswith('a score='):
			a_block = {}
			while 1:
				line = f.readline().strip()
				if len(line) == 0:
					break
				if line.startswith('s '):
					raw = line.split()
					nogap_len = len(raw[6])-raw[6].count('-')	
					assert nogap_len == int(raw[3])
					k = raw[1][ :raw[1].find('.') ] # hopefully all species IDs are like hg18.chr21
					a_block[k] = {'id':raw[1],\
							'start':int(raw[2]),\
							'nogap_len':nogap_len,\
							'len':len(raw[6]),\
							'strand':raw[4],\
							'seq':raw[6].upper().replace('U','T')}
			# discard ANY block with no human
			if 'hg18' not in a_block:
				continue
			# filter 1: see if this block has potential
			N = len(a_block)
			if N < 40 and a_block['hg18']['len'] >= 50:
				del a_block['hg18']['seq']
				stats_filter1.append( a_block['hg18'] )
				discard_notOK_blocks(all_passed_blocks, stats_filter1)
			else:
				if N >= 40:
					a_block['OK'] = True
					# find the last OK and see if the not OKs behind it accumulate to >= 50 bp
					j = len(all_passed_blocks) - 1
					while j >= 0:
						if all_passed_blocks[j]['OK']: break
						j -= 1
					if sum([ b['hg18']['len'] for b in all_passed_blocks[(j+1):] ]) >= 50:
						discard_notOK_blocks(all_passed_blocks, stats_filter1)
				else:
					a_block['OK'] = False
				all_passed_blocks.append( a_block )
	f.close()

	if len(all_passed_blocks) > 0 and not all_passed_blocks[-1]['OK']:
		discard_notOK_blocks(all_passed_blocks, stats_filter1)

	return stats_filter1, all_passed_blocks

def discard_notOK_blocks(all_passed_blocks, stats_filter1):
	sys.stdout.write('.')
	while len(all_passed_blocks) > 0:
		b = all_passed_blocks.pop()
		if b['OK']:
			all_passed_blocks.append ( b )
			break
		
		del b['hg18']['seq']
		stats_filter1.append( b['hg18'] )

def read_pickle_n_merge_blocks(obj_or_picklefilename,type):
	merged_blocks = []

	if type=='filename':
		with open(obj_or_picklefilename) as handle:
			tmp = load( handle )
			all_passed_blocks = tmp['all_passed_blocks']
	else:
		all_passed_blocks = obj_or_picklefilename
	
	print >> sys.stderr, "before merging....{0} blocks".format( len(all_passed_blocks) )

	for b in all_passed_blocks:
		del b['OK']

	if len(all_passed_blocks) <= 1:
		return all_passed_blocks

	merged_blocks = [all_passed_blocks.pop(0)]
	while len(all_passed_blocks) > 0:
		cur = all_passed_blocks.pop(0)
		if not check_merge_blocks(merged_blocks[-1], cur):
			merged_blocks.append( cur )

	print >> sys.stderr, "after merging....{0} blocks".format( len(merged_blocks) )

	return merged_blocks

def filter2_on_merged_blocks(merged_blocks):
	"""
	Remove from merged_blocks all blocks that
	(1) have less than 100 columns AND/OR
	(2) have less than 40 rows that each have at least 80 non-gap columns
	(3) does not pass the potential_UCE test run
	Returns: stats_filter2, all_passed_blocks
	"""
	stats_filter2 = []

	i = 0
	while i < len(merged_blocks):
		cur = merged_blocks[i]
		print >> sys.stderr, "aloha I'm doing ", i
		if cur['hg18']['len'] < 100 or len(filter(lambda x: x['nogap_len']>=80, cur.itervalues())) < 40 or not potential_UCE(cur):
			print >> sys.stderr, "filtering {0}".format(cur['hg18']['start'])		
			merged_blocks.pop(i)
			del cur['hg18']['seq']
			stats_filter2.append( cur['hg18'] )
		else:
			i += 1

	return stats_filter2, merged_blocks

GAP_CHAR = defaultdict(lambda: '-')
GAP_CHAR['hg18'] = 'N'

def check_merge_blocks(a, b):
	"""
	merge block b into block a
	if human coords are < 50 bp apart
	"""
	human_dist_between = b['hg18']['start'] - a['hg18']['start'] - a['hg18']['nogap_len']

	if human_dist_between >= 50:
		return False

	len_of_a = a['hg18']['len']
	len_of_b = b['hg18']['len']

	print >> sys.stderr, "merging blocks {0} to {1}".format(b['hg18']['start'], a['hg18']['start'])

	for sp in b:
		if sp not in a:
			# sp is only in b
			new_seq = GAP_CHAR[sp] * (len_of_a + human_dist_between)
			new_seq += b[sp]['seq']
			a[sp] = {'nogap_len': b[sp]['nogap_len'],\
					'seq': new_seq,\
					'len': len(new_seq),\
					'start': b[sp]['start'],\
					'id': sp,\
					'strand': b[sp]['strand']}
		else:
			# they both have it, excellent!
			a[sp]['seq'] += (GAP_CHAR[sp] * human_dist_between) + b[sp]['seq']
			a[sp]['nogap_len'] += b[sp]['nogap_len']
			a[sp]['len'] += human_dist_between + b[sp]['len']

	for sp in a:
		if sp not in b:
			# sp is only in a
			a[sp]['seq'] += GAP_CHAR[sp] * (human_dist_between + len_of_b)
			a[sp]['len'] += human_dist_between + len_of_b
	
	return True

def iterative_UCE_finding(merged_block):
	species = merged_block.keys()
	N = len(species)

	remaining_regions = [(0, merged_block['hg18']['len'])] # 0-based starts, 1-based ends
	all_found_bests = defaultdict(lambda: []) # key is r (combo of species), value is list of found bests

	prev_failed_conss = set()

	pair_iden = calc_pair_iden(merged_block)

	for subset_size in xrange(N, 39, -1):
		if len(remaining_regions) == 0:
			break

		# PRIORITIZE combinations that have a higher sum of pair_idens
		combos = []
		for r in combinations(species, subset_size):
			good_combo = True
			_score = 0
			for pr in combinations(r, 2):
				if pair_iden[pr] < 80:
					good_combo = False
					break
				_score += pair_iden[pr]
			if good_combo:
				combos.append( {'r':r, 'score':_score} )
		combos.sort(key=lambda x:x['score'], reverse=True)

		for r in combos: #for r in combinations(species, subset_size):
			r = r['r']
			if len(remaining_regions) == 0:
				break

#			if not good_combo_by_pair_iden(r, pair_iden):
#				print >> sys.stderr, "bad combo"
#				continue

			print >> sys.stderr, "trying size-{0} combination, remaining regions: {1}".format(subset_size, remaining_regions)

			s, cons, not_all_gap = setup_for_UCE(r, merged_block)

			# go through each remaining region
			reg_index = 0
			while reg_index < len(remaining_regions):
				reg_start, reg_end_1 = remaining_regions[reg_index]
				if cons[reg_start:reg_end_1] in prev_failed_conss:
					print >> sys.stderr, "no need to try"
					reg_index += 1
					continue

				prebests = find_UCE(s  = s[reg_start:reg_end_1], \
						not_all_gap = not_all_gap[reg_start:reg_end_1], \
						L           = reg_end_1-reg_start, \
						index_offset= reg_start)
				# make the bests 0-based starts with 1-based ends
				# we also want to trim the head and tail
				bests = []
				for (x,y) in prebests:
					# the region is [x-1, y)
					x -= 1
					while cons[x] != '*': x += 1
					while cons[y-1] != '*': y -= 1
					bests.append( (x,y) )

				if len(bests) == 0:
					prev_failed_conss.add( cons[reg_start:reg_end_1] )
					reg_index += 1
				else:
					print >> sys.stderr, "found bests: {0}".format(bests)
					remove_overlapping_bests(bests)
					print >> sys.stderr, "non-overlapping bests: {0}".format(bests)
					all_found_bests[r] += bests

					for i in xrange(1, len(bests)):
						if bests[i][0] - bests[i-1][1] >= 100 and potential_UCE(merged_block,(bests[i-1][1],bests[i][0])):
							remaining_regions.append( (bests[i-1][1], bests[i][0]) )
					if bests[0][0] - reg_start >= 100 and potential_UCE(merged_block,(reg_start,bests[0][0])):
						remaining_regions.append( (reg_start, bests[0][0]) )
					if reg_end_1 - bests[-1][1] >= 100 and potential_UCE(merged_block,(bests[-1][1],reg_end_1)):
						remaining_regions.append( (bests[-1][1], reg_end_1) )
					remaining_regions.pop(reg_index)
					print >> sys.stderr, "after: remaining regions: {0}".format(remaining_regions)
	return dict(all_found_bests)

def find_longest_UCE_for_subset_size(merged_block, subset_size):
	species = merged_block.keys()
	L = merged_block['hg18']['len']

	optimal_best = {'r':None, 'i':None, 'j_1':None, 'len':-1}
	for r in combinations(species, subset_size): 
		s, cons, not_all_gap = setup_for_UCE(r, merged_block)
		bests = find_UCE(s, not_all_gap, L, 0)
		for i,j_1 in bests:
			l = j_1 - i
			if l > optimal_best['len']:
				optimal_best['r'] = r
				optimal_best['i'] = i
				optimal_best['j_1'] = j_1
				optimal_best['len'] = l

	return optimal_best

def dirty_test1():
	merged_blocks = load(open('results/chr20.maf_filter2_passedNmerged.pickle'))['merged_blocks']
	start_t = time.time()
	num_of_blocks = len(merged_blocks)-1
	for ind,m in enumerate(merged_blocks):
		N = len(m)
		for subset_size in xrange(N, 39, -1):
			print >> sys.stderr, "running {0}/{1} for subsetsize {2}".format(ind, num_of_blocks, subset_size)
			find_longest_UCE_for_subset_size(m, subset_size)
			print >> sys.stderr, "current we've spent {0} secs".format(time.time()-start_t)
	end_t = time.time()
	print >> sys.stderr, "the whole thing took {0} secs".forat(end_t-start_t)
	return end_t-start_t

def calc_pair_iden(block):
	L = block['hg18']['len']
	N = len(block)
	pair_iden = {}

	for r in combinations(block.keys(), 2):
		s1 = block[r[0]]['seq']
		s2 = block[r[1]]['seq']
		iden = map(lambda i: s1[i]==s2[i] and s1[i]!='-', xrange(L)).count(True)
		pair_iden[(r[0],r[1])] = iden
		pair_iden[(r[1],r[0])] = iden

	return pair_iden

def good_combo_by_pair_iden(r, pair_iden):
	for pr in combinations(r, 2):
		if pair_iden[pr] < 80:
			return False
	return True

def test(merged_blocks, f):
	result = {}
	for m_i,m in enumerate(merged_blocks):
		all_bs = iterative_UCE_finding(m)
		result[m_i] = all_bs
		for r, list_of_bs in all_bs.iteritems():
			s, cons, not_all_gap = setup_for_UCE(r, m)
			for i,j_1 in list_of_bs:
				print_UCE(m, r, cons, not_all_gap, i, j_1, f)
				f.flush()
	return result

def remove_overlapping_bests(bests):
	def calc_overlap((a1,a2),(b1,b2)): return min(b2,a2)-b1+1 if a1 <= b1 else min(b2,a2)-a1+1

	i = 1
	while i < len(bests):
		if calc_overlap(bests[i-1], bests[i]) > 0:
			len_diff = (bests[i-1][1]-bests[i-1][0]) - (bests[i][1]-bests[i][0])
			if len_diff < 0 or (len_diff == 0 and bests[i][0] < bests[i-1][0]):
				bests[i-1] = bests[i]
			bests.pop(i)
		else:
			i += 1

NT_INDEX = {'A':0,'T':1,'C':2,'G':3}
def potential_UCE(block, region=None):
	if region is None:
		start, end_1 = 0, block['hg18']['len']
		L = end_1
	else:
		start, end_1 = region
		L = end_1 - start

	# s will be the binary array for UCE algorithm, s[0] always 0 for UCE finding
	# s[i] = +1 if at least 40 species are conserved at column i (1-based)
	# otherwise s[i] = -4
	s = [0] * (L + 1)
	tally = np.zeros( (L, 4) )
	for v in block.itervalues():
		for col,x in enumerate(v['seq'][start:end_1]):
			try:
				tally[col, NT_INDEX[x]] += 1
			except KeyError:
				pass

	for col in xrange(L):
		s[col+1] = +1 if any( tally[col,:] >= 40) else -4
	
	# now we can run the UCE algorithm on s
	L += 1
	r = [0] * L
	X = [0] * L
	Y = [0] * L
	for i in xrange(1, L): r[i] = r[i-1] + s[i]
	for i in xrange(1, L): X[i] = min(X[i-1], r[i])
	Y[-1] = r[-1]
	for i in xrange(L-2, -1, -1): Y[i] = max(Y[i+1], r[i])

	j = 0
	i = 0
	while j < L:
		if j==L-1 or Y[j+1] < X[i]:
			if j - i >= 100:
				return True # has potential UCEs! stop here!
			j += 1
			while j < L and i < L and Y[j] < X[i]: i += 1
#			old_i = i
#			while i < L and X[i]==X[old_i]: i+=1
#			j = i
		else:
			j += 1

#	bests = filter(lambda x: x[1]-x[0]+1>=100, bests)
#	return len(bests)>0
	return False

NT_INDEX_PLUS_GAP = {'A':0,'T':1,'C':2,'G':3,'-':4,'N':5}
def setup_for_UCE(r, block):
	n = len(r)
	L = block['hg18']['len']
	s = [0] * L
	tally = np.zeros( (L, 6), dtype=float )
	not_all_gap = [True] * L
	for sp in r:
		for col,x in enumerate(block[sp]['seq']):
			try:
				tally[col, NT_INDEX_PLUS_GAP[x]] += 1
			except KeyError:
				print >> sys.stderr, "weird character {0} encountered".format(x)
				tally[col, 5] += 1

	cons = ''
	for col in xrange(L):
		tally[col,:] /= n
		if tally[col, NT_INDEX_PLUS_GAP['-']] == 1.0:
			# all gaps at this column
			not_all_gap[col] = False
			s[col] = 0
			cons += '-'
		elif any(tally[col,0:4]==1.0):
			s[col] = +1
			cons += '*'
		else:
			s[col] = -4
			cons += ' '
	
	return s, cons, not_all_gap

def find_UCE(s, not_all_gap, L, index_offset):
	# now we can run the UCE algorithm on s
	s = [0] + s
	L += 1
	r = [0] * L
	X = [0] * L
	Y = [0] * L
	for i in xrange(1, L): r[i] = r[i-1] + s[i]
	for i in xrange(1, L): X[i] = min(X[i-1], r[i])
	Y[-1] = r[-1]
	for i in xrange(L-2, -1, -1): Y[i] = max(Y[i+1], r[i])

	# find i,j with Y[j]>=X[i] maximizing j-i
	bests = []
	j = 0
	i = 0
	while j < L:
		if j==L-1 or Y[j+1] < X[i]:
			if j - i >= 100 and sum(not_all_gap[(i+1):(j+1)]) >= 100:
				bests.append( ( (i+1)+index_offset, j+index_offset ) )
			j += 1
			while j < L and i < L and Y[j] < X[i]: i += 1
#			old_i = i
#			while i < L and X[i]==X[old_i]: i+=1
#			j = i
		else:
			j += 1

	return bests

def print_UCE(block, r, cons, not_all_gap, i, j_1, output_f=sys.stdout):
	# the UCE is from block[i : j_1]
	chr = block['hg18']['id']
	chr = chr[ chr.find('.')+1 : ]

	human_start = block['hg18']['start'] + i - block['hg18']['seq'][:i].count('-')
	human_end   = block['hg18']['start'] + j_1 - 1 - block['hg18']['seq'][:j_1].count('-')

	num_species = len(r)
	uce_length  = j_1 - i - not_all_gap[i:j_1].count(False)
	uce_pid     = cons[i:j_1].count('*') * 1. / uce_length
	unmatched   = list(set(species44).difference(r))
	unmatched.sort()
	names = list(species44)
	for x in unmatched: names.remove( x )

	# don't forget to output the human coords as 1-based
	output_f.write(chr + ':' + str(human_start+1) + '-' + str(human_end+1) + '\n')
	output_f.write(str(uce_length) + '\n')
	output_f.write(str(len(r)) + '\tunmatched: ' + ", ".join(unmatched) + '\n')
	output_f.write("{0:.2f}".format(uce_pid) + '\n')

	for sp in names:
		s = ''
		for col in xrange(i, j_1):
			if not_all_gap[col]: s += block[sp]['seq'][col]
		output_f.write("{0:7s}{1:18s}{2}\n".format(sp, '('+species44_common_name[sp]+')', s))

	new_cons = ''
	for i in xrange(i, j_1):
		if not_all_gap[i]: new_cons += cons[i]

	output_f.write(" "*25 + new_cons + "\n")

if __name__=="__main__":
	start_t = time.time()
	stats_filter1, all_passed_blocks = read_maf_n_filter1(sys.argv[1])
	end_t = time.time()
	print >> sys.stderr, "time spent:", (end_t-start_t)
	with open( os.path.join(OUTPUT_DIR, sys.argv[1]+'_filter1_N40_passed.pickle'),'w') as handle:
		dump({'time_spent':end_t-start_t, 'stats_filter1':stats_filter1, 'all_passed_blocks':all_passed_blocks}, handle)


	start_t = time.time()

	merged_blocks = read_pickle_n_merge_blocks(all_passed_blocks, type='obj')
	stats_filter2, merged_blocks = filter2_on_merged_blocks( merged_blocks )

	end_t = time.time()

	with open( os.path.join(OUTPUT_DIR, sys.argv[1]+'_filter2_passedNmerged.pickle'), 'w') as handle:
		dump({'time_spent':end_t-start_t, 'stats_filter2':stats_filter2, 'merged_blocks':merged_blocks}, handle)

#	merged_blocks = load(open(sys.argv[2]))['merged_blocks']

	start_t = time.time()
	with open( os.path.join(OUTPUT_DIR, sys.argv[1]+'_ece.txt'), 'w') as f:
		result = test(merged_blocks, f)
	end_t = time.time()
	print >> sys.stderr, "time spent on UCE finding: {0}".format(end_t-start_t)

	with open( os.path.join(OUTPUT_DIR, sys.argv[1]+'_ece.pickle'), 'w') as handle:
		dump({'time_spent':end_t-start_t, 'result':result}, handle)
