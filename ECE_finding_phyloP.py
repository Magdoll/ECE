import os, re, sys

import findECE_pure
try:
	import findECE
	ECEfunc = findECE.findECE
except:
	print >> sys.stderr, "Cannot import findECE. Using the pure python version"
	ECEfunc = findECE_pure.find_ECE

remove_overlapping_bests = findECE_pure.remove_overlapping_bests

rex = re.compile('fixedStep chrom=(\S+) start=(\d+) step=(\d+)')

def read_phyloP_n_run(filename, phylop_threshold, min_ece_length, outf):
	"""
	Reads the .wigFix (<filename> stored under <INPUT_DIR>) and discard blocks that are
	(a) < min_ece_length bp

	(no merging or num_species counting at this point)

	If a score is negative, pretend it is 0. 
	Otherwise transform the score' = score - phylop_threshold

	For each block that pass through the filter, run it through ECE.
	"""
	outf.write("Processing {0} with phyloP score cutoff {1} and minimum ECE length {2}\n".format(\
			filename, phylop_threshold, min_ece_length))

	s = None

	f = open(filename)
	for line in f:
			if line.startswith('fixedStep'):
				if s is not None and end_2 - start_1 >= min_ece_length:
					bests = map(lambda (x,y): (x-1,y), ECEfunc([0]+s, len(s)+1, min_ece_length))
					if len(bests) > 0:
						remove_overlapping_bests(bests)
						for x,y in bests:
							outf.write("{chr}\t{start}\t{end}\t{len}\n".format(chr=chr,\
									start = start_1 + x,\
									end   = start_1 + y - 1,\
									len   = y - x))

				# a new block
				m = rex.match(line.strip())
				if m is None:
					raise Exception, "Not a valid fixedStep line:\n{0}\nAbort!!".format(line)
				chr = m.group(1)
				start_1 = int(m.group(2)) # 1-based
				step = int(m.group(3))
				assert step == 1
				end_2 = start_1
				s = []
			else: # continue reading the same phyloP block into s
				end_2 += 1
				score = float(line.strip())
				score = max(0, score) - phylop_threshold
				s.append(score)
	f.close()

	if s is not None and end_2 - start_1 >= min_ece_length:
		bests = map(lambda (x,y): (x-1,y), ECEfunc([0]+s, len(s)+1, min_ece_length))
		if len(bests) > 0:
			remove_overlapping_bests(bests)
			for x,y in bests:
				outf.write("{chr}\t{start}\t{end}\t{len}\n".format(chr=chr,\
						start = start_1 + x,\
						end   = start_1 + y - 1,\
						len   = y - x))

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

def find_ECE(s, L, min_ece_length):
	"""
	s           --- the vector containing the -1 and +4s
	L           --- just the length of s

	NOTE: s[0] should be 0, the extra 0 we need for the ECE alg
	      so the real data should start from s[1], s[2]....

	Returns all valid ECEs (>= min_ece_length) in s, which could be overlapping
	Call remove_overlapping_bests afterwards to process overlapping ECEs
	"""
	assert len(s) == L and s[0] == 0
	# now we can run the ECE algorithm on s
	r = [0] * L # r[i] = sum( s[0], s[1], ... s[i] )
	X = [0] * L # X[i] = min( r[0], r[1], ... r[i] )
	Y = [0] * L # Y[i] = max( r[i], r[i+1] ...     )
	for i in xrange(1, L): r[i] = r[i-1] + s[i]
	for i in xrange(1, L): X[i] = min(X[i-1], r[i])
	Y[-1] = r[-1]
	for i in xrange(L-2, -1, -1): Y[i] = max(Y[i+1], r[i])

	# find i,j with Y[j] >= X[i] maximizing j-i
	bests = []
	j = 0
	i = 0
	while j < L:
		if j == L-1 or Y[j+1] < X[i]:
			if j - i >= min_ece_length:
				bests.append((i+1, j))
			j += 1
			while j < L and i < L and Y[j] < X[i]: i += 1
		else:
			j += 1

	return bests

if __name__=="__main__":
	from optparse import OptionParser, make_option

	parser = OptionParser(option_list=[\
			make_option("-f", "--filename", dest="filename", help=".wigFix filename"),\
			make_option("-p", "--phylo", default=5.1, type=float, dest="phylo", help="phyloP score threshold (default 5.1)"),\
			make_option("-m", "--min-ece-length", default=100, type=int, dest="min_ece_length", help="minimum ece length (default 100)"),\
			make_option("-o", "--output", dest="output", help="output filename (default to stdout)")])
	
	options, args = parser.parse_args()

	if options.filename is None:
		raise ValueError, "input filename must be specified!"
	if options.min_ece_length <= 0:
		raise ValueError, "minimum ece length should be > 0 bp!"
	if options.phylo <= 0:
		raise ValueError, "phyloP score cutoff should be > 0!"

	f = sys.stdout if options.output is None else open(options.output, 'w')
	read_phyloP_n_run(options.filename, options.phylo, options.min_ece_length, f)
	f.close()
