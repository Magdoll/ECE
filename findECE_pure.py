import os, re, sys

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
