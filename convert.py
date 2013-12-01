#!/usr/bin/env python

import sys

f = open(sys.argv[1],'r')

for line in f:
	if line[0] == "%":
		continue
	else:
		W, H, P = [int(x) for x in line.split()]
		break

print(W,H,P)

M = [[] for x in xrange(W)]

for line in f:
	c, r, x = line.split()
	r = int(r)-1
	c = int(c)-1
	M[r].append("%s %s" % (str(c+1), x))
	if r != c:
		M[c].append("%s %s" % (str(r+1), x))

fo = open("matrix.txt",'w')

fo.write("%s %s\n" % (str(W), str(H)))

#for r in xrange(H):
#	for c in xrange(W):
#		if M[r][c] != '0':
#			fo.write("%s %s \t" % (str(c), M[r][c]))
#	fo.write("0\n")

for r in xrange(H):
	for x in M[r]:
		fo.write("%s\t" % x)
	fo.write("0\n")

fo.close()

fo = open("vector.txt",'w')

for x in xrange(H):
	fo.write("%s\t" % str((x-H/2)/float(H)))

fo.close()
	