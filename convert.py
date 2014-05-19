#!/usr/bin/env python

import sys

f = open(sys.argv[1],'r')

sim = bool(int(sys.argv[2]))

for line in f:
	if line[0] == "%":
		continue
	else:
		W, H, P = [int(x) for x in line.split()]
		break

print(W,H,P)

M = [[] for x in xrange(W)]

zap = 0

for line in f:
	c, r, x = line.split()
	r = int(r)-1
	c = int(c)-1
	M[r].append("%s %s" % (str(c+1), x))
	zap += 1
	if r != c and sim:
		M[c].append("%s %s" % (str(r+1), x))
		zap += 1

print(zap, float(zap)/(float(H*W))*100.)

fo = open("matrix.txt",'w')

fo.write("%s %s\n" % (str(W), str(H)))

for r in xrange(H):
	for x in M[r]:
		fo.write("%s\t" % x)
	fo.write("0\n")

fo.close()

fo = open("vector.txt",'w')

for x in xrange(H):
	fo.write("%s\t" % str((x-H/2)/float(H)))
fo.close()

res = open("res.csv","a")
#res = open("dec.csv","a")	
res.write("%s; %s; %.2f%%; " % (sys.argv[1].split('/')[2].split('.')[0], str(W), float(zap)/(float(H*W))*100.) )
res.close();
