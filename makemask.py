#!/usr/bin/env python
# coding=UTF-8

import sys
from PIL import Image

f = open(sys.argv[1],'r')

W, H = [int(x) for x in f.readline().split()]

SIZE = min(2000,W,H)

wk = float(SIZE)/float(W)
hk = float(SIZE)/float(H)

image = Image.new('RGB',(SIZE,SIZE),'white')

image.load()

r = 1
for line in f:
	q = 0
	for c in line.split():
	    q += 1
	    if q % 2 == 1 and c != '0':
	        y,x = int(float(r-1)*hk), int(float(int(c)-1)*wk)
	        x = min(SIZE-1,x)
	        y = min(SIZE-1,y)
	        #print(r,int(c),x,y)
	        image.putpixel((x, y),(0,0,255))
	r += 1

image.save("matrix.png", "PNG")



