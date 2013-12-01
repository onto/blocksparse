#!/usr/bin/env python
# coding=UTF-8

import sys
from PIL import Image

if __name__ == '__main__':

    file = sys.argv[1]
    
    permut = sys.argv[2]

    f = open(file,'r')

    W, H = [int(x) for x in f.readline().split()]
    
    fp = open("permut.txt",'r')
    
    if permut == "0":
        P = [int(x) for x in xrange(W+1)]
    else:
        P = [int(x) for x in fp.readline().split()]
    
    #print(W, H)

    if permut == "1" and W != len(P)-1:
        print("Не совпадают размеры")
        print(W, len(P))

    SIZE = min(2000,W,H)

    wk = float(SIZE)/float(W)
    hk = float(SIZE)/float(H)

    #print (wk, hk)

    image = Image.new('RGB',(SIZE,SIZE),'white')

    image.load()

    r = 1
    for line in f:
        q = 0
        for c in line.split():
            q += 1
            if q % 2 == 1 and c != '0':
                x,y = int(float(P[r]-1)*hk), int(float(P[int(c)]-1)*wk)
                x = min(SIZE-1,x)
                y = min(SIZE-1,y)
                #print(r,int(c),x,y)
                image.putpixel((x, y),(0,0,255))
        r += 1

    image.save("matrix.png", "PNG")



