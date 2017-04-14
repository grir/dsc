#!/usr/bin/env python
import random
def Gnp(n,p):
   gr = []
   for i in range(n):
     gr.append([0]*n)
   
   for i in range(n):
     for j in range(i+1,n):
        br = 0
        if random.random() < p:
          br = 1
        gr[i][j] = br
        gr[j][i] = br 
   return gr

def getDegs(gr):
   return [sum(rw) for rw in gr]

gr = Gnp(2**3,1.0/2)  
print getDegs(gr)
