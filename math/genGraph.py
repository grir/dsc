#!/usr/bin/env python
# Gnp
# Autorius: I.Grinis
  
import random
import subprocess
import os
import sets

class Klaida(Exception):
   def __init__(self, reiksme):
      self.value = reiksme
   def __str__(self):
      return repr(self.value)

###############################################

def fact(n):
   z = 1L
   for i in range(2,n+1):
     z = z * i
   return z

###############################################

def permutation(x):
   xx = [] + x
   random.shuffle(xx)
   return xx

###############################################

def Gnp(n,p):
   gr = [[0] * n for i in range(n)]
   for i in range(n):
      for j in range(i+1,n):
         if random.random() < p:
            gr[i][j] = 1
            gr[j][i] = 1     
   return gr

###############################################

def GnM(n, M):
   gr = [[0] * n for i in range(n)]
   nv = []
   for i in range(n):
      for j in range(i+1,n)
         nv.append([i,j])
   random.shuffle(nv)
   for v in range(M):
      gr[v[0]][v[1]] = 1
      gr[v[1]][v[0]] = 1
   return gr


###############################################
def Cube(dim):
   n = 2 ** dim
   gr = [[0]*n for i in range(n)] 
   for i in range(n):
      for j in range(i+1,n):
         z = sum(eval(x) for x in list(bin(i ^ j)[2:]))
         if(z==1):
            gr[i][j] = 1
            gr[j][i] = 1     
   return gr

###############################################

def isSubgraph(gr1, gr2): # Is gr1 a subgraph of gr2
    n = len(gr1)
    for i in range(n):
      for j in range(i+1,n):
        if gr1[i][j] == 1:
           if  gr2[i][j] != 1:
              return False
    return True

###############################################

def isSubgraphPerm(gr1, gr2, perm): # Is gr1 a subgraph of gr2
    n = len(gr1)
    for i in range(n):
      for j in range(i+1,n):
        if gr1[i][j] == 1:
           if  gr2[ perm[i] ][ perm[j] ] != 1:
              return False
    return True


###############################################

def degs(graph): # 
   return [sum(ss) for ss in graph]

###############################################
# main:
z = 0
d = 4
n = 2 ** d
dd = fact(d) * n 
gr1 = Gnp(n, 0.5)
gr2 = Cube(3)
num = fact(n) / dd
print num
for i in range(num):
   perm = permutation(range(n))
   #print perm
   #print gr2
   if isSubgraphPerm(gr2, gr1, perm):
      z = z + 1
print z


