#!/usr/bin/env python
# Gnp
# Autoriai: I.Grinis, V.Gruslys
  
import random
import subprocess
import os
import sets
import math
import numpy
import copy
import matplotlib.pyplot as plt

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
      for j in range(i+1,n):
         nv.append([i,j])
   random.shuffle(nv)
   nnv = nv[0:M]
   for v in nnv:
      gr[v[0]][v[1]] = 1
      gr[v[1]][v[0]] = 1
   return gr

###############################################
def Wn(n):
   gr = [[0] * n for i in range(n)]
   nv = []
   for i in range(1,n):
       nv.append([0,i])
   for i in range(1,n-1):
       nv.append([i,i+1])
   nv.append([n-1,1])
   rv = random.randint(2,n-3)
   nv.append([1,rv])
   for v in nv:
      gr[v[0]][v[1]] = 1
      gr[v[1]][v[0]] = 1
   return [gr, rv]


###############################################
def subgraph(gr, S):
   subgr = [[0] * S for i in range(S)]
   idx = [[0] * S for i in range(S)]
   n = len(gr) 
   rv = range(n)
   random.shuffle(rv)
   rv = rv[0:S]
   grx = gr
   for i in range(n):
      if not (i in rv):
         for j in range(n):
            grx[i][j] = -1
      else:
         for j in range(n):  
           if not j in rv:
              grx[i][j] = -1
   #print grx
   ix = 0 
   for i in range(len(grx)):
      jx = 0
      for j in range(len(grx)):
         if (grx[i][j] != -1):
            subgr[ix][jx] = grx[i][j]
            jx = jx + 1 
      if jx > 0:
         ix = ix +1

   return subgr

###############################################
def subgraph1(gr, idxToRemove):
# remove one vertex
   n = len(gr) 
   S = n - 1
   subgr = [[0] * S for i in range(S)]
   rv = range(n)
   rv = rv[0:idxToRemove] + rv[idxToRemove+1:]
   grx = gr
   for i in range(n):
      if not (i in rv):
         for j in range(n):
            grx[i][j] = -1
      else:
         for j in range(n):  
           if not j in rv:
              grx[i][j] = -1
   #print grx
   ix = 0 
   for i in range(len(grx)):
      jx = 0
      for j in range(len(grx)):
         if (grx[i][j] != -1):
            subgr[ix][jx] = grx[i][j]
            jx = jx + 1 
      if jx > 0:
         ix = ix +1

   return subgr



###############################################
def subgraph2(gr, listToRemove):
   n = len(gr) 
   S = n - len(listToRemove)
   subgr = [[0] * S for i in range(S)] 
   rv = list(set(range(n)) - set(listToRemove))
   random.shuffle(rv)
   grx = gr
   for i in range(n):
      if not (i in rv):
         for j in range(n):
            grx[i][j] = -1
      else:
         for j in range(n):  
           if not j in rv:
              grx[i][j] = -1
   #print grx
   ix = 0 
   for i in range(len(grx)):
      jx = 0
      for j in range(len(grx)):
         if (grx[i][j] != -1):
            subgr[ix][jx] = grx[i][j]
            jx = jx + 1 
      if jx > 0:
         ix = ix +1

   return subgr


###############################################

def degs(graph): # 
   return [sum(ss) for ss in graph]

###############################################

def avgDegs(graph): # 
   return sum([sum(ss) for ss in graph])*1.0/len(graph)

###############################################

def degs3(graph): # 
   zz = [sum(ss)<3 for ss in graph]
   return sum(zz)

###############################################
def getMinVert(graph):
   deg = min(degs(graph))
   mins = []
   for i in range(len(graph)):
      if (sum(graph[i]) - deg)<1:
         mins.append(i)
   random.shuffle(mins)
   return mins

###############################################
def getMinEdge(graph):
   if len(graph)<2:
      return [0,0]
   deg = degs(graph)
   mn = 2 * max(deg) 
   res = [0,1]
   for i1 in range(len(graph)):
      for i2 in range(i1+1,len(graph)):
         if (deg[i1]+deg[i2]) < mn:
             res = [i1,i2]
   if deg[res[0]]>deg[res[1]]:
      res = [res[1], res[0]]
   return res

###############################################

def getVert(graph, deg): # 
# returns a list of v. whose degree is < deg
   v = [i  for i in range(len(graph)) if sum(graph[i]) < deg]
   return v

###############################################
def writeGraph(graph, fn):
   f = open(fn, 'w')
   n = len(graph)
   f.write(str(n)+'\n')
   for i in range(n):
      for j in range(n):
         f.write(str(graph[i][j])+' ')
      f.write(str('\n'))
   f.close() 
###############################################
# main:
n = 10	
avgs = []
ng = 800000
epps = []
maxlast = 0.0
for i in range(ng):
  gr = GnM(n, 2 * n - 1)
  #print gr
  #gr, rv = Wn(n)
  #print n, ' ', sum(degs(gr))/2
#gr = GnM(n, int(2.5 * n))
#gr = Wn(n)
#print [max(degs(gr)), min(degs(gr))] 
#print degs3(gr)
#print sum(degs(gr))/1000.0
  eps = 1.0/math.exp(1.0)
  ssb = (1 - eps)
#print n - math.sqrt(n*1.0/48)
  #print ssb
  mns = []
  deg2get = 3 
# to make proper subgraph:
#vv = getVert(gr, n-1)

  mins = []
  maksai = []
  first = True
  last = -2.0
  grs = copy.deepcopy(gr) 
#  while len(gr)>0:
  while len(gr)>0:

#      avg.append(avgDegs(gr))
      vn = len(gr) 
      #if vn > 10: 
      #  dd = min(degs(gr)[0:vn/2])  
      #else:
      dd = min(degs(gr))  
      if (dd >= 3) and (len(gr)<n) and first:
         #epps.append(len(gr)*1.0/n)
         #first = False
         last = len(gr)*1.0/n
      mins.append(dd)
      #maksai.append(max(degs(gr)))
      gr = subgraph2(gr,getMinVert(gr)[0:1])
      #dds =  degs(gr)
      #dds.sort()
      #gr = subgraph2(gr,dds[0:1])
  epps.append(last)
  if (last>maxlast):
     maxlast = last
  avgs.append(mins)
  avgs.append(maksai)
  if (maxlast>((n-1)*1.0/n-0.00000001)):
     writeGraph(grs, 'gr'+str(n)+'.dat')
     break

print maxlast
print grs

#cl = 'r' 
#for i in range(ng):
#  cl = cl + '-'
#epps.sort()
#plt.plot(epps, 'r')
#plt.hist(epps, bins='auto')
#plt.plot(avgs[1], 'g')


#plt.xlabel(' eps ')
#plt.show()

#for i in range(n-1):
#  print 'Iteracija: ' + str(i)+' ---- '+str(len(gr))
#  vv = getVert(gr, deg2get)
#  gr = subgraph2(gr,vv)
#  print 'max deg ' + str(max(degs(gr))) + ' min deg '+ str(min(degs(gr))) 



