#!/usr/bin/env python
def defA(n):
  return eval('1'*n+'L')

def defB(m):
  return eval('7'*m+'L')

def rez(m,n):
  return defA(n)/defB(m) 

def rez2(m,n):
  return (defA(n)/defA(m))/7 

