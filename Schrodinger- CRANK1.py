#!/usr/bin/env python
# coding: utf-8

# In[26]:


import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.integrate as si
import numpy.linalg as npla


# In[35]:


#this is my attempt at crank with the 1D heat equation

#variables
L=20
t0=0
tmax=L
x0=0
xmax=L
isteps=20
nsteps=20
dt=tmax/nsteps
dx=xmax/isteps
k=L**2/tmax
sigma=L/4
r=(k*dt)/(2*(dx)**2)
a=-r
b=(1+2*r)
c=-r
x=np.zeros(isteps)
t=np.zeros(nsteps)

#tridiagonal matrix
def tridiag(A, B, C, k1=-1, k2=0, k3=1):
    return np.diag(A, k1) + np.diag(B, k2) + np.diag(C, k3)
#put in crank 
#A =[a]*(nsteps-1); B = [b]*(nsteps); C = [c]*(nsteps-1)
#d = tridiag(A, B, C)

#T for t=0
def T(x,t0):
    T=np.exp((-1/2)*(x**2)/(sigma)**2)
    return T

#T=np.zeros(isteps)

def CrankNicolson(T,t,x,dt,dx):#NEED TO FINISH WRITING THIS put in the inverse matrix
    xi1=xi+dx
    tn1=tn+dt
    r=(k*dt)/(2*(dx)**2)
    Tn=r*T[xi1]+(1-2*r)*T[xi]+r*T[xi-dx]
    #inverse matrix goes here of all i
    Tn=-r*Tn1[xi1]+(1+2*r)*Tn1[xi]-r*Tn1[xi-dx]
    return Tn, tn1, xi

#for n in range(nsteps):
    #run crank for all t and plot for each
    

#plt.plot(T,x)
#plt.show()


# In[28]:


#this is me figuring out how to make a matrix...
nsteps=3
A=np.zeros((nsteps,nsteps))
A[0:3]=[2]*1

#print(A)

a=np.array([1,2,3,4])
b=np.diag(a)

#print(b)
a=1
b=2
c=3

def tridiag(A, B, C, k1=-1, k2=0, k3=1):
    return np.diag(A, k1) + np.diag(B, k2) + np.diag(C, k3)

A =[a]*(nsteps-1); B = [b]*(nsteps); C = [c]*(nsteps-1)
d = tridiag(A, B, C)

print(d)

D=npla.inv(d)
print(D)


# In[ ]:




