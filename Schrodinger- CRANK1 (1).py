#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.integrate as si
import numpy.linalg as npla


# In[8]:


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
#

#T for t=0
def T(x,t0):
    xi1=xi+dt
    Ti=np.exp((-1/2)*(xi**2)/(sigma)**2)
    Ti1=np.exp((-1/2)*(xi1**2)/(sigma)**2)
    return Ti1,xi1

T=np.zeros(isteps)

def CrankNicolson(T,t,x,dt,dx):#NEED TO FINISH WRITING THIS put in the inverse matrix
    xi1=xi+dx
    tn1=tn+dt
    r=(k*dt)/(2*(dx)**2)
    a=-r
    b=(1+2*r)
    c=-r
    for i in range(isteps):
        t0=t0
        xi, t0= T(xi,t0)
        x[i+1]=xi
        T[i+1]=Ti
    Tn=r*T[i+1]+(1-2*r)*T[i]+r*T[i-1]    
    D=np.matrix[Tn]
    def tridiag(A, B, C, k1=-1, k2=0, k3=1):
        return np.diag(A, k1) + np.diag(B, k2) + np.diag(C, k3)
    A =[a]*(isteps-1); B = [b]*(isteps); C = [c]*(isteps-1)
    d = tridiag(A, B, C)
    I=npla.inv(d)
    Tn1=np.zeros((isteps,1))
    Tn1=I*D
    #inverse matrix goes here of all i
    Tn1=-r*Tn1[i+1]+(1+2*r)*Tn1[i]-r*Tn1[i-1]
    return Tn1, tn1, xi1

#for n in range(nsteps):
    #run crank for all t and plot for each
    

#plt.plot(T,x)
#plt.show()


# In[7]:


#this is me figuring out how to make a matrix...
nsteps=3
A=np.zeros((nsteps,1))
A[0:3]=[2]*1

print(A)

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

#print(d)

D=npla.inv(d)
#print(D)


# In[ ]:




