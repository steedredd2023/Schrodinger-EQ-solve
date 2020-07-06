#!/usr/bin/env python
# coding: utf-8

# In[10]:


import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.integrate as si
import numpy.linalg as npla


# In[11]:


#this is my 2nd attempt at crank with the 1D heat equation

#variables
L=20
t0=0
tmax=2
x0=0
xmax=L
isteps=4*L
nsteps=10
dt=0.001
dx=xmax/isteps
k=L**2
r=(k*dt)/(2*(dx)**2)
sigma=L/12
x=np.zeros(isteps+1)

for i in range(isteps):
    xi=x0+(i+1)*dx
    x[i+1]=xi
    
T=np.zeros(isteps)

for i in range(1,isteps):#gives first array of T vs x
    xi=x0+(i)*dx
    Ti=np.exp((-1/2)*((xi-(L/2))**2)/(sigma)**2)
    T[i]=Ti

#make array of arrays of T
Tn=[]
Tn.append(T)

def CrankNicolson(T):#NEED TO FINISH WRITING THIS put in the inverse matrix
    a=-r
    b=(1+2*r)
    c=-r#these are also correct for the LHS
    RHS=np.zeros(isteps)
    for i in range(1,isteps-1):
        RHS[i]=r*T[i+1]+(1-2*r)*T[i]+r*T[i-1]#this is correct
    RHS[0]=r*T[1]
    RHS[isteps-1]=r*T[isteps-2]
    D=np.matrix(RHS).T#made matrix and made it vertical
    def tridiag(A, B, C, k1=-1, k2=0, k3=1):#tridiagonal matrix
        return np.diag(A, k1) + np.diag(B, k2) + np.diag(C, k3)
    A=[a]*(isteps-1); B = [b]*(isteps); C = [c]*(isteps-1)
    d=tridiag(A, B, C)
    I=npla.inv(d)#inverse matrix... this might be the issue..?
    T1=I*D
    T1=np.array(T1.T)[0]
    T1[0]=0
    T1[-1]=0
    return T1

plt.plot(x[0:isteps],T)

for n in range(nsteps):#run crank for all t and plot for each
    T=CrankNicolson(T)
    Tn.append(T)
    plt.plot(x[0:isteps],T)


plt.show()


# In[12]:


#this is me figuring out how to make a matrix...
nsteps=5
r=nsteps
isteps=10
A=np.zeros((nsteps,nsteps))
A[0:1]=[2]*1

print(A)

a=-r
b=(1+2*r)
c=-r
e=np.zeros((isteps-1,isteps-1))
e[7:9]=[2]
e[0:7]=[0]
e=np.matrix(e).T
e[8:9]=[2]
e[1:8]=[0]
#print(e)
def tridiag(A, B, C, k1=-1, k2=0, k3=1):#tridiagonal matrix
       return np.diag(A, k1) + np.diag(B, k2) + np.diag(C, k3)
A=[a]*(isteps-1); B = [b]*(isteps); C = [c]*(isteps-1)
d=tridiag(A, B, C)
d[0][isteps-1]=r
d[0][isteps-2]=r
d[isteps-1][0]=r
d[isteps-1][1]=r
print(d)


# In[13]:


#this is my 2nd attempt at crank with the 1D heat equation

#variables
L=20
t0=0
tmax=2
x0=0
xmax=L
isteps=4*L
nsteps=20
dt=0.2
dx=xmax/isteps
k=1j
r=(dt)/(4*k*(dx)**2)
sigma=L/12
x=np.zeros(isteps+1)

for i in range(isteps):
    xi=x0+(i+1)*dx
    x[i+1]=xi
    
T=np.zeros(isteps,np.complex)

for i in range(1,isteps):#gives first array of T vs x
    xi=x0+(i)*dx
    Ti=np.exp((-1)*((xi-(L/2))**2)/(sigma)**2)
    T[i]=Ti


#make array of arrays of T
Tn=[]
Tn.append(T)


def CrankNicolson(T):#NEED TO FINISH WRITING THIS put in the inverse matrix
    a=-r
    b=(1+2*r)
    c=-r
    RHS=np.zeros(isteps,np.complex)
    for i in range(1,isteps-1):
        RHS[i]=r*T[i+1]+(1-2*r)*T[i]+r*T[i-1]
    RHS[0]=r*T[1]
    RHS[isteps-1]=r*T[isteps-2]
    D=np.matrix(RHS).T#made matrix for the RHS
    def tridiag(A, B, C, k1=-1, k2=0, k3=1):#tridiagonal matrix
        return np.diag(A, k1) + np.diag(B, k2) + np.diag(C, k3)
    A=[a]*(isteps-1); B = [b]*(isteps); C = [c]*(isteps-1)
    d=tridiag(A, B, C)
    d[0][isteps-1]=a
    d[0][isteps-2]=c
    d[isteps-1][0]=a
    d[isteps-1][1]=c
    I=npla.inv(d)#inverse matrix... this might be the issue..?
    T1=I*D
    T1=np.array(T1.T)[0]
    T1[0]=0
    T1[-1]=0
    return T1

plt.plot(x[0:isteps],np.real(T))# plots initial condition

for n in range(nsteps):#run crank for all t and plot for each plots real parts
    T=CrankNicolson(T)
    Tn.append(T)
    plt.plot(x[0:isteps],np.real(T))
    

plt.show()

for T in Tn:#Plots imaginary parts
    plt.plot(x[0:isteps],np.imag(T))

plt.show()

for T in Tn:#plots overall average
    plt.plot(x[0:isteps],np.abs(T)**2)
    
plt.show()


# In[14]:


#attempt with periodic boundary conditions

#variables
L=20
t0=0
tmax=2
x0=0
xmax=L
isteps=4*L
nsteps=20
dt=0.2
dx=xmax/isteps
k=1j
r=(dt)/(4*k*(dx)**2)
sigma=L/12
x=np.zeros(isteps+1)

for i in range(isteps):
    xi=x0+(i+1)*dx
    x[i+1]=xi
    
T=np.zeros(isteps,np.complex)

for i in range(1,isteps):#gives first array of T vs x
    xi=x0+(i)*dx
    Ti=np.exp(-4*k*xi+(-1)*((xi-(L/2))**2)/(sigma)**2)
    T[i]=Ti


#make array of arrays of T
Tn=[]
Tn.append(T)


def CrankNicolson(T):#NEED TO FINISH WRITING THIS put in the inverse matrix
    a=-r
    b=(1+2*r)
    c=-r
    RHS=np.zeros(isteps,np.complex)
    for i in range(1,isteps-1):
        RHS[i]=r*T[i+1]+(1-2*r)*T[i]+r*T[i-1]
    RHS[0]=r*T[1]+(1-2*r)*T[0]+r*T[-1]
    RHS[isteps-1]=r*T[0]+(1-2*r)*T[-1]+r*T[-2]
    D=np.matrix(RHS).T#made matrix for the RHS
    def tridiag(A, B, C, k1=-1, k2=0, k3=1):#tridiagonal matrix
        return np.diag(A, k1) + np.diag(B, k2) + np.diag(C, k3)
    A=[a]*(isteps-1); B = [b]*(isteps); C = [c]*(isteps-1)
    d=tridiag(A, B, C)
    d[0][isteps-1]=a#added for periodic boundary conditions, extra terms in upper right and lower left
    d[isteps-1][0]=a
    I=npla.inv(d)
    T1=I*D
    T1=np.array(T1.T)[0]
    return T1

plt.plot(x[0:isteps],np.real(T))# plots initial condition

for n in range(nsteps):#run crank for all t and plot for each plots real parts
    T=CrankNicolson(T)
    Tn.append(T)
    plt.plot(x[0:isteps],np.real(T))
    

plt.show() 

for T in Tn:#Plots imaginary parts
    plt.plot(x[0:isteps],np.imag(T))

plt.show()

for T in Tn:#plots overall average
    plt.plot(x[0:isteps],np.abs(T)**2)
    print(np.sum(np.abs(T)**2))#checking for unitarity
plt.show()


# In[15]:


# trying at particle on ring system with periodic boundary cause ring and vector potential
# THIS ONE ISNT WORKING BUT I AM USING IT AS A TEMPLATE FOR THE NEXT ONE

#variables
L=20
Omeg0=1
Omeg1=np.pi
eps=1
t0=0
tmax=2
x0=0
xmax=L
isteps=4*L
nsteps=20
dt=0.2
dx=xmax/isteps
k=1j
r=(dt)/(4*k*(dx)**2)
sigma=L/12
x=np.zeros(isteps+1)
drv=np.zeros(nsteps+1)
t=np.zeros(nsteps+1)

for i in range(isteps):#x array
    xi=x0+(i+1)*dx
    x[i+1]=xi

for n in range(nsteps):#t array
    tn=t0+(n+1)*dt
    t[n+1]=tn

for n in range(nsteps):#drv array, drv is drive force
    tn=t0+(n+1)*dt
    drvn=(eps/Omeg1)*np.cos(Omeg1*tn)
    drv[n]=drvn
    
#A=drv[0]
#A1=drv[1]


print(drv)    

T=np.zeros(isteps,np.complex)

for i in range(1,isteps):#gives first array of T vs x
    xi=x0+(i)*dx
    Ti=np.exp(-4*k*xi+(-1)*((xi-(L/2))**2)/(sigma)**2)
    T[i]=Ti


#make array of arrays of T
Tn=[]
Tn.append(T)





def CrankNicolson(T):#NEED TO FINISH WRITING THIS put in the inverse matrix
    A=drv[n]
    A1=drv[n+1]
    print(A,A1)
    a=((4*dx-4*(dx**2)*k*A1)/(16*(dx**3)))# do i need to loop this for A also? need 2 variables for this in the loop for the time steps
    b=((-8*(dx**2)+dt-(A1)**2+2*(Omeg0**2*np.cos(x[1])))/(16*dt*(dx**2)))
    c=((4*dx-4*(dx**2)*k*A1)/(16*(dx**3)))
    RHS=np.zeros(isteps,np.complex)
    for i in range(1,isteps-1):
        RHS[i]=((4*dx-4*(dx**2)*k*A)/(16*(dx**3)))*T[i+1]+((-8*(dx**2)+dt+(A)**2+2*(Omeg0**2*np.cos(x[1])))/(16*dt*(dx**2)))*T[i]+((4*dx-4*(dx**2)*k*A)/(16*(dx**3)))*T[i-1]
    RHS[0]=((4*dx-4*(dx**2)*k*A)/(16*(dx**3)))*T[1]+((-8*(dx**2)+dt+(A)**2+2*(Omeg0**2*np.cos(x[1])))/(16*dt*(dx**2)))*T[0]+((4*dx-4*(dx**2)*k*A)/(16*(dx**3)))*T[-1]
    RHS[isteps-1]=((4*dx-4*(dx**2)*k*A)/(16*(dx**3)))*T[0]+((-8*(dx**2)+(A)**2+2*(Omeg0**2*np.cos(x[1])))/(16*dt*(dx**2)))*T[-1]+((4*dx-4*(dx**2)*k*A)/(16*(dx**3)))*T[-2]
    D=np.matrix(RHS).T#made matrix for the RHS
    def tridiag(A, B, C, k1=-1, k2=0, k3=1):#tridiagonal matrix
        return np.diag(A, k1) + np.diag(B, k2) + np.diag(C, k3)
    A=[-a]*(isteps-1); B = [b]*(isteps); C = [-c]*(isteps-1)
    d=tridiag(A, B, C)
    d[0][isteps-1]=a#added for periodic boundary conditions, extra terms in upper right and lower left
    d[isteps-1][0]=a
    I=npla.inv(d)
    T1=I*D
    T1=np.array(T1.T)[0]
    return T1

plt.plot(x[0:isteps],np.real(T))# plots initial condition

for n in range(nsteps):#run crank for all t and plot for each plots real parts
    A=drv[n]
    A1=drv[n+1]
    T=CrankNicolson(T)
    Tn.append(T)
    plt.plot(x[0:isteps],np.real(T))
    

plt.show() 

for T in Tn:#Plots imaginary parts
    plt.plot(x[0:isteps],np.imag(T))

plt.show()

for T in Tn:#plots overall average
    plt.plot(x[0:isteps],np.abs(T)**2)
    print(np.sum(np.abs(T)**2))#checking for unitarity
plt.show()


# In[18]:


# trying at particle on ring system with periodic boundary cause ring and vector potential


#variables
L=20
Omeg0=1
Omeg1=np.pi
eps=1
ntild=1
t0=0
tmax=2
x0=0
xmax=L
isteps=4*L
nsteps=20
dt=0.2
dx=xmax/isteps
k=1j
r=(dt)/(4*k*(dx)**2)
sigma=L/12
x=np.zeros(isteps+1)
drv=np.zeros(nsteps+1)
t=np.zeros(nsteps+1)

for i in range(isteps):#x array
    xi=x0+(i+1)*dx
    x[i+1]=xi

for n in range(nsteps):#t array
    tn=t0+(n+1)*dt
    t[n+1]=tn

for n in range(nsteps):#drv array, drv is drive force
    tn=t0+(n+1)*dt
    drvn=np.cos(Omeg1*tn)
    drv[n]=drvn
    
#A=drv[0]
#A1=drv[1]


print(drv)    

T=np.zeros(isteps,np.complex)

for i in range(1,isteps):#gives first array of T vs x
    xi=x0+(i)*dx
    Ti=np.exp(-4*k*xi+(-1)*((xi-(L/2))**2)/(sigma)**2)
    T[i]=Ti


#make array of arrays of T
Tn=[]
Tn.append(T)


def CrankNicolson(T):#NEED TO FINISH WRITING THIS put in the inverse matrix
    A=drv[0]
    A1=drv[1]
    print(A,A1)
    a=((1/(ntild*4*(dx)**2))-((k*eps*(Omeg0/Omeg1)*A1)/(dx)))
    b=((k/(Omeg0*dt))-(1/(ntild*2*(dx)**2))-(t**2+(Omeg0**2)/(Omeg1**2)*ntild*(A1)**2/2)+(ntild*np.cos(x[1])/2))
    c=((1/(ntild*4*(dx)**2))-((k*eps*(Omeg0/Omeg1)*A1)/(dx)))
    RHS=np.zeros(isteps,np.complex)
    for i in range(1,isteps-1):
        RHS[i]=-((1/(ntild*4*(dx)**2))-((k*eps*(Omeg0/Omeg1)*A)/(dx)))*T[i+1]+((k/(Omeg0*dt))-(1/(ntild*2*(dx)**2))+(t**2+(Omeg0**2)/(Omeg1**2)*ntild*(A)**2/2)-(ntild*np.cos(x[1])/2))*T[i]+((1/(ntild*4*(dx)**2))-((k*eps*(Omeg0/Omeg1)*A)/(dx)))*T[i-1]
    RHS[0]=-((1/(ntild*4*(dx)**2))-((k*eps*(Omeg0/Omeg1)*A)/(dx)))*T[1]+((-8*(dx**2)+dt+(A)**2+((k/(Omeg0*dt))-(1/(ntild*2*(dx)**2))+(t**2+(Omeg0**2)/(Omeg1**2)*ntild*(A)**2/2)-(ntild*np.cos(x[1])/2))*T[0]+((1/(ntild*4*(dx)**2))-((k*eps*(Omeg0/Omeg1)*A)/(dx)))*T[-1]
    RHS[isteps-1]=-((1/(ntild*4*(dx)**2))-((k*eps*(Omeg0/Omeg1)*A)/(dx)))*T[0]+((k/(Omeg0*dt))-(1/(ntild*2*(dx)**2))+(t**2+(Omeg0**2)/(Omeg1**2)*ntild*(A)**2/2)-(ntild*np.cos(x[1])/2))*T[-1]+((1/(ntild*4*(dx)**2))-((k*eps*(Omeg0/Omeg1)*A)/(dx)))*T[-2]
    D=np.matrix(RHS).T#made matrix for the RHS ^there is a syntax error here and i dont know why
    def tridiag(A, B, C, k1=-1, k2=0, k3=1):#tridiagonal matrix
        return np.diag(A, k1) + np.diag(B, k2) + np.diag(C, k3)
    A=[a]*(isteps-1); B = [b]*(isteps); C = [-c]*(isteps-1)
    d=tridiag(A, B, C)
    d[0][isteps-1]=a#added for periodic boundary conditions, extra terms in upper right and lower left
    d[isteps-1][0]=a
    I=npla.inv(d)
    T1=I*D
    T1=np.array(T1.T)[0]
    return T1

plt.plot(x[0:isteps],np.real(T))# plots initial condition

for n in range(nsteps):#run crank for all t and plot for each plots real parts
    A=drv[n]
    A1=drv[n+1]
    T=CrankNicolson(T)
    Tn.append(T)
    plt.plot(x[0:isteps],np.real(T))
    

plt.show() 

for T in Tn:#Plots imaginary parts
    plt.plot(x[0:isteps],np.imag(T))

plt.show()

for T in Tn:#plots overall average
    plt.plot(x[0:isteps],np.abs(T)**2)
    print(np.sum(np.abs(T)**2))#checking for unitarity
plt.show()


# In[ ]:




