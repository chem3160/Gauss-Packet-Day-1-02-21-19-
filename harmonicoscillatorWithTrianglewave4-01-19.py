from numpy.polynomial.hermite import *
import numpy as np
import math
from matplotlib import pyplot as plt

fig = plt.figure()
ax = plt.axes(xlim= (-10, 10), ylim = (-1,1))
line, =ax.plot([], [], lw=2)


hbar=1
m=1
pi=np.pi
K=1
ci=0+1j
r = np.linspace(-10, 10, 200)



def init():
    line.set_data([], [])
    return line,
 
    
def HO_Time (K, m, n, t):
    ci = 0 + 1j
    En = HO_En(K, m, n)
    return np.exp(ci * En * t/ hbar)
def HO_En (K, m, n):
    return hbar * np.sqrt(K/m) * (n + 1/2 )

def HO_Func(K, m,  n, r):
    
  w = np.sqrt(K/m)
  psi = []
  herm_coeff = []

  for i in range(n):
      herm_coeff.append(0)
  herm_coeff.append(1)

  for x in r:
    psi.append(math.exp(-m*w*x**2/(2*hbar)) * hermval((m*w/hbar)**0.5 * x, herm_coeff))
  # normalization factor for the wavefunction:
  psi = np.multiply(psi, 1 / (math.pow(2, n) * math.factorial(n))**0.5 * (m*w/(pi*hbar))**0.25)

      
  return psi

x = np.linspace(-5,5,200)
psi = HO_Func(1, 1, 2, r)

plt.plot(x, psi)

def animate(i):
    line.set_ydata( HO_Func(K,m,n, r)* HO_Time(K,m,n,i/10) )
    line.set_xdata( r )
    return line, 

def Triangle_Wave(xaxis):
    tw = np.zeros(len(xaxis))
    for i in range (0, len(xaxis)):
        xval = xaxis[i]
        if xval <= 2:
            tw[i] = 0
        elif xval < 3:
            tw[i] = xval - 2
        elif xval < 4:
            tw[i] = -xval + 4
        else: 
            tw[i] = 0
    return tw
    
def fourier_analysis (tw, n, r):
    psi = HO_Func(K,m,n,r)
    psistar = np.conj(psi)
    integrand = psistar * tw
    ####rectangle rule width is called w 
    w = r[1] - r[0]
    ###rectangle rule sum will be called rsum
    rsum = 0.+0j
    for i in range(0,len(tw)):
        rsum = rsum +integrand[i] * w
        
    return rsum
    


#anim = animation.FuncAnimation( fig, animate, init_func=init,
  #                             frames= 10000, interval=20, blit=True)

trianglewave = Triangle_Wave(r)


narray = np.linspace(-10,10,21)
c_array = np.zeros(len(narray), dtype=complex) 
psi_exp= np.zeros(len(trianglewave),dtype=complex)
    
    
for i in range(0,len(narray)):
    c_array[i] = fourier_analysis(trianglewave,narray[i], r)
    psi_exp = psi_exp + c_array[i]*HO_Func(K,m,narray[i],r)

trianglewave= Triangle_Wave(r)
psi_1_times_tw = trianglewave * HO_Func(K,m,n,r)



plt.plot(r, trianglewave)
plt.show()
