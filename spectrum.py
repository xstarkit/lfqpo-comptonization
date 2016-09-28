#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

data=np.genfromtxt('mcomp000000.dat')
datax=data[:,0]
ne=50
nangle=10
nphi=4
f=np.zeros((ne,nangle,nphi))
for i in range(ne):
	id=0
        for j in range(nangle):
                for k in range(nphi):
                        f[i][j][k]=data[i,id+1]   #  don't include e
                        id=id+1

datay=datax*f[:,0,0]
#datay=f[:,0,0]
plt.plot(datax, datay, '-', color = 'red', linewidth = 2, label='edge-on' )

datay=datax*f[:,4,0]
#datay=f[:,4,0]
plt.plot(datax, datay, '-', color = 'green', linewidth = 2, label='middle' )

datay=datax*f[:,9,0]
#datay=f[:,9,0]
plt.plot(datax, datay, '-', color = 'blue', linewidth = 2, label='face-on'  )

#===================================================================================

#-----------------------------------------------------------------------------------
plt.xscale('log')
plt.yscale('log')

axis_font={'fontname':'Arial', 'size':'20'}
plt.xlabel(r"E[KeV]", **axis_font)
plt.ylabel(r"$\rm E\cdot F_{\rm E}$ [$\rm KeV/\rm cm^{2}/\rm s$]", **axis_font)
plt.legend(loc='upper left')
plt.minorticks_on()
#plt.tick_params(axis='both', which='major', labelsize=10)
#plt.tick_params(axis='both', which='minor', labelsize=8)
plt.tick_params(which='both',width=2)
plt.tick_params(which='major',length=5)
plt.tick_params(which='minor',length=3)


plt.savefig("spectrum.eps")

plt.show()

