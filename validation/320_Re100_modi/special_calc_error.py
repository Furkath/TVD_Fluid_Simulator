import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import signal 

data = np.loadtxt('u.txt')
u_com=0.5*( data[1:len(data)-1,0:len(data[0])-1] + data[1:len(data)-1,1:len(data[0])] )
u_ref=np.loadtxt('u_ref.txt')
u_err=(u_com-u_ref)#/u_ref
Dx=1.0/(len(data)-2)

data = np.loadtxt('v.txt')
v_com=0.5*( data[0:len(data)-1,1:len(data[0])-1] + data[1:len(data),1:len(data[0])-1] )
v_ref=np.loadtxt('v_ref.txt')
v_err=(v_com-v_ref)#/v_ref
Dy=1.0/(len(data[0])-2)

x,y=np.mgrid[0.5*Dx:1:Dx,0.5*Dy:1:Dy]

aa=3
kernel=np.ones((aa,aa))/(aa**2)
u_modi=signal.convolve2d(u_err,kernel,mode='same')
v_modi=signal.convolve2d(v_err,kernel,mode='same')

fig=plt.figure()
ax=fig.subplots()
cs=ax.contourf(x,y,u_modi*0.24647426205/np.max(u_modi),origin='lower',cmap='nipy_spectral',levels=60)#np.linspace(-30,30,100),extend='both')
cbar=fig.colorbar(cs)
plt.savefig("u_error.png")
plt.close(fig)

fig=plt.figure()
ax=fig.subplots()
cs=ax.contourf(x,y,v_modi*0.1796566369/np.max(v_modi),origin='lower',cmap='nipy_spectral',levels=60)#np.linspace(-30,30,100),extend='both')
cbar=fig.colorbar(cs)
plt.savefig("v_error.png")
plt.close(fig)

print(np.nanmean(abs(u_err)))
print(np.nanmean(abs(v_err)))

print(np.max(abs(u_err)))
print(np.max(v_err))
print(np.min(v_err))

