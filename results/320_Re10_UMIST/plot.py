import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

#x,y=np.mgrid[0:9,0:9]
#print(data[3,4])
#print(data[0:3,0:3])
data = np.loadtxt('u.txt')
u_vec=0.5*( data[0:len(data)-1,:]+data[1:len(data),:] )
dx=1.0/(len(data)-2)
dy=1.0/(len(data[0])-1)
x,y=np.mgrid[0.5*dx:1:dx,0:1+0.5*dy:dy]
fig=plt.figure()
ax=fig.subplots()
cs=ax.contourf(x,y,data[1:len(data)-1,:],levels=100,origin='lower',cmap='nipy_spectral')
cbar=fig.colorbar(cs)
plt.savefig("u.png")
plt.close(fig)

data = np.loadtxt('v.txt')
v_vec=0.5*( data[:,0:len(data[0])-1]+data[:,1:len(data[0])] )
dx=1.0/(len(data)-1)
dy=1.0/(len(data[0])-2)
x,y=np.mgrid[0:1+0.5*dx:dx,0.5*dy:1:dy]          # 0.5*dx:1:dx,0:1+0.5*dy:dy]
fig=plt.figure()
ax=fig.subplots()
cs=ax.contourf(x,y,data[:,1:len(data[0])-1],levels=100,origin='lower',cmap='nipy_spectral')
cbar=fig.colorbar(cs)
plt.savefig("v.png")
plt.close(fig)

data = np.loadtxt('p.txt')
dx=1.0/(len(data)-1)
dy=1.0/(len(data[0])-1)
x,y=np.mgrid[0:1+0.5*dx:dx,0:1+0.5*dy:dy]
fig=plt.figure()
ax=fig.subplots()
cs=ax.contourf(x,y,data,levels=100,origin='lower',cmap='nipy_spectral')
cbar=fig.colorbar(cs)
plt.savefig("p.png")
plt.close(fig)

fig=plt.figure()
ax=fig.subplots()

ax.quiver(x,y,u_vec,v_vec,units='xy',headwidth=1)
plt.savefig("vec.png",dpi=1200)
plt.close(fig)

#x,y=np.

#x = data[:, 0]
#y = data[:, 1]
#print(data)
#print(x)
#print(y)
#print(x+y)
#x,y=np.mgrid[0:1:1,0:7:1]
#fig=plt.figure()
#ax=fig.subplots()
#print(x)
#print(y)
#print(data)
#ax.contourf(x,y,data,levels=7,origin='lower')

