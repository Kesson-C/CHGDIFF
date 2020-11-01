import numpy as np
from vaspy.electro import ChgCar
from scipy.interpolate import interp2d
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
fb = open("/Users/mne081069/Desktop/data/chg/101Tib/CHGCAR",'r',encoding="utf-8")
fs = open("/Users/mne081069/Desktop/data/chg/101Tis/CHGCAR",'r',encoding="utf-8")
foutd = open("CHGCAR_N",'w')
line = fb.readlines()
lines = fs.readlines()
lc = float(line[1].strip('\n'))
lcs = float(lines[1].strip('\n')) 
v1 = np.asarray(list(map(float,line[2].split())))
v2 = np.asarray(list(map(float,line[3].split())))
v3 = np.asarray(list(map(float,line[4].split())))
v1s = np.asarray(list(map(float,lines[2].split())))
v2s = np.asarray(list(map(float,lines[3].split())))
v3s = np.asarray(list(map(float,lines[4].split())))
ctmp = list(line[6].split())
coef = [0]*len(ctmp)
ctmps = list(lines[6].split())
coefs = [0]*len(ctmps)
for i in range(len(ctmp)):
    coef[i]=int(ctmp[i])

for i in range(len(ctmps)):
    coefs[i]=int(ctmps[i])
    
    
pos = np.zeros([np.sum(coef),3])
poss = np.zeros([np.sum(coefs),3])
for i in range(len(pos)):
    str1 = line[i+8].split()
    pos[i]=list(map(float,str1))
    

for i in range(len(poss)):
    str1 = lines[i+8].split()
    poss[i]=list(map(float,str1))

nx = int(line[9+np.sum(coef)].strip('\n').split()[0])
ny = int(line[9+np.sum(coef)].strip('\n').split()[1])
nz = int(line[9+np.sum(coef)].strip('\n').split()[2])
nxs = int(lines[9+np.sum(coefs)].strip('\n').split()[0])
nys = int(lines[9+np.sum(coefs)].strip('\n').split()[1])
nzs = int(lines[9+np.sum(coefs)].strip('\n').split()[2])

tmp1 = np.zeros(nx*ny*nz)
k=0
for i in range(10+np.sum(coef),10+np.sum(coef)+int(np.ceil(nx*ny*nz/5))):
    temp = np.asarray(list(map(float,line[i].split())))
    for j in range(len(temp)):
        tmp1[k]=temp[j]
        k=k+1
 
tmp1=tmp1/(nx*ny*nz)
    
tmp2 = np.zeros(nxs*nys*nzs)
k=0
for i in range(10+np.sum(coefs),10+np.sum(coefs)+int(np.ceil(nxs*nys*nzs/5))):
    temp = np.asarray(list(map(float,lines[i].split())))
    for j in range(len(temp)):
        tmp2[k]=temp[j]
        k=k+1

tmp2=tmp2/(nx*ny*nz)
tmp= tmp2-tmp1
ct=0
for i in range(0,10+np.sum(coefs)):   
    foutd.write(lines[i])

for i in range(10+np.sum(coefs),10+np.sum(coefs)+int(np.ceil(nx*ny*nz/5))):   
        lines[i]=' '
        k=0
        while(k<5 and ct<len(tmp)):
            lines[i] =lines[i] + ''.join('%10.8e'%tmp[ct]) + ' '
            k=k+1
            ct=ct+1
        lines[i]=lines[i]+'\n'
        foutd.write(lines[i])
foutd.close()

c = ChgCar('CHGCAR_N')
z=c.elf_data[:,:,253]
y=np.arange(0,len(z))
x=np.arange(0,len(z[0]))
intef=interp2d(x,y,z)
nex = np.linspace(0, nx, 600)
ney = np.linspace(0, ny, 600)
newz=intef(ney,nex)
#nmap=np.zeros([1200,1200])
#nmap[0:600,0:600]=newz
#nmap[600:1200,0:600]=newz
#nmap[0:600,600:1200]=newz
#nmap[600:1200,600:1200]=newz
nmap=np.transpose(newz)
ms = np.s_[0:nx:600j, 0:ny:600j]
newmx, newmy = np.mgrid[ms]  
#plot 2d contour map
color=plt.cm.rainbow
fig2d_1= plt.figure()
ax1 = fig2d_1.add_subplot(1, 1, 1)
extent = [np.min(nex), np.max(nex), np.min(ney), np.max(ney)]
img = ax1.imshow(nmap, interpolation='spline36', origin='lower', extent=extent,cmap=color,vmin=0,vmax=2e-4)#,vmin=-0.0005,vmax=0.00065
divi=make_axes_locatable(ax1)
pt1 = nmap[np.int(np.round((poss[96,0])*600)),np.int(np.round((poss[96,1])*600))]
pt2 = nmap[np.int(np.round((poss[97,0])*600)),np.int(np.round((poss[97,1])*600))]
pt3 = nmap[np.int(np.round((poss[98,0])*600)),np.int(np.round((poss[98,1])*600))]
pt4 = nmap[np.int(np.round((poss[99,0])*600)),np.int(np.round((poss[99,1])*600))]
plt.scatter((poss[96,0])*120,(poss[96,1])*80,color='b',label=np.str(pt1))
plt.scatter((poss[97,0])*120,poss[97,1]*80,color='g',label=np.str(pt2))
plt.scatter((poss[98,0])*120,(poss[98,1])*80,color='k',label=np.str(pt3))
plt.scatter((poss[99,0])*120,poss[99,1]*80,color='w',label=np.str(pt4))
plt.legend()
#plt.scatter((poss[68,0]-0)*80,poss[68,1]*80,color='g')
#plt.scatter((poss[69,0]-0)*80,poss[69,1]*80,color='g')
#plt.scatter((poss[70,0]-0)*80,poss[70,1]*80,color='g')
#plt.scatter((poss[71,0]-0)*80,poss[71,1]*80,color='g')
#plt.scatter(poss[8,0]*80,poss[8,1]*80,color='k')
#plt.scatter(poss[9,0]*80,poss[9,1]*80,color='k')
#plt.scatter(poss[10,0]*80,poss[10,1]*80,color='k')
#plt.scatter(poss[11,0]*80,poss[11,1]*80,color='k')
#plt.scatter(poss[0,0]*80,poss[0,1]*80,color='g')
#plt.scatter(poss[1,0]*80,poss[1,1]*80,color='g')
#plt.scatter(poss[2,0]*80,poss[2,1]*80,color='g')
#plt.scatter(poss[3,0]*80,poss[3,1]*80,color='g')
#xxx=np.where(nmap>np.max(z)-0.00001)
#plt.scatter(xxx[1]/600*80,xxx[0]/600*80,color='y')
#plt.scatter(poss[96,0]*80,poss[96,1]*80,color='b')
#plt.scatter(poss[72,0]*80,poss[72,1]*80)
cax= divi.append_axes("right",size="5%",pad=0.05)
cbar = ax1.figure.colorbar(img,cax=cax)
fig2d_1.savefig('surface2d.png', dpi=500)
print(np.max(z))
print(np.min(z))
##3d plot
for i in range(30,40,10):
    for j in range(40,60,10):
        fig3d = plt.figure(figsize=(12, 8))
        ax3d = fig3d.add_subplot(111, projection='3d')
        ax3d.view_init(i, j)
        ax3d.set_xlim3d(np.min(nex), np.max(nex))
        ax3d.set_ylim3d(np.min(ney), np.max(ney))
        ax3d.set_zlim3d(-0.006, 0.0001)
        ax3d.plot_surface(newmx, newmy, nmap, cmap=color,vmin=0,vmax=2e-4)#, vmin=-0.0005,vmax=0.00065       
        ax3d.contourf(newmx, newmy, nmap, 100 , cmap=color,extent=extent, zdir='z', offset=-0.0061,vmin=0,vmax=2e-4)#,vmin=-0.0005,vmax=0.00065
        plt.colorbar(mappable=img,cmap=color,boundaries=np.linspace(0,2e-4,100))#,boundaries=np.linspace(-0.0005,0.00065,100)
        fig3d.savefig('surface3d'+str(i)+'_'+str(j)+'.png', dpi=500)
        #plt.close(fig3d)
#
#for angle in range(95, 180, 3):
#    ax3d.set_zlabel("Angle: " + str(angle))
#    ax3d.view_init(10, angle)
#    filename = "./" + str(angle) + ".png"
#    fig3d.savefig(filename,dpi=500)
#    print("Save " + filename + " finish")


fb.close()
fs.close()