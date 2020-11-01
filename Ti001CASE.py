import numpy as np
from vaspy.electro import ChgCar
from scipy.interpolate import interp2d
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
fb = open("/Users/mne081069/Desktop/data/chg/001Tib/CHGCAR",'r',encoding="utf-8")
fs = open("/Users/mne081069/Desktop/data/chg/001Tis/CHGCAR",'r',encoding="utf-8")
ff = open("/Users/mne081069/Desktop/data/chg/001F/CHGCAR",'r',encoding="utf-8")
foutd = open("CHGCAR_N",'w')
line = fb.readlines()
lines = fs.readlines()
linef= ff.readlines()
lc = float(line[1].strip('\n'))
lcs = float(lines[1].strip('\n')) 
lcf = float(linef[1].strip('\n')) 
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
ctmpf = list(linef[6].split())
coeff = [0]*len(ctmps)
for i in range(len(ctmp)):
    coef[i]=int(ctmp[i])

for i in range(len(ctmps)):
    coefs[i]=int(ctmps[i])

for i in range(len(ctmpf)):
    coeff[i]=int(ctmpf[i])  
    
pos = np.zeros([np.sum(coef),3])
poss = np.zeros([np.sum(coefs),3])
posf = np.zeros([np.sum(coeff),3])
for i in range(len(pos)):
    str1 = line[i+8].split()
    pos[i]=list(map(float,str1))
    

for i in range(len(poss)):
    str1 = lines[i+8].split()
    poss[i]=list(map(float,str1))

for i in range(len(posf)):
    str1 = linef[i+8].split()
    posf[i]=list(map(float,str1))
    
nx = int(line[9+np.sum(coef)].strip('\n').split()[0])
ny = int(line[9+np.sum(coef)].strip('\n').split()[1])
nz = int(line[9+np.sum(coef)].strip('\n').split()[2])
nxs = int(lines[9+np.sum(coefs)].strip('\n').split()[0])
nys = int(lines[9+np.sum(coefs)].strip('\n').split()[1])
nzs = int(lines[9+np.sum(coefs)].strip('\n').split()[2])
nxf = int(linef[9+np.sum(coeff)].strip('\n').split()[0])
nyf = int(linef[9+np.sum(coeff)].strip('\n').split()[1])
nzf = int(linef[9+np.sum(coeff)].strip('\n').split()[2])
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

tmp3 = np.zeros(nxf*nyf*nzf)
k=0
for i in range(10+np.sum(coeff),10+np.sum(coeff)+int(np.ceil(nxf*nyf*nzf/5))):
    temp = np.asarray(list(map(float,linef[i].split())))
    for j in range(len(temp)):
        tmp3[k]=temp[j]
        k=k+1

tmp3=tmp3/(nx*ny*nz)
tmp= tmp2-tmp1
tmp=np.array(tmp).reshape((nx, ny, nz), order='F')
tmp[:,:,0:225]=0
tmp[:,:,250:287]=0
tmp[tmp<0]/30
for i in range(8):
    cen=np.round([poss[i,0]*80,poss[i,1]*80,poss[i,2]*288])
    r=14
    for i3 in range(int(cen[2])-r,int(cen[2])+r):
        for i2 in range(int(cen[1])-r,int(cen[1])+r):
            for i1 in range(int(cen[0])-r,int(cen[0])+r):
                cond1 = 1
                cond2 = 1
                if i1>=cen[0]:
                    cond1= tmp[np.mod(i1+1,80),np.mod(i2,80),i3]<=tmp[np.mod(i1,80),np.mod(i2,80),i3]
                elif i1<cen[0]:                        
                    cond1= tmp[np.mod(i1-1,80),np.mod(i2,80),i3]<=tmp[np.mod(i1,80),np.mod(i2,80),i3]
                if i2>=cen[1]:                     
                    cond2= tmp[np.mod(i1,80),np.mod(i2+1,80),i3]<=tmp[np.mod(i1,80),np.mod(i2,80),i3]
                elif i2<cen[1]:                     
                    cond2= tmp[np.mod(i1,80),np.mod(i2-1,80),i3]<=tmp[np.mod(i1,80),np.mod(i2,80),i3]
                if cond1 and cond2:
                    tmp[np.mod(i1,80),np.mod(i2,80),i3]=tmp[np.mod(i1,80),np.mod(i2,80),i3]/30


tmp=np.array(tmp).reshape((nx*ny*nz), order='F')
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
z=c.elf_data[:,:,235]
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
#plot 2d contour map
color=plt.cm.rainbow
fig2d_1= plt.figure()
ax1 = fig2d_1.add_subplot(1, 1, 1)
extent = [np.min(nex), np.max(nex), np.min(ney), np.max(ney)]
img = ax1.imshow(nmap, interpolation='spline36', origin='lower', extent=extent,cmap=color,vmin=0,vmax=2e-4)#,vmin=-0.0005,vmax=0.00065
divi=make_axes_locatable(ax1)
plt.scatter(poss[4,0]*80,poss[4,1]*80,color='b')
plt.scatter(poss[5,0]*80,poss[5,1]*80,color='g')
plt.scatter(poss[6,0]*80,poss[6,1]*80,color='k')
plt.scatter(poss[7,0]*80,poss[7,1]*80,color='w')
#plt.scatter((poss[68,0]-0)*80,poss[68,1]*80,color='g')
#plt.scatter((poss[69,0]-0)*80,poss[69,1]*80,color='g')
#plt.scatter((poss[70,0]-0)*80,poss[70,1]*80,color='g')
#plt.scatter((poss[71,0]-0)*80,poss[71,1]*80,color='g')
#plt.scatter(poss[8,0]*80,poss[8,1]*80,color='k')
#plt.scatter(poss[9,0]*80,poss[9,1]*80,color='k')
#plt.scatter(poss[10,0]*80,poss[10,1]*80,color='k')
#plt.scatter(poss[11,0]*80,poss[11,1]*80,color='k')
plt.scatter(poss[0,0]*80,poss[0,1]*80,color='g')
plt.scatter(poss[1,0]*80,poss[1,1]*80,color='g')
plt.scatter(poss[2,0]*80,poss[2,1]*80,color='g')
plt.scatter(poss[3,0]*80,poss[3,1]*80,color='g')
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
#for i in range(30,40,10):
#    for j in range(40,60,10):
#        fig3d = plt.figure(figsize=(12, 8))
#        ax3d = fig3d.add_subplot(111, projection='3d')
#        ax3d.view_init(40, 40)
#        ax3d.set_xlim3d(np.min(nex), np.max(nex))
#        ax3d.set_ylim3d(np.min(ney), np.max(ney))
#        ax3d.set_zlim3d(-0.0001, 0.0001)
#        ax3d.plot_surface(newmx, newmy, nmap, cmap=color,vmin=0,vmax=2e-4)#, vmin=-0.0005,vmax=0.00065       
#        ax3d.contourf(newmx, newmy, nmap, 100 , cmap=color,extent=extent, zdir='z', offset=-0.00015,vmin=0,vmax=2e-4)#,vmin=-0.0005,vmax=0.00065
#        plt.colorbar(mappable=img,cmap=color,boundaries=np.linspace(0,2e-4,100))#,boundaries=np.linspace(-0.0005,0.00065,100)
#        fig3d.savefig('surface3d'+str(i)+'_'+str(j)+'.png', dpi=500)
#        plt.close(fig3d)
##
#for angle in range(95, 180, 3):
#    ax3d.set_zlabel("Angle: " + str(angle))
#    ax3d.view_init(10, angle)
#    filename = "./" + str(angle) + ".png"
#    fig3d.savefig(filename,dpi=500)
#    print("Save " + filename + " finish")


fb.close()
fs.close()