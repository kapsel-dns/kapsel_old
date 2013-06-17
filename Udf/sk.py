import math
from numpy import *
from UDFManager import *
uobj=UDFManager("./output.udf")
Ns=len(uobj.get("object_type.spherical_particle.Particle_spec[]"))
Nc=uobj.get("object_type.spherical_particle.Particle_spec[].Particle_number")
N1=0
for i in range(Ns):
	N1=N1+Nc[i]
print Ns, Nc, N1
nx=uobj.get("mesh.NPX")
ny=uobj.get("mesh.NPY")
nz=uobj.get("mesh.NPZ")
lx=2**nx
ly=2**ny
lz=2**nz
pi2=3.14159265358979*2.0
dia=uobj.get("A")*2.0
print "#",N1,lx,ly,lz,dia
f_sk = open('sk.dat','w')
mxk=100000
k01=arange(0,mxk)
k02=arange(0,mxk)
k03=arange(0,mxk)
n1    = 10
n2    = 10
delk  = pi2/float(50)
nnk = 0
iz = 0
iy = 0
for ix in range(1, n1):
	nnk = nnk + 1
	k01[nnk] = ix
	k02[nnk] = iy
	k03[nnk] = iz
for iy in range(1, n1):
	for ix in range(-n1, n1):
		nnk = nnk + 1
		k01[nnk] = ix
		k02[nnk] = iy
		k03[nnk] = iz
for iz in range(1, n2):
	for iy in range(-n1, n1):
		for ix in range(-n1, n1):
			nnk = nnk + 1
			k01[nnk] = ix
			k02[nnk] = iy
			k03[nnk] = iz
print nnk,mxk
rr=zeros(1001,float64)
hist=zeros(1001,int32)
for kkk in range(1, nnk):
	k1      = float(k01[kkk]) * pi2 / float(lx)
	k2      = float(k02[kkk]) * pi2 / float(ly)
	k3      = float(k03[kkk]) * pi2 / float(lz)
	kl      = sqrt(k1*k1+k2*k2+k3*k3)
	bin     = int(kl/delk)
	rr[bin]   = rr[bin] + kl*dia
	hist[bin] = hist[bin] + 1
for ij in range(0, 1000):
	if hist[ij]>0:
		rr[ij]=rr[ij]/float(hist[ij])
		print ij,rr[ij],hist[ij]
nt=uobj.totalRecord()
#print nt
suma=zeros(1001,float64)
for n in range(0,nt):
	uobj.jump(n)
	x = uobj.get("Particles[].R.x")
	y = uobj.get("Particles[].R.y")
	z = uobj.get("Particles[].R.z")
	for kkk in range(1,nnk):
		k1      = float(k01[kkk]) * pi2 / float(lx)
		k2      = float(k02[kkk]) * pi2 / float(ly)
		k3      = float(k03[kkk]) * pi2 / float(lz)
		kl      = sqrt(k1*k1+k2*k2+k3*k3)
		bin     = int(kl/delk)
            	sa1     = 0.0
            	sa2     = 0.0
      	  	for ii in range(N1):
               		r1  = k1*x[ii]+k2*y[ii]+k3*z[ii]
               		sa1 = sa1 + cos(r1)
               		sa2 = sa2 + sin(r1)
            	sa = (sa1*sa1+sa2*sa2)/float(N1)
            	suma[bin] = suma[bin] + sa
	time=uobj.get("t")
	print n,time
for ij in range(0,1000):
	if hist[ij]>0:
		suma[ij]=suma[ij]/float(hist[ij]*(nt-1))
		f_sk.write(str(rr[ij])+'  '+str(suma[ij])+'  '+str(n)+'  '+str(hist[ij])+'\n')
f_sk.close()
