from numpy import*
type=$constitutive_eq.type
print type
if type == "Navier_Stokes" :
	dx=$constitutive_eq.Navier_Stokes.DX
elif type == "Shear_Navier_Stokes" :
	dx=$constitutive_eq.Shear_Navier_Stokes.DX
elif type == "Electrolyte" :
	dx=$constitutive_eq.Electrolyte.DX
NX=2**$mesh.NPX
NY=2**$mesh.NPY
NZ=2**$mesh.NPZ
LX=dx*NX
LY=dx*NY
LZ=dx*NZ
objType=$object_type.type
if objType=="spherical_particle":
	Ns=getArray($object_type.spherical_particle.Particle_spec[])
elif objType=="chain":
	Ns=getArray($object_type.chain.Chain_spec[])
size_Ns=len(Ns)
print size_Ns
RAD=($A*dx)*1.
for i in range(size_Ns):
	print Ns[i][0],
print LX,LY,LZ
cells=[[0,0,0],[0,LY,0],[LX,LY,0],[LX,0,0],[0,0,0],[0,0,LZ],[LX,0,LZ],[LX,LY,LZ],[0,LY,LZ],[0,0,LZ]]
polyline(cells,1)
line([ 0,LY, 0],[ 0,LY,LZ],1)
line([LX, 0, 0],[LX, 0,LZ],1)
line([LX,LY, 0],[LX,LY,LZ],1)
#cells2=[[-0.5,-0.5,LZ/2],[-0.5,LY-0.5,LZ/2],[LX-0.5,LY-0.5,LZ/2],[LX-0.5,-0.5,LZ/2]]
cells2=[[0,0,LZ/2],[0,LY,LZ/2],[LX,LY,LZ/2],[LX,0,LZ/2]]
#polyline(cells2,1)
spat=[
	[1.0,0,0,1.0,RAD],
	[0,1.0,0,1.0,RAD],
	[0,0,1.0,1.0,RAD],
	[1,1,0.0,1.0,RAD],
	[0,1,1.0,1.0,RAD],
	[1,0,1.0,1.0,RAD],
	[1,1,1.0,1.0,RAD]
	]
n_offset = 0
vmax=0
for i in range(size_Ns):
	for n in range(Ns[i][0]):
		v1=$Particles[n_offset+n].v
		if (v1[0]*v1[0]> vmax):
			vmax = v1[0]*v1[0] 
n_offset =0
for i in range(size_Ns):
	for n in range(Ns[i][0]):
		r=$Particles[n_offset+n].R
		v0=$Particles[n_offset+n].v
		if v0[0]< 0:
			im = 1
			cv=1
			if(vmax>0):
				cv = v0[0]*v0[0]/vmax
		else:
			im = -1
			cv=1
			if(vmax>0):
				cv = -v0[0]*v0[0]/vmax
		if im == 1:
			sphere(r,[1,cv,0,1.0,RAD])
		elif im == -1:
			sphere(r,[1,0,-cv,1.0,RAD])
			#elif im == 0:
			#	sphere(r,[0,1,fabs(cv),1.0,RAD])
			#else:
			#	sphere(r,[0,1-fabs(cv),1,1.0,RAD])
	else:
		h=1
			#sphere(r,spat[i%len(spat)])
	n_offset += Ns[i][0]
#n_offset = 0sphere(r,spat[i%len(spat)])
#for i in range(size_Ns):
#	for n in range(Ns[i][0]):
#		r=$Particles[n_offset+n].R
#		if(n<4):
#			r[2] = LZ/2+0.2
#			disk(r,[0.6,0.6,0.6,0.8,RAD,0,0,1])
#n_offset += Ns[i][0]
#sphere([0,0,0],[1,0,0,1,0])
