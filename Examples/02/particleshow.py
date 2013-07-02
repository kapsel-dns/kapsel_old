showRigidCOM=0
type=$constitutive_eq.type
print type
if type == "Navier_Stokes" :
	dx=$constitutive_eq.Navier_Stokes.DX
elif type == "Shear_Navier_Stokes" :
	dx=$constitutive_eq.Shear_Navier_Stokes.DX
elif type == "Shear_Navier_Stokes_Lees_Edwards" :
	dx=$constitutive_eq.Shear_Navier_Stokes_Lees_Edwards.DX
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
elif objType == "rigid":
	Ns=getArray($object_type.rigid.Rigid_spec[])
size_Ns=len(Ns)
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
if objType=="rigid" or objType=="chain":
        if showRigidCOM == 0:
                #species
                for i in range(size_Ns):
                        #chains
                        for m in range(Ns[i][1]):
                                #beads
                                for n in range(Ns[i][0]):
                                        r=$Particles[n_offset+n].R
                                        sphere(r,spat[i%len(spat)])
                                n_offset += Ns[i][0]
        elif showRigidCOM == 1:
                #species
                for i in range(size_Ns):
                        #chains
                        for m in range(Ns[i][1]):
                                r=$RigidParticles[n_offset+m].R
                                sphere(r,spat[i%len(spat)])
                        n_offset+=Ns[i][1]
else:
        #species
        for i in range(size_Ns):
                #particles
                for n in range(Ns[i][0]):
                        r=$Particles[n_offset+n].R
                        sphere(r,spat[i%len(spat)])
                n_offset+=Ns[i][0]
