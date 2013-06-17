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
shear_type=$constitutive_eq.Shear_Navier_Stokes.External_field.type
if shear_type == "AC":
	shear_rate_v=0.25*LY*$constitutive_eq.Shear_Navier_Stokes.External_field.AC.Shear_rate
elif shear_type == "DC":
	shear_rate_v=0.25*LY*$constitutive_eq.Shear_Navier_Stokes.External_field.DC.Shear_rate
if objType=="spherical_particle":
	Ns=getArray($object_type.spherical_particle.Particle_spec[])
elif objType=="chain":
	Ns=getArray($object_type.chain.Chain_spec[])
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
for i in range(size_Ns):
	for n in range(Ns[i][0]):
		r=$Particles[n_offset+n].R
		sphere(r,spat[i%len(spat)])
	n_offset += Ns[i][0]
#n_offset = 0
#for i in range(size_Ns):
#	for n in range(Ns[i][0]):
#		r=$Particles[n_offset+n].R
#		if(n<4):
#			r[2] = LZ/2+0.2
#			disk(r,[0.6,0.6,0.6,0.8,RAD,0,0,1])
#n_offset += Ns[i][0]
#sphere([0,0,0],[1,0,0,1,0])
################################################################
# Show Field
################################################################
from numpy import *
from array import *
from math import *
from string import *
import os
import colorsys
SHOW_VELOCITY = 1
SHOW_VEC_VOLUME = 0
### Pressure
SHOW_PRESSURE = 0
SHOW_VOLUME = 0
### Velocity
coord_skip = 1
color_start = 0.66666666666
z_factor = 1
arrow_factor = 4
### Common
msize=NX*NY*NZ
coN = [NX,NY,NZ]
#slice = [LX/2,LY/2,LZ/2]
slice = [LX/2,LY/2,LZ/2]
view_axes = 2
ii = [[1,0,0],
      [0,1,0],
      [0,0,1]]
vn = ii[view_axes]
def setAVS():
	avs_ts = str(int($t))
	avs_file_path = join(os.path.join(os.path.split(_udf_.udfFile())[:-1])) + "/" + $output.ON.Out_dir + "/avs"
	avs_file_name = join([$output.ON.Out_name,"_",avs_ts,".dat"],"")
	print os.path.join(avs_file_path,avs_file_name)
	return os.path.join(avs_file_path,avs_file_name)
def setData():
	global	u
	u = array('f')
	fin = open(setAVS(), "rb")	
#############################################
	u.fromfile(fin,msize*7)
#############################################
	fin.close()
def showVelocity():
	minmax = [0.,0.]
	offset=0
	for k in range(NZ):
		for j in range(NY):
			for i in range(NX):
				iijk = [i,j,k]
				if((((i%coord_skip)==0)
				    and ((j%coord_skip)==0)
				    and ((k%coord_skip)==0))
				   and (iijk[view_axes]==slice[view_axes])):
					[ux,uy,uz]=[u[offset],u[offset+msize],u[offset+msize*2]]
					usize = dot([ux,uy,uz],[ux,uy,uz])
					if usize > shear_rate_v*shear_rate_v:
						usize = shear_rate_v*shear_rate_v
					if(i==0 and j==0 and k==0):
						minmax = [usize,usize]
					if(minmax[0] > usize):
						minmax[0] = usize
					if(minmax[1] < usize):
						minmax[1] = usize
				offset+=1
	if(minmax[1]):
		minmax = [math.sqrt(minmax[0]),math.sqrt(minmax[1])]
		minmax = [minmax[0],minmax[1]]
		norm = 1/minmax[1]
	else :
		minmax = [math.sqrt(minmax[0]),math.sqrt(minmax[1])]
		norm = 0.0
	offset=0
	for k in range(NZ):
		for j in range(NY):
			for i in range(NX):
				iijk = [i,j,k]
				if((((i%coord_skip)==0)
				    and ((j%coord_skip)==0)
				    and ((k%coord_skip)==0))
				   and (iijk[view_axes]==slice[view_axes])):
					[ux,uy,uz]=[u[offset],u[offset+msize],u[offset+msize*2]]
					usize = dot([ux,uy,uz],[ux,uy,uz])
					if usize > shear_rate_v*shear_rate_v:
						usize = shear_rate_v*shear_rate_v
					if(norm):
						std_usize = norm*math.sqrt(usize)
						hue = color_start*(1 - std_usize)
						if(hue < (color_start - 1)):
							print 'hue over'
						while (hue < 0):
							hue += 1
						li_rgb = colorsys.hsv_to_rgb(hue,1.0,1.0)
						scale = 40.
						dmyvn=[1,1,1]
						dmyvn[view_axes]=z_factor
						ux *= scale*dmyvn[0]
						uy *= scale*dmyvn[1]
						uz *= scale*dmyvn[2]
						l_length = std_usize*coord_skip*arrow_factor
						a_length = l_length/3.*2.
						a_radius = a_length/3.
						#arrow([i,j,k],
						#      [i+ux,j+uy,k+uz],
						#      [li_rgb[0],li_rgb[1],li_rgb[2],1,a_radius,a_length,l_length]
						#      )
						cone ([i,j,k],
						      [i+ux,j+uy,k+uz],
						      [li_rgb[0],li_rgb[1],li_rgb[2],1.0,0.2]
						      )
				offset+=1
def showPressure():
	field = meshfield("regular",[[0,LX],[0,LY],[0,LZ]],[NX,NY,NZ])
	offset=0
	for k in range(NZ):
		for j in range(NY):
			for i in range(NX):
				[ux,uy,uz,phi,s_ch,rho,e_p]=[u[offset],u[offset+msize],u[offset+msize*2],u[offset+msize*3],u[offset+msize*4],u[offset+msize*5],u[offset+msize*6]]
#				[ux,uy,uz,phi,conc]=[u[offset],u[offset+msize],u[offset+msize*2],u[offset+msize*3],u[offset+msize*4]]
				field.set([i,j,k],e_p)
				field.setSubValue([i,j,k],e_p)
				offset+=1
	if(SHOW_PRESSURE):
		wid = 1.0
		dmy = [slice[0]+wid*vn[0],
		       slice[1]+wid*vn[1],
		       slice[2]+wid*vn[2],
		       slice[0],
		       slice[1],
		       slice[2],
		       slice[0]-wid*vn[0],
		       slice[1]-wid*vn[1],
		       slice[2]-wid*vn[2]
		       ]
#		field.cplane(dmy[0],dmy[1],dmy[2],vn[0],vn[1],vn[2])
		field.cplane(dmy[3],dmy[4],dmy[5],vn[0],vn[1],vn[2])
#		field.cplane(dmy[6],dmy[7],dmy[8],vn[0],vn[1],vn[2])
		field.ccolor([0,0,1,1,-240])
		field.draw(frame = 1,subdivision = 4)
#		field.draw(subdivision = 4)
		field.delete_clevel()
		field.delete_cplane()
	if(SHOW_VOLUME):
		field.drawVolume([0,0,1,1,0,1,0,1],intensity=2,slope_factor=4)
		field.draw(frame = 1,subdivision = 4)
def showVecVolume():
	field = meshfield("regular",[[0,LX],[0,LY],[0,LZ]],[NX,NY,NZ])
	offset=0
	for k in range(NZ):
		for j in range(NY):
			for i in range(NX):
				[ux,uy,uz,phi,pressure]=[u[offset],u[offset+msize],u[offset+msize*2],u[offset+msize*3],u[offset+msize*4]]
				usize = dot([ux,uy,uz],[ux,uy,uz])
				msusize = math.sqrt(usize)
				field.set([i,j,k],msusize)
				field.setSubValue([i,j,k],msusize)
				offset+=1
	if(SHOW_VEC_VOLUME):
		field.drawVolume([0,0,1,1,0,1,1,1],intensity=2,slope_factor=4)
setData()
if(SHOW_VEC_VOLUME):
	showVecVolume()
elif(SHOW_VELOCITY):
	showVelocity()
else:
	showPressure()
