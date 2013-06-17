#
# $Id: plot.py,v 1.1 2006/06/27 18:31:22 nakayama Exp $
#
import math
locationList = ['t','Particles[].R','Particles[].v']
timeList = ['timestep']
particleList = ['rx','ry','rz','vx','vy','vz']
dataList = []
N1=$object_type.spherical_particle.Particle_spec[0].Particle_number
nLocation = 1 + N1*len(particleList)
for n in range(0,nLocation):
	dataList.append([])
for n in range(0,totalRecord()):
	jump(n)
	tmpData = []
	time = get(locationList[0])
	rx = $Particles[0].R.x
	ry = $Particles[0].R.y
	rz = $Particles[0].R.z
	vx = $Particles[0].v.x
	vy = $Particles[0].v.y
	vz = $Particles[0].v.z
#	print n,rx,ry,rz,vx,vy,vz,time
	tmpData.append(time)
	for i in range(0,N1):
		tmpData.append(rx)
		tmpData.append(ry)
		tmpData.append(rz)
		tmpData.append(vx)
		tmpData.append(vy)
		tmpData.append(vz)
	for n in range(0, nLocation):
		dataList[n].append(tmpData[n])
#print dataList
for i in range(0,N1):
	for j in range(0,len(particleList)):
		timeList.append(particleList[j] + `i`)
tagList = timeList
#print tagList
#print len(tagList)
for n in range(0, nLocation):
	createSheetCol(n,tagList[n])
	setSheetCol(n,dataList[n])
