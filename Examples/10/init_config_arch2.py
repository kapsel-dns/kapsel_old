import math
from numpy import *

n=13
d=4.0
g=3.0
r=18.0
th0=0.025
x0=32.0
y0=64.0
z0=32.0
nn=(n-1)/2
ii=0
cord=[0.,0.,0.]
cord[0]=x0+r*sin(th0)
cord[1]=y0-r*cos(th0)
cord[2]=z0
print ii,cord
for i in range(nn):
	ii=ii+1
	th=th0+g*float(i+1)/r
	cord[0]=x0+r*sin(th)
	cord[1]=y0-r*cos(th)
	cord[2]=z0
	print ii,cord
	ii=ii+1
	th=th0-g*float(i+1)/r
	cord[0]=x0+r*sin(th)
	cord[1]=y0-r*cos(th)
	cord[2]=z0
	print ii,cord
