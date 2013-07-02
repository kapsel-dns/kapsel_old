import math
from numpy import *
n=13
d=4.0
g=3.0
r=18.0
x0=16.0
y0=32.0
z0=16.0
cord=[0.,0.,0.]
for i in range(n):
	th=g*i/r
	cord[0]=x0+r-r*cos(th)
	cord[1]=y0-r*sin(th)
	cord[2]=z0
	print i,cord
