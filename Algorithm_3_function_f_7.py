"""
Copyright © Institute of Mathematical and Computational Sciences (IMACS),
Ho Chi Minh City University of Technology (HCMUT).

Contact: Prof. Phan Thanh An <thanhan@hcmut.edu.vn>

#Title: Finding the robust index of the quasiconvex function f_7 on D_7
"""

from Algorithm_1 import *
from Algorithm_2 import *

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.widgets import Slider, Button

#Get the line segments [u,v] with u,v on the boundary of D_7=[0,1]x[-1,1]
def line_segments(D,m):
	SEGs=[]
	P = [(D[0][0], D[1][0]), (D[0][1], D[1][0]), (D[0][1], D[1][1]), (D[0][0],D[1][1]),(D[0][0], D[1][0])]

	a = 0
	b = len(P)-1

	def xC(t):
		if t<b:
			return (1-t+np.floor(t))*P[int(np.floor(t))][0] + (t-np.floor(t))*P[int(np.floor(t)+1)][0]
		else:
			return P[0][0]
	def yC(t):
		if t<b:
			return (1-t+np.floor(t))*P[int(np.floor(t))][1] + (t-np.floor(t))*P[int(np.floor(t)+1)][1]
		else:
			return P[0][1]
	
	del_D = [(xC(a+(b-a)/2**m*i),yC(a+(b-a)/2**m*i)) for i in range(2**m)]

	for i in range(2**m):
		u = del_D[i]
		for j in range(i+1,2**m):
			v = del_D[j]
			SEGs = SEGs + [(u,v)]
	return SEGs

#Algorithm 3: Finding the robust index of the quasiconvex function f_7 on D_7
x,y = symbols('x y')
D = [[1,2],[1,2]]
f = log(x**2+2*y**2)

m=3
S_m = []

SEGs=line_segments(D,m)
for seg in SEGs:
	u = seg[0]
	v = seg[1]
	t = symbols('t')
	
	uv = np.sqrt((u[0]-v[0])**2+(u[1]-v[1])**2)
	a=0
	b=uv
	I = Interval(a, b)
	xt = u[0] + t*(v[0]-u[0])/uv
	yt = u[1] + t*(v[1]-u[1])/uv

	g = f.subs([(x,xt),(y,yt)])

	sg = Algorithm_1(g,t,a,b)
	S_m = S_m + [sg]
	if sg<0:
		break

print('sf('+str(m)+') = '+str(min(S_m)))

#Illustrative for S_m
fig = plt.figure(figsize=(5, 5))
i=0
for seg in SEGs:
	u = seg[0]
	v = seg[1]

	RGB = np.random.rand(1,3)
	plt.plot([u[0],v[0]], [u[1],v[1]], color = RGB[0])
	plt.text(u[0]+0.2*(v[0]-u[0]), (u[1]+0.2*(v[1]-u[1])), str(round(S_m[i],3)),ha='center',color=RGB[0])

	i=i+1
plt.show()
