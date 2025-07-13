"""
Copyright © Institute of Mathematical and Computational Sciences (IMACS),
Ho Chi Minh City University of Technology (HCMUT).

Contact: Prof. Phan Thanh An <thanhan@hcmut.edu.vn>

#Title: Finding the robust index of the quasiconvex function f_8 on D_8
"""

from Algorithm_1 import *
from Algorithm_2 import *

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

#Get the points on the boundary of the unit sphere, m=(m_1,m_2)
def partial_D(m):
    def C(t):
        r=1
        h=t[0]
        theta=t[1]
        x = r * np.sin(np.arccos(h/r)) * np.cos(theta)
        y = r * np.sin(np.arccos(h/r)) * np.sin(theta)
        z = h
        return (x,y,z)
    partial_D_m = []
    r=1
    for i in range(2**m[0]+1):
        h=-r+2*r*i/2**m[0]
        for j in range(2**m[1]):
            if h in [-r,r]:
                theta=0
                t=[h,theta]
                partial_D_m = partial_D_m + [C(t)]
                break
            else:
                theta = 2*np.pi*j/2**m[1]

                t=[h,theta]
                partial_D_m = partial_D_m + [C(t)]
    return partial_D_m

#Algorithm 3: Finding the robust index of the quasiconvex function f_8
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x,y,z = symbols('x y z')
f = log(1*(x**2+2*y**2+3*z**2+1))
m=(1,2)

partial_D_m=partial_D(m)
S_m=[]
for i in range(len(partial_D_m)):
    u=partial_D_m[i]
    ax.scatter(u[0], u[1], u[2], c='r', s=10)
    for j in range(i+1,len(partial_D_m)):
        v=partial_D_m[j]

        t = symbols('t')
        uv = np.sqrt((u[0]-v[0])**2+(u[1]-v[1])**2+(u[2]-v[2])**2)

        xt = u[0] + t*(v[0]-u[0])/uv
        yt = u[1] + t*(v[1]-u[1])/uv
        zt = u[2] + t*(v[2]-u[2])/uv
        g = f.subs([(x,xt),(y,yt),(z,zt)])

        s_g = algorithm_2(g,t,0,uv)
        S_m = S_m+[s_g]
        
        ax.plot([u[0], v[0]], [u[1], v[1]], [u[2], v[2]], c='b',linewidth=0.5)

print('sf('+str(m[0])+','+str(m[1])+') = '+str(min(S_m)))

r = 1
u, v = np.mgrid[0:2*np.pi:30j, 0:np.pi:30j]
x = r*np.cos(u)*np.sin(v)
y = r*np.sin(u)*np.sin(v)
z = r*np.cos(v)

surf = ax.plot_surface(x, y, z, cmap='viridis')
surf.set_alpha(0.2)
ax.grid(False)
ax.set_axis_off()
plt.show()
