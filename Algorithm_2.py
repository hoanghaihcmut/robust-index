"""
This code is copyrighted by Institute of Mathematical and Computational Sciences - IMACS, 
	Ho Chi Minh City University of Technology (HCMUT).  
Contact: Prof. Phan Thanh An thanhan@hcmut.edu.vn
"""

#Title: Computing the Robust Indices of Quasiconvex Functions

from sympy import symbols, solveset, diff, Interval, log, oo, is_convex
import numpy as np

#Algorithm_2: Finding the robust index of quasiconvex functions on [a,b]
def Algorithm_2(f,x,a,b,gamma=10**-2,zmin = 10**-323,alpha=0):
	x = symbols(str(x))
	I = Interval(a,b)

	def L_abs_df(alpha):
		return solveset(abs(diff(f,x))<=alpha,x,domain=I)

	if is_convex(f,x,domain=I) == True:
		sf = oo
	else:
		flag = 0
		alpha = alpha + gamma
		while is_convex(f,x,domain=L_abs_df(alpha)) == True:
			flag = 1
			alpha = alpha + gamma
		
		if flag==1:
			sf = alpha - gamma
		else:
			while is_convex(f,x,domain=L_abs_df(alpha)) == False:
				alpha = alpha - gamma
				if alpha <= zmin:
					sf = -oo
					break
				elif is_convex(f,x,domain=L_abs_df(alpha)) == True:
					sf = alpha
					break

	return sf

#Examples
if __name__ == "__main__":
	x = symbols('x')

	a=0
	b=1
	f = 1/3*x**3-2*x**2+4*x

	sf = Algorithm_2(f,x,a,b,gamma=10**-1)

	print('f(x) = '+str(f)+', D = ['+str(a)+','+str(b)+']')
	print('sf = '+str(sf))
