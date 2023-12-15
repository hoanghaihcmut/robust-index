"""
This code is copyrighted by Institute of Mathematical and Computational Sciences - IMACS, 
	Ho Chi Minh City University of Technology (HCMUT).  
Contact: Prof. Phan Thanh An thanhan@hcmut.edu.vn
"""

#Title: Computing the Robust Indices of Quasiconvex Functions

from sympy import symbols, solveset, diff, Interval, log, oo,EmptySet, is_decreasing, is_increasing, is_convex, is_monotonic
import numpy as np

#Checking for the quasiconvex function on [a,b]
def is_quasiconvex(f,x,a,b):
	x = symbols(str(x))
	I = Interval(a,b)

	isQC = False
	
	df = diff(f,x)

	if is_decreasing(f,I) or is_increasing(f,I):
		isQC = True
	else:
		X = list(solveset(df,x,domain=I))
		x_min=X[0]
		for i in range(1,len(X)):
			if f.subs(x,X[i])<f.subs(x,x_min):
				x_min=X[i]

		Ia = Interval(I.start, x_min)
		Ib = Interval(x_min, I.end)

		if is_decreasing(f,Ia) and is_increasing(f,Ib):
			isQC = True
	return isQC

#Algorithm_1: Finding the robust index of quasiconvex functions on [a,b]
def Algorithm_1(f,x,a,b):
    x = symbols(str(x))

    I = Interval(a,b)
    sf = -oo
    if is_convex(f,x,domain=I):
        sf = oo
    else:
	    df = f.diff()
	    dff = df.diff()
	    sol_dff = solveset(dff < 0, x, domain=Interval.open(a, b))

	    S=[]
	    sol_df_1 = solveset(df > 0, x, domain=sol_dff)
	    if sol_df_1!=EmptySet:
	    	S=S+[df.subs(x, sol_df_1.sup)]
	    sol_df_2 = solveset(-df > 0, x, domain=sol_dff)
	    if sol_df_2!=EmptySet:
	    	S=S+[-df.subs(x, sol_df_2.inf)]

	    sf_diamond = min(S)
	    if sf_diamond==0:
	    	sf=-oo
	    else:
	    	sf = sf_diamond
	    	for t in list(solveset(dff,x,domain=Interval.open(a,b))):
		    	df_t = df.subs(df,t)
		    	if is_quasiconvex(f-df_t*x,x,a,b)==False:
		    		sf = min(sf,abs(df_t))
    return sf

#Examples
if __name__ == "__main__":
	x = symbols('x')

	a=0
	b=1
	f = 1/3*x**3-2*x**2+2*x

	sf = Algorithm_1(f,x,a,b)

	print('f(x) = '+str(f)+', D = ['+str(a)+','+str(b)+']')
	print('sf = '+str(sf))