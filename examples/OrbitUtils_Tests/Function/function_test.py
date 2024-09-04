#------------------------------------------------------
# This is an example of Function and SplineCH
# They are containers for (x,y) (or (x,y,err) ) points
# Function provides a linear interpolation, and SplineCH
# uses 3-rd order polynomials. SplineCH can be used for 
# derivatives calculations.
#-------------------------------------------------------

import sys
import math

from orbit.core.orbit_utils import Function
from orbit.core.orbit_utils import SplineCH

f = Function()

def FF(x):
	return math.sin(x)

def FFP(x):
	return math.cos(x)

#------------------------------------------
# Let's create PyORBIT Function = sin(x)
#------------------------------------------
n = 50
step = 2*math.pi/n
for i in range(n):
	x = step*i+0.1*((1.0*i)/n)**2;
	y = FF(x)
	f.add(x,y)

"""
#---- Test of creating of x,y,err lists
(x_arr,y_arr,err_arr) = f.getXYErrLists()
print ("debug x_arr  =",x_arr)
print ("debug y_arr  =",y_arr)
print ("debug err_arr=",err_arr)
sys.exit(0)
"""

"""
#----- Test of initialization form PyLists
x_arr   = [0.,1.,2.,3.]
y_arr   = [1.,2.,3.,4.]
err_arr = [0.1,0.2,0.3,0.4]
f_test = Function()
#---------------------------------------
f_test.initFromLists(x_arr,y_arr)
f_test.dump()
#-----------------------------------------
f_test.initFromLists(x_arr,y_arr,err_arr)
f_test.dump()
sys.exit(0)
"""

#--------------------------------
# This will print out Function
#--------------------------------
#f.dump()
#f.dump("funcion_x_y_err.dat")

#--------------------------------------
# This will create spline from Function
#--------------------------------------
spline = SplineCH()
spline.compile(f)

#--------------------------------
# This will print out Spline
#--------------------------------	
#spline.dump()
#spline.dump("spline_x_y_err.dat")

print ("====================================")
n = 30
step = 0.9*(f.getMaxX() - f.getMinX())/(n-1)
y_dev_max = 0.
yp_dev_max = 0.
for j in range(n):
	x = f.getMinX() + j*step
	y = FF(x)
	yp = FFP(x)
	ys = spline.getY(x)
	yps = spline.getYP(x)
	ys = spline.getY(x)
	yps = spline.getYP(x)
	dy = abs(ys - y)
	dyp = abs(yps - yp)
	if(y_dev_max < dy): y_dev_max = dy
	if(yp_dev_max < dyp): yp_dev_max = dyp
	#------------------------
	st  = " x= %+8.6f "%x
	st += " y = %+8.6f ys = %+8.6f "%(y,ys)
	st += " err = %+8.6f   "%dy
	st += " yp = %+8.6f yps = %+8.6f "%(yp,yps)
	st += " err = %+8.6f   "%dyp
	print (st)
	
print ("====================================")

print (" max deviation y  =",y_dev_max)
print (" max deviation yp =",yp_dev_max)

print ("Stop.")

#sys.exit(0)

#----------------------------------
# Test removePoint and updatePoint 
#----------------------------------
print ("====================================")
f = Function()
f.add(1.,1.)
f.add(2.,2.)
f.add(3.,3.)

f_test = Function()

iteration = 0

while(1 < 2):
	#f.dump()
	#print ("====================================")
	f.removePoint(0)
	#f.dump()
	#print ("====================================")
	f.add(1.,1.)
	#f.dump()
	#print ("====================================")
	f.removePoint(1)
	#f.dump()
	#print ("====================================")
	f.add(2.,2.)
	#f.dump()
	#print ("====================================")
	f.removePoint(2)
	#f.dump()
	#print ("====================================")
	f.add(3.,3.)
	#f.dump()
	#print ("====================================")
	f.updatePoint(0,1.)
	f.updatePoint(1,2.)
	f.updatePoint(2,2.)
	#print ("====================================")
	(x_arr,y_arr,err_arr) = f.getXYErrLists()
	#print ("====================================")
	f_test.initFromLists(x_arr,y_arr)
	#-----------------------------------------
	f_test.initFromLists(x_arr,y_arr,err_arr)
	#print ("====================================")
	if(iteration%100000 == 0):
		print ("iter=",iteration)
	iteration += 1

print ("====================================")
print ("====================================")
print ("====================================")

f.updatePoint(0,5.)
f.dump()
print ("====================================")
f.updatePoint(2,7.,3.)
f.dump()
print ("====================================")