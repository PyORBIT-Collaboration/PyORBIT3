#-------------------------------------------------------------
# This is an example of Polynomial class
# The methods "derivativeTo","derivative", "copyTo" and "value" are included.
#-------------------------------------------------------------

import sys
import math

from orbit.core.orbit_utils import Polynomial

poly = Polynomial(4)
print ("initial order = ",poly.order())

poly.coefficient(2,2.0)

print ("========Polynomial=======")
order = poly.order()
for i in range(order+1):
	print ("i=",i," coef=",poly.coefficient(i))
	
x = 10.
y = poly.value(x)
print ("x=",x," y=",y)

poly1 = Polynomial()
poly.derivativeTo(poly1)

print ("========Derivative polynomial=======")
order = poly1.order()
for i in range(order+1):
	print ("i=",i," coef=",poly1.coefficient(i))

x = 10.
y = poly1.value(x)
yp = poly.derivative(x)
print ("x=",x," y_prime=",y,"  yp =",yp)

print ("========Copy of initial polynomial=======")
poly2 = Polynomial()
poly.copyTo(poly2)
order = poly2.order()
for i in range(order+1):
	print ("i=",i," coef=",poly2.coefficient(i))

x = 10.
y = poly2.value(x)
yp = poly2.derivative(x)
print ("x=",x," y=",y," y_prime=",yp)

"""
# Memory leak test
count = 1
while(1 < 2):
 poly1 = Polynomial()
 poly2 = Polynomial()
 poly.derivativeTo(poly1)
 poly.copyTo(poly2)
 count += 1
 if(count % 100000 == 0):
	 print ("count=",count)
"""

print ("Stop.")
sys.exit(0)

