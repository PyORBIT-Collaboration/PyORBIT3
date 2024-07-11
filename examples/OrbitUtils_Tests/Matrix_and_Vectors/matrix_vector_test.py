#============================================================
# This script is a test for Matrix and PhaseVector operations
#============================================================

import sys

from orbit.core.orbit_utils import Matrix, PhaseVector

def printM(m):
	print ("----matrix--- size=",m.size())
	for i in range(m.size()[0]):
		for j in range(m.size()[1]):
			print ("m(" + str(i) + "," + str(j)+")="+str(m.get(i,j)) + " ",end="")
		print ("")	

def printV(v):
	print ("----vector--- size=",v.size())
	for i in range(v.size()):
		print ("v(" + str(i) + ")="+str(v.get(i)) + " ", end="")
	print ("")

def setM(m):
	for i in range(m.size()[0]):
		for j in range(m.size()[1]):
			if(i == j): m.set(i,j,i+1.0)
			if(i > j): m.set(i,j,i+j)

def setV(v):
	for i in range(v.size()):
		v.set(i,i+1.0)
		
print ("Start.")

			
#check vector
print ("-----------------Set up v and v1 --------------")
v = PhaseVector(5)
v1 = PhaseVector(5)
setV(v)
setV(v1)
printV(v)
printV(v1)
print ("-----------------v3 = v + 1 as scalar ---------")
v3 = v.add(1)
printV(v3)
print ("------------------v3 = v + v1 -----------------")
v3 = v.add(v1)
printV(v3)
print ("------------------print v*v1 ------------------")
print ("v*v1 =", v.mult(v1))
print ("-----------------------------------------------")

#check matrix
print ("----mult -------------")
m = Matrix(5,4)
m1 = Matrix(4,5)
print ("----------------Print m and m1 ----------------")
printM(m)
printM(m1)
print ("----------------Print m3 = m*m1 ---------------")
setM(m)
setM(m1)
m3 = m.mult(m1)
printM(m3)

print ("----add --------------")
m = Matrix(4,4)
m1 = Matrix(4,4)
setM(m)
setM(m1)
m3 = m.add(m1)
printM(m3)

print ("----invert Result should be unit matrix -------")
m4 = m3.invert()
if(m4 == None): print ("Cannot invert matrix:"); printM(m3); sys.exit(1)
printM(m4.mult(m3))

print ("-----------------------------------------------")

#check vector and matrix
print ("---- v.mult(m) -------------------------")
v = PhaseVector(5)
m = Matrix(5,4)
setV(v)
setM(m)
printM(m)
printV(v)
v1 = v.mult(m)
printV(v1)

print ("---- m.mult(v) -------------------------")
v = PhaseVector(5)
m = Matrix(4,5)
setV(v)
setM(m)
printM(m)
printV(v)
v1 = m.mult(v)
printV(v1)

#----------------------------------------------------------
#---- Comment sys.exit(0) if you want to check memeory leak
#---- in Matrix and Vector operations.
#----------------------------------------------------------
sys.exit(0)

count = 0
while(2 > 1):
	v = PhaseVector(5)
	v1 = PhaseVector(5)
	v_tmp = v.copy()
	setV(v)
	setV(v1)
	v3 = v.add(1)
	v3 = v.add(v1)
	v.mult(v1)
	#------------------------
	m = Matrix(5,4)
	m1 = Matrix(4,5)
	setM(m)
	setM(m1)
	m3 = m.mult(m1)
	#------------------
	m = Matrix(4,4)
	m1 = Matrix(4,4)
	setM(m)
	setM(m1)
	m3 = m.add(m1)
	m4 = m3.invert()
	#-----------------
	v = PhaseVector(5)
	m = Matrix(5,4)
	setV(v)
	setM(m)
	v1 = v.mult(m)
	#---------------
	v = PhaseVector(5)
	m = Matrix(4,5)
	v1 = m.mult(v)
	v_tmp = v.copy()
	count += 1
	if(count%1000 == 0):
		print ("i=",count)

