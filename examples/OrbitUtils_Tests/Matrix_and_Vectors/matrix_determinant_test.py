#==================================
# This script is a test for Matrix 
# 1. Determinant
# 2. Transpose
# 3. Inversion
# calculations.
#==================================

import sys
import time

from orbit.core.orbit_utils import Matrix

def printM(m):
	print ("----matrix--- size=",m.size())
	for i in range(m.size()[0]):
		for j in range(m.size()[1]):
			print ("m(" + str(i) + "," + str(j)+")="+str(m.get(i,j)) + " ",end="")
		print("")

def setM(arr):
	n_row = len(arr)
	n_col = len(arr[0])
	m = Matrix(n_row,n_col)
	for i in range(n_row):
		for j in range(n_col):
			m.set(i,j,arr[i][j])
	return m
		
print ("Start.")

"""

arr = [[],[]]
arr[0] = [ 1., 0.]
arr[1] = [ 0., 3.]

arr = [[],[],[]]
arr[0] = [ 1., 0., 0.]
arr[1] = [ 0., 2., 0.]
arr[2] = [ 0., 0., 3.]

arr = [[],[],[],[]]
arr[0] = [ 1., 0., 0., 0.]
arr[1] = [ 0., 2., 0., 0.]
arr[2] = [ 0., 0., 3., 0.]
arr[3] = [ 0., 0., 0., 4.]

"""

arr = [[],[],[],[],[]]
arr[0] = [ 1., 1., 2., 3., 4.]
arr[1] = [ 0., 2., 3., 4., 5.]
arr[2] = [ 0., 0., 3., 6., 7.]
arr[3] = [ 0., 0., 0., 4., 8.]
arr[4] = [ 0., 0., 0., 0., 5.]         
             
mtrx = setM(arr)
print ("===================================")
printM(mtrx)
print ("===================================")
print ("Det(mtrx)=",mtrx.det())
print ("===========Inverted================")
mtrx_invert = mtrx.invert()
printM(mtrx_invert)
print ("===========Caclulate M*M^-1 =======")
mtrx_res = mtrx_invert.mult(mtrx)
printM(mtrx_res)
print ("===================================")

arr = [[],[],[]]
arr[0] = [ 1., 1., 2., 3., 4., 5.]
arr[1] = [ 0., 2., 3., 4., 5., 6.]
arr[2] = [ 0., 0., 3., 6., 7., 8.]
mtrx_new = setM(arr)
printM(mtrx_new)
print ("===============Transpose===========")
mtrx_new_transp = mtrx_new.transpose()
printM(mtrx_new_transp)
print ("===============Initial=============")
printM(mtrx_new)
print ("===================================")

print ("Stop.")

#----------------------------------------------------------
#---- Comment sys.exit(0) if you want to check memeory leak
#---- in Matrix operations.
#----------------------------------------------------------
sys.exit(0)

count = 1
start_time = time.time()
while(1<2):
	det = mtrx.det()
	mtrx_invert = mtrx.invert()
	mtrx_res = mtrx_invert.mult(mtrx)
	mtrx_new_transp = mtrx_new.transpose()
	mtrx_invert = mtrx.invert()
	if(count % 100000 == 0):
		print ("count =",count," speed [calc/sec] = ",count/(time.time()-start_time))
	count += 1
	
