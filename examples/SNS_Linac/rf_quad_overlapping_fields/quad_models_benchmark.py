#---------------------------------------------------------------
# This script will read the dB/dr of the quad magnetic field distribution 
# anong the longitudinal axis and will peform Twiss parameters tracking
# along this lattice. The user can compare these results with the PyORBIT 
# tracking results without space charge. They should be be very close.

import sys
import math
import orbit.core

from orbit_utils import Matrix
from orbit_utils import PhaseVector
from orbit_utils import Function

#---- Print matrix function
def printM(m):
	print("----matrix--- size=",m.size())
	for i in range(m.size()[0]):
		st = ""
		for j in range(m.size()[1]):
			st += "m(%1d,%1d)="%(i,j)+" %12.5g  "%(m.get(i,j))
			#print ("m(" + str(i) + "," + str(j)+")="+str(m.get(i,j)) + " "),
		#print ""
		print(st)
		
#---- Print vector function
def printV(v):
	print("----vector--- size=",v.size())
	for i in range(v.size()):
		print(("v(" + str(i) + ")="+str(v.get(i)) + " "), end=' ')
	print("")

#---- Let's read the dB/dr distribution. It could be equals to zero, 
#---- and in this case we will treat this as drift
file_name = "MEBT_fields_rf_and_quads.dat"
fl_in = open(file_name,"r")
lns = fl_in.readlines()
fl_in.close()

#---- first line in the input files is a headers line
g_func = Function()
for ind in range(1,len(lns)):
	ln = lns[ind].strip()
	res_arr = ln.split()
	if(len(res_arr) > 5):
		#---- here we parse the input file: 1st position - z, 4th - dB/dr(z
		z = float(res_arr[0])
		g = float(res_arr[3])
		g_func.add(z,g)

#---- let's set energy and momentum
mass = 0.938272 + 2*0.000511
eKin = 0.0025
momentum = math.sqrt(eKin*(eKin + 2*mass))
Brho = 3.33564*momentum  # [T*m]
charge = -1.0

#---- now let's prepare array of 4x4 (transport) 
#---- and 6x6 (Twiss params transport) matrices
matrix_arr = []
for ind in range(1,g_func.getSize()):
	z0 = g_func.x(ind-1)
	z1 = g_func.x(ind)
	g0 = g_func.y(ind-1)
	g1 = g_func.y(ind)
	z_avg = (z1+z0)/2
	g_avg = (g0+g1)/2
	g_avg *= charge
	dL = z1-z0
	#print "debug dL=",dL, " y= %12.3f "%y_avg
	#---- this is a transport 4x4 matrix for coordinates x,x',y,y'
	trMtrx = Matrix(4,4)
	k2 = abs(g_avg/Brho)
	kq = math.sqrt(k2)
	sin_kL = math.sin(kq*dL)
	cos_kL = math.cos(kq*dL)
	sinh_kL = math.sinh(kq*dL)
	cosh_kL = math.cosh(kq*dL)
	if(g_avg == 0.):
		trMtrx.set(0,0,1.)
		trMtrx.set(1,1,1.)
		trMtrx.set(2,2,1.)
		trMtrx.set(3,3,1.)
		trMtrx.set(0,1,dL)
		trMtrx.set(2,3,dL)
	else:
		if(g_avg > 0.):
			trMtrx.set(0,0,cos_kL)
			trMtrx.set(1,1,cos_kL)
			trMtrx.set(0,1,sin_kL/kq)	
			trMtrx.set(1,0,-kq*sin_kL)
			#-------------------------
			trMtrx.set(2,2,cosh_kL)
			trMtrx.set(3,3,cosh_kL)
			trMtrx.set(2,3,sinh_kL/kq)
			trMtrx.set(3,2,kq*sinh_kL)
		else:
			trMtrx.set(2,2,cos_kL)
			trMtrx.set(3,3,cos_kL)
			trMtrx.set(2,3,sin_kL/kq)	
			trMtrx.set(3,2,-kq*sin_kL)
			#-------------------------
			trMtrx.set(0,0,cosh_kL)
			trMtrx.set(1,1,cosh_kL)
			trMtrx.set(0,1,sinh_kL/kq)
			trMtrx.set(1,0,kq*sinh_kL)
	#------------------- this is 6x6 matrix for Twiss parameters for x and y axis
	trTwissMtrx = Matrix(6,6)
	#---- dir_ind = 0 for x and 1 for y-directions
	for dir_ind in range(2):
		#---- index shift
		ish = 3*dir_ind
		ind1 = 0 + 2*dir_ind
		ind2 = 1 + 2*dir_ind
		trTwissMtrx.set(0+ish,0+ish,trMtrx.get(ind1,ind1)*trMtrx.get(ind2,ind2) + trMtrx.get(ind1,ind2)*trMtrx.get(ind2,ind1))
		trTwissMtrx.set(0+ish,1+ish,- trMtrx.get(ind1,ind1)*trMtrx.get(ind2,ind1))
		trTwissMtrx.set(0+ish,2+ish,- trMtrx.get(ind1,ind2)*trMtrx.get(ind2,ind2))
		trTwissMtrx.set(1+ish,0+ish,- 2*trMtrx.get(ind1,ind1)*trMtrx.get(ind1,ind2))
		trTwissMtrx.set(1+ish,1+ish,trMtrx.get(ind1,ind1)*trMtrx.get(ind1,ind1))
		trTwissMtrx.set(1+ish,2+ish,trMtrx.get(ind1,ind2)*trMtrx.get(ind1,ind2))
		trTwissMtrx.set(2+ish,0+ish,- 2*trMtrx.get(ind2,ind1)*trMtrx.get(ind2,ind2))
		trTwissMtrx.set(2+ish,1+ish,trMtrx.get(ind2,ind1)*trMtrx.get(ind2,ind1))
		trTwissMtrx.set(2+ish,2+ish,trMtrx.get(ind2,ind2)*trMtrx.get(ind2,ind2))		
	matrix_arr.append([z_avg,trMtrx,trTwissMtrx])
	#print "debug trMtrx ===="
	#printM(trMtrx)
	#print "debug trTwissMtrx ===="
	#printM(trTwissMtrx)
	#sys.exit()

#----------------------------------------------
(alphaX,betaX,emittX) = (-1.9620, 0.1831, 0.21)
(alphaY,betaY,emittY) = ( 1.7681, 0.1620, 0.21)

gammaX = (1+alphaX**2)/betaX
gammaY = (1+alphaY**2)/betaY


twiss_v_ini = PhaseVector(6)
twiss_v_ini.set(0,alphaX)
twiss_v_ini.set(1,betaX)
twiss_v_ini.set(2,gammaX)
twiss_v_ini.set(3,alphaY)
twiss_v_ini.set(4,betaY)
twiss_v_ini.set(5,gammaY)

print("===========================================")
print("Initial Twiss parameters:")
printV(twiss_v_ini)
print("===========================================")

fl_out = open("mebt_twiss_matr_distributed_quads.dat","w")
st = " ind x[m]  alphaX  betaX    alphaY  betaY "
fl_out.write(st+"\n")

v_in = twiss_v_ini
for ind in range(len(matrix_arr)):
	[z_avg,trMtrx,trTwissMtrx] = matrix_arr[ind]
	v_out = trTwissMtrx.mult(v_in)
	alpha_x=v_out.get(0)
	beta_x= v_out.get(1)
	gamma_x=v_out.get(2)
	alpha_y=v_out.get(3)
	beta_y= v_out.get(4)
	gamma_y=v_out.get(5)
	v_in = v_out
	st = " %4d %9.3f   %12.6f  %12.6f    %12.6f  %12.6f "%(ind,z_avg,alpha_x,beta_x,alpha_y,beta_y)
	if(ind % 100 == 0): print(st)
	fl_out.write(st+"\n")
	
fl_out.close()	

  

		
	
	
