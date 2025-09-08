"""
This is a pure Python subclass of the C++ Bunch class from orbit.core.bunch.
Users can extend the C++ Bunch with pure Python methods.
"""

import sys

from orbit.core.bunch import Bunch as cppBunch

class Bunch(cppBunch):
    def __init__(self):
        cppBunch.__init__(self)
        
    def getNumpyVersion(self):
        import numpy as np
        return np.__version__

if __name__ == '__main__':
    from orbit.py_linac.lattice import Drift
    
    py_bunch = Bunch()
    eKin = 1.4
    py_bunch.getSyncParticle().kinEnergy(eKin)
    
    x = 1.0; y = 2.0; z = 3.0
    xp = -1.0; yp = -2.0; zp = -3.0
    py_bunch.addParticle(x*1.e-3,xp*1.e-3,y*1.e-3,yp*1.e-3,z*1.e-3,zp*1.e-3)
    
    print ("============== Initial coordinates =============")
    py_bunch.dumpBunch() 
    
    drift = Drift("drift")
    drift.setLength(1.0)
    drift.trackBunch(py_bunch)
    print ("============== New coordinates =================")
    py_bunch.dumpBunch() 
    
    print ("Stop.")
    sys.exit(0)