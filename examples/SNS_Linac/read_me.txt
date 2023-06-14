This a set of test scripts for the new edition of the PyORBIT linac package
which is located in the directory:
py/orbit/py_linac  

This new edition was committed to the repository in August, 2015. 
The lattice factory and the parser for the new edition is based on 
the XmlDataAdaptor (package orbit.utils.xml). This class is a functional 
copy of the OpenXAL class for XML parsing. 

The linac lattice factory was also clean up and improved.
The RF gap models in the new edition are based on

A. Shishlo, J. Holmes, 
"Physical Models for Particle Tracking Simulations in the RF Gap", 
ORNL Tech. Note ORNL/TM-2015/247, June 2015

The old linac package will be removed after all functionality transferred to the 
new package.
The location of the old linac package is
py/orbit/sns_linac/sns_linac
