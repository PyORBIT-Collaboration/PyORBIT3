. .po3/bin/activate

# Simple test to check if MPI was linked at all
mpirun -np 2 python examples/MPI_Tests/mpi_initialization_test.py | tee test.txt
calculated_sha1=$(sha1sum test.txt | awk '{ print $1 }')
[ "$calculated_sha1" != "d94e03044bf06c7a42d07505b50fa58b4b30e49a" ] && exit 1

# Uniformly charged sphere
echo "Run with Ellipsoid space charge"
mpirun -np 2 python examples/SpaceCharge/sc3D/sc_3d_drift_latt_uniform_sphere_bunch.py --SC ellipsoid| tee results.txt
python examples/SpaceCharge/sc3D/read_sc_numbers.py
[ $? != 0 ] && exit 1

echo "Run with FFT3D space charge"
mpirun -np 2 python examples/SpaceCharge/sc3D/sc_3d_drift_latt_uniform_sphere_bunch.py --SC fft3d| tee results.txt
python examples/SpaceCharge/sc3D/read_sc_numbers.py
[ $? != 0 ] && exit 1

# this should fail as no space charge will give incorrect sphere size
echo "Run with no space charge"
mpirun -np 2 python examples/SpaceCharge/sc3D/sc_3d_drift_latt_uniform_sphere_bunch.py --SC none| tee results.txt
python examples/SpaceCharge/sc3D/read_sc_numbers.py
[ $? == 0 ] && exit 1


# This tests that  five nodes give the same result as one node
echo "Run with 1 node."
mpirun -np 1 python examples/SpaceCharge/sc3D/sc_3d_drift_latt_uniform_sphere_bunch.py | tee results1.txt
calculated_sha1=$(sha1sum results1.txt | awk '{ print $1 }')

echo "Run with 5 nodes."
mpirun -np 5 python examples/SpaceCharge/sc3D/sc_3d_drift_latt_uniform_sphere_bunch.py | tee results5.txt
calculated_sha5=$(sha1sum results5.txt | awk '{ print $1 }')

echo "Diff between two runs."
echo "$calculated_sha1" "$calculated_sha5"
diff results1.txt results5.txt

[ "$calculated_sha1" != "$calculated_sha5" ] && exit 1

exit 0

