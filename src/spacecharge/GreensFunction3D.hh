#ifndef SC_GREENS_FUNCTION_3D_H
#define SC_GREENS_FUNCTION_3D_H

namespace GreensFunction3D {

void fillPointKernel(double* kernel, int x_size, int y_size, int z_size,
			double dx, double dy, double dz, int n_bunches, double lambda);

void fillIntegratedKernel(double* kernel, int x_size, int y_size, int z_size,
			double dx, double dy, double dz, int n_bunches, double lambda);

} // namespace GreensFunction3D

#endif
