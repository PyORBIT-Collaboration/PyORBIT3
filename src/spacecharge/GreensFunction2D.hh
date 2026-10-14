#ifndef SC_GREENS_FUNCTION_2D_H
#define SC_GREENS_FUNCTION_2D_H

namespace GreensFunction2D {

void fillPointKernel(double* kernel, int x_size, int y_size, double dx, double dy);
void fillIntegratedKernel(double* kernel, int x_size, int y_size, double dx, double dy);

} // namespace GreensFunction2D

#endif
