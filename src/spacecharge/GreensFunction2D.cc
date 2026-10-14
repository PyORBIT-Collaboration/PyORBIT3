#include "GreensFunction2D.hh"

#include <cmath>
#include <vector>

namespace {

int kernelIndex(int x, int y, int y_size)
{
	return y + y_size * x;
}

void reflectKernel(double* kernel, int x_size, int y_size)
{
	const int x_half = x_size / 2;
	const int y_half = y_size / 2;
	for (int y = 0; y <= y_half; ++y) {
		for (int x = x_half + 1; x < x_size; ++x) {
			kernel[kernelIndex(x, y, y_size)] = kernel[kernelIndex(x_size - x, y, y_size)];
		}
	}
	for (int y = y_half + 1; y < y_size; ++y) {
		for (int x = 0; x < x_size; ++x) {
			kernel[kernelIndex(x, y, y_size)] = kernel[kernelIndex(x, y_size - y, y_size)];
		}
	}
}

// An antiderivative of -log(hypot(x, y)).  The integration table uses
// half-cell coordinates, so neither coordinate is zero in normal use.
double integratedGreenAntiderivative(double x, double y)
{
	if (x == 0.0 || y == 0.0) return 0.0;
	const double radius = std::hypot(x, y);
	return -x * y * std::log(radius) + 1.5 * x * y
		- 0.5 * x * x * std::atan(y / x)
		- 0.5 * y * y * std::atan(x / y);
}

} // namespace

namespace GreensFunction2D {

void fillPointKernel(double* kernel, int x_size, int y_size, double dx, double dy)
{
	const int x_half = x_size / 2;
	const int y_half = y_size / 2;
	for (int y = 0; y <= y_half; ++y) {
		const double y_offset = y * dy;
		for (int x = 0; x <= x_half; ++x) {
			const double x_offset = x * dx;
			const double radius_squared = x_offset * x_offset + y_offset * y_offset;
			kernel[kernelIndex(x, y, y_size)] = radius_squared == 0.0 ? 0.0 : -std::log(radius_squared) / 2.0;
		}
	}
	reflectKernel(kernel, x_size, y_size);
}

void fillIntegratedKernel(double* kernel, int x_size, int y_size, double dx, double dy)
{
	const int x_half = x_size / 2;
	const int y_half = y_size / 2;
	const int table_x = x_half + 2;
	const int table_y = y_half + 2;

	std::vector<double> corners(table_x * table_y);
	for (int x = 0; x < table_x; ++x) {
		const double x_corner = (x - 0.5) * dx;
		for (int y = 0; y < table_y; ++y) {
			const double y_corner = (y - 0.5) * dy;
			corners[kernelIndex(x, y, table_y)] = integratedGreenAntiderivative(x_corner, y_corner);
		}
	}

	for (int x = 0; x <= x_half; ++x) {
		for (int y = 0; y <= y_half; ++y) {
			const double integral =
				corners[kernelIndex(x + 1, y + 1, table_y)] - corners[kernelIndex(x, y + 1, table_y)]
				- corners[kernelIndex(x + 1, y, table_y)] + corners[kernelIndex(x, y, table_y)];
			kernel[kernelIndex(x, y, y_size)] = integral / (dx * dy);
		}
	}
	reflectKernel(kernel, x_size, y_size);
}

} // namespace GreensFunction2D
