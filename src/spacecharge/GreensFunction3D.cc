#include "GreensFunction3D.hh"

#include <cmath>
#include <vector>

namespace {

int kernelIndex(int x, int y, int z, int y_size, int z_size)
{
	return z + z_size * (y + y_size * x);
}

void reflectKernel(double* kernel, int x_size, int y_size, int z_size)
{
	const int x_half = x_size / 2;
	const int y_half = y_size / 2;
	const int z_half = z_size / 2;

	for (int z = 0; z <= z_half; ++z) {
		for (int y = 0; y <= y_half; ++y) {
			for (int x = x_half + 1; x < x_size; ++x) {
				kernel[kernelIndex(x, y, z, y_size, z_size)] = kernel[kernelIndex(x_size - x, y, z, y_size, z_size)];
			}
		}
		for (int y = y_half + 1; y < y_size; ++y) {
			for (int x = 0; x < x_size; ++x) {
				kernel[kernelIndex(x, y, z, y_size, z_size)] = kernel[kernelIndex(x, y_size - y, z, y_size, z_size)];
			}
		}
	}

	for (int z = z_half + 1; z < z_size; ++z) {
		for (int x = 0; x < x_size; ++x) {
			for (int y = 0; y < y_size; ++y) {
				kernel[kernelIndex(x, y, z, y_size, z_size)] = kernel[kernelIndex(x, y, z_size - z, y_size, z_size)];
			}
		}
	}
}

double stableLogArgument(double coordinate, double other_a, double other_b, double radius)
{
	if (coordinate >= 0.0) {
		return std::log(coordinate + radius);
	}
	return std::log((other_a * other_a + other_b * other_b) / (radius - coordinate));
}

// Equation (8) of Qiang et al.  The zero-prefactor checks make the
// antiderivative well-defined on coordinate planes used by image bunches.
double integratedGreenAntiderivative(double x, double y, double z)
{
	const double radius = std::hypot(std::hypot(x, y), z);
	double value = 0.0;

	if (y != 0.0 && z != 0.0) {
		value += y * z * stableLogArgument(x, y, z, radius);
	}
	if (x != 0.0 && z != 0.0) {
		value += x * z * stableLogArgument(y, x, z, radius);
	}
	if (x != 0.0 && y != 0.0) {
		value += x * y * stableLogArgument(z, x, y, radius);
	}
	if (z != 0.0) {
		value -= 0.5 * z * z * std::atan((x * y) / (z * radius));
	}
	if (y != 0.0) {
		value -= 0.5 * y * y * std::atan((x * z) / (y * radius));
	}
	if (x != 0.0) {
		value -= 0.5 * x * x * std::atan((y * z) / (x * radius));
	}
	return value;
}

} // namespace

namespace GreensFunction3D {

void fillPointKernel(double* kernel, int x_size, int y_size, int z_size,
		double dx, double dy, double dz, int n_bunches, double lambda)
{
	const int x_half = x_size / 2;
	const int y_half = y_size / 2;
	const int z_half = z_size / 2;
	for (int z = 0; z <= z_half; ++z) {
		const double z_offset = z * dz;
		for (int y = 0; y <= y_half; ++y) {
			const double y_offset = y * dy;
			for (int x = 0; x <= x_half; ++x) {
				const double x_offset = x * dx;
				const double radius = std::sqrt(x_offset * x_offset + y_offset * y_offset + z_offset * z_offset);
				double external_phi = 0.0;
				for (int bunch = -n_bunches / 2; bunch <= n_bunches / 2; ++bunch) {
					if (bunch == 0) continue;
					const double image_z = z_offset + bunch * lambda;
					const double image_radius = std::sqrt(x_offset * x_offset + y_offset * y_offset + image_z * image_z);
					if (image_radius > 1.0e-10) external_phi += 1.0 / image_radius;
				}
				kernel[kernelIndex(x, y, z, y_size, z_size)] = radius == 0.0 ? external_phi : 1.0 / radius + external_phi;
			}
		}
	}
	reflectKernel(kernel, x_size, y_size, z_size);
}

void fillIntegratedKernel(double* kernel, int x_size, int y_size, int z_size,
		double dx, double dy, double dz, int n_bunches, double lambda)
{
	const int x_half = x_size / 2;
	const int y_half = y_size / 2;
	const int z_half = z_size / 2;
	const int table_x = x_half + 2;
	const int table_y = y_half + 2;
	const int table_z = z_half + 2;
	const double cell_volume = dx * dy * dz;
	const auto index = [table_y, table_z](int x, int y, int z) {
		return z + table_z * (y + table_y * x);
	};

	for (int x = 0; x <= x_half; ++x) {
		for (int y = 0; y <= y_half; ++y) {
			for (int z = 0; z <= z_half; ++z) {
				kernel[kernelIndex(x, y, z, y_size, z_size)] = 0.0;
			}
		}
	}

	for (int bunch = -n_bunches / 2; bunch <= n_bunches / 2; ++bunch) {
		std::vector<double> corners(table_x * table_y * table_z);
		for (int x = 0; x < table_x; ++x) {
			const double x_corner = (x - 0.5) * dx;
			for (int y = 0; y < table_y; ++y) {
				const double y_corner = (y - 0.5) * dy;
				for (int z = 0; z < table_z; ++z) {
					const double z_corner = (z - 0.5) * dz + bunch * lambda;
					corners[index(x, y, z)] = integratedGreenAntiderivative(x_corner, y_corner, z_corner);
				}
			}
		}

		for (int x = 0; x <= x_half; ++x) {
			for (int y = 0; y <= y_half; ++y) {
				for (int z = 0; z <= z_half; ++z) {
					const double integral =
						corners[index(x + 1, y + 1, z + 1)] - corners[index(x, y + 1, z + 1)]
						- corners[index(x + 1, y, z + 1)] + corners[index(x, y, z + 1)]
						- corners[index(x + 1, y + 1, z)] + corners[index(x, y + 1, z)]
						+ corners[index(x + 1, y, z)] - corners[index(x, y, z)];
					kernel[kernelIndex(x, y, z, y_size, z_size)] += integral / cell_volume;
				}
			}
		}
	}
	reflectKernel(kernel, x_size, y_size, z_size);
}

} // namespace GreensFunction3D
