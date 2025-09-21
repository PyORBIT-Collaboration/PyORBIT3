#include "DanilovEnvelopeTracker.hh"

DanilovEnvelopeTracker::DanilovEnvelopeTracker(double perveance) : CppPyWrapper(NULL) {
    Q = perveance;
}

void DanilovEnvelopeTracker::setPerveance(double perveance) { 
    Q = perveance;
}

double DanilovEnvelopeTracker::getPerveance() { 
    return Q;
}

void DanilovEnvelopeTracker::trackBunch(Bunch *bunch, double length) {
    // Compute ellipse size and orientation.
    // -----------------------------------------------------------------------------------

    // Get envelope parameters from first two bunch particles.
    double a = bunch->x(0);
    double b = bunch->x(1);
    double e = bunch->y(0);
    double f = bunch->y(1);

    // Compute covariance matrix in x-y plane
    double cov_xx = a * a + b * b; // 4 * <xx>
    double cov_yy = e * e + f * f; // 4 * <yy>
    double cov_xy = a * e + b * f; // 4 * <xy>

    // Compute tilt angle phi in x-y plane (below x axis).
    double phi = -0.5 * atan2(2.0 * cov_xy, cov_xx - cov_yy);

    // Store sin(phi) and cos(phi).
    double cs = cos(phi);
    double sn = sin(phi);
    double cs2 = cs * cs;
    double sn2 = sn * sn;
    double sn_cs = sn * cs;

    // Compute ellipse radii cu and cv in upright frame.
    double cu = sqrt(abs(cov_xx * cs2 + cov_yy * sn2 - 2.0 * cov_xy * sn_cs));
    double cv = sqrt(abs(cov_xx * sn2 + cov_yy * cs2 + 2.0 * cov_xy * sn_cs));
    double cu2 = cu * cu;
    double cv2 = cv * cv;

    // Kick envelope
    // -----------------------------------------------------------------------------------

    double sc_term = (2.0 * Q / (cu + cv));
    if (cu > 0.0) {
        bunch->xp(0) += length * sc_term * (a * cs2 - e * sn_cs) / cu;
        bunch->xp(1) += length * sc_term * (b * cs2 - f * sn_cs) / cu;
        bunch->yp(0) += length * sc_term * (e * sn2 - a * sn_cs) / cu;
        bunch->yp(1) += length * sc_term * (f * sn2 - b * sn_cs) / cu;
    }
    if (cv > 0.0) {
        bunch->xp(0) += length * sc_term * (a * sn2 + e * sn_cs) / cv;
        bunch->xp(1) += length * sc_term * (b * sn2 + f * sn_cs) / cv;
        bunch->yp(0) += length * sc_term * (e * cs2 + a * sn_cs) / cv;
        bunch->yp(1) += length * sc_term * (f * cs2 + b * sn_cs) / cv;
    }

    // Kick particles
    // -----------------------------------------------------------------------------------

    for (int i = 2; i < bunch->getSize(); i++) {
        // Collect particle coordinates.
        double x = bunch->x(i);
        double y = bunch->y(i);
        double u = x * cs - y * sn;
        double v = x * sn + y * cs;

        // Check if particle is inside beam ellipse.
        double u2 = u * u;
        double v2 = v * v;
        bool in_ellipse = ((u2 / cu2) + (v2 / cv2)) <= 1.0;

        // Update momentum.
        double delta_up = 0.0;
        double delta_vp = 0.0;
        
        if (in_ellipse) {
            if (cu > 0.0) {
                delta_up = length * sc_term * u / cu;
            }
            if (cv > 0.0) {
                delta_vp = length * sc_term * v / cv;
            }
        } 
        else {
            // https://arxiv.org/abs/physics/0108040
            double B = u2 + v2 - cu2 - cv2;
            double C = u2 * cv2 + v2 * cu2 - cu2 * cv2;
            double t1 = pow(0.25 * B * B + C, 0.5) + 0.5 * B;
            double Dx = pow(cu2 + t1, 0.5);
            double Dy = pow(cv2 + t1, 0.5);
            delta_up = length * 2.0 * Q * u / (Dx * (Dx + Dy));
            delta_vp = length * 2.0 * Q * v / (Dy * (Dx + Dy));
        }
        
        bunch->xp(i) += (+delta_up * cs + delta_vp * sn);
        bunch->yp(i) += (-delta_up * sn + delta_vp * cs);
    }
}