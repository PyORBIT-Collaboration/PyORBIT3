#include "KVEnvelopeTracker.hh"

KVEnvelopeTracker::KVEnvelopeTracker(double perveance, double emittance_x, double emittance_y): CppPyWrapper(NULL) {
    Q = perveance;
    eps_x = emittance_x;
    eps_y = emittance_y;
}

void KVEnvelopeTracker::setPerveance(double perveance) { 
    Q = perveance; 
}

void KVEnvelopeTracker::setEmittanceX(double emittance) { 
    eps_x = emittance; 
}

void KVEnvelopeTracker::setEmittanceY(double emittance) { 
    eps_y = emittance; 
}

double KVEnvelopeTracker::getPerveance() { 
    return Q; 
}

double KVEnvelopeTracker::getEmittanceX() { 
    return eps_x; 
}

double KVEnvelopeTracker::getEmittanceY() { 
    return eps_y; 
}

void KVEnvelopeTracker::trackBunch(Bunch *bunch, double length) {
    // Get envelope parameters
    double cx = bunch->x(0);
    double cy = bunch->y(0);
    double cx2 = cx * cx;
    double cy2 = cy * cy;

    // Kick envelope
    double eps_x_term = (eps_x * eps_x) / (cx * cx * cx);
    double eps_y_term = (eps_y * eps_y) / (cy * cy * cy);
    double sc_term = 2.0 * Q / (cx + cy);
    bunch->xp(0) += length * (sc_term + eps_x_term);
    bunch->yp(0) += length * (sc_term + eps_y_term);

    // Kick particles
    for (int i = 1; i < bunch->getSize(); i++) {
        // Check if particle is inside beam ellipse.
        double x = bunch->x(i);
        double y = bunch->y(i);
        double x2 = x * x;
        double y2 = y * y;
        bool in_ellipse = ((x2 / cx2) + (y2 / cy2)) <= 1.0;

        // Update momentum coordinates
        if (in_ellipse) {
            bunch->xp(i) += length * sc_term * (x / cx);
            bunch->yp(i) += length * sc_term * (y / cy);
        } 
        else {
            // https://arxiv.org/abs/physics/0108040
            double B = x2 + y2 - cx2 - cy2;
            double C = x2 * cy2 + y2 * cx2 - cx2 * cy2;
            double t1 = pow(0.25 * B * B + C, 0.5) + 0.5 * B;
            double Dx = pow(cx2 + t1, 0.5);
            double Dy = pow(cy2 + t1, 0.5);
            bunch->xp(i) += length * (2.0 * Q * x / (Dx * (Dx + Dy)));
            bunch->yp(i) += length * (2.0 * Q * y / (Dy * (Dx + Dy)));
        }
    }
}