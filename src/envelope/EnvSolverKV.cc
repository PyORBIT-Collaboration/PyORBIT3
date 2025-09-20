#include "EnvSolverKV.hh"

EnvSolverKV::EnvSolverKV(double perveance, double emittance_x, double emittance_y) : CppPyWrapper(NULL) {
  Q = perveance;
  eps_x = emittance_x;
  eps_y = emittance_y;
}

void EnvSolverKV::setPerveance(double perveance) {
  Q = perveance;
}

void EnvSolverKV::setEmittanceX(double emittance) {
  eps_x = emittance;
}

void EnvSolverKV::setEmittanceY(double emittance) {
  eps_y = emittance;
}

double EnvSolverKV::getPerveance() {
  return Q;
}

double EnvSolverKV::getEmittanceX() {
  return eps_x;
}

double EnvSolverKV::getEmittanceY() {
  return eps_y;
}

void EnvSolverKV::trackBunch(Bunch *bunch, double length) {
    // Track envelope
    double cx = bunch->x(0);
    double cy = bunch->y(0);
    double sc_term = 2.0 * Q / (cx + cy);
    double eps_x_term = (eps_x * eps_x) / (cx * cx * cx);
    double eps_y_term = (eps_y * eps_y) / (cy * cy * cy);
    bunch->xp(0) += length * (sc_term + eps_x_term);
    bunch->yp(0) += length * (sc_term + eps_y_term);

    // Track test particles
    double cx2 = cx * cx;
    double cy2 = cy * cy;

    double x;
    double y;
    double x2;
    double y2;

    double B;
    double C;
    double Dx;
    double Dy;
    double t1;

    double delta_xp;
    double delta_yp;
    bool inside;

    for (int i = 1; i < bunch->getSize(); i++) {
        x = bunch->x(i);
        y = bunch->y(i);

        x2 = x * x;
        y2 = y * y;

        inside = ((x2 / cx2) + (y2 / cy2)) <= 1.0;

        if (inside) {
            delta_xp = sc_term * x / cx;
            delta_yp = sc_term * y / cy;
        } 
        else {
            // https://arxiv.org/abs/physics/0108040
            // UNTESTED!
            B = x2 + y2 - cx2 - cy2;
            C = x2 * cy2 + y2 * cx2 - cx2 * cy2;
            t1 = pow(0.25 * B * B + C, 0.5) + 0.5 * B;
            Dx = pow(cx2 + t1, 0.5);
            Dy = pow(cy2 + t1, 0.5);
            delta_xp = 2.0 * Q * x / (Dx * (Dx + Dy));
            delta_yp = 2.0 * Q * y / (Dy * (Dx + Dy));
        }
        bunch->xp(i) += delta_xp;
        bunch->yp(i) += delta_yp;
    }
  }