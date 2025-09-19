#include "Danilov20EnvelopeSolver.hh"

Danilov20EnvelopeSolver::Danilov20EnvelopeSolver(double perveance, double emittanceX, double emittanceY) : CppPyWrapper(NULL) {
  Q = perveance;
  epsX = emittanceX;
  epsY = emittanceY;
}

void Danilov20EnvelopeSolver::setPerveance(double perveance) {
  Q = perveance;
}

void Danilov20EnvelopeSolver::setEmittanceX(double emittanceX) {
  epsX = emittanceX;
}

void Danilov20EnvelopeSolver::setEmittanceY(double emittanceY) {
  epsY = emittanceY;
}

double Danilov20EnvelopeSolver::getPerveance() {
  return Q;
}

double Danilov20EnvelopeSolver::getEmittanceX() {
  return epsX;
}

double Danilov20EnvelopeSolver::getEmittanceY() {
  return epsY;
}

void Danilov20EnvelopeSolver::trackBunch(Bunch *bunch, double length) {
    // Track envelope
    double cx = bunch->x(0);
    double cy = bunch->y(0);
    double sc_term = 2.0 * Q / (cx + cy);
    double emit_term_x = (epsX * epsX) / (cx * cx * cx);
    double emit_term_y = (epsY * epsY) / (cy * cy * cy);
    bunch->xp(0) += length * (sc_term + emit_term_x);
    bunch->yp(0) += length * (sc_term + emit_term_y);

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

    double deltaXP;
    double deltaYP;
    bool inside;

    for (int i = 1; i < bunch->getSize(); i++) {
        x = bunch->x(i);
        y = bunch->y(i);

        x2 = x * x;
        y2 = y * y;

        inside = ((x2 / cx2) + (y2 / cy2)) <= 1.0;

        if (inside) {
            deltaXP = sc_term * x / cx;
            deltaYP = sc_term * y / cy;
        } 
        else {
            // https://arxiv.org/abs/physics/0108040
            B = x2 + y2 - cx2 - cy2;
            C = x2 * cy2 + y2 * cx2 - cx2 * cy2;
            t1 = pow(0.25 * B * B + C, 0.5) + 0.5 * B;
            Dx = pow(cx2 + t1, 0.5);
            Dy = pow(cy2 + t1, 0.5);
            deltaXP = 2.0 * Q * x / (Dx * (Dx + Dy));
            deltaYP = 2.0 * Q * y / (Dy * (Dx + Dy));
        }
        bunch->xp(i) += deltaXP;
        bunch->yp(i) += deltaYP;
    }
  }
