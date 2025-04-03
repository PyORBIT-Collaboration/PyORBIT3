#include "DanilovEnvelopeSolver20.hh"

DanilovEnvelopeSolver20::DanilovEnvelopeSolver20(double perveance, double eps_x, double eps_y) : CppPyWrapper(NULL) {
  _perveance = perveance;
  _eps_x = eps_x;
  _eps_y = eps_y;
}

void DanilovEnvelopeSolver20::setPerveance(double perveance) {
  _perveance = perveance;
}

void DanilovEnvelopeSolver20::setEmittanceX(double eps_x) {
  _eps_x = eps_x;
}

void DanilovEnvelopeSolver20::setEmittanceY(double eps_y) {
  _eps_y = eps_y;
}

double DanilovEnvelopeSolver20::getPerveance() {
  return _perveance;
}

double DanilovEnvelopeSolver20::getEmittanceX() {
  return _eps_x;
}

double DanilovEnvelopeSolver20::getEmittanceY() {
  return _eps_y;
}

void DanilovEnvelopeSolver20::trackBunch(Bunch *bunch, double length) {
  // Track envelope
  double cx = bunch->x(0);
  double cy = bunch->y(0);
  double sc_term = 2.0 * _perveance / (cx + cy);
  double emit_term_x = (_eps_x * _eps_x) / (cx * cx * cx);
  double emit_term_y = (_eps_y * _eps_y) / (cy * cy * cy);
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

  double delta_xp;
  double delta_yp;
  bool in_ellipse;

  for (int i = 1; i < bunch->getSize(); i++) {
    x = bunch->x(i);
    y = bunch->y(i);

    x2 = x * x;
    y2 = y * y;

    in_ellipse = ((x2 / cx2) + (y2 / cy2)) <= 1.0;

    if (in_ellipse) {
      delta_xp = sc_term * x / cx;
      delta_yp = sc_term * y / cy;
    } else {
      // https://arxiv.org/abs/physics/0108040
      B = x2 + y2 - cx2 - cy2;
      C = x2 * cy2 + y2 * cx2 - cx2 * cy2;
      t1 = pow(0.25 * B * B + C, 0.5) + 0.5 * B;
      Dx = pow(cx2 + t1, 0.5);
      Dy = pow(cy2 + t1, 0.5);
      delta_xp = 2.0 * _perveance * x / (Dx * (Dx + Dy));
      delta_yp = 2.0 * _perveance * y / (Dy * (Dx + Dy));
    }
    bunch->xp(i) += delta_xp;
    bunch->yp(i) += delta_yp;
  }
}
