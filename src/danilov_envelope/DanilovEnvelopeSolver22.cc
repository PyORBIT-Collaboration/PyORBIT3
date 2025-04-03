#include "DanilovEnvelopeSolver22.hh"

DanilovEnvelopeSolver22::DanilovEnvelopeSolver22(double perveance) : CppPyWrapper(NULL) {
  Q = perveance;
}

void DanilovEnvelopeSolver22::setPerveance(double perveance) {
  Q = perveance; 
}

double DanilovEnvelopeSolver22::getPerveance() {
  return Q;
}

void DanilovEnvelopeSolver22::trackBunch(Bunch *bunch, double length) {
  // Compute ellipse size and orientation.
  double a = bunch->x(0);
  double b = bunch->x(1);
  double e = bunch->y(0);
  double f = bunch->y(1);

  double cov_xx = a * a + b * b; // 4 * <x^2>
  double cov_yy = e * e + f * f; // 4 * <y^2>
  double cov_xy = a * e + b * f; // 4 * <xy>

  double phi = -0.5 * atan2(2.0 * cov_xy, cov_xx - cov_yy);

  double _cos = cos(phi);
  double _sin = sin(phi);
  double cos2 = _cos * _cos;
  double sin2 = _sin * _sin;
  double sin_cos = _sin * _cos;

  double cxn = sqrt(abs(cov_xx * cos2 + cov_yy * sin2 - 2.0 * cov_xy * sin_cos));
  double cyn = sqrt(abs(cov_xx * sin2 + cov_yy * cos2 + 2.0 * cov_xy * sin_cos));
  double factor = length * (2.0 * Q / (cxn + cyn));

  // Track envelope
  if (cxn > 0.0) {
    bunch->xp(0) += factor * (a * cos2 - e * sin_cos) / cxn;
    bunch->xp(1) += factor * (b * cos2 - f * sin_cos) / cxn;
    bunch->yp(0) += factor * (e * sin2 - a * sin_cos) / cxn;
    bunch->yp(1) += factor * (f * sin2 - b * sin_cos) / cxn;
  }
  if (cyn > 0.0) {
    bunch->xp(0) += factor * (a * sin2 + e * sin_cos) / cyn;
    bunch->xp(1) += factor * (b * sin2 + f * sin_cos) / cyn;
    bunch->yp(0) += factor * (e * cos2 + a * sin_cos) / cyn;
    bunch->yp(1) += factor * (f * cos2 + b * sin_cos) / cyn;
  }

  // Track test particles
  double cxn2 = cxn * cxn;
  double cyn2 = cyn * cyn;
  double x;
  double y;
  double x2;
  double y2;

  double xn;
  double yn;
  double xn2;
  double yn2;

  double t1;
  double B;
  double C;
  double Dx;
  double Dy;

  double delta_xp;
  double delta_yp;
  double delta_xpn;
  double delta_ypn;
  bool in_ellipse;

  for (int i = 2; i < bunch->getSize(); i++) {
    x = bunch->x(i);
    y = bunch->y(i);

    x2 = x * x;
    y2 = y * y;

    xn = x * _cos - y * _sin;
    yn = x * _sin + y * _cos;

    xn2 = xn * xn;
    yn2 = yn * yn;

    in_ellipse = ((xn2 / cxn2) + (yn2 / cyn2)) <= 1.0;

    if (in_ellipse) {
      if (cxn > 0.0) {
        delta_xpn = factor * xn / cxn;
      }
      if (cyn > 0.0) {
        delta_ypn = factor * yn / cyn;
      }
    } else {
      // https://arxiv.org/abs/physics/0108040
      B = xn2 + yn2 - cxn2 - cyn2;
      C = xn2 * cyn2 + yn2 * cxn2 - cxn2 * cyn2;
      t1 = pow(0.25 * B * B + C, 0.5) + 0.5 * B;
      Dx = pow(cxn2 + t1, 0.5);
      Dy = pow(cyn2 + t1, 0.5);
      delta_xpn = 2.0 * Q * xn / (Dx * (Dx + Dy));
      delta_ypn = 2.0 * Q * yn / (Dy * (Dx + Dy));
    }
    delta_xp = +delta_xpn * _cos + delta_ypn * _sin;
    delta_yp = -delta_xpn * _sin + delta_ypn * _cos;
    bunch->xp(i) += delta_xp;
    bunch->yp(i) += delta_yp;
  }
}
