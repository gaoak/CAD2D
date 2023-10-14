#include "airfoil.h"
#include "util.h"
#include <cmath>
#include <iostream>
#include <vector>

NACAmpxx::NACAmpxx(std::vector<double> &params) {
  if(params.size()<3) {
    params.push_back(0.);
    params.push_back(0.);
    params.push_back(0.);
  }
  m_m = params[0];
  m_p = params[1];
  if (m_m < 0.001) {
    m_m = 0.;
    m_p = 100.;
  }
  m_t = params[2];
  calculateArcTable();
  m_rRoundTrailing = -1.;
}

NACAmpxx::NACAmpxx(double m, double p, double t) {
  m_m = m;
  m_p = p;
  if (m_m < 0.001) {
    m_m = 0.;
    m_p = 100.;
  }
  m_t = t;
  calculateArcTable();
  m_rRoundTrailing = -1.;
}

NACAmpxx::NACAmpxx(std::string name) {
  if (name.size() < 4) {
    std::cout << "incorrect NACA 4-digit format" << std::endl;
    exit(-1);
  }
  m_m = (name[0] - '0') * 0.01;
  m_p = (name[1] - '0') * 0.1;
  m_t = (name[2] - '0') * 0.1 + (name[3] - '0') * 0.01;
  if (m_m < 0.001) {
    m_m = 0.;
    m_p = 100.;
  }
  calculateArcTable();
  m_rRoundTrailing = -1.;
}

double NACAmpxx::findx(double s, int surf) {
  if (surf > 0) {
    return findx(s, m_arcu);
  } else {
    return findx(s, m_arcd);
  }
}

double NACAmpxx::findx(double s, std::vector<std::vector<double>> &arc) {
  s = fabs(s);
  if (s >= arc[arc.size() - 1][1])
    return 1.;
  int i = 1;
  while (arc[i][1] < s)
    ++i;
  return arc[i - 1][0] + (s - arc[i - 1][1]) * (arc[i][0] - arc[i - 1][0]) /
                             (arc[i][1] - arc[i - 1][1]);
}

double NACAmpxx::finds(double x, int surf) {
  if (surf > 0) {
    return finds(x, m_arcu);
  } else {
    return finds(x, m_arcd);
  }
}

double NACAmpxx::finds(double x, std::vector<std::vector<double>> &arc) {
  x = fabs(x);
  if (x >= 1.)
    return arc[arc.size() - 1][1];
  int i = 1;
  while (arc[i][0] < x)
    ++i;
  return arc[i - 1][1] + (x - arc[i - 1][0]) * (arc[i][1] - arc[i - 1][1]) /
                             (arc[i][0] - arc[i - 1][0]);
}

double NACAmpxx::halft(double x) {
  x = fabs(x);
  std::vector<double> xs(5, 1.);
  for (int i = 1; i < 5; ++i)
    xs[i] = xs[i - 1] * x;
  return 5. * m_t *
         (0.2969 * sqrt(x) - 0.1260 * xs[1] - 0.3516 * xs[2] + 0.2843 * xs[3] -
          0.1015 * xs[4]);
}

double NACAmpxx::halfSt(double x) {
  x = fabs(x);
  std::vector<double> xs(5, 1.);
  for (int i = 1; i < 5; ++i)
    xs[i] = xs[i - 1] * x;
  return 5. * m_t *
         (0.2969 * sqrt(x) * 2./3. - 0.1260 * xs[1] / 2. - 0.3516 * xs[2] / 3. + 0.2843 * xs[3] / 4. -
          0.1015 * xs[4] / 5.) * x;
}

double NACAmpxx::halfdt(double x) {
  x = fabs(x);
  std::vector<double> xs(5, 1.);
  for (int i = 1; i < 5; ++i)
    xs[i] = xs[i - 1] * x;
  return 5. * m_t *
         (0.2969 * 0.5 / sqrt(x) - 0.1260 - 0.3516 * 2. * xs[1] +
          0.2843 * 3. * xs[2] - 0.1015 * 4. * xs[3]);
}

double NACAmpxx::roundTrailingSize() {
  double x;
  calculateTrailingRadius(x);
  return 1. - x;
}

/***
 * remapping the points to a round trailing if the point lies on the airfoil
 * return the original point if it is not on the wall
***/
std::vector<double> NACAmpxx::roundTrailingEdge(std::vector<double> &p0,
                                                double eps) {
  if (m_rRoundTrailing < 0.) {
    m_rRoundTrailing = calculateTrailingRadius(m_xRoundTrailing);
  }
  std::vector<double> p1;
  bool notexist = fabs(halft(p0[0]) - fabs(p0[1])) > eps;
  if (fabs(1. - p0[0]) < eps && fabs(p0[1]) <= halft(1) + eps)
    notexist = false;
  if (m_rRoundTrailing < 0. || p0[0] <= m_xRoundTrailing || notexist) {
    p1.push_back(p0[0]);
    p1.push_back(p0[1]);
  } else {
    double r0 = p0[0] - 1. + m_rRoundTrailing;
    double r1 = p0[1];
    double r = sqrt(r0 * r0 + r1 * r1);
    p1.push_back(1. + m_rRoundTrailing * (r0 / r - 1.));
    p1.push_back(m_rRoundTrailing * r1 / r);
  }
  return p1;
}

double NACAmpxx::calculateTrailingRadius(double &xtmp) {
  if (m_m > 1.E-6) {
    return -1.;
  }
  double eps = 1.E-14, rtmp, ftmp;
  double x[2];
  double f[2];
  x[0] = 0.9;
  x[1] = 0.9999;
  xtmp = x[1];
  f[0] = testRadius(x[0], rtmp);
  f[1] = testRadius(xtmp, rtmp);
  while (fabs(f[1]) > eps) {
    xtmp = (x[0] * f[1] - x[1] * f[0]) / (f[1] - f[0]);
    ftmp = testRadius(xtmp, rtmp);
    x[0] = x[1];
    f[0] = f[1];
    x[1] = xtmp;
    f[1] = ftmp;
  }
  return rtmp;
}

double NACAmpxx::testRadius(double x, double &r) {
  double y = halft(x);
  double t = halfdt(x);
  double xc = x + t * y;
  r = 1. - xc;
  return sqrt((x - xc) * (x - xc) + y * y) - r;
}

double NACAmpxx::chamber(double x) {
  double yc;
  if (x <= m_p) {
    yc = m_m / (m_p * m_p) * (2. * m_p * x - x * x);
  } else {
    yc = m_m / (1. - m_p) / (1. - m_p) *
         ((1. - 2. * m_p) + 2. * m_p * x - x * x);
  }
  return yc;
}

std::vector<double> NACAmpxx::up(double x) {
  std::vector<double> res(2);
  if (x > 1.)
    x = 1.;
  if (x < 0.)
    x = 0.;
  double yt = halft(x);
  double yc = chamber(x);
  double the = theta(x);
  res[0] = x - yt * sin(the);
  res[1] = yc + yt * cos(the);
  return res;
}

std::vector<double> NACAmpxx::down(double x) {
  std::vector<double> res(2);
  if (x > 1.)
    x = 1.;
  if (x < 0.)
    x = 0.;
  double yt = halft(x);
  double yc = chamber(x);
  double the = theta(x);
  res[0] = x + yt * sin(the);
  res[1] = yc - yt * cos(the);
  return res;
}

double NACAmpxx::theta(double x) {
  double dyc;
  if (x <= m_p) {
    dyc = 2. * m_m / (m_p * m_p) * (m_p - x);
  } else {
    dyc = 2. * m_m / (1. - m_p) / (1. - m_p) * (m_p - x);
  }
  return atan(dyc);
}

/***
 * Construct a table of (x, s) for both upper and lower surfaces
 * s starts from 0 (leading edge) to the trailing edge
***/
void NACAmpxx::calculateArcTable() {
  double xmid = 0.1;
  int N0 = 1000;
  int N1 = 900;
  std::vector<double> tmp(2, 0.);
  std::vector<double> p0(2, 0.);
  m_arcu.push_back(tmp);
  m_arcd.push_back(tmp);
  //////upper surface
  double x, dx = xmid / N0 / 10., arcl = 0.;
  for (int i = 1; i <= 10 * N0; ++i) {
    x = dx * i;
    std::vector<double> p1 = up(x);
    arcl += distance(p0, p1);
    if (i % 10 == 0) {
      std::vector<double> tmp1(2);
      tmp1[0] = x;
      tmp1[1] = arcl;
      m_arcu.push_back(tmp1);
    }
    p0 = p1;
  }
  dx = (1. - xmid) / N1;
  for (int i = 1; i <= N1; ++i) {
    x = dx * i + xmid;
    std::vector<double> p1 = up(x);
    arcl += distance(p0, p1);
    if (1) {
      std::vector<double> tmp1(2);
      tmp1[0] = x;
      tmp1[1] = arcl;
      m_arcu.push_back(tmp1);
    }
    p0 = p1;
  }
  //////lower surface
  dx = xmid / N0 / 10.;
  arcl = 0.;
  p0[0] = 0.;
  p0[1] = 0.;
  for (int i = 1; i <= 10 * N0; ++i) {
    x = dx * i;
    std::vector<double> p1 = down(x);
    arcl += distance(p0, p1);
    if (i % 10 == 0) {
      std::vector<double> tmp1(2);
      tmp1[0] = x;
      tmp1[1] = arcl;
      m_arcd.push_back(tmp1);
    }
    p0 = p1;
  }
  dx = (1. - xmid) / N1;
  for (int i = 1; i <= N1; ++i) {
    x = dx * i + xmid;
    std::vector<double> p1 = down(x);
    arcl += distance(p0, p1);
    if (1) {
      std::vector<double> tmp1(2);
      tmp1[0] = x;
      tmp1[1] = arcl;
      m_arcd.push_back(tmp1);
    }
    p0 = p1;
  }
}

double NACAmpxx::Area(bool roundTrailingEdge) {
  double res = 0.;
  if(roundTrailingEdge) {
    double xinterc;
    double rt = calculateTrailingRadius(xinterc);
    double theta = acos((rt - 1. + xinterc)/rt);
    res = halfSt(xinterc) * 2. + rt * rt * theta;
  } else {
    res = halfSt(1.) * 2.;
  }
  return res;
}

WedgeFoil::WedgeFoil(std::vector<double> &params) {
  m_LEDiameter = params[0];
  m_TEThich = params[1];
  double det = m_TEThich * m_TEThich + 4. * 1. * (1. - m_LEDiameter);
  m_theta = 2. * atan(0.5*(m_TEThich + sqrt(det)));
  m_LETangencyX = 0.5 * m_LEDiameter*(1. + cos(m_theta));
}

WedgeFoil::WedgeFoil(double D0, double D1, double D2) {
  m_LEDiameter = D0;
  m_TEThich = D1;
  double det = m_TEThich * m_TEThich + 4. * 1. * (1. - m_LEDiameter);
  m_theta = 2. * atan(0.5*(m_TEThich + sqrt(det)));
  m_LETangencyX = 0.5 * m_LEDiameter*(1. + cos(m_theta));
}

std::vector<double> WedgeFoil::up(double x) {
  double radius = 0.5 * m_LEDiameter;
  std::vector<double> res(2);
  res[0] = x;
  if(x<m_LETangencyX) {
    res[1] = sqrt(radius * radius - (x - radius)*(x - radius));
  } else {
    res[1] = 0.5 * m_TEThich - (x - 1.)/tan(m_theta);
  }
  return res;
}

std::vector<double> WedgeFoil::down(double x) {
  double radius = 0.5 * m_LEDiameter;
  std::vector<double> res(2);
  res[0] = x;
  if(x<m_LETangencyX) {
    res[1] = -sqrt(radius * radius - (x - radius)*(x - radius));
  } else {
    res[1] = -0.5 * m_TEThich + (x - 1.)/tan(m_theta);
  }
  return res;
}

double WedgeFoil::findx(double s, int up) {
  double radius = 0.5 * m_LEDiameter;
  double x, sinter = radius * (M_PI - m_theta);
  if(s < sinter) {
    x = radius * (1 - cos(s / radius));
  } else {
    x = m_LETangencyX + (s-sinter) * sin(m_theta);
  }
  return x;
}

double WedgeFoil::finds(double x, int up) {
  double radius = 0.5 * m_LEDiameter;
  double s;
  if(x<m_LETangencyX) {
    s = radius * (M_PI - acos((x - radius)/radius));
  } else {
    s = radius * m_theta + (x - m_LETangencyX)/sin(m_theta);
  }
  return s;
}

/***
 * remapping the points to a round trailing if the point lies on the airfoil
 * return the original point if it is not on the wall
***/
std::vector<double> WedgeFoil::roundTrailingEdge(std::vector<double> &p0,
                                        double eps) {
  double radius = roundTrailingSize();
  double xintercept = 1. - radius * (1 - cos(m_theta)), xcenter = 1. - radius;
  std::vector<double> res = p0;
  if(p0[0] > xintercept && p0[0] <= 1.+eps) {
    bool onwall = fabs(up(p0[0])[1] - fabs(p0[1])) <= eps;
    if (fabs(1. - p0[0]) < eps && fabs(p0[1]) <= up(1.)[1] + eps) onwall = true;
    if(onwall) {
      double theta  = atan(p0[1] / (p0[0] - xcenter));
      res[0] = xcenter + radius * cos(theta);
      res[1] = radius * sin(theta);
    }
  }
  return res;
}

double WedgeFoil::roundTrailingSize() {
  return 0.5 * m_TEThich / tan(0.5 * m_theta);
}

double WedgeFoil::Area(bool roundTrailingEdge) {
  double rl = 0.5 * m_LEDiameter;
  double rt = roundTrailingSize();
  double res = (M_PI - m_theta) * rl * rl;
  if(roundTrailingEdge) {
    res += (rl+rt)*(1.-rl-rt)*sin(m_theta) + m_theta * rt * rt;
  } else {
    res += rl * (1. - m_LETangencyX) / sin(m_theta) + 0.5 * m_TEThich * (1 - rl);
  }
  return res;
}

/*int main() {
    NACAmpxx airf("0012");
    double x;
    printf("%20.12f, %20.12f, %20.12f, %20.12f\n", x,
airf.calculateTrailingRadius(x), 1.-airf.calculateTrailingRadius(x), x);
}*/

/*int main() {
    WedgeFoil airf(0.1, 0.01, 0.);
    for(double x = 0.; x<=1.; x+=0.01) {
      std::vector<double> p = airf.up(x);
      printf("%20.12f %20.12f\n", p[0], p[1]);
    }
    printf("Area %20.12f\n", airf.Area(true));
}*/