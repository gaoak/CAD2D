#include "airfoil.h"
#include "util.h"
#include <cmath>
#include <iostream>
#include <vector>

NACAmpxx::NACAmpxx(std::map<std::string, double> params) {
  m_m = params["Chamber"];
  m_p = params["MaxChamberX"];
  if (m_m < 0.001) {
    m_m = 0.;
    m_p = 100.;
  }
  m_t = params["Thickness"];
  m_isRoundTrailing = params["RoundTE"] >= 0.5;
  InitAirfoil();
}

NACAmpxx::NACAmpxx(double m, double p, double t, bool roundtrailing) {
  m_m = m;
  m_p = p;
  if (m_m < 0.001) {
    m_m = 0.;
    m_p = 100.;
  }
  m_t = t;
  m_isRoundTrailing = roundtrailing;
  InitAirfoil();
}

NACAmpxx::NACAmpxx(std::string name) {
  if (name.size() < 4) {
    std::cout << "incorrect NACA 4-digit format" << std::endl;
    exit(-1);
  }
  m_m = (name[0] - '0') * 0.01;
  m_p = (name[1] - '0') * 0.1;
  m_t = (name[2] - '0') * 0.1 + (name[3] - '0') * 0.01;
  if (name.size() >= 4 && name[3] == 'r') {
    m_isRoundTrailing = true;
  } else {
    m_isRoundTrailing = false;
  }
  if (m_m < 0.001) {
    m_m = 0.;
    m_p = 100.;
  }
  InitAirfoil();
}

void NACAmpxx::InitAirfoil() {
  m_TERadius = calculateTrailingRadius(m_TETangencyX);
  calculateArcTable();
  printf("NACA airfoil with Area and AreaMomentx (round trailing edge %d) "
         "%26.18f, %26.18f\n",
         m_isRoundTrailing, Area(), AreaMoment());
}

double NACAmpxx::Findx(double s, int surf) {
  if (surf > 0) {
    return Findx(s, m_arcu);
  } else {
    return Findx(s, m_arcd);
  }
}

double NACAmpxx::Findx(double s, std::vector<std::vector<double>> &arc) {
  s = fabs(s);
  if (s >= arc[arc.size() - 1][1])
    return 1.;
  int i = 1;
  while (arc[i][1] < s)
    ++i;
  return arc[i - 1][0] + (s - arc[i - 1][1]) * (arc[i][0] - arc[i - 1][0]) /
                             (arc[i][1] - arc[i - 1][1]);
}

double NACAmpxx::Finds(double x, int surf) {
  if (surf > 0) {
    return Finds(x, m_arcu);
  } else {
    return Finds(x, m_arcd);
  }
}

double NACAmpxx::Finds(double x, std::vector<std::vector<double>> &arc) {
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
         (0.2969 * sqrt(x) * 2. / 3. - 0.1260 * xs[1] / 2. -
          0.3516 * xs[2] / 3. + 0.2843 * xs[3] / 4. - 0.1015 * xs[4] / 5.) *
         x;
}

double NACAmpxx::halfSxt(double x) {
  x = fabs(x);
  std::vector<double> xs(5, 1.);
  for (int i = 1; i < 5; ++i)
    xs[i] = xs[i - 1] * x;
  return 5. * m_t *
         (0.2969 * sqrt(x) * 2. / 5. - 0.1260 * xs[1] / 3. -
          0.3516 * xs[2] / 4. + 0.2843 * xs[3] / 5. - 0.1015 * xs[4] / 6.) *
         xs[2];
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

double NACAmpxx::roundTrailingSize() { return m_TERadius; }

void NACAmpxx::GetInfo(std::map<std::string, double> &p) {
  p["LERadius"] = 0.5 * (5. * m_t * 0.2969) * (5. * m_t * 0.2969);
  p["LEInterX"] = p["LERadius"];
  if (m_isRoundTrailing) {
    p["TERadius"] = m_TERadius;
    p["TEInterX"] = m_TETangencyX;
  } else {
    p["TEThickness"] = halft(1.) * 2.;
  }
  p["Thickness"] = m_t;
  p["Area"] = Area();
  p["Areax"] = AreaMoment();
}

/***
 * remapping the points to a round trailing if the point lies on the airfoil
 * return the original point if it is not on the wall
 ***/
std::vector<double> NACAmpxx::roundTrailingEdge(std::vector<double> &p0,
                                                double eps) {
  std::vector<double> p1;
  bool notexist = fabs(halft(p0[0]) - fabs(p0[1])) > eps;
  if (fabs(1. - p0[0]) < eps && fabs(p0[1]) <= halft(1) + eps)
    notexist = false;
  if (p0[0] <= m_TETangencyX || notexist) {
    p1.push_back(p0[0]);
    p1.push_back(p0[1]);
  } else {
    double r0 = p0[0] - 1. + m_TERadius;
    double r1 = p0[1];
    double r = sqrt(r0 * r0 + r1 * r1);
    p1.push_back(1. + m_TERadius * (r0 / r - 1.));
    p1.push_back(m_TERadius * r1 / r);
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

std::vector<double> NACAmpxx::Upper(double x) {
  std::vector<double> res(2);
  if (x > 1.)
    x = 1.;
  if (x < 0.)
    x = 0.;
  double yt = halft(x);
  double yc = chamber(x);
  double the = theta(x);
  res[0] = x - yt * sin(the);
  if (m_isRoundTrailing && x > m_TETangencyX) {
    res[1] = sqrt((1. - x) * (x - 1. + 2. * m_TERadius));
  } else {
    res[1] = yc + yt * cos(the);
  }
  return res;
}

std::vector<double> NACAmpxx::Lower(double x) {
  std::vector<double> res(2);
  if (x > 1.)
    x = 1.;
  if (x < 0.)
    x = 0.;
  double yt = halft(x);
  double yc = chamber(x);
  double the = theta(x);
  res[0] = x + yt * sin(the);
  if (m_isRoundTrailing && x > m_TETangencyX) {
    res[1] = -sqrt((1. - x) * (x - 1. + 2. * m_TERadius));
  } else {
    res[1] = yc - yt * cos(the);
  }
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
    std::vector<double> p1 = Upper(x);
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
    std::vector<double> p1 = Upper(x);
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
    std::vector<double> p1 = Lower(x);
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
    std::vector<double> p1 = Lower(x);
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

double NACAmpxx::Area() {
  double res = 0.;
  if (m_isRoundTrailing) {
    double xinterc;
    double rt = calculateTrailingRadius(xinterc);
    double theta = acos((rt - 1. + xinterc) / rt);
    res = halfSt(xinterc) * 2. + rt * rt * (theta - cos(theta) * sin(theta));
  } else {
    res = halfSt(1.) * 2.;
  }
  return res;
}

double NACAmpxx::AreaMoment() {
  double res = 0.;
  if (m_isRoundTrailing) {
    double xinterc;
    double rt = calculateTrailingRadius(xinterc);
    double theta = acos((rt - 1. + xinterc) / rt);
    res = halfSxt(xinterc) * 2. +
          rt * rt * (1. - rt) * (theta - cos(theta) * sin(theta)) +
          2. / 3. * pow(rt * sin(theta), 3);
  } else {
    res = halfSxt(1.) * 2.;
  }
  return res;
}

WedgeFoil::WedgeFoil(std::map<std::string, double> params) {
  m_LEDiameter = params["Thickness"];
  m_TEThich = params["TEThickness"];
  m_isRoundTrailing = params["RoundTE"] >= 0.5;
  InitAirfoil();
}

WedgeFoil::WedgeFoil(double D0, double D1, double D2, bool roundtrailing) {
  m_LEDiameter = D0;
  m_TEThich = D1;
  m_isRoundTrailing = roundtrailing;
  InitAirfoil();
}

void WedgeFoil::InitAirfoil() {
  double det = m_TEThich * m_TEThich + 4. * 1. * (1. - m_LEDiameter);
  m_theta = 2. * atan(0.5 * (m_TEThich + sqrt(det)));
  m_LETangencyX = 0.5 * m_LEDiameter * (1. + cos(m_theta));
  m_TERadius = roundTrailingSize();
  m_TETangencyX = 1. - m_TERadius * (1 - cos(m_theta));
  printf("Wedge-shape airfoil with Area (round trailing edge %d) %26.18f\n",
         m_isRoundTrailing, Area());
}

void WedgeFoil::GetInfo(std::map<std::string, double> &p) {
  p["LERadius"] = m_LEDiameter * 0.5;
  p["LEInterX"] = m_LETangencyX;
  if (m_isRoundTrailing) {
    p["TERadius"] = m_TERadius;
    p["TEInterX"] = m_TETangencyX;
  } else {
    p["TEThickness"] = m_TEThich;
  }
  p["Thickness"] = m_LEDiameter;
  p["Area"] = Area();
}

std::vector<double> WedgeFoil::Upper(double x) {
  double radius = 0.5 * m_LEDiameter;
  std::vector<double> res(2);
  res[0] = x;
  if (x < m_LETangencyX) {
    res[1] = sqrt(radius * radius - (x - radius) * (x - radius));
  } else {
    if (m_isRoundTrailing && x > m_TETangencyX) {
      res[1] = sqrt((1. - x) * (x - 1. + 2. * m_TERadius));
    } else {
      res[1] = 0.5 * m_TEThich - (x - 1.) / tan(m_theta);
    }
  }
  return res;
}

std::vector<double> WedgeFoil::Lower(double x) {
  double radius = 0.5 * m_LEDiameter;
  std::vector<double> res(2);
  res[0] = x;
  if (x < m_LETangencyX) {
    res[1] = -sqrt(radius * radius - (x - radius) * (x - radius));
  } else {
    if (m_isRoundTrailing && x > m_TETangencyX) {
      res[1] = -sqrt((1. - x) * (x - 1. + 2. * m_TERadius));
    } else {
      res[1] = -0.5 * m_TEThich + (x - 1.) / tan(m_theta);
    }
  }
  return res;
}

double WedgeFoil::Findx(double s, int up) {
  double radius = 0.5 * m_LEDiameter;
  double x, sinter = radius * (M_PI - m_theta);
  if (s < sinter) {
    x = radius * (1 - cos(s / radius));
  } else {
    x = m_LETangencyX + (s - sinter) * sin(m_theta);
  }
  return x;
}

double WedgeFoil::Finds(double x, int up) {
  double radius = 0.5 * m_LEDiameter;
  double s;
  if (x < m_LETangencyX) {
    s = radius * (M_PI - acos((x - radius) / radius));
  } else {
    s = radius * m_theta + (x - m_LETangencyX) / sin(m_theta);
  }
  return s;
}

/***
 * remapping the points to a round trailing if the point lies on the airfoil
 * return the original point if it is not on the wall
 ***/
std::vector<double> WedgeFoil::roundTrailingEdge(std::vector<double> &p0,
                                                 double eps) {
  std::vector<double> res = p0;
  if (p0[0] > m_TETangencyX && p0[0] <= 1. + eps) {
    bool onwall = fabs(Upper(p0[0])[1] - fabs(p0[1])) <= eps;
    if (fabs(1. - p0[0]) < eps && fabs(p0[1]) <= Upper(1.)[1] + eps)
      onwall = true;
    if (onwall) {
      double theta = atan(p0[1] / (p0[0] - 1. + m_TERadius));
      res[0] = 1. - m_TERadius + m_TERadius * cos(theta);
      res[1] = m_TERadius * sin(theta);
    }
  }
  return res;
}

double WedgeFoil::roundTrailingSize() {
  return 0.5 * m_TEThich / tan(0.5 * m_theta);
}

double WedgeFoil::Area() {
  double rl = 0.5 * m_LEDiameter;
  double rt = roundTrailingSize();
  double res = (M_PI - m_theta) * rl * rl;
  if (m_isRoundTrailing) {
    res += (rl + rt) * (1. - rl - rt) * sin(m_theta) + m_theta * rt * rt;
  } else {
    res +=
        rl * (1. - m_LETangencyX) / sin(m_theta) + 0.5 * m_TEThich * (1 - rl);
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