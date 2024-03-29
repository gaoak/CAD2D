#ifndef CYLINDER_H
#define CYLINDER_H
#include <string>
#include <vector>
class Cylinder {
public:
  Cylinder(double radius, double thetaStart, double thetaEnd);
  std::vector<double> GetPoint(double s);

protected:
  double m_radius;
  double m_thetaStart;
  double m_thetaEnd;
};

#endif
