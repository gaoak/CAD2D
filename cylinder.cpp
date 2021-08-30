#include<cmath>
#include<vector>
#include<iostream>
#include"cylinder.h"

Cylinder::Cylinder(double radius, double thetaStart, double thetaEnd)
    :m_radius(radius),m_thetaStart(thetaStart),m_thetaEnd(thetaEnd) {
}

std::vector<double> Cylinder::GetPoint(double s) {
    double t = 0.5*(-s+1.)*m_thetaStart + 0.5*(s+1.)*m_thetaEnd;
    std::vector<double> res(2);
    res[0] = m_radius * cos(t);
    res[1] = m_radius * sin(t);
    return res;
}