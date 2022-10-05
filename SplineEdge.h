#ifndef SPLINEEDGE_H
#define SPLINEEDGE_H
#include <vector>
#include <cmath>
#include <map>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>

class SplineEdge {
public:
    SplineEdge(std::string filename);
    std::vector<double> Evaluate(double x, int d = 0);
    std::vector<double> findx(double s);
    double finds(double x, int d);
    std::vector<double> finds(double x, int d, std::vector<int> &index);
    double m_ArcLength;
private:
    std::vector<double> EvaluateTau(double tau);
    void calculateArcTable();
    void appendArcLength(std::vector<double> p, double s, double tau);
    std::vector<std::vector<double>> m_arc;
    std::vector<std::vector<double> > m_pts;
    std::vector<boost::math::interpolators::cardinal_cubic_b_spline<double> > m_spline;
};

#endif
