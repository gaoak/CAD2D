#ifndef SPLINEEDGE_H
#define SPLINEEDGE_H
#include <vector>
#include <cmath>
#include <map>
#include <boost/math/interpolators/cubic_b_spline.hpp>

class SplineEdge {
public:
    SplineEdge(std::string filename, std::map<int, double> params, int dim);
    std::vector<double> EvaluateTau(double s);
    std::vector<double> findx(double s);
    double finds(std::vector<double> x, int d);
    double m_ArcLength;
private:
    void calculateArcTable();
    void appendArcLength(double s, std::vector<double> p);
    std::vector<std::vector<double>> m_arc;
    std::vector<std::vector<double> > m_pts;
    std::vector<boost::math::cubic_b_spline<double> > m_spline;
};

#endif
