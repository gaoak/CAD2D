#ifndef SPLINEEDGE_H
#define SPLINEEDGE_H
#include <vector>
#include <cmath>
#include <map>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>

enum ArchIndex {
    eX,
    eY,
    eS,
    eTau
};

class SplineEdge {
public:
    SplineEdge(std::string filename);
    std::vector<double> Evaluate(double s);
    std::vector<double> Evaluatex(double x, int d = 0);
    std::vector<double> Evaluates(double s);
    double LookFors(double x, int d);
    double m_ArcLength;
private:
    std::vector<double> LookFor(double x, int d, std::vector<int> &index);
    std::vector<double> EvaluateTau(double tau);
    void calculateArcTable();
    void appendArcLength(std::vector<double> p, double s, double tau);
    std::vector<std::vector<double>> m_arc;
    std::vector<std::vector<double> > m_pts;
    std::vector<boost::math::interpolators::cardinal_cubic_b_spline<double> > m_spline;
};

#endif
