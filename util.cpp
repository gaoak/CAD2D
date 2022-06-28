#include "util.h"
#include <cmath>

double distance(std::vector<double> p0, std::vector<double> p1) {
    double res = 0.;
    for(size_t i=0; i<p0.size() && i<p1.size(); ++i) {
        res += (p0[i] - p1[i]) * (p0[i] - p1[i]);
    }
    return sqrt(res);
}

double distance(std::vector<std::vector<double> > p) {
    double res = 0.;
    for(size_t i=1; i<p.size(); ++i) {
        res += distance(p[i-1], p[i]);
    }
    return res;
}