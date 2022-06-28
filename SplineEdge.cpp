#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "SplineEdge.h"
#include "util.h"

SplineEdge::SplineEdge(std::string filename,
    std::map<int, double> params, int dim) {
    // clean variables
    m_pts.clear();
    m_spline.clear();
    m_pts.resize(2);
    // working space
    char buffer[1000];
    double p0, p1;
    // read points from file
    std::ifstream infile(filename.c_str());
    infile.getline(buffer, sizeof(buffer));
    while(!infile.eof()) {
        sscanf(buffer, "%lf%lf", &p0, &p1);
        m_pts[0].push_back(p0);
        m_pts[1].push_back(p1);
        infile.getline(buffer, sizeof(buffer));
    }
    infile.close();
    int np = m_pts[0].size();
    if (np<2) {
        std::cout << "Error: the number of points smaller than 2 in Spline " << filename << std::endl;
    }
    // create spline
    double step = 1./(np-1.);
    // boundary slope
    for (int d=0; d<dim; ++d) {
        if (params.find(2*d)==params.end()) {
            params[2*d] = (m_pts[d][1] - m_pts[d][0])/step;
        }
        if (params.find(2*d+1)==params.end()) {
            params[2*d+1] = (m_pts[d][np-1] - m_pts[d][np-2])/step;
        }
        m_spline.push_back(boost::math::cubic_b_spline<double>(
            m_pts[d].begin(), m_pts[d].end(),
            0, step,
            params[2*d], params[2*d+1]));
    }
    // calculate arclength on a refined grid
    calculateArcTable();
}

std::vector<double> SplineEdge::EvaluateTau(double tau) {
    std::vector<double> res(m_spline.size());
    for(size_t d=0; d<m_spline.size(); ++d) {
        res[d] = m_spline[d](tau);
    }
    return res;
}

void SplineEdge::appendArcLength(std::vector<double> p, double s, double tau) {
    for(size_t i=0; i<m_pts.size(); ++i) {
        m_arc[i].push_back(p[i]);
    }
    m_arc[m_pts.size()].push_back(s);
    m_arc[m_pts.size()+1].push_back(tau);
}

// arc[i]: x, y, z, s, tau
void SplineEdge::calculateArcTable() {
    m_arc.clear();
    m_arc.resize(2+m_pts.size());
    int Nelem = m_pts[0].size() - 1;
    int Nrefine = Nelem * 20;
    double tau = 0., dt = 1./Nrefine, arcl = 0.;
    std::vector<double> p0 = EvaluateTau(tau);
    appendArcLength(p0, 0., 0.);
    for(int i=1; i<=Nrefine; ++i) {
        tau = dt*i;
        std::vector<double> p1 = EvaluateTau(tau);
        arcl += distance(p0, p1);
        appendArcLength(p1, arcl, tau);
        p0 = p1;
    }
    m_ArcLength = arcl;
}

std::vector<double> SplineEdge::findx(double s) {
    int is = m_pts.size();
    int it = is + 1;
    auto lower = std::lower_bound(m_arc[is].begin(), m_arc[is].end(), s);
    int l = lower - m_arc[is].begin();
    if(l>m_arc[is].size()-2) {
        l = m_arc[is].size() - 2;
    }
    double tau = m_arc[it][l] + (m_arc[it][l+1]-m_arc[it][l])*(s - m_arc[is][l])/(m_arc[is][l+1]-m_arc[is][l]);
    //std::cout << "(" << l << "," << m_arc[is].size() << ":tau:" << tau << std::endl;
    std::vector<double> p = EvaluateTau(tau);
    return p;
}

double SplineEdge::finds(double x, int d) {
    size_t i = 0;
    if((m_arc[d][1]-m_arc[d][0])*(m_arc[d][0]-x)>=0.) {
        i = 0;
    } else {
        for(i=0; i<m_arc[d].size()-2; ++i) {
            if((m_arc[d][i+1]-x)*(x-m_arc[d][i]) >= 0.) {
                break;
            }
        }
    }
    int is = m_pts.size();
    return m_arc[is][i] + (m_arc[is][i+1]-m_arc[is][i])*(x - m_arc[d][i])/(m_arc[d][i+1]-m_arc[d][i]);
}

/*
double SplineEdge::finds(std::vector<double> x, int d) {
    if((m_arc[d][1]-m_arc[d][0])*(m_arc[d][0]-x[d])>=0.) {
        std::vector<double> p{m_arc[0][0], m_arc[1][0], m_arc[2][0]};
        return - distance(x, p);
    } else {
        size_t i = 1;
        for(; i<m_arc[d].size(); ++i) {
            if((m_arc[d][i]-x[d])*(x[d]-m_arc[d][i-1]) >= 0.) {
                break;
            }
        }
        int is = m_pts.size();
        std::vector<double> p{m_arc[0][i-1], m_arc[1][i-1], m_arc[2][i-1]};
        return m_arc[is][i-1] + distance(x, p);
    }
}*/

int main() {
    std::map<int, double> params;
    params[1] = 0.;
    SplineEdge spline("clarky_up.dat", params, 2);


    std::ofstream outfile("interped.dat");
    for(double x=0.; x<=0.1; x+=0.000020) {
        double s = spline.finds(x, 0);
        //std::cout << "x:" << x[0] << ", s:" << s << std::endl;
        std::vector<double> p = spline.findx(s);
        outfile << p[0] << " " << p[1] << std::endl;
    }
    outfile.close();
}