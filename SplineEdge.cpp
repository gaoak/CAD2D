#include "SplineEdge.h"
#include "util.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

/***
 * find the point from the tau
 * return a point on the spline
 ***/
std::vector<double> SplineEdge::EvaluateTau(double tau) {
  std::vector<double> res(m_spline.size());
  for (size_t d = 0; d < m_spline.size(); ++d) {
    res[d] = m_spline[d](tau);
  }
  return res;
}

/***
 * find the point from a monotonically varying x = m_arc[d]
 * first find tau by linear interpolation
 * then EvaluateTau
 * return a point on the spline
 ***/
std::vector<double> SplineEdge::Evaluatex(double x, int d) {
  std::vector<int> target(1, ArchIndex::eTau);
  return EvaluateTau(LookFor(x, d, target)[0]);
}

/***
 * find the point from the arclength
 * first find tau by linear interpolation
 * then EvaluateTau
 * return a point on the spline
 ***/
std::vector<double> SplineEdge::Evaluates(double s) {
  std::vector<int> target(1, ArchIndex::eTau);
  return EvaluateTau(LookFor(s, ArchIndex::eS, target)[0]);
}

/***
 * find the point from the scaled arclength [-1, 1]
 * first find tau by linear interpolation
 * then EvaluateTau
 * return a point on the spline
 ***/
std::vector<double> SplineEdge::Evaluate(double s) {
  s = 0.5 * (s + 1.0) * m_ArcLength;
  return Evaluates(s);
}

double SplineEdge::LookFors(double x, int d) {
  std::vector<int> target(1, ArchIndex::eS);
  return LookFor(x, d, target)[0];
}

/***
 * find approximately values from given x[monotone]
 ***/
std::vector<double> SplineEdge::LookFor(double x, int d,
                                        std::vector<int> &target) {
  size_t i = 0;
  if ((m_arc[d][1] - m_arc[d][0]) * (m_arc[d][0] - x) >= 0.) {
    i = 0;
  } else {
    for (i = 0; i < m_arc[d].size() - 2; ++i) {
      if ((m_arc[d][i + 1] - x) * (x - m_arc[d][i]) >= 0.) {
        break;
      }
    }
  }
  std::vector<double> res(target.size());
  for (size_t k = 0; k < target.size(); ++k) {
    int is = target[k];
    res[k] = m_arc[is][i] + (m_arc[is][i + 1] - m_arc[is][i]) *
                                (x - m_arc[d][i]) /
                                (m_arc[d][i + 1] - m_arc[d][i]);
  }
  return res;
}

SplineEdge::SplineEdge(std::string filename) {
  // clean variables
  m_pts.clear();
  m_spline.clear();
  m_pts.resize(2);
  // working space
  char buffer[1000];
  double p0, p1;
  std::ifstream infile(filename.c_str());
  if (!infile.is_open()) {
    std::cout << "Error: unable to open file " << filename << std::endl;
    exit(-1);
  }
  // read header
  // params format is (d, slope)
  // slope on boundary d%2 in dimension [d/2]
  // if slope is not given, one-side finite difference is used
  std::map<int, double> params;
  int dim;
  infile.getline(buffer, sizeof(buffer));
  sscanf(buffer, "%d", &dim);
  infile.getline(buffer, sizeof(buffer));
  while (!NanString(buffer, sizeof(buffer)) && !infile.eof()) {
    int i1;
    double v1;
    sscanf(buffer, "%d%lf", &i1, &v1);
    params[i1] = v1;
    infile.getline(buffer, sizeof(buffer));
  }
  // read points from file
  infile.getline(buffer, sizeof(buffer));
  while (!infile.eof()) {
    sscanf(buffer, "%lf%lf", &p0, &p1);
    m_pts[0].push_back(p0);
    m_pts[1].push_back(p1);
    infile.getline(buffer, sizeof(buffer));
  }
  infile.close();
  int np = m_pts[0].size();
  if (np < 2) {
    std::cout << "Error: the number of points smaller than 2 in Spline "
              << filename << std::endl;
  }
  // create spline
  double step = 1. / (np - 1.);
  // boundary slope
  for (int d = 0; d < dim; ++d) {
    if (params.find(2 * d) == params.end()) {
      params[2 * d] = (m_pts[d][1] - m_pts[d][0]) / step;
    }
    if (params.find(2 * d + 1) == params.end()) {
      params[2 * d + 1] = (m_pts[d][np - 1] - m_pts[d][np - 2]) / step;
    }
    m_spline.push_back(
        boost::math::interpolators::cardinal_cubic_b_spline<double>(
            m_pts[d].begin(), m_pts[d].end(), 0, step, params[2 * d],
            params[2 * d + 1]));
  }
  // calculate arclength on a refined grid
  calculateArcTable();
}

void SplineEdge::appendArcLength(std::vector<double> p, double s, double tau) {
  for (size_t i = 0; i < m_pts.size(); ++i) {
    m_arc[i].push_back(p[i]);
  }
  m_arc[m_pts.size()].push_back(s);
  m_arc[m_pts.size() + 1].push_back(tau);
}

// arc[i]: x, y, s, tau
void SplineEdge::calculateArcTable() {
  m_arc.clear();
  m_arc.resize(2 + m_pts.size());
  int Nelem = m_pts[0].size() - 1;
  int Nrefine = Nelem * 20;
  double tau = 0., dt = 1. / Nrefine, arcl = 0.;
  std::vector<double> p0 = EvaluateTau(tau);
  appendArcLength(p0, 0., 0.);
  for (int i = 1; i <= Nrefine; ++i) {
    tau = dt * i;
    std::vector<double> p1 = EvaluateTau(tau);
    arcl += distance(p0, p1);
    appendArcLength(p1, arcl, tau);
    p0 = p1;
  }
  m_ArcLength = arcl;
}

/*
double SplineEdge::LookFors(std::vector<double> x, int d) {
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

/*
int main() {
    std::map<int, double> params;
    params[0] = 0.;
    SplineEdge spline("clarky_low.dat", params, 2);


    std::ofstream outfile("interped.dat");
    for(double x=0.; x<=1.001; x+=0.0020) {
        double s = spline.LookFors(x, 0);
        //std::cout << "x:" << x[0] << ", s:" << s << std::endl;
        std::vector<double> p = spline.Evaluates(s);
        outfile << p[0] << " " << p[1] << std::endl;
    }
    outfile.close();
}*/