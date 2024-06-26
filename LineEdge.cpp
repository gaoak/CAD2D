#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
using namespace std;

#include "LineEdge.h"

static double alphaExpStretch(double h) {
  double alpha = 1.2;
  for (int i = 0; i < 50; ++i) {
    alpha = 0.5 * log(1. + 2. / h * alpha);
  }
  return alpha;
}

static double alphaSinStretch(double h) {
  double alpha = 1.4;
  for (int i = 0; i < 50; ++i) {
    alpha = 0.5 * acos(h * sin(2. * alpha) / (2. * alpha));
  }
  return alpha;
}

static double alphaSin2Stretch(double h0, double h1, double &alpha,
                               double &beta) {
  alpha = 1.4;
  beta = 0.;
  double tatb = (h0 - h1) / (h0 + h1);
  for (int i = 0; i < 50; ++i) {
    alpha = beta -
            acos(h0 * (sin(alpha + beta) - sin(-alpha + beta)) / (2. * alpha));
    beta = atan(tatb / tan(alpha));
  }
  return alpha;
}

double ExpStretch0(double s, double h0) {
  double alpha = alphaExpStretch(h0);
  return -1. + (exp(s * alpha) - exp(-alpha)) / (exp(alpha) - exp(-alpha)) * 2.;
}

double ExpStretch1(double s, double h1) {
  double alpha = alphaExpStretch(h1);
  return -1. +
         (-exp(-s * alpha) + exp(alpha)) / (-exp(-alpha) + exp(alpha)) * 2.;
}

double QuadStretch0(double s, double h0) {
  double alpha = 0.5 * (1. - h0);
  return alpha * (s * s - 1.) + s;
}

double QuadStretch1(double s, double h1) {
  double alpha = -0.5 * (1. - h1);
  return alpha * (s * s - 1.) + s;
}

double SinStretch0(double s, double h0) {
  double alpha = alphaSinStretch(h0);
  double beta = -alpha;
  return -1. + 2. * (sin(alpha * s + beta) - sin(-alpha + beta)) /
                   (sin(alpha + beta) - sin(-alpha + beta));
}

double SinStretch1(double s, double h1) {
  double alpha = alphaSinStretch(h1);
  double beta = alpha;
  return -1. + 2. * (sin(alpha * s + beta) - sin(-alpha + beta)) /
                   (sin(alpha + beta) - sin(-alpha + beta));
}

double SinStretch2(double s, double h0, double h1) {
  double alpha, beta;
  alphaSin2Stretch(h0, h1, alpha, beta);
  return -1. + 2. * (sin(alpha * s + beta) - sin(-alpha + beta)) /
                   (sin(alpha + beta) - sin(-alpha + beta));
}

double LineEdge::BuildDiscretes(double ds0, double ds1) {
  double ds;
  if (m_refineType == BOUNDARYLAYER0) {
    ds = ds0;
    m_discretes = std::vector<double>(m_N + 1, -1.);
    for (int i = 1; i <= m_nBLayers0; ++i) {
      if (ds > (1. - m_discretes[i - 1]) / (m_N - i + 1)) {
        m_nBLayers0 = i - 1;
        break;
      } else {
        m_discretes[i] = m_discretes[i - 1] + ds;
        ds *= m_q0;
      }
    }
    ds = (1. - m_discretes[m_nBLayers0]) / (m_N - m_nBLayers0);
    for (int i = m_nBLayers0 + 1; i <= m_N; ++i) {
      m_discretes[i] = m_discretes[i - 1] + ds;
    }
    m_discretes[m_N] = 1.;
  } else if (m_refineType == BOUNDARYLAYER1) {
    ds = ds1;
    m_discretes = std::vector<double>(m_N + 1, 1.);
    for (int i = 1; i <= m_nBLayers1; ++i) {
      if (ds > (1. + m_discretes[m_N - i + 1]) / (m_N - i + 1)) {
        m_nBLayers1 = i - 1;
        break;
      } else {
        m_discretes[m_N - i] = m_discretes[m_N - i + 1] - ds;
        ds *= m_q1;
      }
    }
    ds = (1. + m_discretes[m_N - m_nBLayers1]) / (m_N - m_nBLayers1);
    for (int i = m_nBLayers1 + 1; i <= m_N; ++i) {
      m_discretes[m_N - i] = m_discretes[m_N - i + 1] - ds;
    }
    m_discretes[0] = -1.;
  } else if (m_refineType == BOUNDARYLAYER2) {
    m_discretes = std::vector<double>(m_N + 1, -1.);
    m_discretes[m_N] = 1.;
    ds = ds0;
    for (int i = 1; i <= m_nBLayers0; ++i) {
      if (ds > (1. - m_discretes[i - 1]) / (m_N - i + 1)) {
        m_nBLayers0 = i - 1;
        break;
      } else {
        m_discretes[i] = m_discretes[i - 1] + ds;
        ds *= m_q0;
      }
    }
    ds = ds1;
    for (int i = 1; i <= m_nBLayers1; ++i) {
      if (ds > (m_discretes[m_N - i + 1] - m_discretes[m_nBLayers0]) /
                   (m_N - i + 1)) {
        m_nBLayers1 = i - 1;
        break;
      } else {
        m_discretes[m_N - i] = m_discretes[m_N - i + 1] - ds;
        ds *= m_q1;
      }
    }
    ds = (m_discretes[m_N - m_nBLayers1] - m_discretes[m_nBLayers0]) /
         (m_N - m_nBLayers0 - m_nBLayers1);
    if (ds <= 0) {
      std::cout << "error: unsupported BL refine type for edge number [r] "
                << m_N << std::endl;
    }
    for (int i = m_nBLayers0 + 1; i < m_N - m_nBLayers1; ++i) {
      m_discretes[i] = m_discretes[i - 1] + ds;
    }
  } else if (m_refineType == BOUNDARYLAYER4) {
    m_discretes = std::vector<double>(m_N + 1, -1.);
    m_discretes[m_N] = 1.;
    ds = ds1;
    for (int i = 1; i <= m_nBLayers1; ++i) {
      m_discretes[m_N - i] = m_discretes[m_N - i + 1] - ds;
      ds *= m_q1;
      if (ds > (m_discretes[m_N - i] + 1.) / (m_N - i)) {
        m_nBLayers1 = i;
        break;
      }
    }
    ds = ds0;
    for (int i = 1; i <= m_nBLayers0; ++i) {
      m_discretes[i] = m_discretes[i - 1] + ds;
      ds *= m_q0;
      if (ds > (m_discretes[m_N - m_nBLayers1] - m_discretes[i]) / (m_N - i)) {
        m_nBLayers0 = i;
        break;
      }
    }
    ds = (m_discretes[m_N - m_nBLayers1] - m_discretes[m_nBLayers0]) /
         (m_N - m_nBLayers0 - m_nBLayers1);
    if (ds <= 0) {
      std::cout << "error: unsupported BL refine type for edge number [r] "
                << m_N << std::endl;
    }
    for (int i = m_nBLayers0 + 1; i < m_N - m_nBLayers1; ++i) {
      m_discretes[i] = m_discretes[i - 1] + ds;
    }
  } else if (m_refineType == SMALLSTRETCH0) {
    int Nu = 0;
    m_discretes.clear();
    double stmp = -1.;
    m_discretes.push_back(stmp);
    double dstmp = ds0;
    while (true) {
      Nu = (1. - stmp) / ds1 + 0.5;
      double dux = (1. - stmp) / Nu;
      if (dstmp < dux) {
        stmp += dstmp;
        m_discretes.push_back(stmp);
        dstmp *= m_q0;
      } else {
        break;
      }
    }
    double dux = (1. - stmp) / Nu;
    double tmpratio = dux / ds1;
    if (tmpratio < 1.) {
      tmpratio = 1. / tmpratio;
    }
    std::cout << "Segment [" << m_p0[0] << ", " << m_p1[0] << "](" << m_q1
              << "), right end " << ds1 << "->" << dux << "[" << tmpratio << "]"
              << std::endl;
    for (int i = Nu - 1; i >= 0; --i) {
      stmp = 1. - dux * i;
      m_discretes.push_back(stmp);
    }
    m_N = m_discretes.size() - 1;
  } else if (m_refineType == SMALLSTRETCH1) {
    int Nu = 0;
    m_discretes.clear();
    double stmp = 1.;
    m_discretes.push_back(stmp);
    double dstmp = ds1;
    while (true) {
      Nu = (1. + stmp) / ds0 + 0.5;
      double dux = (1. + stmp) / Nu;
      if (dstmp < dux) {
        stmp -= dstmp;
        m_discretes.push_back(stmp);
        dstmp *= m_q1;
      } else {
        break;
      }
    }
    double dux = (1. + stmp) / Nu;
    double tmpratio = dux / ds0;
    if (tmpratio < 1.) {
      tmpratio = 1. / tmpratio;
    }
    std::cout << "Segment [" << m_p0[0] << ", " << m_p1[0] << "](" << m_q1
              << "), left end " << ds0 << "->" << dux << "[" << tmpratio << "]"
              << std::endl;
    for (int i = Nu - 1; i >= 0; --i) {
      stmp = -1. + dux * i;
      m_discretes.push_back(stmp);
    }
    std::reverse(m_discretes.begin(), m_discretes.end());
    m_N = m_discretes.size() - 1;
  }
  return 0.;
}

double LineEdge::DiscreteStretch(double s, double h0, double h1) {
  if (m_discretes.size() < m_N + 1)
    BuildDiscretes(h0, h1);
  double x = (s + 1.) / 2. * m_N;
  if (x < 0.)
    x = 0.;
  int ilower = floor(x);
  if (ilower == m_N)
    return 1.;
  return (x - ilower) * (m_discretes[ilower + 1] - m_discretes[ilower]) +
         m_discretes[ilower];
}

LineEdge::LineEdge(double *p0, double *p1, int N, int refineType, double h0,
                   double h1) {
  m_p0 = p0;
  m_p1 = p1;
  m_N = N;
  m_refineType = refineType;
  m_h0 = h0 * N;
  m_h1 = h1 * N;
}

// s = [-1., 1.]
std::vector<double> LineEdge::Evaluate(double s) {
  double Len = sqrt((m_p0[0] - m_p1[0]) * (m_p0[0] - m_p1[0]) +
                    (m_p0[1] - m_p1[1]) * (m_p0[1] - m_p1[1]));
  m_g0 = m_h0 / Len;
  m_g1 = m_h1 / Len;
  std::vector<double> res(2, 0.);
  if (m_refineType == EXPREFINE0) {
    s = ExpStretch0(s, m_g0);
  }
  if (m_refineType == EXPREFINE1) {
    s = ExpStretch1(s, m_g1);
  }
  if (m_refineType == SINREFINE0) {
    s = SinStretch0(s, m_g0);
  }
  if (m_refineType == SINREFINE1) {
    s = SinStretch1(s, m_g1);
  }
  if (m_refineType == SINREFINE2) {
    s = SinStretch2(s, m_g0, m_g1);
  }
  if (m_refineType == QUDREFINE0) {
    s = QuadStretch0(s, m_g0);
  }
  if (m_refineType == QUDREFINE1) {
    s = QuadStretch1(s, m_g1);
  }
  if (m_refineType == BOUNDARYLAYER0 || m_refineType == BOUNDARYLAYER1 ||
      m_refineType == BOUNDARYLAYER2 || m_refineType == BOUNDARYLAYER4 ||
      m_refineType == SMALLSTRETCH0 || m_refineType == SMALLSTRETCH1) {
    s = DiscreteStretch(s, m_g0 * 2., m_g1 * 2.);
  }
  for (int i = 0; i < 2; ++i) {
    res[i] = 0.5 * (1. - s) * m_p0[i] + 0.5 * (s + 1.) * m_p1[i];
  }
  return res;
}

LineEdge::LineEdge(double *p0, double *p1, int N, int refineType, double h0,
                   double q0, int NBlayers0, double h1, double q1,
                   int NBlayers1) {
  m_p0 = p0;
  m_p1 = p1;
  m_N = N;
  m_refineType = refineType;
  if ((BOUNDARYLAYER0 == refineType && (q0 < 1. || N <= NBlayers0)) ||
      (BOUNDARYLAYER1 == refineType && (q1 < 1. || N <= NBlayers1)) ||
      ((BOUNDARYLAYER2 == refineType || BOUNDARYLAYER4 == refineType) &&
       (q0 < 1. || q1 < 1. || N <= NBlayers0 + NBlayers1))) {
    std::cout << "error: unsupported BL refine type for edge number " << N
              << std::endl;
  }
  m_h0 = h0;
  m_h1 = h1;
  m_q0 = q0;
  m_q1 = q1;
  m_nBLayers0 = NBlayers0;
  m_nBLayers1 = NBlayers1;
}

LineEdge::LineEdge(double *p0, double *p1, int refineType, double h0, double h1,
                   double q) {
  m_p0 = p0;
  m_p1 = p1;
  m_refineType = refineType;
  if (refineType == SMALLSTRETCH0 && h0 > h1 ||
      refineType == SMALLSTRETCH1 && h0 < h1 || q < 1.) {
    std::cout
        << "error: incorrect definition for SMALLSTRETCH0 or SMALLSTRETCH1, "
        << h0 << ", " << h1 << std::endl;
    exit(-1);
  }
  m_h0 = h0;
  m_h1 = h1;
  m_q0 = q;
  m_q1 = q;
  Evaluate(0.);
}

// int main()
// {
//     double pts[2][2];
//     pts[0][0] = -5.;
//     pts[0][1] = 0.;
//     pts[1][0] = -1.;
//     pts[1][1] = 0.;
//     LineEdge Cedge2(pts[0], pts[1], SMALLSTRETCH1, 0.05, 0.01, 1.05);
//     double ds = 2./Cedge2.m_N;
//     for(int i=0; i<=Cedge2.m_N; ++i) {
//         printf("%lf %lf\n", Cedge2.Evaluate(ds*i-1.)[0],
//         Cedge2.Evaluate(ds*i-1.)[1]);
//     }
//     printf("\n");
// }

// int main()
// {
//     double pts[2][2];
//     pts[0][0] = -1.;
//     pts[0][1] = 0.;
//     pts[1][0] = 1.;
//     pts[1][1] = 0.;
//     int N=10;
//     double ds = 2./N;
//     LineEdge Cedge2(pts[0], pts[1], N , BOUNDARYLAYER2, 0.1, 2. , 4, 0.2, 1.,
//     2); LineEdge Cedge3(pts[1], pts[0], N , BOUNDARYLAYER4, 0.2, 1., 2,
//     0.1, 2. , 4); for(int i=0; i<=N; ++i) {
//         printf("%lf ", Cedge2.Evaluate(ds*i-1.)[0]);
//     }
//     printf("\n");
//     for(int i=N; i>=0; --i) {
//         printf("%lf ", Cedge3.Evaluate(ds*i-1.)[0]);
//     }
//     printf("\n");
// }