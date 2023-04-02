#ifndef COMPOSITEDGE_H
#define COMPOSITEDGE_H
#include "LineEdge.h"
#include <cmath>
#include <vector>

class CompositEdge {
public:
  CompositEdge();
  std::vector<double> Evaluate(double s);
  std::vector<LineEdge> m_edges;
  std::vector<void *> m_functions;
  void addEdge(LineEdge edge, void *edgefunction);
  int m_N;

private:
  std::vector<int> m_Ns;
  int m_Nedges;
};

#endif
