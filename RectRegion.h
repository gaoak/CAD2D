#ifndef RECTREGION_H
#define RECTREGION_H
#include "MeshRegion.h"
#include <set>
#include <vector>

enum EdgeType { eContinuous, eDiscrete };

enum MeshType {
  eIsoparametric,
  eBoundaryLayer0, // discrete norm
  eBoundaryLayer1, // continuous norm
  eStraightLine
};

class RectRegion : public MeshRegion {
public:
  RectRegion(std::vector<void *> edges, std::string name,
             bool connectivityCheck = true, double tolerance = 1.E-6,
             int edgeAttractMesh = 0,
             double attractEps = 1.); // for boundary layer mesh,
                                      // connectivityCheck needs to be false
  RectRegion(std::vector<std::vector<std::vector<double>>> edges,
             std::string name, bool connectivityCheck = true,
             double tolerance = 1.E-6, int edgeAttractMesh = 0,
             double attractEps = 1.);
  std::vector<double> EvaluateEdgePts(int i, double s);
  void PrintEdge(int eID, int N = 11);
  int MeshGen(int M, int N, MeshType method = eStraightLine,
              bool trimWakeSlope = false, double offset0 = 0.,
              double offset1 = 0.);
  void Tec360Pts(std::string fileName);
  int m_M; // number of elements
  int m_N; // number of elements
  std::vector<void *> m_edgesFun;
  std::vector<std::vector<std::vector<double>>> m_edgesPts;
  std::vector<double> m_edgesDirec;
  std::vector<std::vector<double>> m_vertex;
  std::vector<double> getVertex(int i);
  std::vector<double> getVertexOffset(int i);
  void SetEdgesDirec(std::vector<int> &dirs);

private:
  EdgeType m_edgeType;
  int CheckConnectivity();
  int EdgeAttractMesh(double x0, double y0, double &x, double &y);
  int m_edgeAttractMesh;
  double m_attractEps;
  int ptsByIsoParametric();
  int ptsByBoundaryLayer(bool trimWakeSlope, MeshType method,
                         double offset0 = 0., double offset1 = 0.);
  int ptsByStraightLine();
  std::vector<double> EvaluateEdgePtsDerivOneSide(int i, double s);
  std::vector<std::vector<double>> m_offsetPts;
};
#endif // RECTREGION_H
