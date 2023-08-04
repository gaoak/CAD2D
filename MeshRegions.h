#ifndef MESHREGIONS_H
#define MESHREGIONS_H
#include "MeshRegion.h"
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <vector>
// regions manipulation
// related with physics
class MeshRegions : public MeshRegion {
public:
  MeshRegions(std::string name, double tolerance);
  int AddRegion(const MeshRegion &region);
  int AddRegion(const MeshRegion &region, std::set<int> &excludePts);
  int outXml(std::string filename);
  int defineBoundary(void *edgeFun, int N, int bndID, int Ncurve = 2,
                     double AoA = 0., int direction = 1,
                     void *mapFun = nullptr);
  int defineBoundary(std::vector<void *> edgeFuns,
                     double angle = 10. / 180. * M_PI);
  int outCOMPO(std::string filename, std::vector<int> comps);
  void outOuterRegion(std::string filename,
                      std::vector<std::vector<double>> box,
                      std::vector<double> center, double radius, bool exclude);
  void outInnerRegion(std::string filename,
                      std::vector<std::vector<double>> breakpts,
                      std::vector<double> center, double radius);
  int RemapPts(void *mapFun);
  int omeshBoundaryMapping(std::string filename, std::vector<double> center,
                           double radius);

private:
  std::map<int, std::vector<int>> m_boundary;
  void findAllBoundaryEdges();
  std::set<int> m_allBoundaryEdges;
  std::vector<int> m_curvedEdges;
  std::vector<std::vector<double>> m_curvedPoints;
  void outGeo(std::string filename, std::vector<std::vector<double>> box,
              std::vector<std::vector<int>> unSharedPts);
  int ExcludePts(std::set<int> &excludePts);
  std::map<int, int> m_boundarydefinitionStats;
};
#endif
