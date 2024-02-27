#ifndef BLAIRFOIL_H
#define BLAIRFOIL_H
#include "BLMeshModule.h"
#include "LineEdge.h"
#include "MeshRegions.h"
#include "MeshTool.h"
#include <memory>
#include <vector>

class BLAirfoil : public BLMeshModule {
public:
  BLAirfoil(std::map<std::string, double> &doubleparams,
            std::map<std::string, int> &intparams);
  int MeshGen(MeshRegions &combinedReg, std::vector<void *> &BLedge) override;
  int DefineBCs(MeshRegions &combinedReg, int offset,
                std::vector<void *> &BLedge) override;
  std::vector<double> edge0(double s) override;
  std::vector<double> edge1(double s) override;
  std::vector<double> edge2(double s) override;
  std::vector<double> edge3(double s) override;
  std::vector<double> edge4(double s) override;
  std::vector<double> edge5(double s) override;
  std::vector<double> edge6(double s) override;
  std::vector<double> edge7(double s) override;
  std::vector<double> edge8(double s) override;
  double pts[20][2];
  double virtualpts[20][2];

  LineEdge Cedge0;
  LineEdge Cedge1;
  LineEdge Cedge2;
  LineEdge Cedge3;
  LineEdge Cedge4;
  LineEdge Cedge5;
  LineEdge Cedge6;
  LineEdge Cedge7;
  LineEdge Cedge8;
};

#endif
