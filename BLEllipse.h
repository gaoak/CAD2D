#ifndef BLELLIPSE_H
#define BLELLIPSE_H
#include "BLMeshModule.h"
#include "LineEdge.h"
#include "MeshRegions.h"
#include "MeshTool.h"
#include <memory>
#include <vector>

class BLEllipse : public BLMeshModule {
public:
  BLEllipse() {}
  BLEllipse(std::map<std::string, double> &doubleparams,
            std::map<std::string, int> &intparams);
  void Initialise() override;
  int MeshGen(MeshRegions &combinedReg, std::vector<void *> &BLedge) override;
  int DefineBCs(MeshRegions &combinedReg, int offset,
                std::vector<void *> &BLedge) override;
  std::vector<double> edge0(double s) override;
  double g_ptsC[20][2];
};

#endif
