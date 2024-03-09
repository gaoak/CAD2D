#ifndef BLRECTANGLE_H
#define BLRECTANGLE_H
#include "BLFlatPlate.h"
#include "LineEdge.h"
#include "MeshRegions.h"
#include "MeshTool.h"
#include <memory>
#include <vector>

class BLRectangle : public BLFlatPlate {
public:
  BLRectangle() {}
  BLRectangle(std::map<std::string, double> &doubleparams,
              std::map<std::string, int> &intparams);
  void Initialise() override;
  int MeshGen(MeshRegions &combinedReg, std::vector<void *> &BLedge) override;
  std::vector<double> edge0(double s) override;
  std::vector<double> edge1(double s) override;
  std::vector<double> edge2(double s) override;
  std::vector<double> edge3(double s) override;
  int MeshGenLEdge(MeshRegions &combinedReg,
                   std::vector<void *> &BLedge) override;
  int MeshGenTEdge(MeshRegions &combinedReg,
                   std::vector<void *> &BLedge) override;
};

#endif
