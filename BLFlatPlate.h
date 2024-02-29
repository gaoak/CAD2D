#ifndef BLFLATPLATE_H
#define BLFLATPLATE_H
#include "BLMeshModule.h"
#include "LineEdge.h"
#include "MeshRegions.h"
#include "MeshTool.h"
#include <memory>
#include <vector>

class BLFlatPlate : public BLMeshModule {
public:
  BLFlatPlate() {}
  BLFlatPlate(std::map<std::string, double> &doubleparams,
              std::map<std::string, int> &intparams);
  void Initialise() override;
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
  std::vector<double> edge9(double s) override;
  std::vector<double> edge10(double s) override;
  std::vector<double> edge11(double s) override;
  double g_ptsA[20][2];
  double g_thetaA[20][2];

  LineEdge Cedge1011;
  LineEdge Cedge110;
  LineEdge Cedge45;
  LineEdge Cedge56;
  LineEdge Cedge01;
  LineEdge Cedge12;
  LineEdge Cedge23;
  LineEdge Cedge34;
  LineEdge Cedge67;
  LineEdge Cedge78;
  LineEdge Cedge89;
  LineEdge Cedge910;
};

#endif
