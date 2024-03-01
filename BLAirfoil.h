#ifndef BLAIRFOIL_H
#define BLAIRFOIL_H
#include "BLFlatPlate.h"
#include "LineEdge.h"
#include "MeshRegions.h"
#include "MeshTool.h"
#include "airfoil.h"
#include <memory>
#include <vector>

class BLAirfoil : public BLFlatPlate {
public:
  BLAirfoil() {}
  BLAirfoil(std::map<std::string, double> &doubleparams,
            std::map<std::string, int> &intparams);
  void Initialise() override;
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

private:
  AirFoilShPtr m_foil;
  double m_LER;
  double m_TER;
};

#endif
