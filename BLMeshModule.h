#ifndef BLMESHMODULE_H
#define BLMESHMODULE_H
#include "MeshRegions.h"
#include <memory>
#include <vector>
class BLMeshModule;

typedef std::shared_ptr<BLMeshModule> BLMeshModuleShPtr;

class BLMeshModule {
public:
  BLMeshModule(std::map<std::string, double> &doubleparams,
               std::map<std::string, int> &intparams) {
    p = doubleparams;
    q = intparams;
  }
  virtual int MeshGen(MeshRegions &combinedReg, std::vector<void *> &BLedge) {
    return 0;
  }
  virtual int DefineBCs(MeshRegions &combinedReg, int offset,
                        std::vector<void *> &BLedge) {
    return 0;
  }
  virtual std::vector<double> edge0(double s) {
    return std::vector<double>(3, 0.);
  }
  virtual std::vector<double> edge1(double s) {
    return std::vector<double>(3, 0.);
  }
  virtual std::vector<double> edge2(double s) {
    return std::vector<double>(3, 0.);
  }
  virtual std::vector<double> edge3(double s) {
    return std::vector<double>(3, 0.);
  }
  virtual std::vector<double> edge4(double s) {
    return std::vector<double>(3, 0.);
  }
  virtual std::vector<double> edge5(double s) {
    return std::vector<double>(3, 0.);
  }
  virtual std::vector<double> edge6(double s) {
    return std::vector<double>(3, 0.);
  }
  virtual std::vector<double> edge7(double s) {
    return std::vector<double>(3, 0.);
  }
  virtual std::vector<double> edge8(double s) {
    return std::vector<double>(3, 0.);
  }
  virtual std::vector<double> edge9(double s) {
    return std::vector<double>(3, 0.);
  }
  virtual std::vector<double> edge10(double s) {
    return std::vector<double>(3, 0.);
  }
  virtual std::vector<double> edge11(double s) {
    return std::vector<double>(3, 0.);
  }
  virtual std::vector<double> edge12(double s) {
    return std::vector<double>(3, 0.);
  }
  virtual std::vector<double> MapFunc(std::vector<double> p) { return p; }
  std::map<std::string, double> p;
  std::map<std::string, int> q;
  std::vector<double> Transform(std::vector<double> pt) {
    double x = pt[0], y = pt[1];
    pt[0] = x * cos(p["AoA"]) + y * sin(p["AoA"]) + p["Tx0"];
    pt[1] = -x * sin(p["AoA"]) + y * cos(p["AoA"]) + p["Ty0"];
    return pt;
  }
  void TransformMesh(MeshRegions &combinedReg) {
    combinedReg.transformation(p["AoA"], p["Tx0"], p["Ty0"]);
  }
};
#endif
