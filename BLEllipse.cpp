#include "BLEllipse.h"
#include "RectRegion.h"
#include "util.h"
using namespace std;

BLEllipse::BLEllipse(std::map<std::string, double> &doubleparams,
                     std::map<std::string, int> &intparams)
    : BLMeshModule(doubleparams, intparams) {}

void BLEllipse::Initialise() {
  // cylinder
  g_ptsC[0][0] = 0.;
  g_ptsC[0][1] = 0.;
}

int BLEllipse::MeshGen(MeshRegions &combinedReg, std::vector<void *> &BLedge) {
  double hFirstLayer = p["hFirstLayer"];
  double progress = p["progress"];
  double maxLayerh = p["maxLayerh"];
  int nBLayers =
      findNlayers(hFirstLayer, progress, p["wallBLThickness0"], maxLayerh);
  setRadiusMesh(hFirstLayer, progress, maxLayerh);
  /////////near body region////////////////
  std::vector<RectRegion> Rects;
  // boundary layer region 0
  std::vector<void *> edges0;
  void *edgetmp;
  edges0.push_back((void *)BLedge[0]);
  edges0.push_back((void *)radiusEdge);
  edges0.push_back(edgetmp);
  edges0.push_back(edgetmp);
  Rects.push_back(RectRegion(edges0, "cyl0", false));
  setRadiusLayers(nBLayers);
  Rects[Rects.size() - 1].MeshGen(q["Ncylinder"], nBLayers, eBoundaryLayer1);
  Rects[Rects.size() - 1].Tec360Pts("cyl0.dat");
  ///////////// combine the near field mesh
  for (unsigned int i = 0; i < Rects.size(); ++i) {
    combinedReg.AddRegion(Rects[i]);
  }
  return 0;
}

int BLEllipse::DefineBCs(MeshRegions &combinedReg, int offset,
                         std::vector<void *> &BLedge) {
  int curvedpts = q["curvedpts"];
  combinedReg.defineBoundary(BLedge[0], q["Ncylinder"], 0 + offset, curvedpts);
  return 1 + offset;
}

std::vector<double> BLEllipse::edge0(double s) {
  double t = -M_PI * s;
  double a = p["ChordLen"] * 0.5;
  double b = p["Thickness"] * 0.5;
  std::vector<double> res(2, 0.);
  res[0] = g_ptsC[0][0] + a * cos(t);
  res[1] = g_ptsC[0][1] + b * sin(t);
  return Transform(res);
}