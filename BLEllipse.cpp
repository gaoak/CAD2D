#include "BLEllipse.h"
#include "RectRegion.h"
#include "util.h"
using namespace std;

BLEllipse::BLEllipse(std::map<std::string, double> &doubleparams,
                     std::map<std::string, int> &intparams)
    : BLMeshModule(doubleparams, intparams) {}

void BLEllipse::Initialise() {
  // cylinder
  g_thetaA[0][0] = p["Theta0"];
  g_thetaA[0][1] = 0.;
  g_thetaA[1][0] = p["Theta1"];
  g_thetaA[1][1] = 0.;
  g_thetaA[2][0] = p["Theta2"];
  g_thetaA[2][1] = 0.;
  g_thetaA[3][0] = p["Theta3"];
  g_thetaA[3][1] = 0.;
  g_thetaA[4][0] = p["Theta4"];
  g_thetaA[4][1] = 0.;

  Cedge0 = LineEdge(g_thetaA[0], g_thetaA[1], q["nLE"], UNIFORM, 0., 0.);
  Cedge2 = LineEdge(g_thetaA[2], g_thetaA[3], q["nTE"], UNIFORM, 0., 0.);
  Cedge1 = LineEdge(g_thetaA[1], g_thetaA[2], q["nUp"], BOUNDARYLAYER2,
                    fabs(g_thetaA[1][0] - g_thetaA[0][0]) / Cedge0.m_N, 1.5, 4,
                    fabs(g_thetaA[3][0] - g_thetaA[2][0]) / Cedge2.m_N, 1.5, 4);
  Cedge3 = LineEdge(g_thetaA[3], g_thetaA[4], q["nLow"], BOUNDARYLAYER2,
                    fabs(g_thetaA[3][0] - g_thetaA[2][0]) / Cedge2.m_N, 1.5, 4,
                    fabs(g_thetaA[1][0] - g_thetaA[0][0]) / Cedge0.m_N, 1.5, 4);

  a = p["ChordLen"] * 0.5;
  b = p["Thickness"] * 0.5;
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
  Rects.push_back(RectRegion(edges0, "LE", false));
  setRadiusLayers(nBLayers);
  Rects[Rects.size() - 1].MeshGen(q["nLE"], nBLayers, eBoundaryLayer1);
  Rects[Rects.size() - 1].Tec360Pts("LE.dat");
  // edge 01
  edges0[0] = (void *)BLedge[1];
  Rects.push_back(RectRegion(edges0, "Up", false));
  Rects[Rects.size() - 1].MeshGen(q["nUp"], nBLayers, eBoundaryLayer1);
  Rects[Rects.size() - 1].Tec360Pts("Up.dat");
  // edge 01
  edges0[0] = (void *)BLedge[2];
  Rects.push_back(RectRegion(edges0, "TE", false));
  Rects[Rects.size() - 1].MeshGen(q["nTE"], nBLayers, eBoundaryLayer1);
  Rects[Rects.size() - 1].Tec360Pts("TE.dat");
  // edge 01
  edges0[0] = (void *)BLedge[3];
  Rects.push_back(RectRegion(edges0, "Low", false));
  Rects[Rects.size() - 1].MeshGen(q["nLow"], nBLayers, eBoundaryLayer1);
  Rects[Rects.size() - 1].Tec360Pts("Low.dat");
  ///////////// combine the near field mesh
  for (unsigned int i = 0; i < Rects.size(); ++i) {
    combinedReg.AddRegion(Rects[i]);
  }
  return 0;
}

int BLEllipse::DefineBCs(MeshRegions &combinedReg, int offset,
                         std::vector<void *> &BLedge) {
  int curvedpts = q["curvedpts"];
  combinedReg.defineBoundary(BLedge[0], Cedge0.m_N, 0 + offset, curvedpts);
  combinedReg.defineBoundary(BLedge[1], Cedge1.m_N, 0 + offset, curvedpts);
  combinedReg.defineBoundary(BLedge[2], Cedge2.m_N, 0 + offset, curvedpts);
  combinedReg.defineBoundary(BLedge[3], Cedge3.m_N, 0 + offset, curvedpts);
  return 1 + offset;
}

std::vector<double> BLEllipse::edge0(double s) {
  double t = Cedge0.Evaluate(s)[0];
  std::vector<double> res(2, 0.);
  res[0] = a * cos(t);
  res[1] = b * sin(t);
  return Transform(res);
}

std::vector<double> BLEllipse::edge1(double s) {
  double t = Cedge1.Evaluate(s)[0];
  std::vector<double> res(2, 0.);
  res[0] = a * cos(t);
  res[1] = b * sin(t);
  return Transform(res);
}

std::vector<double> BLEllipse::edge2(double s) {
  double t = Cedge2.Evaluate(s)[0];
  std::vector<double> res(2, 0.);
  res[0] = a * cos(t);
  res[1] = b * sin(t);
  return Transform(res);
}

std::vector<double> BLEllipse::edge3(double s) {
  double t = Cedge3.Evaluate(s)[0];
  std::vector<double> res(2, 0.);
  res[0] = a * cos(t);
  res[1] = b * sin(t);
  return Transform(res);
}