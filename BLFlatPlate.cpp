#include "BLFlatPlate.h"
#include "RectRegion.h"
#include "util.h"
using namespace std;

BLFlatPlate::BLFlatPlate(std::map<std::string, double> &doubleparams,
                         std::map<std::string, int> &intparams)
    : BLMeshModule(doubleparams, intparams) {
  // define other quantities
  setRadiusMesh(p["hFirstLayer"], p["progress"], p["maxLayerh"]);
  // airfoil
  g_thetaA[0][0] = 0.5 * M_PI;
  g_thetaA[4][0] = 0.5 * M_PI;
  g_thetaA[5][0] = 0.;
  g_thetaA[6][0] = -0.5 * M_PI;
  g_thetaA[10][0] = 1.5 * M_PI;
  g_thetaA[11][0] = M_PI;

  g_ptsA[0][0] = 0.5 * p["Thickness"];
  g_ptsA[1][0] = p["upperx1"];
  g_ptsA[2][0] = p["upperx2"];
  g_ptsA[3][0] = p["upperx3"];
  g_ptsA[4][0] = p["ChordLen"] - 0.5 * p["Thickness"];
  g_ptsA[0][1] = 0.5 * p["Thickness"];
  g_ptsA[1][1] = 0.5 * p["Thickness"];
  g_ptsA[2][1] = 0.5 * p["Thickness"];
  g_ptsA[3][1] = 0.5 * p["Thickness"];
  g_ptsA[4][1] = 0.5 * p["Thickness"];

  g_ptsA[6][0] = p["ChordLen"] - 0.5 * p["Thickness"];
  g_ptsA[7][0] = p["lowerx3"];
  g_ptsA[8][0] = p["lowerx2"];
  g_ptsA[9][0] = p["lowerx1"];
  g_ptsA[10][0] = 0.5 * p["Thickness"];
  g_ptsA[6][1] = -0.5 * p["Thickness"];
  g_ptsA[7][1] = -0.5 * p["Thickness"];
  g_ptsA[8][1] = -0.5 * p["Thickness"];
  g_ptsA[9][1] = -0.5 * p["Thickness"];
  g_ptsA[10][1] = -0.5 * p["Thickness"];
  Cedge1011 = LineEdge(g_thetaA[10], g_thetaA[11], q["nLow0"], UNIFORM, 0., 0.);
  Cedge110 = LineEdge(g_thetaA[11], g_thetaA[0], q["nUp0"], UNIFORM, 0., 0.);
  Cedge45 = LineEdge(g_thetaA[4], g_thetaA[5], q["nUp5"], UNIFORM, 0., 0.);
  Cedge56 = LineEdge(g_thetaA[5], g_thetaA[6], q["nLow5"], UNIFORM, 0., 0.);
  Cedge01 = LineEdge(g_ptsA[0], g_ptsA[1], q["nUp1"], BOUNDARYLAYER0,
                     p["hFirstLayer"], 1.6, 5, 0., 0., 0);
  Cedge12 = LineEdge(g_ptsA[1], g_ptsA[2], q["nUp2"], UNIFORM, 0., 0.);
  Cedge23 = LineEdge(g_ptsA[2], g_ptsA[3], q["nUp3"], UNIFORM, 0., 0.);
  Cedge34 = LineEdge(g_ptsA[3], g_ptsA[4], q["nUp4"], BOUNDARYLAYER1, 0., 0., 0,
                     p["hFirstLayer"], 1.6, 5);
  Cedge67 = LineEdge(g_ptsA[6], g_ptsA[7], q["nLow4"], BOUNDARYLAYER0,
                     p["hFirstLayer"], 1.6, 5, 0., 0., 0);
  Cedge78 = LineEdge(g_ptsA[7], g_ptsA[8], q["nLow3"], UNIFORM, 0., 0.);
  Cedge89 = LineEdge(g_ptsA[8], g_ptsA[9], q["nLow2"], UNIFORM, 0., 0.);
  Cedge910 = LineEdge(g_ptsA[9], g_ptsA[10], q["nLow1"], BOUNDARYLAYER1, 0., 0.,
                      0, p["hFirstLayer"], 1.6, 5);
}

int BLFlatPlate::MeshGen(MeshRegions &combinedReg,
                         std::vector<void *> &BLedge) {
  double hFirstLayer = p["hFirstLayer"];
  double progress = p["progress"];
  double maxLayerh = p["maxLayerh"];
  int nLayersU0 = findNlayers(hFirstLayer, progress, p["upperBL0"], maxLayerh);
  int nLayersU1 = findNlayers(hFirstLayer, progress, p["upperBL1"], maxLayerh);
  int nLayersU2 = findNlayers(hFirstLayer, progress, p["upperBL2"], maxLayerh);
  int nLayersU3 = findNlayers(hFirstLayer, progress, p["upperBL3"], maxLayerh);
  int nLayersU4 = findNlayers(hFirstLayer, progress, p["upperBL4"], maxLayerh);
  int nLayersU5 = findNlayers(hFirstLayer, progress, p["upperBL5"], maxLayerh);

  int nLayersL0 = findNlayers(hFirstLayer, progress, p["lowerBL0"], maxLayerh);
  int nLayersL1 = findNlayers(hFirstLayer, progress, p["lowerBL1"], maxLayerh);
  int nLayersL2 = findNlayers(hFirstLayer, progress, p["lowerBL2"], maxLayerh);
  int nLayersL3 = findNlayers(hFirstLayer, progress, p["lowerBL3"], maxLayerh);
  int nLayersL4 = findNlayers(hFirstLayer, progress, p["lowerBL4"], maxLayerh);
  int nLayersL5 = findNlayers(hFirstLayer, progress, p["lowerBL5"], maxLayerh);
  /////////near body region////////////////
  std::vector<RectRegion> Rects;
  // boundary layer region 0
  std::vector<void *> edges0;
  void *edgetmp;
  // edge 11-0
  edges0.push_back((void *)BLedge[1]);
  edges0.push_back((void *)radiusEdge);
  edges0.push_back(edgetmp);
  edges0.push_back(edgetmp);
  Rects.push_back(RectRegion(edges0, "up0", false));
  setRadiusLayers(nLayersU0);
  Rects[Rects.size() - 1].MeshGen(Cedge110.m_N, nLayersU0, eBoundaryLayer1);
  Rects[Rects.size() - 1].Tec360Pts("up0.dat");
  // edge 01
  edges0[0] = (void *)BLedge[4];
  Rects.push_back(RectRegion(edges0, "up1", false));
  setRadiusLayers(nLayersU1);
  Rects[Rects.size() - 1].MeshGen(Cedge01.m_N, nLayersU1, eBoundaryLayer1);
  Rects[Rects.size() - 1].Tec360Pts("up1.dat");
  // edge 12
  edges0[0] = (void *)BLedge[5];
  Rects.push_back(RectRegion(edges0, "up2", false));
  setRadiusLayers(nLayersU2);
  Rects[Rects.size() - 1].MeshGen(Cedge12.m_N, nLayersU2, eBoundaryLayer1);
  Rects[Rects.size() - 1].Tec360Pts("up2.dat");
  // edge 23
  edges0[0] = (void *)BLedge[6];
  Rects.push_back(RectRegion(edges0, "up3", false));
  setRadiusLayers(nLayersU3);
  Rects[Rects.size() - 1].MeshGen(Cedge23.m_N, nLayersU3, eBoundaryLayer1);
  Rects[Rects.size() - 1].Tec360Pts("up3.dat");
  // edge 34
  edges0[0] = (void *)BLedge[7];
  Rects.push_back(RectRegion(edges0, "up4", false));
  setRadiusLayers(nLayersU4);
  Rects[Rects.size() - 1].MeshGen(Cedge34.m_N, nLayersU4, eBoundaryLayer1);
  Rects[Rects.size() - 1].Tec360Pts("up4.dat");
  // edge 45
  edges0[0] = (void *)BLedge[2];
  Rects.push_back(RectRegion(edges0, "up5", false));
  setRadiusLayers(nLayersU5);
  Rects[Rects.size() - 1].MeshGen(Cedge45.m_N, nLayersU5, eBoundaryLayer1);
  Rects[Rects.size() - 1].Tec360Pts("up5.dat");

  // edge 5-6
  edges0[0] = (void *)BLedge[3];
  Rects.push_back(RectRegion(edges0, "low5", false));
  setRadiusLayers(nLayersL5);
  Rects[Rects.size() - 1].MeshGen(Cedge56.m_N, nLayersL5, eBoundaryLayer1);
  Rects[Rects.size() - 1].Tec360Pts("low5.dat");
  // edge 6-7
  edges0[0] = (void *)BLedge[8];
  Rects.push_back(RectRegion(edges0, "low4", false));
  setRadiusLayers(nLayersL4);
  Rects[Rects.size() - 1].MeshGen(Cedge67.m_N, nLayersL4, eBoundaryLayer1);
  Rects[Rects.size() - 1].Tec360Pts("low4.dat");
  // edge 7-8
  edges0[0] = (void *)BLedge[9];
  Rects.push_back(RectRegion(edges0, "low3", false));
  setRadiusLayers(nLayersL3);
  Rects[Rects.size() - 1].MeshGen(Cedge78.m_N, nLayersL3, eBoundaryLayer1);
  Rects[Rects.size() - 1].Tec360Pts("low3.dat");
  // edge 8-9
  edges0[0] = (void *)BLedge[10];
  Rects.push_back(RectRegion(edges0, "low2", false));
  setRadiusLayers(nLayersL2);
  Rects[Rects.size() - 1].MeshGen(Cedge89.m_N, nLayersL2, eBoundaryLayer1);
  Rects[Rects.size() - 1].Tec360Pts("low2.dat");
  // edge 9-10
  edges0[0] = (void *)BLedge[11];
  Rects.push_back(RectRegion(edges0, "low1", false));
  setRadiusLayers(nLayersL1);
  Rects[Rects.size() - 1].MeshGen(Cedge910.m_N, nLayersL1, eBoundaryLayer1);
  Rects[Rects.size() - 1].Tec360Pts("low1.dat");
  // edge 10-11
  edges0[0] = (void *)BLedge[0];
  Rects.push_back(RectRegion(edges0, "low0", false));
  setRadiusLayers(nLayersL0);
  Rects[Rects.size() - 1].MeshGen(Cedge1011.m_N, nLayersL0, eBoundaryLayer1);
  Rects[Rects.size() - 1].Tec360Pts("low0.dat");

  ///////////// combine the near field mesh
  for (unsigned int i = 0; i < Rects.size(); ++i) {
    combinedReg.AddRegion(Rects[i]);
  }
  return 0;
}

int BLFlatPlate::DefineBCs(MeshRegions &combinedReg, int offset,
                           std::vector<void *> &BLedge) {
  int curvedpts = q["curvedpts"];
  combinedReg.defineBoundary(BLedge[0], Cedge1011.m_N, 0 + offset, curvedpts);
  combinedReg.defineBoundary(BLedge[1], Cedge110.m_N, 0 + offset, curvedpts);
  combinedReg.defineBoundary(BLedge[2], Cedge45.m_N, 1 + offset, curvedpts);
  combinedReg.defineBoundary(BLedge[3], Cedge56.m_N, 1 + offset, curvedpts);
  combinedReg.defineBoundary(BLedge[4], Cedge01.m_N, 2 + offset);
  combinedReg.defineBoundary(BLedge[5], Cedge12.m_N, 2 + offset);
  combinedReg.defineBoundary(BLedge[6], Cedge23.m_N, 2 + offset);
  combinedReg.defineBoundary(BLedge[7], Cedge34.m_N, 2 + offset);
  combinedReg.defineBoundary(BLedge[8], Cedge67.m_N, 3 + offset);
  combinedReg.defineBoundary(BLedge[9], Cedge78.m_N, 3 + offset);
  combinedReg.defineBoundary(BLedge[10], Cedge89.m_N, 3 + offset);
  combinedReg.defineBoundary(BLedge[11], Cedge910.m_N, 3 + offset);
  return 4 + offset;
}

// leading edge half circule
std::vector<double> BLFlatPlate::edge0(double s) {
  double x0 = 0.5 * p["Thickness"], radius = 0.5 * p["Thickness"];
  double t = Cedge1011.Evaluate(s)[0];
  std::vector<double> res(2, 0.);
  res[0] = x0 + radius * cos(t);
  res[1] = radius * sin(t);
  return Transform(res);
}
std::vector<double> BLFlatPlate::edge1(double s) {
  double x0 = 0.5 * p["Thickness"], radius = 0.5 * p["Thickness"];
  double t = Cedge110.Evaluate(s)[0];
  std::vector<double> res(2, 0.);
  res[0] = x0 + radius * cos(t);
  res[1] = radius * sin(t);
  return Transform(res);
}
// trailing edge half circule
std::vector<double> BLFlatPlate::edge2(double s) {
  double x0 = p["ChordLen"] - 0.5 * p["Thickness"],
         radius = 0.5 * p["Thickness"];
  double t = Cedge45.Evaluate(s)[0];
  std::vector<double> res(2, 0.);
  res[0] = x0 + radius * cos(t);
  res[1] = radius * sin(t);
  return Transform(res);
}
std::vector<double> BLFlatPlate::edge3(double s) {
  double x0 = p["ChordLen"] - 0.5 * p["Thickness"],
         radius = 0.5 * p["Thickness"];
  double t = Cedge56.Evaluate(s)[0];
  std::vector<double> res(2, 0.);
  res[0] = x0 + radius * cos(t);
  res[1] = radius * sin(t);
  return Transform(res);
}
// upper surface
std::vector<double> BLFlatPlate::edge4(double s) {
  return Transform(Cedge01.Evaluate(s));
}
std::vector<double> BLFlatPlate::edge5(double s) {
  return Transform(Cedge12.Evaluate(s));
}
std::vector<double> BLFlatPlate::edge6(double s) {
  return Transform(Cedge23.Evaluate(s));
}
std::vector<double> BLFlatPlate::edge7(double s) {
  return Transform(Cedge34.Evaluate(s));
}
// lower surface
std::vector<double> BLFlatPlate::edge8(double s) {
  return Transform(Cedge67.Evaluate(s));
}
std::vector<double> BLFlatPlate::edge9(double s) {
  return Transform(Cedge78.Evaluate(s));
}
std::vector<double> BLFlatPlate::edge10(double s) {
  return Transform(Cedge89.Evaluate(s));
}
std::vector<double> BLFlatPlate::edge11(double s) {
  return Transform(Cedge910.Evaluate(s));
}