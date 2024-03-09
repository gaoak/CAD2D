#include "BLRectangle.h"
#include "RectRegion.h"
#include "util.h"
using namespace std;

BLRectangle::BLRectangle(std::map<std::string, double> &doubleparams,
                         std::map<std::string, int> &intparams)
    : BLFlatPlate(doubleparams, intparams) {}

void BLRectangle::Initialise() {
  BLFlatPlate::Initialise();
  g_thetaA[0][0] = 0.;
  g_thetaA[10][0] = 0.;
  g_thetaA[11][0] = 0.;
  g_thetaA[0][1] = 0.5 * p["Thickness"];
  g_thetaA[10][1] = -0.5 * p["Thickness"];
  g_thetaA[11][1] = 0.;

  g_thetaA[4][0] = p["ChordLen"];
  g_thetaA[5][0] = p["ChordLen"];
  g_thetaA[6][0] = p["ChordLen"];
  g_thetaA[4][1] = 0.5 * p["Thickness"];
  g_thetaA[5][1] = 0.;
  g_thetaA[6][1] = -0.5 * p["Thickness"];

  g_ptsA[0][0] = 0;
  g_ptsA[4][0] = p["ChordLen"];
  g_ptsA[6][0] = p["ChordLen"];
  g_ptsA[10][0] = 0.;

  Cedge1011 = LineEdge(g_thetaA[10], g_thetaA[11], q["nLow0"], BOUNDARYLAYER0,
                       p["hFirstLayer"], 1.2, 4, 0., 0., 0);
  Cedge110 = LineEdge(g_thetaA[11], g_thetaA[0], q["nUp0"], BOUNDARYLAYER1, 0.,
                      0., 0, p["hFirstLayer"], 1.2, 4);
  Cedge45 = LineEdge(g_thetaA[4], g_thetaA[5], q["nUp5"], BOUNDARYLAYER0,
                     p["hFirstLayer"], 1.2, 4, 0., 0., 0);
  Cedge56 = LineEdge(g_thetaA[5], g_thetaA[6], q["nLow5"], BOUNDARYLAYER1, 0.,
                     0., 0, p["hFirstLayer"], 1.2, 4);
}

int BLRectangle::MeshGenLEdge(MeshRegions &combinedReg,
                              std::vector<void *> &BLedge) {
  double hFirstLayer = p["hFirstLayer"];
  double progress = p["progress"];
  double maxLayerh = p["maxLayerh"];
  int nLayersU0 = findNlayers(hFirstLayer, progress, p["upperBL0"], maxLayerh);
  int nLayersU1 = findNlayers(hFirstLayer, progress, p["upperBL1"], maxLayerh);
  int nLayersL0 = findNlayers(hFirstLayer, progress, p["lowerBL0"], maxLayerh);
  int nLayersL1 = findNlayers(hFirstLayer, progress, p["lowerBL1"], maxLayerh);
  /////////near body region////////////////
  std::vector<RectRegion> Rects;
  // boundary layer region 0
  std::vector<void *> edges0;
  void *edgetmp;
  // edge 0-1
  edges0.push_back((void *)BLedge[4]);
  edges0.push_back((void *)radiusEdge);
  edges0.push_back(edgetmp);
  edges0.push_back(edgetmp);
  Rects.push_back(RectRegion(edges0, "up1", false));
  setRadiusLayers(nLayersU1);
  Rects[Rects.size() - 1].MeshGen(Cedge01.m_N, nLayersU1, eBoundaryLayer1,
                                  false, hFirstLayer, 0.);
  Rects[Rects.size() - 1].Tec360Pts("up1.dat");
  // edge 11-0
  edges0[0] = (void *)BLedge[1];
  Rects.push_back(RectRegion(edges0, "up0", false));
  setRadiusLayers(nLayersU0);
  Rects[Rects.size() - 1].MeshGen(Cedge110.m_N, nLayersU0, eBoundaryLayer1,
                                  false, 0., hFirstLayer);
  Rects[Rects.size() - 1].Tec360Pts("up0.dat");
  // edge 10-11
  edges0[0] = (void *)BLedge[0];
  Rects.push_back(RectRegion(edges0, "low0", false));
  setRadiusLayers(nLayersL0);
  Rects[Rects.size() - 1].MeshGen(Cedge1011.m_N, nLayersL0, eBoundaryLayer1,
                                  false, hFirstLayer, 0.);
  Rects[Rects.size() - 1].Tec360Pts("low0.dat");
  // edge 9-10
  edges0[0] = (void *)BLedge[11];
  Rects.push_back(RectRegion(edges0, "low1", false));
  setRadiusLayers(nLayersL1);
  Rects[Rects.size() - 1].MeshGen(Cedge910.m_N, nLayersL1, eBoundaryLayer1,
                                  false, 0., hFirstLayer);
  Rects[Rects.size() - 1].Tec360Pts("low1.dat");

  ///////////// combine the near field mesh
  for (unsigned int i = 0; i < Rects.size(); ++i) {
    combinedReg.AddRegion(Rects[i]);
  }
  return 0;
}
int BLRectangle::MeshGenTEdge(MeshRegions &combinedReg,
                              std::vector<void *> &BLedge) {
  double hFirstLayer = p["hFirstLayer"];
  double progress = p["progress"];
  double maxLayerh = p["maxLayerh"];
  int nLayersU4 = findNlayers(hFirstLayer, progress, p["upperBL4"], maxLayerh);
  int nLayersU5 = findNlayers(hFirstLayer, progress, p["upperBL5"], maxLayerh);
  int nLayersL4 = findNlayers(hFirstLayer, progress, p["lowerBL4"], maxLayerh);
  int nLayersL5 = findNlayers(hFirstLayer, progress, p["lowerBL5"], maxLayerh);
  /////////near body region////////////////
  std::vector<RectRegion> Rects;
  // boundary layer region 0
  std::vector<void *> edges0;
  void *edgetmp;
  // edge 6-7
  edges0.push_back((void *)BLedge[8]);
  edges0.push_back((void *)radiusEdge);
  edges0.push_back(edgetmp);
  edges0.push_back(edgetmp);
  Rects.push_back(RectRegion(edges0, "low4", false));
  setRadiusLayers(nLayersL4);
  Rects[Rects.size() - 1].MeshGen(Cedge67.m_N, nLayersL4, eBoundaryLayer1,
                                  false, hFirstLayer, 0.);
  Rects[Rects.size() - 1].Tec360Pts("low4.dat");
  // edge 5-6
  edges0[0] = (void *)BLedge[3];
  Rects.push_back(RectRegion(edges0, "low5", false));
  setRadiusLayers(nLayersL5);
  Rects[Rects.size() - 1].MeshGen(Cedge56.m_N, nLayersL5, eBoundaryLayer1,
                                  false, 0., hFirstLayer);
  Rects[Rects.size() - 1].Tec360Pts("low5.dat");
  // edge 4-5
  edges0[0] = (void *)BLedge[2];
  Rects.push_back(RectRegion(edges0, "up5", false));
  setRadiusLayers(nLayersU5);
  Rects[Rects.size() - 1].MeshGen(Cedge45.m_N, nLayersU5, eBoundaryLayer1,
                                  false, hFirstLayer, 0.);
  Rects[Rects.size() - 1].Tec360Pts("up5.dat");
  // edge 3-4
  edges0[0] = (void *)BLedge[7];
  Rects.push_back(RectRegion(edges0, "up4", false));
  setRadiusLayers(nLayersU4);
  Rects[Rects.size() - 1].MeshGen(Cedge34.m_N, nLayersU4, eBoundaryLayer1,
                                  false, 0., hFirstLayer);
  Rects[Rects.size() - 1].Tec360Pts("up4.dat");

  ///////////// combine the near field mesh
  for (unsigned int i = 0; i < Rects.size(); ++i) {
    combinedReg.AddRegion(Rects[i]);
  }
  return 0;
}

int BLRectangle::MeshGen(MeshRegions &combinedReg,
                         std::vector<void *> &BLedge) {
  MeshGenLEdge(combinedReg, BLedge);
  MeshGenUpper(combinedReg, BLedge);
  MeshGenTEdge(combinedReg, BLedge);
  MeshGenLower(combinedReg, BLedge);
  return 0;
}

// leading edge half circule
std::vector<double> BLRectangle::edge0(double s) {
  return Transform(Cedge1011.Evaluate(s));
}
std::vector<double> BLRectangle::edge1(double s) {
  return Transform(Cedge110.Evaluate(s));
}
// trailing edge half circule
std::vector<double> BLRectangle::edge2(double s) {
  return Transform(Cedge45.Evaluate(s));
}
std::vector<double> BLRectangle::edge3(double s) {
  return Transform(Cedge56.Evaluate(s));
}