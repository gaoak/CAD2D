#include "BLAirfoil.h"
#include "airfoil.h"
#include "util.h"
using namespace std;

BLAirfoil::BLAirfoil(std::map<std::string, double> &doubleparams,
                     std::map<std::string, int> &intparams)
    : BLMeshModule(doubleparams, intparams) {
  // define other quantities
  double chordLen = p["chordLen"];
  double hFirstLayer = p["hFirstLayer"];
  double wakeLen = p["wakeLen"];
  double progress = p["progress"];
  double maxLayerh = p["maxLayerh"];
  int nLayers1 =
      findNlayers(hFirstLayer, progress, p["rBoundaryLayer1"], maxLayerh);
  int nLayers2 =
      findNlayers(hFirstLayer, progress, p["rBoundaryLayer2"], maxLayerh);
  int nLayers3 =
      findNlayers(hFirstLayer, progress, p["rBoundaryLayer3"], maxLayerh);
  int nLayers4 =
      findNlayers(hFirstLayer, progress, p["rBoundaryLayer4"], maxLayerh);
  int nLayers5 =
      findNlayers(hFirstLayer, progress, p["rBoundaryLayer5"], maxLayerh);
  int nLayers6 =
      findNlayers(hFirstLayer, progress, p["rBoundaryLayer6"], maxLayerh);
  int nLayers7 =
      findNlayers(hFirstLayer, progress, p["rBoundaryLayer7"], maxLayerh);
  AirFoilShPtr foil;
  if (q["NACAFOIL"]) {
    foil =
        std::make_shared<NACAmpxx>(p["chamber"], p["chamberp"], p["Thickness"]);
  } else if (q["WEDGEFOIL"]) {
    foil = std::make_shared<WedgeFoil>(p["D0"], p["D1"], 0.);
  }
  double sFrontUp = foil->finds(xmidUp2, 1);
  double sFrontLow = foil->finds(xmidLow2, -1);
  double sMidUp = foil->finds(xmidUp1, 1) - foil->finds(xmidUp2, 1);
  double sMidLow = foil->finds(xmidLow1, -1) - foil->finds(xmidLow2, -1);
  // airfoil
  setRadiusMesh(hFirstLayer, progress, maxLayerh);
  pts[0][0] = chordLen + wakeLen;
  pts[0][1] = p["wakeyUp"];

  pts[1][0] = chordLen + wakeLen;
  pts[1][1] = p["wakeDown"];

  pts[2][0] = foil->down(chordLen)[0];
  pts[2][1] = foil->down(chordLen)[1];

  pts[3][0] = foil->down(p["xmidLow1"])[0];
  pts[3][1] = foil->down(p["xmidLow1"])[1];

  pts[4][0] = foil->down(p["xmidLow2"])[0];
  pts[4][1] = foil->down(p["xmidLow2"])[1];

  pts[5][0] = foil->up(p["xmidUp2"])[0];
  pts[5][1] = foil->up(p["xmidUp2"])[1];

  pts[6][0] = foil->up(p["xmidUp1"])[0];
  pts[6][1] = foil->up(p["xmidUp1"])[1];

  pts[7][0] = foil->up(chordLen)[0];
  pts[7][1] = foil->up(chordLen)[1];
  setRadiusLayers(nLayers1);
  pts[13][0] = chordLen + wakeLen;
  pts[13][1] = -radiusEdge(1.)[0] - wakeLen * tan(p["nearWakeDiffuseAngle"]);

  setRadiusLayers(nLayers7);
  pts[14][0] = chordLen + wakeLen;
  pts[14][1] = radiusEdge(1.)[0] + wakeLen * tan(p["nearWakeDiffuseAngle"]);

  virtualpts[2][0] = pts[2][0];
  virtualpts[2][1] = pts[2][1];
  virtualpts[7][0] = pts[7][0];
  virtualpts[7][1] = pts[7][1];

  virtualpts[3][0] = p["xmidLow1"];
  virtualpts[3][1] = foil->up(p["xmidLow1"])[1];

  virtualpts[4][0] = p["xmidLow2"];
  virtualpts[4][1] = foil->up(p["xmidLow2"])[1];

  virtualpts[5][0] = p["xmidUp2"];
  virtualpts[5][1] = foil->up(p["xmidUp2"])[1];

  virtualpts[6][0] = p["xmidUp1"];
  virtualpts[6][1] = foil->up(p["xmidUp1"])[1];
  double hTrailingEdge = foil->roundTrailingSize();
  Cedge2 =
      LineEdge(virtualpts[2], virtualpts[3], q["nLow1"], BOUNDARYLAYER2,
               hTrailingEdge, (hTrailingEdge + hFirstLayer) / hTrailingEdge, 5,
               sMidLow / q["nLow2"], 2, 3);
  Cedge3 = LineEdge(virtualpts[3], virtualpts[4], q["nLow2"], BOUNDARYLAYER1,
                    0., 0., 0, (sFrontUp + sFrontLow) / q["nFront"],
                    p["growthrateLow2"], std::min(10, q["nLow2"] - 1));
  Cedge4 = LineEdge(virtualpts[4], virtualpts[5], q["nFront"], UNIFORM, 0., 0.);
  Cedge5 = LineEdge(virtualpts[5], virtualpts[6], q["nUp2"], BOUNDARYLAYER0,
                    (sFrontUp + sFrontLow) / q["nFront"], p["growthrateUp2"],
                    std::min(10, q["nUp2"] - 1), 0., 0., 1);
  Cedge6 = LineEdge(virtualpts[6], virtualpts[7], q["nUp1"], BOUNDARYLAYER2,
                    sMidUp / q["nUp2"], 2, 3, hTrailingEdge,
                    (hTrailingEdge + hFirstLayer) / hTrailingEdge, 5);
  Cedge0 = LineEdge(pts[0], pts[1], q["nTrailingEdge"], UNIFORM, 0., 0.);
  Cedge1 = LineEdge(pts[1], pts[2], q["nWake"], BOUNDARYLAYER1, 0., 0., 0,
                    hTrailingEdge + hFirstLayer, 2., 2);
  Cedge7 = LineEdge(pts[7], pts[0], q["nWake"], BOUNDARYLAYER0,
                    hTrailingEdge + hFirstLayer, 2., 2, 0., 0., 0);
  Cedge8 = LineEdge(pts[2], pts[7], q["nTrailingEdge"], UNIFORM, 0., 0.);
}

int BLAirfoil::MeshGen(MeshRegions &combinedReg, std::vector<void *> &BLedge) {
  double hFirstLayer = p["hFirstLayer"];
  double progress = p["progress"];
  double maxLayerh = p["maxLayerh"];
  int nLayers1 =
      findNlayers(hFirstLayer, progress, p["rBoundaryLayer1"], maxLayerh);
  int nLayers2 =
      findNlayers(hFirstLayer, progress, p["rBoundaryLayer2"], maxLayerh);
  int nLayers3 =
      findNlayers(hFirstLayer, progress, p["rBoundaryLayer3"], maxLayerh);
  int nLayers4 =
      findNlayers(hFirstLayer, progress, p["rBoundaryLayer4"], maxLayerh);
  int nLayers5 =
      findNlayers(hFirstLayer, progress, p["rBoundaryLayer5"], maxLayerh);
  int nLayers6 =
      findNlayers(hFirstLayer, progress, p["rBoundaryLayer6"], maxLayerh);
  int nLayers7 =
      findNlayers(hFirstLayer, progress, p["rBoundaryLayer7"], maxLayerh);
  /////////near body region////////////////
  std::vector<RectRegion> Rects;
  // boundary layer region 0
  std::vector<void *> edges0;
  void *edgetmp;
  // edge 2
  edges0.push_back((void *)edge2);
  edges0.push_back((void *)radiusEdge);
  edges0.push_back((void *)edgetmp);
  edges0.push_back((void *)edgetmp);
  Rects.push_back(RectRegion(edges0, "Rwall2", false));
  setRadiusLayers(nLayers2);
  Rects[Rects.size() - 1].MeshGen(Cedge2.m_N, nLayers2, eBoundaryLayer1, false,
                                  hFirstLayer, 0.);
  Rects[Rects.size() - 1].Tec360Pts("test0.dat");
  // edge 3
  edges0[0] = (void *)edge3;
  Rects.push_back(RectRegion(edges0, "Rwall3", false));
  setRadiusLayers(nLayers3);
  Rects[Rects.size() - 1].MeshGen(Cedge3.m_N, nLayers3, eBoundaryLayer1, false);
  Rects[Rects.size() - 1].Tec360Pts("test1.dat");
  // edge 4
  edges0[0] = (void *)edge4;
  Rects.push_back(RectRegion(edges0, "Rwall4", false));
  setRadiusLayers(nLayers4);
  Rects[Rects.size() - 1].MeshGen(Cedge4.m_N, nLayers4, eBoundaryLayer1, false);
  Rects[Rects.size() - 1].Tec360Pts("test3.dat");
  // edge 5
  edges0[0] = (void *)edge5;
  Rects.push_back(RectRegion(edges0, "Rwall5", false));
  setRadiusLayers(nLayers5);
  Rects[Rects.size() - 1].MeshGen(Cedge5.m_N, nLayers5, eBoundaryLayer1, false);
  Rects[Rects.size() - 1].Tec360Pts("test4.dat");
  // edge 6
  edges0[0] = (void *)edge6;
  Rects.push_back(RectRegion(edges0, "Rwall6", false));
  setRadiusLayers(nLayers6);
  Rects[Rects.size() - 1].MeshGen(Cedge6.m_N, nLayers6, eBoundaryLayer1, false,
                                  0., hFirstLayer);
  Rects[Rects.size() - 1].Tec360Pts("test5.dat");

  // give the normal direction, and points 12, 15
  vector<double> p1 = Rects[Rects.size() - 1].getVertexOffset(1);
  vector<double> n12 = Rects[Rects.size() - 1].getVertex(2);
  pts[7][0] = p1[0];
  pts[7][1] = p1[1];
  pts[15][0] = n12[0];
  pts[15][1] = n12[1];
  n12[0] = n12[0] - p1[0];
  n12[1] = n12[1] - p1[1];
  norm18[0] = n12[0] / sqrt(n12[0] * n12[0] + n12[1] * n12[1]);
  norm18[1] = n12[1] / sqrt(n12[0] * n12[0] + n12[1] * n12[1]);
  setRadiusLayers(nLayers7);
  pts[14][0] = pts[0][0] + norm18[0] * radiusEdge(1.)[0];
  pts[14][1] = 0. + norm18[1] * radiusEdge(1.)[0];
  p1 = Rects[0].getVertexOffset(0);
  n12 = Rects[0].getVertex(3);
  pts[2][0] = p1[0];
  pts[2][1] = p1[1];
  pts[12][0] = n12[0];
  pts[12][1] = n12[1];
  n12[0] = n12[0] - p1[0];
  n12[1] = n12[1] - p1[1];
  norm13[0] = n12[0] / sqrt(n12[0] * n12[0] + n12[1] * n12[1]);
  norm13[1] = n12[1] / sqrt(n12[0] * n12[0] + n12[1] * n12[1]);
  setRadiusLayers(nLayers1);
  pts[13][0] = pts[1][0] + norm13[0] * radiusEdge(1.)[0];
  pts[13][1] = 0. + norm13[1] * radiusEdge(1.)[0];
  // region 1
  std::vector<void *> edges1;
  edges1.push_back((void *)edge1);
  edges1.push_back((void *)edge0);
  edges1.push_back((void *)edge7);
  edges1.push_back((void *)edge8);
  Rects.push_back(RectRegion(edges1, "RTrailing"));
  Rects[Rects.size() - 1].MeshGen(Cedge1.m_N, Cedge0.m_N);
  Rects[Rects.size() - 1].Tec360Pts("testwake.dat");
  // region 2
  std::vector<void *> edges2;
  edges2.push_back((void *)edge14);
  edges2.push_back((void *)edge15);
  edges2.push_back((void *)edge1);
  edges2.push_back((void *)edge13);
  Rects.push_back(RectRegion(edges2, "RTrailDown"));
  Rects[Rects.size() - 1].MeshGen(Cedge14.m_N, Cedge15.m_N);
  Rects[Rects.size() - 1].Tec360Pts("test6.dat");
  // region 3
  std::vector<void *> edges3;
  edges3.push_back((void *)edge7);
  edges3.push_back((void *)edge16);
  edges3.push_back((void *)edge17);
  edges3.push_back((void *)edge18);
  Rects.push_back(RectRegion(edges3, "RTrailUp"));
  Rects[Rects.size() - 1].MeshGen(Cedge7.m_N, Cedge16.m_N);
  Rects[Rects.size() - 1].Tec360Pts("test7.dat");
  // regin trailing edge
  std::vector<void *> edgest;
  edgest.push_back((void *)edge8);
  edgest.push_back((void *)edgeb7_7);
  edgest.push_back((void *)edgeb2_b7);
  edgest.push_back((void *)edgeb2_2);
  Rects.push_back(RectRegion(edgest, "R_wake_Up"));
  Rects[Rects.size() - 1].MeshGen(Cedge8.m_N, Cedgeb7_7.m_N);
  Rects[Rects.size() - 1].Tec360Pts("test8.dat");

  ///////////// combine the near field mesh
  for (unsigned int i = 0; i < Rects.size(); ++i) {
    combinedReg.AddRegion(Rects[i]);
  }
  combinedReg.RemapPts((void *)roundTrailingEdge);
  return 0;
}

int BLAirfoil::DefineBCs(MeshRegions &combinedReg, int offset,
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
std::vector<double> BLAirfoil::edge0(double s) {
  double x0 = 0.5 * p["Thickness"], radius = 0.5 * p["Thickness"];
  double t = Cedge1011.Evaluate(s)[0];
  std::vector<double> res(2, 0.);
  res[0] = x0 + radius * cos(t);
  res[1] = radius * sin(t);
  return Transform(res);
}
std::vector<double> BLAirfoil::edge1(double s) {
  double x0 = 0.5 * p["Thickness"], radius = 0.5 * p["Thickness"];
  double t = Cedge110.Evaluate(s)[0];
  std::vector<double> res(2, 0.);
  res[0] = x0 + radius * cos(t);
  res[1] = radius * sin(t);
  return Transform(res);
}
// trailing edge half circule
std::vector<double> BLAirfoil::edge2(double s) {
  double x0 = p["ChordLen"] - 0.5 * p["Thickness"],
         radius = 0.5 * p["Thickness"];
  double t = Cedge45.Evaluate(s)[0];
  std::vector<double> res(2, 0.);
  res[0] = x0 + radius * cos(t);
  res[1] = radius * sin(t);
  return Transform(res);
}
std::vector<double> BLAirfoil::edge3(double s) {
  double x0 = p["ChordLen"] - 0.5 * p["Thickness"],
         radius = 0.5 * p["Thickness"];
  double t = Cedge56.Evaluate(s)[0];
  std::vector<double> res(2, 0.);
  res[0] = x0 + radius * cos(t);
  res[1] = radius * sin(t);
  return Transform(res);
}
// upper surface
std::vector<double> BLAirfoil::edge4(double s) {
  return Transform(Cedge01.Evaluate(s));
}
std::vector<double> BLAirfoil::edge5(double s) {
  return Transform(Cedge12.Evaluate(s));
}
std::vector<double> BLAirfoil::edge6(double s) {
  return Transform(Cedge23.Evaluate(s));
}
std::vector<double> BLAirfoil::edge7(double s) {
  return Transform(Cedge34.Evaluate(s));
}
// lower surface
std::vector<double> BLAirfoil::edge8(double s) {
  return Transform(Cedge67.Evaluate(s));
}
std::vector<double> BLAirfoil::edge9(double s) {
  return Transform(Cedge78.Evaluate(s));
}
std::vector<double> BLAirfoil::edge10(double s) {
  return Transform(Cedge89.Evaluate(s));
}
std::vector<double> BLAirfoil::edge11(double s) {
  return Transform(Cedge910.Evaluate(s));
}