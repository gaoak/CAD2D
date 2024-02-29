#include "BLAirfoil.h"
#include "RectRegion.h"
#include "util.h"
using namespace std;

BLAirfoil::BLAirfoil(std::map<std::string, double> &doubleparams,
                     std::map<std::string, int> &intparams)
    : BLFlatPlate(doubleparams, intparams) {}
void BLAirfoil::Initialise() { // define other quantities
  double upx0, lowx0, upx4, lowx4;
  if (q["NACAFOIL"]) {
    m_foil = std::make_shared<NACAmpxx>(p["chamber"], p["chamberp"],
                                        p["Thickness"], true);
  } else if (q["WEDGEFOIL"]) {
    m_foil =
        std::make_shared<WedgeFoil>(p["Thickness"], p["TEThickness"], 0., true);
  }
  std::map<std::string, double> foilparams;
  m_foil->GetInfo(foilparams);
  m_LER = foilparams["LERadius"];
  if (q["CutFore"]) {
    upx0 = m_LER;
    lowx0 = m_LER;
  } else {
    upx0 = foilparams["LEInterX"];
    lowx0 = foilparams["LEInterX"];
  }
  m_TER = foilparams["TERadius"];
  upx4 = foilparams["TEInterX"];
  lowx4 = foilparams["TEInterX"];
  //
  setRadiusMesh(p["hFirstLayer"], p["progress"], p["maxLayerh"]);
  // airfoil
  g_thetaA[0][0] = acos((upx0 - m_LER) / m_LER);
  g_thetaA[4][0] = acos((upx4 - 1. + m_TER) / m_TER);
  g_thetaA[5][0] = 0.;
  g_thetaA[6][0] = -acos((lowx4 - 1. + m_TER) / m_TER);
  g_thetaA[10][0] = 2. * M_PI - acos((lowx0 - m_LER) / m_LER);
  g_thetaA[11][0] = M_PI;

  g_ptsA[0][0] = upx0;
  g_ptsA[1][0] = p["upperx1"];
  g_ptsA[2][0] = p["upperx2"];
  g_ptsA[3][0] = p["upperx3"];
  g_ptsA[4][0] = upx4;
  g_ptsA[0][1] = m_foil->Upper(g_ptsA[0][0])[1];
  g_ptsA[1][1] = m_foil->Upper(g_ptsA[1][0])[1];
  g_ptsA[2][1] = m_foil->Upper(g_ptsA[2][0])[1];
  g_ptsA[3][1] = m_foil->Upper(g_ptsA[3][0])[1];
  g_ptsA[4][1] = m_foil->Upper(g_ptsA[4][0])[1];

  g_ptsA[6][0] = lowx4;
  g_ptsA[7][0] = p["lowerx3"];
  g_ptsA[8][0] = p["lowerx2"];
  g_ptsA[9][0] = p["lowerx1"];
  g_ptsA[10][0] = lowx0;
  g_ptsA[6][1] = m_foil->Lower(g_ptsA[6][0])[1];
  g_ptsA[7][1] = m_foil->Lower(g_ptsA[7][0])[1];
  g_ptsA[8][1] = m_foil->Lower(g_ptsA[8][0])[1];
  g_ptsA[9][1] = m_foil->Lower(g_ptsA[9][0])[1];
  g_ptsA[10][1] = m_foil->Lower(g_ptsA[10][0])[1];
  // straight edges
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

// leading edge half circule
std::vector<double> BLAirfoil::edge0(double s) {
  double t = Cedge1011.Evaluate(s)[0];
  double x = m_LER + m_LER * cos(t);
  return Transform(m_foil->Lower(x));
}
std::vector<double> BLAirfoil::edge1(double s) {
  double t = Cedge110.Evaluate(s)[0];
  double x = m_LER + m_LER * cos(t);
  return Transform(m_foil->Upper(x));
}
// trailing edge half circule
std::vector<double> BLAirfoil::edge2(double s) {
  double t = Cedge45.Evaluate(s)[0];
  double x = 1. - m_TER + m_TER * cos(t);
  return Transform(m_foil->Upper(x));
}
std::vector<double> BLAirfoil::edge3(double s) {
  double t = Cedge56.Evaluate(s)[0];
  double x = 1. - m_TER + m_TER * cos(t);
  return Transform(m_foil->Lower(x));
}
// upper surface
std::vector<double> BLAirfoil::edge4(double s) {
  return Transform(m_foil->Upper(Cedge01.Evaluate(s)[0]));
}
std::vector<double> BLAirfoil::edge5(double s) {
  return Transform(m_foil->Upper(Cedge12.Evaluate(s)[0]));
}
std::vector<double> BLAirfoil::edge6(double s) {
  return Transform(m_foil->Upper(Cedge23.Evaluate(s)[0]));
}
std::vector<double> BLAirfoil::edge7(double s) {
  return Transform(m_foil->Upper(Cedge34.Evaluate(s)[0]));
}
// lower surface
std::vector<double> BLAirfoil::edge8(double s) {
  return Transform(m_foil->Lower(Cedge67.Evaluate(s)[0]));
}
std::vector<double> BLAirfoil::edge9(double s) {
  return Transform(m_foil->Lower(Cedge78.Evaluate(s)[0]));
}
std::vector<double> BLAirfoil::edge10(double s) {
  return Transform(m_foil->Lower(Cedge89.Evaluate(s)[0]));
}
std::vector<double> BLAirfoil::edge11(double s) {
  return Transform(m_foil->Lower(Cedge910.Evaluate(s)[0]));
}