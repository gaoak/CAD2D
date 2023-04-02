#include "MeshRegions.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <vector>

int MeshRegions::outXml(std::string filename) {
  std::ofstream outxml(filename.c_str());
  // xml head
  outxml << "<?xml version=\"1.0\" encoding=\"utf-8\" ?>"
         << "\n";
  outxml << "<NEKTAR>"
         << "\n";
  outxml << "    <GEOMETRY DIM=\"2\" SPACE=\"2\">"
         << "\n";
  outxml << "        <VERTEX>"
         << "\n";
  //
  // vertex
  //
  for (unsigned int i = 0; i < m_pts.size(); i++) {
    outxml << "            <V ID=\"" << i << "\">";
    outxml << std::scientific << std::setprecision(17);
    outxml << std::setw(26) << m_pts[i][0] << std::setw(26) << m_pts[i][1]
           << std::setw(26) << 0.;
    outxml << "</V>"
           << "\n";
  }
  outxml << "        </VERTEX>"
         << "\n";
  outxml << "        <EDGE>"
         << "\n";
  //
  // edge
  //
  for (unsigned int i = 0; i < m_edges.size(); i++) {
    outxml << "            <E ID=\"" << i << "\">";
    outxml << std::setw(7) << m_edges[i][0] << std::setw(7) << m_edges[i][1]
           << "</E>"
           << "\n";
  }
  outxml << "        </EDGE>"
         << "\n";
  //
  // elements
  //
  sortCellFromQtoT();
  outxml << "        <ELEMENT>"
         << "\n";
  for (unsigned int i = 0; i < m_cells.size(); i++) {
    if (m_cells[i].size() == 4)
      outxml << "            <Q ID=\"" << i << "\">" << std::setw(7)
             << m_cells[i][0] << std::setw(7) << m_cells[i][1] << std::setw(7)
             << m_cells[i][2] << std::setw(7) << m_cells[i][3] << " </Q>"
             << "\n";
    if (m_cells[i].size() == 3)
      outxml << "            <T ID=\"" << i << "\">" << std::setw(7)
             << m_cells[i][0] << std::setw(7) << m_cells[i][1] << std::setw(7)
             << m_cells[i][2] << " </T>"
             << "\n";
  }
  outxml << "        </ELEMENT>"
         << "\n";
  //
  // curved edges
  //
  outxml << "        <CURVED>"
         << "\n";
  for (unsigned int i = 0; i < m_curvedEdges.size(); i++) {
    int curvedpoints = m_curvedPoints[i].size() / 3;
    outxml << "            <E ID=\"" << i << "\" EDGEID=\"" << m_curvedEdges[i]
           << "\" NUMPOINTS=\"" << curvedpoints
           << "\" TYPE=\"PolyEvenlySpaced\">";
    for (int k = 0; k < curvedpoints * 3; k++) {
      outxml << std::setw(26) << m_curvedPoints[i][k];
    }
    outxml << "</E>"
           << "\n";
  }
  outxml << "        </CURVED>"
         << "\n";
  std::cout << "mesh infomation in region " << m_name << std::endl;
  std::cout << "points number: " << m_pts.size() << std::endl;
  std::cout << "edges number : " << m_edges.size() << std::endl;
  std::cout << "cells number : " << m_cells.size() << std::endl;
  return 0;
}
int MeshRegions::outCOMPO(std::string filename, std::vector<int> compsi) {
  std::ofstream outxml(filename.c_str(), std::ofstream::app);
  int NUMMODES = 4, NUMPOINT = 5;
  // sort
  sortCellFromQtoT();
  std::vector<int> comps;
  std::set<int> compsSet;
  for (int i = 0; i < compsi.size(); ++i) {
    if (compsSet.find(compsi[i]) == compsSet.end()) {
      comps.push_back(compsi[i]);
      compsSet.insert(compsi[i]);
    }
  }
  if (comps[0] != 0)
    comps.push_back(0);
  if (comps[comps.size() - 1] != m_cells.size())
    comps.push_back(m_cells.size());
  int NTs = 0;
  for (NTs = 0; NTs < m_cells.size(); ++NTs) {
    if (m_cells[NTs].size() == 3)
      break;
  }
  bool toAdd = true;
  for (int i = 0; i < comps.size(); ++i) {
    if (comps[i] == NTs)
      toAdd = false;
  }
  if (toAdd) {
    comps.push_back(NTs);
    std::sort(comps.begin(), comps.end());
  }
  //
  // composite
  //
  outxml << "        <COMPOSITE>"
         << "\n";
  // check boundary definition
  findAllBoundaryEdges();
  std::set<int> allbndEdges = m_allBoundaryEdges;
  // output boundary
  int bndEdgeNumber = 0;
  int bndIDMax = 0;
  for (int i = 0; i < 99; ++i) {
    if (m_boundary.find(i) != m_boundary.end()) {
      bndEdgeNumber += m_boundary[i].size();
      for (unsigned int j = 0; j < m_boundary[i].size(); j++) {
        if (m_allBoundaryEdges.find(m_boundary[i][j]) ==
            m_allBoundaryEdges.end()) {
          std::cout << "Error: in region " << m_name << ", Edge "
                    << m_boundary[i][j] << " should not be a boundary "
                    << std::endl;
        } else {
          allbndEdges.erase(m_boundary[i][j]);
        }
      }
      outxml << "            <C ID=\"" << bndIDMax << "\"> E[";
      for (unsigned int j = 0; j < m_boundary[i].size() - 1; j++) {
        outxml << m_boundary[i][j] << ",";
      }
      outxml << m_boundary[i][m_boundary[i].size() - 1] << "] </C>"
             << "\n";
      ++bndIDMax;
    }
  }
  if (bndEdgeNumber != m_allBoundaryEdges.size()) {
    std::cout << "warning: in region " << m_name
              << ", boudnary condition incomplete C[" << bndIDMax << "]"
              << std::endl;
    outxml << "            <C ID=\"" << bndIDMax << "\"> E[";
    auto it = allbndEdges.begin();
    outxml << (*it);
    for (++it; it != allbndEdges.end(); ++it) {
      outxml << "," << (*it);
    }
    outxml << "] </C>"
           << "\n";
    ++bndIDMax;
  }
  std::vector<int> quadComp;
  std::vector<int> trigComp;
  for (int i = 0; i < comps.size() - 1; ++i) {
    if (comps[i] < NTs) {
      outxml << "            <C ID=\"" << bndIDMax << "\"> Q[" << comps[i]
             << "-" << comps[i + 1] - 1 << "] </C>"
             << "\n";
      quadComp.push_back(bndIDMax);
    } else {
      outxml << "            <C ID=\"" << bndIDMax << "\"> T[" << comps[i]
             << "-" << comps[i + 1] - 1 << "] </C>"
             << "\n";
      trigComp.push_back(bndIDMax);
    }
    ++bndIDMax;
  }
  outxml << "        </COMPOSITE>"
         << "\n";

  outxml << "        <DOMAIN> C["
         << bndIDMax - (quadComp.size() + trigComp.size());
  if (quadComp.size() + trigComp.size() > 1)
    outxml << "-" << bndIDMax - 1;
  outxml << "] </DOMAIN>"
         << "\n";
  outxml << "    </GEOMETRY>"
         << "\n";
  // expansion
#ifdef OUTPUTEXP
  outxml << "    <EXPANSIONS>"
         << "\n";
  // outxml << "        <E COMPOSITE=\"C[5]\" NUMMODES=\"11\" TYPE=\"MODIFIED\"
  // FIELDS=\"u,v,p\" />" << "\n"
  for (int i = 0; i < quadComp.size(); ++i) {
    outxml << "        <E COMPOSITE=\"C[" << quadComp[i] << "]\" NUMMODES=\""
           << NUMMODES << ", " << NUMMODES
           << "\" BASISTYPE=\"Modified_A,Modified_A\" "
              "POINTSTYPE=\"GaussLobattoLegendre,GaussLobattoLegendre\" "
              "NUMPOINTS=\""
           << NUMPOINT << "," << NUMPOINT << "\" FIELDS=\"u,v,p\" />"
           << "\n";
  }
  for (int i = 0; i < trigComp.size(); ++i) {
    outxml << "        <E COMPOSITE=\"C[" << trigComp[i] << "]\" NUMMODES=\""
           << NUMMODES << ", " << NUMMODES
           << "\" BASISTYPE=\"Modified_A,Modified_B\" "
              "POINTSTYPE=\"GaussLobattoLegendre,GaussRadauMAlpha1Beta0\" "
              "NUMPOINTS=\""
           << NUMPOINT << "," << NUMPOINT - 1 << "\" FIELDS=\"u,v,p\" />"
           << "\n";
  }
  outxml << "    </EXPANSIONS>"
         << "\n";
#endif
  outxml << "</NEKTAR>"
         << "\n";
  std::cout << "Output file " << filename << std::endl;
  return 0;
}

MeshRegions::MeshRegions(std::string name, double tolerance)
    : MeshRegion(name, tolerance) {}

int MeshRegions::AddRegion(const MeshRegion &region,
                           std::set<int> &excludePts) {
  MeshRegions rtmp("tmpreg", m_tolerance);
  rtmp.AddRegion(region);
  rtmp.ExcludePts(excludePts);
  AddRegion(rtmp);
  return 0;
}

int MeshRegions::ExcludePts(std::set<int> &excludePts) {
  // rewrite pts
  std::vector<std::vector<double>> pts;
  std::map<int, int> ptsmap;
  for (int i = 0; i < m_pts.size(); ++i) {
    if (excludePts.find(i) == excludePts.end()) {
      pts.push_back(m_pts[i]);
      ptsmap[i] = pts.size() - 1;
    }
  }
  m_pts = pts;

  std::vector<std::vector<int>> edges;
  std::set<int> excludeEdge;
  std::map<int, int> edgemap;
  for (int i = 0; i < m_edges.size(); ++i) {
    if (excludePts.find(m_edges[i][0]) != excludePts.end() ||
        excludePts.find(m_edges[i][1]) != excludePts.end()) {
      excludeEdge.insert(i);
    } else {
      std::vector<int> e;
      e.push_back(ptsmap[m_edges[i][0]]);
      e.push_back(ptsmap[m_edges[i][1]]);
      edges.push_back(e);
      edgemap[i] = edges.size() - 1;
    }
  }
  m_edges = edges;

  std::vector<std::vector<int>> cells;
  for (int i = 0; i < m_cells.size(); ++i) {
    bool toremove = false;
    std::vector<int> c;
    for (int j = 0; j < m_cells[i].size(); ++j) {
      c.push_back(edgemap[m_cells[i][j]]);
      if (excludeEdge.find(m_cells[i][j]) != excludeEdge.end())
        toremove = true;
    }
    if (!toremove) {
      cells.push_back(c);
    }
  }
  m_cells = cells;
  ResetBndPts();
  return 0;
}

int MeshRegions::AddRegion(const MeshRegion &region) {
  rebuildEdgesIndex();
  // merge points
  std::set<int> ptsExistSet;
  std::map<int, int> ptsMap;
  for (unsigned int i = 0; i < region.m_pts.size(); ++i) {
    int pId;
    if (region.m_bndPts.find(i) != region.m_bndPts.end() &&
        pointIsExist(region.m_pts[i], pId)) {
      ptsExistSet.insert(i);
      ptsMap[i] = pId;
    } else {
      m_pts.push_back(region.m_pts[i]);
      ptsMap[i] = m_pts.size() - 1;
      if (region.m_bndPts.find(i) != region.m_bndPts.end()) {
        m_bndPts.insert(m_pts.size() - 1);
      }
    }
  }
  // merge edges
  std::map<int, int> edgeMap;
  for (unsigned int i = 0; i < region.m_edges.size(); ++i) {
    std::vector<int> e;
    e.push_back(ptsMap[region.m_edges[i][0]]);
    e.push_back(ptsMap[region.m_edges[i][1]]);
    std::set<int> es;
    es.insert(e[0]);
    es.insert(e[1]);
    if (m_edgesIndex.find(es) == m_edgesIndex.end()) {
      m_edges.push_back(e);
      edgeMap[i] = m_edges.size() - 1;
      m_edgesIndex[es] = m_edges.size() - 1;
    } else {
      edgeMap[i] = m_edgesIndex[es];
    }
  }
  // merge cells
  for (unsigned int i = 0; i < region.m_cells.size(); ++i) {
    std::vector<int> c;
    for (int k = 0; k < region.m_cells[i].size(); ++k)
      c.push_back(edgeMap[region.m_cells[i][k]]);
    m_cells.push_back(c);
  }
  m_name += ":";
  m_name += region.m_name;
  return 0;
}

int MeshRegions::RemapPts(void *mapFun) {
  if (mapFun == nullptr)
    return 0;
  std::vector<double> (*mapFunction)(std::vector<double>) =
      (std::vector<double>(*)(std::vector<double>))mapFun;
  for (int i = 0; i < m_pts.size(); ++i) {
    m_pts[i] = mapFunction(m_pts[i]);
  }
  return m_pts.size();
}

int MeshRegions::defineBoundary(void *edgeFun, int N, int bndID, int Ncurve,
                                double AoA, int direction, void *mapFun) {
  std::vector<double> (*edgeFunction)(double) =
      (std::vector<double>(*)(double))edgeFun;
  std::vector<double> (*mapFunction)(std::vector<double>) =
      (std::vector<double>(*)(std::vector<double>))mapFun;
  double ds = 2. / N;
  double dir = 1.;
  if (direction < 0)
    dir = -1.;
  // for(auto it = m_edgesIndex.begin(); it!=m_edgesIndex.end(); ++it) {
  //    std::cout << "EIndex " << *(it->first.begin()) << ", " << (it->second)
  //    << std::endl;
  //}
  for (int i = 0; i < N; ++i) {
    std::vector<double> p0 = edgeFunction(dir * (-1. + i * ds));
    std::vector<double> p1 = edgeFunction(dir * (-1. + (i + 1) * ds));
    if (mapFun) {
      p0 = mapFunction(p0);
      p1 = mapFunction(p1);
    }
    transform(p0, AoA);
    transform(p1, AoA);
    int id0, id1;
    if (!pointIsExist(p0, id0) || !pointIsExist(p1, id1)) {
      std::cout << "error points 0(" << p0[0] << "," << p0[1]
                << ") not found on boundary " << bndID << std::endl;
      std::cout << "error points 1(" << p1[0] << "," << p1[1]
                << ") not found on boundary " << bndID << std::endl;
      m_boundarydefinitionStats[bndID] = -1;
    }
    std::set<int> es;
    es.insert(id0);
    es.insert(id1);
    int eId;
    if (m_edgesIndex.find(es) == m_edgesIndex.end()) {
      std::cout << "error Edge (" << id0 << "---" << id1
                << ") not found on boundary " << bndID << std::endl;
      m_boundarydefinitionStats[bndID] = -1;
    }
    eId = m_edgesIndex[es];
    double edgeDirec = 1.;
    std::vector<double> pe0 = m_pts[m_edges[eId][0]];
    if (0.5 * fabs(pe0[0] - p0[0]) + 0.5 * fabs(pe0[1] - p0[1]) > m_tolerance)
      edgeDirec = -1;
    // std::cout << "insert bnd edge " << id0 << ", " << id1 << ": " << eId <<
    // std::endl;
    m_boundary[bndID].push_back(eId);
    if (Ncurve > 2) {
      m_curvedEdges.push_back(eId);
      std::vector<double> points;
      for (int j = 0; j < Ncurve; ++j) {
        double dds = ds / (Ncurve - 1);
        std::vector<double> pt;
        if (edgeDirec > 0)
          pt = edgeFunction(dir * (-1. + i * ds + j * dds));
        else
          pt = edgeFunction(dir * (-1. + (i + 1) * ds - j * dds));
        if (mapFun) {
          pt = mapFunction(pt);
        }
        transform(pt, AoA);
        points.push_back(pt[0]);
        points.push_back(pt[1]);
        points.push_back(0.);
      }
      m_curvedPoints.push_back(points);
    }
  }
  if (m_boundarydefinitionStats.find(bndID) ==
      m_boundarydefinitionStats.end()) {
    m_boundarydefinitionStats[bndID] = 1;
  }
  return 0;
}

static double getAngle(std::vector<double> p0, std::vector<double> p1, std::vector<double> p2) {
  double x0 = p1[0] - p0[0];
  double y0 = p1[1] - p0[1];
  double x1 = p2[0] - p1[0];
  double y1 = p2[1] - p1[1];
  return asin((x0*y1-y0*x1)/sqrt(x0*x0+y0*y0)/sqrt(x1*x1+y1*y1));
}

std::vector<std::vector<int>> MeshRegions::splitBoundaryPts(std::vector<std::vector<int>> &allbndpts, double angle) {
  std::vector<std::vector<int>> res;
  for(auto &bnd : allbndpts) {
    std::vector<int> split;
    for(size_t i=0; i<bnd.size(); ++i) {
      int im1 = (i-1+bnd.size())%bnd.size();
      int ip1 = (i+1+bnd.size())%bnd.size();
      double a = fabs(getAngle(m_pts[im1], m_pts[i], m_pts[ip1]));
      if(a>angle) {
        split.push_back(i);
      }
    }
    std::cout << "split size " << split.size() << std::endl;
    if(split.size()>1) {
      for(int i=0; i<split.size()-1; ++i) {
        std::vector<int> bn;
        for(int j=split[i]; j<=split[i+1]; ++j) {
          bn.push_back(j);
        }
        res.push_back(bnd);
      }
      std::vector<int> bn;
      for(int j=split[split.size()-1]; j<= split[0] + bnd.size();++j) {
        bn.push_back(j%bnd.size());
      }
      res.push_back(bnd);
    } else {
      res.push_back(bnd);
    }
  }
  return res;
}

int MeshRegions::defineBoundary(std::vector<void*> edgeFuns, double angle) {
  std::vector<std::vector<int>> tmp = extractBoundaryPoints();
  std::vector<std::vector<int>> allbndpts = splitBoundaryPts(tmp, angle);
  std::set<int> defined;
  int bndID = -1;
  for(auto fun : edgeFuns) {
    bool (*condition)(std::vector<double>&) =
      (bool(*)(std::vector<double>&))fun;
    for (int i = 0; i < allbndpts.size(); ++i) {
      bool findedge = false;
      for (int j = 0; j < allbndpts[i].size(); ++j) {
        int j1 = (j+1)%allbndpts[i].size();
        if(condition(m_pts[allbndpts[i][j]]) && condition(m_pts[allbndpts[i][j1]])) {
          std::set<int> es;
          es.insert(allbndpts[i][j]);
          es.insert(allbndpts[i][j1]);
          if (m_edgesIndex.find(es) != m_edgesIndex.end() && defined.find(m_edgesIndex[es])==defined.end()) {
            if(!findedge) {
              findedge = true;
              ++bndID;
              m_boundary[bndID] = std::vector<int>();
            }
            m_boundary[bndID].push_back(m_edgesIndex[es]);
            defined.insert(m_edgesIndex[es]);
          }
        }
      }
    }
  }
  return 0;
}

/***
 * Find all boundary edges and store them in a set m_allBoundaryEdges
***/
void MeshRegions::findAllBoundaryEdges() {
  std::vector<std::vector<int>> allbnd = extractBoundaryPoints();
  m_allBoundaryEdges.clear();
  for (int i = 0; i < allbnd.size(); ++i) {
    for (int j = 0; j < allbnd[i].size(); ++j) {
      std::set<int> es;
      es.insert(allbnd[i][j]);
      if (j < allbnd[i].size() - 1)
        es.insert(allbnd[i][j + 1]);
      else
        es.insert(allbnd[i][0]);
      if (m_edgesIndex.find(es) != m_edgesIndex.end()) {
        int eid = m_edgesIndex[es];
        m_allBoundaryEdges.insert(eid);
      }
    }
  }
}

int MeshRegions::omeshBoundaryMapping(std::string filename,
                                      std::vector<double> center,
                                      double radius) {
  extractBoundaryPoints();
  int j;
  // find outer edges
  std::vector<int> unSharedPts;
  std::set<int> unSharedSet;
  std::map<int, int> mapping;
  for (int i = 0; i < m_unSharedPts.size(); ++i) {
    for (j = 0; j < m_unSharedPts[i].size(); ++j) {
      if (fabs(m_pts[m_unSharedPts[i][j]][0] - center[0]) +
              fabs(m_pts[m_unSharedPts[i][j]][1] - center[1]) <
          radius) {
        unSharedPts = m_unSharedPts[i];
        break;
      }
    }
  }
  for (int i = 0; i < unSharedPts.size(); ++i) {
    unSharedSet.insert(unSharedPts[i]);
  }
  for (int i = 0; i < m_edges.size(); ++i) {
    if (unSharedSet.find(m_edges[i][0]) != unSharedSet.end() &&
        unSharedSet.find(m_edges[i][1]) == unSharedSet.end()) {
      mapping[m_edges[i][0]] = m_edges[i][1];
    } else if (unSharedSet.find(m_edges[i][1]) != unSharedSet.end() &&
               unSharedSet.find(m_edges[i][0]) == unSharedSet.end()) {
      mapping[m_edges[i][1]] = m_edges[i][0];
    }
  }
  std::ofstream ofile(filename.c_str());
  char buffer[1000];
  ofile << mapping.size() << " 2" << std::endl;
  for (auto it = mapping.begin(); it != mapping.end(); ++it) {
    sprintf(buffer, "%25.18lf %25.18lf\n%25.18lf %25.18lf\n",
            m_pts[it->first][0], m_pts[it->first][1], m_pts[it->second][0],
            m_pts[it->second][1]);
    ofile << buffer;
  }
  return 0;
}

void MeshRegions::outOuterRegion(std::string filename,
                                 std::vector<std::vector<double>> box,
                                 std::vector<double> center, double radius,
                                 bool exclude) {
  extractBoundaryPoints();
  int j;
  // find outer edges
  std::vector<std::vector<int>> unSharedPts;
  if (exclude) {
    for (int i = 0; i < m_unSharedPts.size(); ++i) {
      for (j = 0; j < m_unSharedPts[i].size(); ++j) {
        if (fabs(m_pts[m_unSharedPts[i][j]][0] - center[0]) +
                fabs(m_pts[m_unSharedPts[i][j]][1] - center[1]) <
            radius) {
          break;
        }
      }
      if (j < m_unSharedPts[i].size())
        continue;
      unSharedPts.push_back(m_unSharedPts[i]);
    }
  } else {
    for (int i = 0; i < m_unSharedPts.size(); ++i) {
      for (j = 0; j < m_unSharedPts[i].size(); ++j) {
        if (fabs(m_pts[m_unSharedPts[i][j]][0] - center[0]) +
                fabs(m_pts[m_unSharedPts[i][j]][1] - center[1]) <
            radius) {
          unSharedPts.push_back(m_unSharedPts[i]);
          break;
        }
      }
    }
  }
  outGeo(filename, box, unSharedPts);
}

void MeshRegions::outInnerRegion(std::string filename,
                                 std::vector<std::vector<double>> breakpts,
                                 std::vector<double> center, double radius) {
  extractBoundaryPoints();
  // find inner edges
  std::vector<std::vector<int>> unSharedPtsFront;
  std::vector<std::vector<int>> unSharedPtsBack;
  std::vector<int> pts;
  for (int i = 0; i < m_unSharedPts.size(); ++i) {
    for (int j = 0; j < m_unSharedPts[i].size(); ++j) {
      if (fabs(m_pts[m_unSharedPts[i][j]][0] - center[0]) +
              fabs(m_pts[m_unSharedPts[i][j]][1] - center[1]) <
          radius) {
        pts = m_unSharedPts[i];
        break;
      }
    }
  }
  // break inner boundary
  int index0 = -1, index1 = -1;
  for (int i = 0; i < pts.size(); ++i) {
    if (fabs(m_pts[pts[i]][0] - breakpts[0][0]) +
            fabs(m_pts[pts[i]][1] - breakpts[0][1]) <
        2. * m_tolerance) {
      index0 = i;
    }
    if (fabs(m_pts[pts[i]][0] - breakpts[1][0]) +
            fabs(m_pts[pts[i]][1] - breakpts[1][1]) <
        2. * m_tolerance) {
      index1 = i;
    }
  }
  if (index0 == -1 || index1 == -1 || index0 == index1) {
    std::cout << "error failed to break inner boundary." << std::endl;
  }
  if (index0 > index1) {
    int tmp = index1;
    index1 = index0;
    index0 = tmp;
  }
  std::vector<int> pts0;
  std::vector<int> pts1;
  for (int i = 0; i <= index0; ++i) {
    pts0.push_back(pts[i]);
  }
  for (int i = index1; i < pts.size(); ++i) {
    pts0.push_back(pts[i]);
  }
  for (int i = index0; i <= index1; ++i) {
    pts1.push_back(pts[i]);
  }
  if (m_pts[pts[index0 + 1]][0] > m_pts[pts[index0]][0]) {
    unSharedPtsBack.push_back(pts1);
    unSharedPtsFront.push_back(pts0);
  } else {
    unSharedPtsBack.push_back(pts0);
    unSharedPtsFront.push_back(pts1);
  }
  //
  std::string frontname("front_");
  std::string backname("back_");
  frontname += filename;
  backname += filename;
  std::vector<std::vector<double>> box;
  outGeo(frontname, box, unSharedPtsFront);
  outGeo(backname, box, unSharedPtsBack);
}

void MeshRegions::outGeo(std::string filename,
                         std::vector<std::vector<double>> box,
                         std::vector<std::vector<int>> unSharedPts) {
  std::ofstream outgmsh(filename.c_str());
  int index = 1, j;
  // find outer edges
  // points
  std::map<int, int> pts;
  for (int i = 0; i < unSharedPts.size(); ++i) {
    for (j = 0; j < unSharedPts[i].size(); ++j) {
      pts[unSharedPts[i][j]] = index;
      outgmsh << std::scientific << std::setprecision(17);
      outgmsh << "Point(" << index << ") = {" << std::setw(26)
              << m_pts[unSharedPts[i][j]][0] << ", "
              << m_pts[unSharedPts[i][j]][1] << ", 0};\n";
      ++index;
    }
  }
  std::vector<std::vector<int>> line;
  if (box.size() >= 3) {
    int indexs = index;
    for (j = 0; j < box.size() - 1; ++j) {
      outgmsh << std::scientific << std::setprecision(17);
      outgmsh << "Point(" << index << ") = {" << std::setw(26) << box[j][0]
              << ", " << box[j][1] << ", 0};\n";
      std::vector<int> l;
      l.push_back(index);
      l.push_back(index + 1);
      line.push_back(l);
      ++index;
    }
    outgmsh << std::scientific << std::setprecision(17);
    outgmsh << "Point(" << index << ") = {" << std::setw(26) << box[j][0]
            << ", " << box[j][1] << ", 0};\n";
    std::vector<int> l;
    l.push_back(index);
    l.push_back(indexs);
    line.push_back(l);
    ++index;
  }
  // lines
  std::vector<std::vector<int>> lineloop;
  for (int i = 0; i < unSharedPts.size(); ++i) {
    std::vector<int> c;
    lineloop.push_back(c);
    for (j = 0; j < unSharedPts[i].size() - 1; ++j) {
      outgmsh << "Line(" << index << ") = {" << pts[unSharedPts[i][j]] << ", "
              << pts[unSharedPts[i][j + 1]] << "};\n";
      lineloop[i].push_back(index);
      ++index;
    }
    outgmsh << "Line(" << index << ") = {" << pts[unSharedPts[i][j]] << ", "
            << pts[unSharedPts[i][0]] << "};\n";
    lineloop[i].push_back(index);
    ++index;
  }
  std::vector<int> c;
  for (j = 0; j < line.size(); ++j) {
    outgmsh << "Line(" << index << ") = {" << line[j][0] << ", " << line[j][1]
            << "};\n";
    c.push_back(index);
    ++index;
  }
  if (c.size() > 0) {
    lineloop.push_back(c);
  }
  // lineloop
  std::vector<int> surf;
  for (int i = 0; i < lineloop.size(); ++i) {
    outgmsh << "Line Loop(" << index << ") = {";
    for (int j = 0; j < lineloop[i].size() - 1; ++j) {
      outgmsh << lineloop[i][j] << ", ";
    }
    outgmsh << lineloop[i][lineloop[i].size() - 1] << "};\n";
    surf.push_back(index);
    ++index;
  }
  // surf
  outgmsh << "Plane Surface(" << index << ") = {";
  for (int i = 0; i < surf.size() - 1; ++i) {
    outgmsh << surf[i] << ", ";
  }
  outgmsh << surf[surf.size() - 1] << "};\n";
  outgmsh << "Recombine Surface {" << index << "};\n";
  std::cout << "Output file " << filename << std::endl;
}