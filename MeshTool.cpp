#include "MeshTool.h"
#include "RectRegion.h"
#include "util.h"
#include<vector>
static struct {
    int nLayers;
    double hFirstLayer;
    double progress;
    double maxLayerh;
} BLlayer;

void setRadiusMesh(double h, double p, double m) {
    BLlayer.hFirstLayer = h;
    BLlayer.progress = p;
    BLlayer.maxLayerh = m;
}

void setRadiusLayers(int n) {
    BLlayer.nLayers = n;
}

std::vector<double> radiusEdge(double s) {
  int n = round(0.5 * (1. + s) * BLlayer.nLayers);
  static std::vector<std::vector<double>> reses;
  if (reses.size() < BLlayer.nLayers + 1) {
    reses.clear();
    std::vector<double> p0(2, 0.);
    reses.push_back(p0);
    double delta = BLlayer.hFirstLayer;
    for (int n = 1; n <= BLlayer.nLayers; ++n) {
      std::vector<double> p1(2, 0.);
      if (delta >= BLlayer.maxLayerh)
        delta = BLlayer.maxLayerh;
      p1[0] = reses[reses.size() - 1][0] + delta;
      delta *= BLlayer.progress;
      reses.push_back(p1);
    }
  }
  return reses[n];
}

int meshingBoundayLayer(MeshRegions &region, void* thickFunc, void* edge0,
                        std::string name, std::vector<std::vector<double>> trimnorm) {
  double (*ThicknessFunction)(std::vector<double>) =
        (double(*)(std::vector<double>))thickFunc;
    std::vector<void *> edges;
  // edge 2
  setRadiusLayers(1);
  edges.push_back((void *)edge0);
  edges.push_back((void *)radiusEdge);
  edges.push_back((void *)edge0);
  edges.push_back((void *)edge0);
  RectRegion firstlayer(edges, name, false);
  int Nslice = 100;
  firstlayer.MeshGen(Nslice, 1, eBoundaryLayer1);
  firstlayer.Tec360Pts(name+".dat");
  //get boundary points & normals
  std::vector<std::vector<double>> bndpts;
  std::vector<std::vector<double>> norms;
  for(int i=0; i<=Nslice; ++i) {
    bndpts.push_back(firstlayer.m_pts[i]);
    std::vector<double> p2 = firstlayer.m_pts[i+Nslice+1];
    double len = 0.;
    for(size_t j=0; j<p2.size(); ++j) {
      p2[j] -= firstlayer.m_pts[i][j];
      len += p2[j]*p2[j];
    }
    len = sqrt(len);
    for(size_t j=0; j<p2.size(); ++j) {
      p2[j] /= len;
    }
    norms.push_back(p2);
  }
  //trim normals
  if(trimnorm.size()) {
    norms[0] = trimnorm[0];
    norms[Nslice] = trimnorm[1];
  }
  //add boundary layer mesh
  for(int i=0; i<Nslice; ++i) {
    std::vector<double> k = intersect(bndpts[i], norms[i], bndpts[i+1], norms[i+1]);
    std::vector<double> pc(bndpts[i].size());
    for(size_t j=0; j<pc.size(); ++j) {
      pc[j] = 0.5*(bndpts[i][j]+bndpts[i+1][j]);
    }
    double thick = ThicknessFunction(pc);
    if(k[0]>0 &&k[1]>0) {
      thick = std::min(thick, k[0]);
      thick = std::min(thick, k[1]);
    }
    int nLayers = findNlayers(BLlayer.hFirstLayer, BLlayer.progress, thick, BLlayer.maxLayerh);
    std::vector<double> rad;
    double ds = 2./nLayers;
    setRadiusLayers(nLayers);
    std::vector<std::vector<double>> e0, e1;
    for(int j=0; j<nLayers; ++j) {
      rad.push_back(radiusEdge(-1. + j*ds)[0]);
      std::vector<double> p0 = AddVect(1., bndpts[i], rad[j], norms[i]);
      std::vector<double> p1 = AddVect(1., bndpts[i+1], rad[j], norms[i+1]);
      e0.push_back(p0);
      e1.push_back(p1);
    }
    std::vector<std::vector<std::vector<double>>> edges(4);
    edges[0].push_back(bndpts[i]); edges[0].push_back(bndpts[i+1]);
    edges[1] = e1;
    edges[2].push_back(e1[e1.size()-1]);edges[2].push_back(e0[e0.size()-1]);
    edges[3] = e0;
    RectRegion pic(edges, "p");
    pic.MeshGen(1, e0.size()-1);
    region.AddRegion(pic);
  }
  region.outXml(name+".xml");
  std::vector<int> comp2;
  comp2.push_back(0);
  region.outCOMPO(name+".xml", comp2);
  return 0;
}