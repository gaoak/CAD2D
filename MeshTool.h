#ifndef MESHTOOL_H
#define MESHTOOL_H
#include "MeshRegions.h"
#include <vector>
void setRadiusMesh(double h, double p, double m);
std::vector<double> radiusEdge(double s);
void setRadiusLayers(int n);
int meshingBoundayLayer(MeshRegions &region, int Nx, void *thickFunc,
                        void *edges, std::string name,
                        std::vector<std::vector<double>> trimnorm);
int outputGeo(MeshRegions &combinedReg, std::vector<int> OutLevels);
#endif