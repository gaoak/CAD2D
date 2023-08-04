#ifndef MESHTOOL_H
#define MESHTOOL_H
#include<vector>
#include "MeshRegions.h"
void setRadiusMesh(double h, double p, double m);
std::vector<double> radiusEdge(double s);
void setRadiusLayers(int n);
int meshingBoundayLayer(MeshRegions &region, void* thickFunc, void * edges,
                        std::string name, std::vector<std::vector<double>> trimnorm);
#endif