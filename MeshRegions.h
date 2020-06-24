#ifndef MESHREGIONS_H
#define MESHREGIONS_H
#include<vector>
#include<set>
#include<map>
#include<string>
#include"MeshRegion.h"
//regions manipulation
//related with physics
class MeshRegions : public MeshRegion {
public:
    MeshRegions(std::string name, double tolerance);
    int AddRegion(const MeshRegion &region);
    int outXml(std::string filename);
    int defineBoundary(void* edgeFun, int N, int bndID, int Ncurve = 2, double AoA = 0., int direction = 1);
    int outCOMPO(std::string filename, std::vector<int> comps);
    void outOuterRegion(std::string filename, std::vector<std::vector<double>> box, std::vector<double> center,double radius, bool exclude);
private:
    std::map<int, std::vector<int>> m_boundary;
    void findAllBoundaryEdges();
    std::set<int> m_allBoundaryEdges;
    std::vector<int> m_curvedEdges;
    std::vector<std::vector<double>> m_curvedPoints;
};
#endif
