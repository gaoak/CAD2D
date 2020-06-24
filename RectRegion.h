#ifndef RECTREGION_H
#define RECTREGION_H
#include<vector>
#include<set>
#include"MeshRegion.h"

enum EdgeType {
    eContinuous,
    eDiscrete
};

enum MeshType {
    eIsoparametric,
    eBoundaryLayer0,
    eStraightLine
};

class RectRegion : public MeshRegion {
public:
    RectRegion(std::vector<void*> edges, std::string name, bool connectivityCheck = true, double tolerance = 1.E-6, int edgeAttractMesh = 0, double attractEps = 1.);
    RectRegion(std::vector<std::vector<std::vector<double> > > edges, std::string name, bool connectivityCheck = true, double tolerance = 1.E-6, int edgeAttractMesh = 0, double attractEps = 1.);
    std::vector<double> EvaluateEdgePts(int i, double s);
    void PrintEdge(int eID, int N = 11);
    int MeshGen(int M, int N, MeshType method = eStraightLine, bool trimWakeSlope = false);
    void Tec360Pts(std::string fileName);
    int m_M; // number of elements
    int m_N; // number of elements
    std::vector<void*> m_edgesFun;
    std::vector<std::vector<std::vector<double> > > m_edgesPts;
    std::vector<double> m_edgesDirec;
    std::vector<std::vector<double>> m_vertex;
    std::vector<double> getVertex(int i);
private:
    EdgeType m_edgeType;
    int CheckConnectivity();
    int EdgeAttractMesh(double x0, double y0, double &x, double &y);
    int m_edgeAttractMesh;
    double m_attractEps;
    int ptsByIsoParametric();
    int ptsByBoundaryLayer(bool trimWakeSlope);
    int ptsByStraightLine();
    std::vector<double> EvaluateEdgePtsDerivOneSide(int i, double s);
};
#endif // RECTREGION_H
