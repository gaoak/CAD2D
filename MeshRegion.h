#ifndef MESHREGION_H
#define MESHREGION_H
#include<vector>
#include<set>
#include<map>
#include<string>
#define TRIGELEMENT 2
#define QUADELEMENT 3
#define LINEBUFFERSIZE 1000
//only topology structures
//no physics
class MeshRegion {
public:
    MeshRegion(std::string name, double tolerance = 1.E-6);
    MeshRegion(std::string name, std::vector<std::vector<double> > &p, bool sort = true, double tolerance = 1.E-6);
    double m_tolerance;
    std::string m_name;
    std::vector<std::vector<double> > m_pts;
    std::set<int> m_bndPts;
    std::vector<std::vector<int> > m_edges;
    std::vector<std::vector<int> > m_cells;
    int getCellsNumber();
    void transformation(double AoA);
    void sortCellFromQtoT();
    int loadFromMsh(std::string filename, double maxInnerAngle = 2.35619449019234);
    int addTriangle(std::vector<std::vector<double>> pts);
    std::vector<std::vector<int>> extractBoundary();
    int pointIsExist(std::vector<double> p, int &pId);
    void rebuildEdgesIndex();
    std::vector<std::vector<int>> m_unSharedPts;
    std::vector<std::set<int>> m_unSharedPtsSet;
    std::map<std::set<int>, int> m_edgesIndex;
    bool consistancyCheck(MeshRegion m);
    void transform(std::vector<double> &p, double AoA);
    std::vector<double> getLineFromVertex(std::vector<double> A, std::vector<double>B);
    std::vector<double> intersection(std::vector<double> l0, std::vector<double> l1);
private:
    int loadNode(std::ifstream &inxml, int N, char buffer[]);
    int loadElements(std::ifstream &inxml, int N, char buffer[], std::set<int> types, double maxInnerAngle);
    int AddSplitElemens(std::vector<int> pts, double maxInnerAngle);
    void CalculateCosAngle(std::vector<int> pts, std::vector<double> &cangle, int & maxc, int & minabs);
    std::map<int, int> m_pIDf2s;
};
#endif // MESHREGION_H
