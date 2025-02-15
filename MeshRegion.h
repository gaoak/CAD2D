#ifndef MESHREGION_H
#define MESHREGION_H
#include <map>
#include <set>
#include <string>
#include <vector>
#include "tinyxml2.h"
#define TRIGELEMENT 2
#define QUADELEMENT 3
#define LINEBUFFERSIZE 1000
// only topology structures
// no physics
class MeshRegion {
public:
  MeshRegion(std::string name, double tolerance = 1.E-6);
  MeshRegion(std::string name, std::vector<std::vector<double>> &p,
             bool sort = true, double tolerance = 1.E-6);
  double m_tolerance;
  std::string m_name;
  std::vector<std::vector<double>> m_pts;
  std::set<int> m_bndPts;
  std::vector<std::vector<int>> m_edges;
  std::vector<std::vector<int>> m_cells;
  int getCellsNumber();
  void transformation(double AoA, double x0, double y0);
  void sortCellFromQtoT();
  int loadFromMsh(std::string filename,
                  double maxInnerAngle = 2.35619449019234);
  int loadFromXml(std::string filename, std::string bodytype);
  int AddElement(std::vector<std::vector<double>> pts);
  std::vector<std::vector<int>> extractBoundaryPoints();
  int pointIsExist(std::vector<double> p, int &pId);
  void rebuildEdgesIndex();
  std::vector<std::vector<int>> m_unSharedPts;
  std::vector<std::set<int>> m_unSharedPtsSet;
  std::map<std::set<int>, int> m_edgesIndex;
  bool consistancyCheck(MeshRegion m);
  bool consistancyCheck(std::vector<std::vector<double>> &pts);
  void transform(std::vector<double> &p, double AoA, double x0, double y0);
  std::vector<double> getLineFromVertex(std::vector<double> A,
                                        std::vector<double> B);
  std::vector<double> intersection(std::vector<double> l0,
                                   std::vector<double> l1);
  int ResetBndPts();
  void CheckMesh(double angle = 75. / 180. * 3.1415926);
  void FixMesh();
  void GetBoundBox(std::vector<double> &box);
  int RemoveElements(void *ptsFunc);
  std::vector<std::vector<int>> splitBoundaryPts(double angle);

private:
  void XmlLoadModifyPts(tinyxml2::XMLDocument &doc, std::string bodytype, std::map<int, int> &ptsmap);
  void XmlLoadEdge(tinyxml2::XMLDocument &doc, std::map<int, int> &ptsmap, std::map<int, int> &edgemap);
  void XmlLoadElement(tinyxml2::XMLDocument &doc, std::map<int, int> &edgemap);
  int loadNode(std::ifstream &inxml, int N, char buffer[]);
  int loadElements(std::ifstream &inxml, int N, char buffer[],
                   std::set<int> types, double maxInnerAngle);
  int AddSplitElemens(std::vector<int> pts, double maxInnerAngle);
  void CalculateCosAngle(std::vector<int> pts, std::vector<double> &cangle,
                         int &maxc, int &minabs);
  void GetFacePts(int index, std::vector<int> &pts);
  double ElementArea(int index);
  std::vector<std::vector<int>>
  splitBoundaryPts(std::vector<std::vector<int>> &allbndpts, double angle);
  std::map<int, int> m_pIDf2s;
};
#endif // MESHREGION_H
