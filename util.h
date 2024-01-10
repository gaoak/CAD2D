//
#ifndef UTIL_H
#define UTIL_H
#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <vector>
double distance(std::vector<double> p0, std::vector<double> p1);
double distance(std::vector<std::vector<double>> p);
bool NanString(const char *buffer, size_t len);
void OutGeo(std::string filename, std::vector<std::vector<double>> outerbox,
            std::vector<std::vector<std::vector<double>>> innerbnd);
std::vector<double> intersect(std::vector<double> p0, std::vector<double> n0,
                              std::vector<double> p1, std::vector<double> n1);
std::vector<double> AddVect(const double a1, const std::vector<double> &a,
                            const double b1, const std::vector<double> &b);
std::vector<double> GetBoundingBox(std::vector<std::vector<double>> &pts);
int findNlayers(double h, double q, double R, double m);
void FindTreesDepths(int root, int lroot, std::map<int, std::set<int>> &trees,
                     std::vector<int> &levels);
void BuildTopoTree(std::vector<std::vector<std::vector<double>>> &pts,
                   std::map<int, std::set<int>> &trees, std::set<int> &roots);
int FindRelationByBoundBox(
    std::vector<double> &box0,
    std::vector<double>
        &box1); // 1, 0 contains 1; -1, 1 contains 0; 0 overlap; 2, no contact.
#endif