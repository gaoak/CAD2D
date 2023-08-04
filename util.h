//
#ifndef UTIL_H
#define UTIL_H
#include <cstddef>
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
int findNlayers(double h, double q, double R, double m);
#endif