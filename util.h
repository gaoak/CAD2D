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
#endif