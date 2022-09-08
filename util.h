#ifndef UTIL_H
#define UTIL_H
#include <vector>
#include <cstddef>
double distance(std::vector<double> p0, std::vector<double> p1);
double distance(std::vector<std::vector<double> > p);
bool NanString(const char *buffer, size_t len);

#endif