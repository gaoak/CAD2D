#include "util.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>

double distance(std::vector<double> p0, std::vector<double> p1) {
  double res = 0.;
  for (size_t i = 0; i < p0.size() && i < p1.size(); ++i) {
    res += (p0[i] - p1[i]) * (p0[i] - p1[i]);
  }
  return sqrt(res);
}

double distance(std::vector<std::vector<double>> p) {
  double res = 0.;
  for (size_t i = 1; i < p.size(); ++i) {
    res += distance(p[i - 1], p[i]);
  }
  return res;
}

bool NanString(const char *buffer, size_t N) {
  for (size_t i = 0; i < N && buffer[i]; ++i) {
    if (buffer[i] != ' ') {
      return (buffer[i] < '0' || buffer[i] > '9');
    }
  }
  return true;
}

void OutGeo(std::string filename, std::vector<std::vector<double>> outerbox,
            std::vector<std::vector<std::vector<double>>> innerbnd) {
  std::ofstream outgmsh(filename.c_str());
  int index = 1, j;
  // inner boundary points
  std::vector<std::vector<int>> line;
  for (int i = 0; i < innerbnd.size(); ++i) {
    int indexs = index;
    for (j = 0; j < innerbnd[i].size(); ++j) {
      outgmsh << std::scientific << std::setprecision(17);
      outgmsh << "Point(" << index << ") = {" << std::setw(26)
              << innerbnd[i][j][0] << ", " << innerbnd[i][j][1] << ", 0};\n";
      std::vector<int> l;
      l.push_back(index);
      if (j < innerbnd[i].size() - 1) {
        l.push_back(index + 1);
      } else {
        l.push_back(indexs);
      }
      line.push_back(l);
      ++index;
    }
  }
  // outer boundary points
  if (outerbox.size() >= 3) {
    int indexs = index;
    for (j = 0; j < outerbox.size(); ++j) {
      outgmsh << std::scientific << std::setprecision(17);
      outgmsh << "Point(" << index << ") = {" << std::setw(26) << outerbox[j][0]
              << ", " << outerbox[j][1] << ", 0};\n";
      std::vector<int> l;
      l.push_back(index);
      if (j < outerbox.size() - 1) {
        l.push_back(index + 1);
      } else {
        l.push_back(indexs);
      }
      line.push_back(l);
      ++index;
    }
  }
  // lines
  int linestart = index;
  for (j = 0; j < line.size(); ++j) {
    outgmsh << "Line(" << index << ") = {" << line[j][0] << ", " << line[j][1]
            << "};\n";
    ++index;
  }
  // line loops
  int count = 0;
  std::vector<std::vector<int>> lineloop;
  for (int i = 0; i < innerbnd.size(); ++i) {
    std::vector<int> c;
    lineloop.push_back(c);
    for (j = 0; j < innerbnd[i].size(); ++j) {
      lineloop[i].push_back(count + linestart);
      ++count;
    }
  }
  std::vector<int> c;
  for (j = 0; j < outerbox.size(); ++j) {
    c.push_back(count + linestart);
    ++count;
  }
  if (c.size() > 0) {
    lineloop.push_back(c);
  }
  // lineloop
  std::vector<int> surf;
  for (int i = 0; i < lineloop.size(); ++i) {
    outgmsh << "Line Loop(" << index << ") = {";
    for (int j = 0; j < lineloop[i].size() - 1; ++j) {
      outgmsh << lineloop[i][j] << ", ";
    }
    outgmsh << lineloop[i][lineloop[i].size() - 1] << "};\n";
    surf.push_back(index);
    ++index;
  }
  // surf
  outgmsh << "Plane Surface(" << index << ") = {";
  for (int i = 0; i < surf.size() - 1; ++i) {
    outgmsh << surf[i] << ", ";
  }
  outgmsh << surf[surf.size() - 1] << "};\n";
  outgmsh << "Recombine Surface {" << index << "};\n";
  std::cout << "Output file " << filename << std::endl;
}