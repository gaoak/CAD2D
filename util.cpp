#include "util.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
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
      double zcoord = 0.;
      if (outerbox[j].size() > 2) {
        zcoord = outerbox[j][2];
      }
      outgmsh << std::scientific << std::setprecision(17);
      outgmsh << "Point(" << index << ") = {" << std::setw(26) << outerbox[j][0]
              << ", " << outerbox[j][1] << ", " << zcoord << "};\n";
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

std::vector<double> intersect(std::vector<double> p0, std::vector<double> n0,
                              std::vector<double> p1, std::vector<double> n1) {
  std::vector<double> res(2, 1E38);
  double jac = -n0[0] * n1[1] + n1[0] * n0[1];
  if (fabs(jac) < 1E-38) {
    return res;
  }
  double b0 = p1[0] - p0[0], b1 = p1[1] - p0[1];
  jac = 1. / jac;
  res[0] = jac * (-n1[1] * b0 + n1[0] * b1);
  res[1] = jac * (-n0[1] * b0 + n0[0] * b1);
  return res;
}

std::vector<double> AddVect(const double a1, const std::vector<double> &a,
                            const double b1, const std::vector<double> &b) {
  size_t n = std::min(a.size(), b.size());
  std::vector<double> res(n);
  for (size_t i = 0; i < n; ++i) {
    res[i] = a1 * a[i] + b1 * b[i];
  }
  return res;
}

int findNlayers(double h, double q, double R, double m) {
  int n = 0;
  double len = 0;
  double delta = h;
  for (n = 1; n <= 1000000; ++n) {
    if (delta >= m)
      delta = m;
    len += delta;
    if (len >= R)
      return n;
    delta *= q;
  }
  return n;
}

// relation :A in C in D; B in C in D
// [A; empty]
// [B; empty]
// [C; A, B];
// [D, C];
void BuildTopoTree(std::vector<std::vector<std::vector<double>>> &pts,
                   std::map<int, std::set<int>> &trees, std::set<int> &roots) {
  trees.clear();
  roots.clear();
  std::vector<std::vector<double>> boxes;
  for (size_t i = 0; i < pts.size(); ++i) {
    boxes.push_back(GetBoundingBox(pts[i]));
  }
  // buid the relation tree
  for (size_t p = 0; p < pts.size(); ++p) {
    trees[p] = std::set<int>();
  }
  std::set<int> sons;
  for (size_t p0 = 0; p0 + 1 < pts.size(); ++p0) {
    for (size_t p1 = p0 + 1; p1 < pts.size(); ++p1) {
      int r = FindRelationByBoundBox(boxes[p0], boxes[p1]);
      if (r == 1) {
        trees[p0].insert(p1);
        sons.insert(p1);
      } else if (r == -1) {
        trees[p1].insert(p0);
        sons.insert(p0);
      }
    }
  }
  // trim the tree
  std::map<int, std::set<int>> toremove;
  for (size_t p = 0; p < trees.size(); ++p) {
    std::vector<int> list(trees[p].begin(), trees[p].end());
    std::set<int> rm;
    for (size_t i1 = 0; i1 + 1 < list.size(); ++i1) {
      for (size_t i2 = i1 + 1; i2 < list.size(); ++i2) {
        if (trees[list[i2]].find(list[i1]) != trees[list[i2]].end()) {
          rm.insert(list[i1]);
        } else if (trees[list[i1]].find(list[i2]) != trees[list[i1]].end()) {
          rm.insert(list[i2]);
        }
      }
    }
    toremove[p] = rm;
  }
  for (size_t p = 0; p < trees.size(); ++p) {
    if (sons.find(p) == sons.end()) {
      roots.insert(p);
    }
    for (auto p1 : toremove[p]) {
      trees[p].erase(p1);
    }
  }
}

// 1, 0 contains 1; -1, 1 contains 0; 0 overlap; 2, no contact.
int FindRelationByBoundBox(std::vector<double> &box0,
                           std::vector<double> &box1) {
  int dim = (int)box0.size() / 2;
  int sum = 0;
  for (int i = 0; i < dim; ++i) {
    if (box0[2 * i] < box1[2 * i] && box1[2 * i + 1] < box0[2 * i + 1]) {
      sum += 1;
    } else if (box1[2 * i] < box0[2 * i] && box0[2 * i + 1] < box1[2 * i + 1]) {
      sum -= 1;
    } else if (box0[2 * i + 1] < box1[2 * i] || box1[2 * i + 1] < box0[2 * i]) {
      return 2;
    }
  }
  if (sum == dim) {
    return 1;
  } else if (sum == -dim) {
    return -1;
  } else {
    return 0; // to do improve this using Argz check
  }
}

std::vector<double> GetBoundingBox(std::vector<std::vector<double>> &pts) {
  std::vector<double> res;
  if (pts.size() == 0) {
    return res;
  }
  size_t dim = pts[0].size();
  res.resize(2 * dim);
  for (size_t i = 0; i < dim; ++i) {
    res[2 * i] = std::numeric_limits<double>::max();
    res[2 * i + 1] = std::numeric_limits<double>::lowest();
  }
  for (auto &p : pts) {
    for (size_t i = 0; i < dim; ++i) {
      if (res[2 * i] > p[i]) {
        res[2 * i] = p[i];
      }
      if (res[2 * i + 1] < p[i]) {
        res[2 * i + 1] = p[i];
      }
    }
  }
  return res;
}

void FindTreesDepths(int root, int lroot, std::map<int, std::set<int>> &trees,
                     std::vector<int> &levels) {
  levels[root] = lroot;
  for (auto p : trees[root]) {
    FindTreesDepths(p, 1 + lroot, trees, levels);
  }
}

void parserDouble(const char *cstr, std::vector<double> &value) {
  value.clear();
  std::vector<int> digs;
  std::vector<int> dige;
  int i = 0;
  int flag = 0; // digit chunk
  while (1) {
    if ((cstr[i] >= '0' && cstr[i] <= '9') || cstr[i] == '.' ||
        cstr[i] == 'e' || cstr[i] == 'E' || cstr[i] == '+' || cstr[i] == '-') {
      if (flag == 0) {
        digs.push_back(i);
      }
      flag = 1;
    } else {
      if (flag == 1) {
        dige.push_back(i);
      }
      flag = 0;
    }
    if (cstr[i] == 0)
      break;
    ++i;
  }
  double k;
  for (int i = 0; i < digs.size(); ++i) {
    std::string cuts(cstr + digs[i], dige[i] - digs[i]);
    if (sscanf(cuts.c_str(), "%lf", &k) < 1) {
      std::cout << "error: parser double " << cuts << std::endl;
    }
    value.push_back(k);
  }
}