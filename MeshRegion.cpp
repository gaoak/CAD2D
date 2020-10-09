#include"MeshRegion.h"
#include<cmath>
#include<fstream>
#include<iostream>
#include<string>
#include<vector>
#include<map>
#include<set>
MeshRegion::MeshRegion(std::string name, double tolerance)
{
    m_name = name;
    m_tolerance = tolerance;
}

MeshRegion::MeshRegion(std::string name, std::vector<std::vector<double> > &p, bool sort, double tolerance)
{
    m_name = name;
    m_tolerance = tolerance;
    int n = p.size();
    if(n<3 || n>4) {
        std::cout << "error: illegal call of MeshRegion using points." << std::endl;
        return;
    }
    std::vector<int> cell, ptsvec;
    std::vector<double> cangle;
    int maxindex, minabsindex;
    for(int i=0; i<n; ++i) {
        m_pts.push_back(p[i]);
        std::vector<int> e;
        e.push_back(i); e.push_back((i+1)%n);
        m_edges.push_back(e);
        ptsvec.push_back(i);
    }
    CalculateCosAngle(ptsvec, cangle, maxindex, minabsindex);
    if(!sort) minabsindex = 0;
    for(int j= minabsindex; j<n+minabsindex; ++j) {
        cell.push_back(j%n);
    }
    m_cells.push_back(cell);
    ResetBndPts();
}

void MeshRegion::transform(std::vector<double> &p, double AoA) {
    double x = p[0], y = p[1];
    p[0] = x*cos(AoA) + y*sin(AoA);
    p[1] =-x*sin(AoA) + y*cos(AoA);
}

void MeshRegion::transformation(double AoA) {
    for(int i=0; i<m_pts.size(); ++i) {
        std::vector<double> p = m_pts[i];
        transform(p, AoA);
        m_pts[i][0] = p[0];
        m_pts[i][1] = p[1];
    }
}

std::vector<double> MeshRegion::getLineFromVertex(std::vector<double> A, std::vector<double>B){
    std::vector<double> l;
    l.push_back(-(A[1] - B[1]));
    l.push_back( (A[0] - B[0]));
    l.push_back( A[0]*B[1] - A[1]*B[0]);
    return l;
}
std::vector<double> MeshRegion::intersection(std::vector<double> l0, std::vector<double> l1){
    //line = (a,b,c), a x + b y = c
    std::vector<double> res;
    double Jac = l0[0]*l1[1] - l0[1]*l1[0];
    if(fabs(Jac)<m_tolerance) {
        std::cout << "Error in region " << m_name << ", lines in parallel have no intersection" << std::endl;
    }
    res.push_back((l0[2]*l1[1] - l0[1]*l1[2])/Jac);
    res.push_back((l0[0]*l1[2] - l0[2]*l1[0])/Jac);
    return res;
}

int MeshRegion::loadFromMsh(std::string filename, double maxInnerAngle) {
    std::ifstream inxml(filename.c_str());
    if(!inxml.is_open()) return -1;
    char buffer[LINEBUFFERSIZE];
    std::string strNode("$Nodes");
    std::string strElem("$Elements");
    std::string strVers("$MeshFormat");
    inxml.getline(buffer, LINEBUFFERSIZE);
    while(!inxml.eof()) {
        if(strNode.compare(buffer) == 0) {
            inxml.getline(buffer, LINEBUFFERSIZE);
            int nNode = 0;
            std::sscanf(buffer, "%d", &nNode);
            loadNode(inxml, nNode, buffer);
        }
        if(strElem.compare(buffer) == 0) {
            inxml.getline(buffer, LINEBUFFERSIZE);
            int nElem = 0;
            std::sscanf(buffer, "%d", &nElem);
            std::set<int> elemtypes;
            elemtypes.insert(QUADELEMENT);
            elemtypes.insert(TRIGELEMENT);
            loadElements(inxml, nElem, buffer, elemtypes, maxInnerAngle);
        }
        if(strVers.compare(buffer) == 0) {
            inxml.getline(buffer, LINEBUFFERSIZE);
            if(buffer[0]>'3') {
                std::cout << "Error: in region " << m_name << ", version too high " << std::endl;
                return -2;
            }
        }
        inxml.getline(buffer, LINEBUFFERSIZE);
    }
    ResetBndPts();
    return 0;
}

int MeshRegion::ResetBndPts() {
    m_bndPts.clear();
    extractBoundary();
    for(int i=0; i<m_unSharedPts.size(); ++i) {
        for(int j=0; j<m_unSharedPts[i].size(); ++j) {
            m_bndPts.insert(m_unSharedPts[i][j]);
        }
    }
    return 0;
}

int MeshRegion::loadNode(std::ifstream &inxml, int N, char buffer[]) {
    for(int i=0; i<N; ++i) {
        inxml.getline(buffer, LINEBUFFERSIZE);
        double x, y, z;
        int ip;
        sscanf(buffer, "%d%lf%lf%lf", &ip, &x, &y, &z);
        m_pIDf2s[ip] = i;
        std::vector<double> p;
        p.push_back(x); p.push_back(y);
        m_pts.push_back(p);
    }
    return 0;
}
int MeshRegion::loadElements(std::ifstream &inxml, int N, char buffer[], std::set<int> types, double maxInnerAngle) {
    for(int i=0; i<N; ++i) {
        inxml.getline(buffer, LINEBUFFERSIZE);
        int elemType = 0, ntag = 0, ie = 0;
        int val[10];
        sscanf(buffer, "%d%d%d", &ie, &elemType, &ntag);
        if(types.find(elemType) != types.end()) {
            sscanf(buffer, "%d%d%d%d%d%d%d%d%d", val, val+1, val+2, val+3, val+4, val+5, val+6, val+7, val+8);
        }
        else {
            continue;
        }
        std::vector<int> pointseries;
        for(int j= ntag+3+elemType; j>=ntag+3; --j) {
            pointseries.push_back(m_pIDf2s[val[j]]);
        }
        AddSplitElemens(pointseries, maxInnerAngle);
    }
    return 0;
}

void MeshRegion::CalculateCosAngle(std::vector<int> pts, std::vector<double> &cangle,
                                  int & maxc, int & minabs) {
    double maxcangle = -1E21;
    double minabscangle = 1E21;
    int n = pts.size();
    cangle = std::vector<double>(n);
    for(int i=0; i<n; ++i) {
        double xm1 = m_pts[pts[(i+n-1)%n]][0];
        double ym1 = m_pts[pts[(i+n-1)%n]][1];
        double x   = m_pts[pts[(i    )%n]][0];
        double y   = m_pts[pts[(i    )%n]][1];
        double xp1 = m_pts[pts[(i+1  )%n]][0];
        double yp1 = m_pts[pts[(i+1  )%n]][1];
        double nmx = (x-xm1)/sqrt((x-xm1)*(x-xm1)+(y-ym1)*(y-ym1));
        double nmy = (y-ym1)/sqrt((x-xm1)*(x-xm1)+(y-ym1)*(y-ym1));
        double npx = (xp1-x)/sqrt((xp1-x)*(xp1-x)+(yp1-y)*(yp1-y));
        double npy = (yp1-y)/sqrt((xp1-x)*(xp1-x)+(yp1-y)*(yp1-y));
        cangle[i] = nmx*npx + nmy*npy;
        if(cangle[i]>maxcangle) {
            maxcangle = cangle[i];
            maxc = i;
        }
        if(fabs(cangle[i])<minabscangle) {
            minabscangle = fabs(cangle[i]);
            minabs = i;
        }
    }
}

int MeshRegion::AddSplitElemens(std::vector<int> pts, double maxInnerAngle) {
    std::vector<double> cangle;
    int maxindex, minabsindex;
    CalculateCosAngle(pts, cangle, maxindex, minabsindex);
    std::vector<std::vector<int> > ptsvec;
    int n = pts.size();
    if(n == 4 && cos(maxInnerAngle)<0 && cangle[maxindex] > -cos(maxInnerAngle)) {
        std::vector<int> p1;
        std::vector<int> p2;
        for(int i=0; i<3; ++i) {
            p1.push_back(pts[(maxindex+i)%n]);
            p2.push_back(pts[(maxindex+i+2)%n]);
        }
        ptsvec.push_back(p1);
        ptsvec.push_back(p2);
    } else {
        ptsvec.push_back(pts);
    }
    for(int i=0; i<ptsvec.size(); ++i) {
        std::vector<int> cell;
        n = ptsvec[i].size();
        CalculateCosAngle(ptsvec[i], cangle, maxindex, minabsindex);
        for(int j= minabsindex; j<n+minabsindex; ++j) {
            std::set<int> p;
            std::vector<int> pv;
            p.insert(    ptsvec[i][j%n]); p.insert(    ptsvec[i][(j+1)%n]);
            pv.push_back(ptsvec[i][j%n]); pv.push_back(ptsvec[i][(j+1)%n]);
            if(m_edgesIndex.find(p) == m_edgesIndex.end()) {
                m_edgesIndex[p] = m_edges.size();
                m_edges.push_back(pv);
            }
            cell.push_back(m_edgesIndex[p]);
        }
        m_cells.push_back(cell);
    }
    return ptsvec.size();
}

bool MeshRegion::consistancyCheck(MeshRegion m) {
    int tmp;
    std::vector<std::vector<int>> pts = m.extractBoundary();
    for(int i=0; i<pts.size(); ++i) {
        bool exist = pointIsExist(m.m_pts[pts[i][0]], tmp);
        for(int j=1; j<pts.size(); ++j)
            if(exist != pointIsExist(m.m_pts[pts[i][j]], tmp)) {
                std::cout << "Error: between regions " << m_name << " and " << m.m_name  << ", mismatch points (" << m.m_pts[pts[i][0]][0] << ", " << m.m_pts[pts[i][0]][1] << "),\n";
                return false;
            }
    }
    return true;
}

void MeshRegion::sortCellFromQtoT(){
    std::vector<std::vector<int> > tmpcells;
    for(int i=0; i<m_cells.size(); ++i) {
        if(m_cells[i].size()==4) tmpcells.push_back(m_cells[i]);
    }
    for(int i=0; i<m_cells.size(); ++i) {
        if(m_cells[i].size()==3) tmpcells.push_back(m_cells[i]);
    }
    m_cells = tmpcells;
}

int MeshRegion::getCellsNumber() {
    return m_cells.size();
}

void MeshRegion::rebuildEdgesIndex(){
    m_edgesIndex.clear();
    for(int i=0; i<m_edges.size(); ++i) {
        std::set<int> p;
        p.insert(m_edges[i][0]);
        p.insert(m_edges[i][1]);
        m_edgesIndex[p] = i;
    }
}

std::vector<std::vector<int>> MeshRegion::extractBoundary() {
    m_unSharedPts.clear();
    m_unSharedPtsSet.clear();
    std::vector<int> edgeCount(m_edges.size(), 0);
    for(int i=0; i<m_cells.size(); ++i) {
        for(int k=0; k<m_cells[i].size(); ++k) {
            edgeCount[m_cells[i][k]] += 1;
        }
    }
    std::map<int, std::vector<int>> ptsToUnSharedEdges;
    std::set<int> pts;
    for(int i=0; i<edgeCount.size(); ++i) {
        if(1 == edgeCount[i]) {
            ptsToUnSharedEdges[m_edges[i][0]].push_back(i);
            ptsToUnSharedEdges[m_edges[i][1]].push_back(i);
            pts.insert(m_edges[i][0]);
            pts.insert(m_edges[i][1]);
        }else {
            if(edgeCount[i]!=2) std::cout << "Warning: in region " << m_name <<  ", edge " << i << ", times used by cells " << edgeCount[i] << std::endl;
        }
    }
    for(auto it=ptsToUnSharedEdges.begin(); it!=ptsToUnSharedEdges.end(); ++it) {
        if(it->second.size()!=2) {
            std::cout << "Warning: in region" << m_name <<  ", pts " << (it->first) << ", = (" << m_pts[(it->first)][0] << ", " << m_pts[(it->first)][1] << "), on edge " << it->second[0] << "only used once." << std::endl;
        }
    }
    //std::cout << "Extract " << pts.size() << " points on the boundary in region " << m_name << std::endl;
    int ptsIndex = -1;
    int nbnds = -1;
    while(pts.size()>0) {
        if(ptsIndex==-1) {
            ptsIndex = *(pts.begin());
            std::vector<int> b0;
            std::set<int>    bs0;
            b0.push_back(ptsIndex);
            bs0.insert(ptsIndex);
            m_unSharedPts.push_back(b0);
            m_unSharedPtsSet.push_back(bs0);
            nbnds++;
            pts.erase(ptsIndex);
        }else {
            int tmpindex;
            if(ptsIndex != m_edges[ptsToUnSharedEdges[ptsIndex][0]][0]) {
                tmpindex = m_edges[ptsToUnSharedEdges[ptsIndex][0]][0];
            }else {
                tmpindex = m_edges[ptsToUnSharedEdges[ptsIndex][0]][1];
            }
            if(m_unSharedPtsSet[nbnds].find(tmpindex)==m_unSharedPtsSet[nbnds].end()) {
                m_unSharedPts[nbnds].push_back(tmpindex);
                m_unSharedPtsSet[nbnds].insert(tmpindex);
                pts.erase(tmpindex);
                ptsIndex = tmpindex;
                continue;
            }
            if(ptsIndex != m_edges[ptsToUnSharedEdges[ptsIndex][1]][0]) {
                tmpindex = m_edges[ptsToUnSharedEdges[ptsIndex][1]][0];
            }else {
                tmpindex = m_edges[ptsToUnSharedEdges[ptsIndex][1]][1];
            }
            if(m_unSharedPtsSet[nbnds].find(tmpindex)==m_unSharedPtsSet[nbnds].end()) {
                m_unSharedPts[nbnds].push_back(tmpindex);
                m_unSharedPtsSet[nbnds].insert(tmpindex);
                pts.erase(tmpindex);
                ptsIndex = tmpindex;
                continue;
            }
            ptsIndex = -1;
        }
    }
    //std::cout << "Extract " << m_unSharedPts.size() << " unconnected boundaries in region " << m_name << std::endl;
    return m_unSharedPts;
}

int MeshRegion::pointIsExist(std::vector<double> p, int &pId) {
    for(auto it = m_bndPts.begin(); it!=m_bndPts.end(); ++it) {
        pId = *it;
        if(0.5*(fabs(m_pts[pId][0] - p[0]) + fabs(m_pts[pId][1] - p[1]))<m_tolerance) {
            return 1;
        }
    }
    return 0;
}
