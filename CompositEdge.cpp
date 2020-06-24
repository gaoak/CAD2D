#include <vector>
#include <cmath>
#include <iostream>

#include "CompositEdge.h"

CompositEdge::CompositEdge() {
    m_Nedges = 0;
    m_N = 0;
}

void CompositEdge::addEdge(LineEdge edge, void* edgefunction){
    m_edges.push_back(edge);
    m_functions.push_back(edgefunction);
    if(m_Nedges==0) {
        m_Ns.push_back(edge.m_N);
    }else {
        m_Ns.push_back(m_Ns[m_Nedges-1] + edge.m_N);
    }
    ++m_Nedges;
    m_N += edge.m_N;
}
std::vector<double> CompositEdge::Evaluate(double s) {
    double index = 0.5*(1.+s)*m_Ns[m_Nedges-1];
	if(index<0.) index = 0.;
	if(index>m_Ns[m_Nedges-1]) index = double(m_Ns[m_Nedges-1])-1.E-15;
	int i;
    double indexCut = index;
    for(i=0; i<m_Nedges; ++i) {
        if(index<=m_Ns[i]) break;
        indexCut = index - m_Ns[i];
    }
    std::vector<double>(*edgeFunction)(double) = (std::vector<double>(*)(double))m_functions[i];
    return edgeFunction(2.*indexCut/double(m_edges[i].m_N)-1.);
}
