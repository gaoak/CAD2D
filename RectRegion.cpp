#include<iostream>
#include<vector>
#include<cmath>
#include<fstream>
#include"RectRegion.h"
RectRegion::RectRegion(std::vector<void*> edges, std::string name, bool connectivityCheck, double tolerance, int edgeAttractMesh, double attractEps) : MeshRegion(name, tolerance)
{
    m_edgeAttractMesh = edgeAttractMesh;
    m_attractEps = attractEps;
    if(m_attractEps>1.) m_attractEps = 1.;
    if(m_edgeAttractMesh<0) m_edgeAttractMesh = 0;
    if(m_edgeAttractMesh>3) m_edgeAttractMesh = 3;
    if(edges.size()<4) {
        std::cout << "error!: not enough edges" << std::endl;
    }
    for(int i=0; i<4; ++i) {
        m_edgesFun.push_back(edges[i]);
    }
    m_edgeType = eContinuous;
    m_edgesDirec = std::vector<double>(4, 1.);
    if(connectivityCheck && CheckConnectivity()!=0) {
        std::cout << "error!: connectivity check fail" << std::endl;
    }
}

RectRegion::RectRegion(std::vector<std::vector<std::vector<double> > > edges, std::string name, bool connectivityCheck, double tolerance, int edgeAttractMesh, double attractEps) : MeshRegion(name, tolerance)
{
    m_edgeAttractMesh = edgeAttractMesh;
    m_attractEps = attractEps;
    if(m_attractEps>1.) m_attractEps = 1.;
    if(m_edgeAttractMesh<0) m_edgeAttractMesh = 0;
    if(m_edgeAttractMesh>3) m_edgeAttractMesh = 3;
    if(edges.size()<4) {
        std::cout << "error!: not enough edges" << std::endl;
    }
    for(int i=0; i<4; ++i) {
        if(edges[i].size()<2) {
            std::cout << "error!: not enough points on edge " << i << std::endl;
        }
        m_edgesPts.push_back(edges[i]);
    }
    m_edgeType = eDiscrete;
    m_edgesDirec = std::vector<double>(4, 1.);
    if(connectivityCheck && CheckConnectivity()!=0) {
        std::cout << "error!: connectivity check fail" << std::endl;
    }
}

std::vector<double> RectRegion::EvaluateEdgePts(int i, double s)
{
    if(s<-1.) {
        s = -1.;
    }
    if(s>1) {
        s = 1.;
    }
    if(m_edgeType == eContinuous) {
        std::vector<double>(*edgeFunction)(double) = (std::vector<double>(*)(double))m_edgesFun[i];
        return edgeFunction(s*m_edgesDirec[i]);
    } else {
        int npts = m_edgesPts[i].size() - 1;
        double x = (s*m_edgesDirec[i] + 1.)*npts/2.;
        int nH = ceil (x);
        int nL = floor(x);
        if(nL<0)    nL = 0;
        if(nH>npts) nH = npts;
        if(nL==nH) {
            return m_edgesPts[i][nL];
        } else {
            double wL = nH - x;
            double wH = x - nL;
            std::vector<double> res;
            for(int j=0; j<m_edgesPts[i][nL].size(); ++j) {
                res.push_back( wH*m_edgesPts[i][nH][j] + wL*m_edgesPts[i][nL][j] );
            }
            return res;
        }
    }
}

int RectRegion::CheckConnectivity()
{
    std::vector<double> s1s = { 1., 1., -1., -1.};
    std::vector<double> p10, p11, p20, p21;
    for(int i=0; i<4; ++i) {
		p11 = EvaluateEdgePts(i,   s1s[i]);
		p20 = EvaluateEdgePts((i+1)%4, -1.);
		p21 = EvaluateEdgePts((i+1)%4, 1.);
		
		if(0.5*fabs(p11[0]-p20[0]) + 0.5*fabs(p11[1]-p20[1])<=m_tolerance ||
		   0.5*fabs(p11[0]-p21[0]) + 0.5*fabs(p11[1]-p21[1])<=m_tolerance) continue;
		
    	p10 = EvaluateEdgePts(i,  -s1s[i]);
		if(0.5*fabs(p10[0]-p20[0]) + 0.5*fabs(p10[1]-p20[1])<=m_tolerance ||
		   0.5*fabs(p10[0]-p21[0]) + 0.5*fabs(p10[1]-p21[1])<=m_tolerance) {
		    m_edgesDirec[i] = -1.;
		} else{
			std::cout << "error: point unconnect at edge (" << p10[0] << ", " << p10[1] << ":" << p11[0] << ", " << p11[1] << "), (" << p20[0] << ", " << p20[1] << ":" << p21[0] << ", " << p21[1] << std::endl;
			return -1;
		}
    }
    return 0;
}

void RectRegion::PrintEdge(int eID, int N)
{
    for(int i=0; i<N; ++i)
    {
        std::vector<double> pts = EvaluateEdgePts(eID, -1. + i*2./(N-1.));
        std::cout << i << ": " << pts[0] << ", " << pts[1] << std::endl;
    }
}

int RectRegion::ptsByIsoParametric() {
    // points generation
    //four vertexes + bounday points
    double dx = 2./m_M, dy = 2./m_N;
    std::vector<std::vector<double> > edge0, edge1, edge2, edge3;
    for(int i=0; i<=m_M; ++i) {
        edge0.push_back(EvaluateEdgePts(0, -1. + i*dx));
        edge2.push_back(EvaluateEdgePts(2, -1. + i*dx));
    }
    for(int j=0; j<=m_N; ++j) {
        edge1.push_back(EvaluateEdgePts(1, -1. + j*dy));
        edge3.push_back(EvaluateEdgePts(3, -1. + j*dy));
    }
    for(int j=0; j<=m_N; ++j) {
        double y0 = -1. + j*dy;
        for(int i=0; i<=m_M; ++i) {
            double x0 = -1. + i*dx;
            double x = 0., y = 0.;
            EdgeAttractMesh(x0, y0, x, y);
            std::vector<double> p(2);
            for(int k=0; k<2; ++k)
                p[k] = edge0[i][k]*(1.-y)*0.5 + edge2[i][k]*(1.+y)*0.5 + edge1[j][k]*(1.+x)*0.5 + edge3[j][k]*(1.-x)*0.5 - 0.25*(1.-x)*(1.-y)*edge0[0][k] - 0.25*(1.+x)*(1.-y)*edge0[m_M][k] - 0.25*(1.-x)*(1.+y)*edge2[0][k] - 0.25*(1.+x)*(1.+y)*edge2[m_M][k];
            m_pts.push_back(p);
            if(i==0 || j==0 || i==m_M || j==m_N) m_bndPts.insert(m_pts.size()-1);
        }
    }
    m_vertex.push_back(EvaluateEdgePts(0, -1.));
    m_vertex.push_back(EvaluateEdgePts(0, 1.));
    m_vertex.push_back(EvaluateEdgePts(2, 1.));
    m_vertex.push_back(EvaluateEdgePts(2, -1.));
    return m_pts.size();
}

int RectRegion::ptsByStraightLine() {
    // points generation
    //four vertexes + bounday points
    double dx = 2./m_M, dy = 2./m_N;
    std::vector<std::vector<double> > edge0, edge1, edge2, edge3;
    for(int i=0; i<=m_M; ++i) {
        edge0.push_back(EvaluateEdgePts(0, -1. + i*dx));
        edge2.push_back(EvaluateEdgePts(2, -1. + i*dx));
    }
    for(int j=0; j<=m_N; ++j) {
        edge1.push_back(EvaluateEdgePts(1, -1. + j*dy));
        edge3.push_back(EvaluateEdgePts(3, -1. + j*dy));
    }
    std::vector<std::vector<double> > linesx, linesy;
    for(int j=0; j<=m_N; ++j) {
        linesx.push_back(getLineFromVertex(edge3[j], edge1[j]));
    }
    for(int i=0; i<=m_M; ++i) {
        linesy.push_back(getLineFromVertex(edge0[i], edge2[i]));
    }
    for(int j=0; j<=m_N; ++j) {
        for(int i=0; i<=m_M; ++i) {
            m_pts.push_back(intersection(linesx[j], linesy[i]));
            if(i==0 || j==0 || i==m_M || j==m_N) m_bndPts.insert(m_pts.size()-1);
        }
    }
    m_vertex.push_back(EvaluateEdgePts(0, -1.));
    m_vertex.push_back(EvaluateEdgePts(0, 1.));
    m_vertex.push_back(EvaluateEdgePts(2, 1.));
    m_vertex.push_back(EvaluateEdgePts(2, -1.));
    return m_pts.size();
}

static double stretch(double eps, double h, double s, double R) {
    return s;
    double res = (R-1.)/eps;
    res = log(res)/log(h);
    std::cout << R << ", res " << res << std::endl;
    return s + eps*h*pow(s, res);
}

std::vector<double> RectRegion::EvaluateEdgePtsDerivOneSide(int i, double s)
{
    double dir = 1.;
    if(s>0.) dir = -1.;
    double h = 0.001*dir;
    std::vector<std::vector<double>> f;
    for(int j=0; j<7; ++j) {
        f.push_back(EvaluateEdgePts(i, s + j*h));
    }
    std::vector<double> res(2);
    for(int k=0; k<res.size(); ++k) res[k] = (-147.*f[0][k]+360.*f[1][k]-450.*f[2][k]+400.*f[3][k]-225.*f[4][k]+72.*f[5][k]-10.*f[6][k])/(60.*h);
    return res;
}

int RectRegion::ptsByBoundaryLayer(bool trimWakeSlope) {
    // points generation
    //four vertexes + bounday points
    double dx = 2./m_M, dy = 2./m_N;
    std::vector<std::vector<double> > edge0, norm0;
    std::vector<double> grows;
    for(int i=0; i<=m_M; ++i) {
        edge0.push_back(EvaluateEdgePts(0, -1. + i*dx));
    }
    double maxSlopx = 0., maxx = 0.;
    for(int i=0; i<=m_M; ++i) {
        double nf0 = 0., nf1 = 0., nb0 = 0., nb1 = 0., n0 = 0., n1 = 0.;
        if(i<m_M && i>0) {
            nf0 = edge0[i+1][0] - edge0[i][0];
            nf1 = edge0[i+1][1] - edge0[i][1];
            nb0 = edge0[i][0] - edge0[i-1][0];
            nb1 = edge0[i][1] - edge0[i-1][1];
        } else {
            std::vector<double> t = EvaluateEdgePtsDerivOneSide(0, -1. + i*dx);
            nf0 = t[0];
            nf1 = t[1];
        }
        n0 = nf0 + nb0;
        n1 = nf1 + nb1;
        std::vector<double> norm(2);
	    norm[0] =-n1/sqrt(n0*n0+n1*n1);
        norm[1] = n0/sqrt(n0*n0+n1*n1);
        if(norm[0]>maxSlopx) maxSlopx = norm[0];
        if(edge0[i][0]>maxx) maxx = edge0[i][0];
	    norm0.push_back(norm);
    }
    int concave1 = 0, concave2=m_N;
    double eps = 0.0001;
    for(int i=1; i<m_M; ++i) {
        if(norm0[i][0]>norm0[i-1][0]+eps && norm0[i][0]>norm0[i-1][0]+eps) {
            concave1 = i;
            break;
        }
    }
    for(int i=m_M-1; i>0; --i) {
        if(norm0[i][0]>norm0[i-1][0]+eps && norm0[i][0]>norm0[i-1][0]+eps) {
            concave2 = i;
            break;
        }
    }
    double smax = EvaluateEdgePts(1, 1.)[0], sfirst = EvaluateEdgePts(1, -1. + 1*dy)[0];
    if(trimWakeSlope) {
        for(int i=0; i<=concave1; ++i) {
            double n0, n1 = -sqrt(1. - maxSlopx*maxSlopx);
            n0 = maxSlopx/(maxx - edge0[concave1][0])*(maxx - edge0[i][0]);
            norm0[i][0] = n0/sqrt(n0*n0+n1*n1);
            norm0[i][1] = n1/sqrt(n0*n0+n1*n1);
        }
        for(int i=concave2; i<=m_M; ++i) {
            double n0, n1 = sqrt(1. - maxSlopx*maxSlopx);
            n0 = maxSlopx/(maxx - edge0[concave1][0])*(maxx - edge0[i][0]);
            norm0[i][0] = n0/sqrt(n0*n0+n1*n1);
            norm0[i][1] = n1/sqrt(n0*n0+n1*n1);
        }
    }
    //calculate vertex
    m_vertex.push_back(EvaluateEdgePts(0, -1.));
    m_vertex.push_back(EvaluateEdgePts(0, 1.));
    std::vector<double> vertex2 = EvaluateEdgePts(0, 1.);
    std::vector<double> vertex3 = EvaluateEdgePts(0, -1.);
    for(int i=0; i<2; ++i) {
        vertex2[i] = vertex2[i] + smax*norm0[m_M][i];
        vertex3[i] = vertex3[i] + smax*norm0[0][i];
    }
    m_vertex.push_back(vertex2);
    m_vertex.push_back(vertex3);
    //
    for(int j=0; j<=m_N; ++j) {
        grows.push_back(EvaluateEdgePts(1, -1. + j*dy)[0]);
    }
    for(int j=0; j<=m_N; ++j) {
        for(int i=0; i<=m_M; ++i) {
            std::vector<double> p(2);
            double s = grows[j];
	        if(trimWakeSlope && edge0[i][0]>1.) {
	            double ratio = fabs(edge0[i][1])/sfirst;
	            s = stretch(ratio*2.+5., sfirst, s, ratio);
	        }
	        for(int k=0; k<2; ++k) p[k] = edge0[i][k] + s*norm0[i][k];
            m_pts.push_back(p);
            if(i==0 || j==0 || i==m_M || j==m_N) m_bndPts.insert(m_pts.size()-1);
        }
    }
    return m_pts.size();
}

std::vector<double> RectRegion::getVertex(int i) {
    return m_vertex[i];
}

int RectRegion::MeshGen(int M, int N, MeshType method, bool trimWakeSlope)
{
    m_M = M;
    m_N = N;
    if(method == eBoundaryLayer0) {
        ptsByBoundaryLayer(trimWakeSlope);
    }else if(method == eIsoparametric){
        ptsByIsoParametric();
    }else {
        ptsByStraightLine();
    }

    //edge generation
    for(int j=0; j<N; ++j) {
        for(int i=0; i<M; ++i) {
            std::vector<int> p;
            p.push_back(i   + j*(M+1));
            p.push_back(i+1 + j*(M+1));
            m_edges.push_back(p);
        }
        for(int i=0; i<=M; ++i) {
            std::vector<int> p;
            p.push_back(i +     j*(M+1));
            p.push_back(i + (j+1)*(M+1));
            m_edges.push_back(p);
        }
    }
    for(int i=0; i<M; ++i) {
        std::vector<int> p;
        p.push_back(i   + N*(M+1));
        p.push_back(i+1 + N*(M+1));
        m_edges.push_back(p);
    }
    rebuildEdgesIndex();
    //cell generation
    for(int j=0; j<N; ++j)
    for(int i=0; i<M; ++i)
    {
        std::vector<int> p;
        p.push_back(i      +     j*(2*M+1));
        p.push_back(i+1+ M +     j*(2*M+1));
        p.push_back(i      + (j+1)*(2*M+1));
        p.push_back(i  + M +     j*(2*M+1));
        m_cells.push_back(p);
    }
    return 0;
}
void RectRegion::Tec360Pts(std::string fileName)
{
    std::ofstream outtec(fileName.c_str());
    outtec << "title = mesh\n";
    outtec << "variables = x, y\n";
    outtec << "zone i = " << m_M+1 << " j = " << m_N+1 << "\n";
    for(unsigned int j=0; j<m_pts.size(); j++)
    {
	    outtec << m_pts[j][0] << " " << m_pts[j][1] << std::endl;
    }
}
int RectRegion::EdgeAttractMesh(double x0, double y0, double &x, double &y)
{
    double a = 0.5*(1.-m_attractEps);
    x = x0;
    y = y0;
    if(m_edgeAttractMesh==1 || m_edgeAttractMesh==2) a = -a;
    if(m_edgeAttractMesh==0 || m_edgeAttractMesh==2) {
        x = x0;
        y = x0*x0*y0 + (1.-x0*x0)*(a*(y0*y0-1.)+y0);
    }else {
        x = y0*y0*x0 + (1.-y0*y0)*(a*(x0*x0-1.)+x0);
        y = y0;
    }
    return 0;
}
