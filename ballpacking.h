
#include"pinpolyhedron.h"
#include"localh.hpp"
#include<vector>
#include"adtree.hpp"
using namespace std;
typedef double jf_point[3];
class BallPacking
{
private:
	LocalH * lochfunc;
	///ȫ�����ߴ�
	double hglob;
	///ȫ����С�ߴ�
	double hgmin;

	double grading;

	PointInPolyhedron *ptpoly;

	//����Ľڵ�����,����Ƭ
	double(*vert)[3];
	int numvert, (*tris)[3], numtri;

	//���еĽڵ�����
	vector<Poi<3>> allpoints;
	vector<double> radis;
	//����Ľڵ�����
	vector<Poi<3>> genpoints;



public:
	BallPacking(double(*vti)[3], int numvi, int(*tris)[3], int numti);
	~BallPacking();
	void  SetGlobalH(double h);
	void  SetMinH(double h);
	double GetGloabalH(){ return hglob; }
	void  SetLocalH(const Poi<3> & pmin, const Poi<3> & pmax, double grading);
	void  RestrictLocalH(const Poi<3> & p, double hloc);
	double  GetH(const Poi<3> & p) const;

	void creatLocalHFromSurface();

	vector<Poi<3>>& generateInnerPoints();

	void outInnerPoints(char *strs);

	void solid_fill();
	void outInnerPoints2(char *strs);
	void outInnernodes(char *strs);

private:
	void getNearNodes(int p, Poi3dTree * pointtree, vector<int>& pis);
	int set_node_3d(int ni, int nj, int nk, jf_point *p);
	int is_near_3d(jf_point p, double r, vector<int>& pis, int nj, int nk);
	int is_near_3d(jf_point p, double r, int nearnodes[3200], int numnear, int nj, int nk);
};