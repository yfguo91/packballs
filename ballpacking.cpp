#include"ballpacking.h"
#include"geomfuncs.hpp"
#include"adtree.hpp"
#include "math.h"
#include<iostream>
using namespace std;

double global_size;
double(*coord_nodes3)[3] = 0, (*coordbak)[3] = 0;
int numnodes3;
int coordMode = 0;
double midCoordX, midCoordY;

const double PI = 4 * atan(1.);

int divide = 1;
double ratio_g = 0.175, finest_g = 0.01, g_maxincresize = (1 + ratio_g) / (1 - ratio_g); //0.1

int nums3n = 0, nums4n = 0, numnodes3b = 0, numnodes3s = 0, numvertnodes = 0, numedgenodes = 0, maxnodesface = 0, maxnodes3 = 0,
status_feat = 0, status_div = 0;
double *nodes3_r = 0, Tmin_size = 0, min_size = 0;

extern int *startaddress;
Kodtree *ktreeofnodes;
void   *mallocate(unsigned int size){

	void  *pointer;

	if (size <= 0){
		cout << " allocaterr ";
		return NULL;
	}
	pointer = malloc(size);
	if (pointer == NULL){
		cout << ("Memory allocation error --- heap overflow. \n");
		return(NULL);
	}
	else return(pointer);
}

void *reallocat(void *block, unsigned int size){

	if (size <= 0) cout << ("allocate errt", 1);
	block = realloc(block, size);
	if (block == NULL){
		cout << ("Memory allocation error --- heap overflow. \n", 1);
		return(NULL);
	}
	else return(block);
}
BallPacking::BallPacking(double(*vti)[3], int numvi, int(*tris)[3], int numti)
{
	Poi<3> point;
	this->numvert = numvi;
	this->vert = (double(*)[3]) new double[3 * numvert];
	coord_nodes3 = (double(*)[3])mallocate(numvert*sizeof(double[3]));
	for (int i = 0; i < numvi; i++)
	{
		point[0] = coord_nodes3[i][0] = vert[i][0] = vti[i][0];
		point[1] = coord_nodes3[i][1] = vert[i][1] = vti[i][1];
		point[2] = coord_nodes3[i][2] = vert[i][2] = vti[i][2];
		allpoints.push_back(point);
	}
	this->numtri = numti;
	this->tris = (int(*)[3]) new int[3 * numtri];
	for (int i = 0; i < numti; i++)
	{
		this->tris[i][0] = tris[i][0];
		this->tris[i][1] = tris[i][1];
		this->tris[i][2] = tris[i][2];
	}
	ptpoly = new PointInPolyhedron(this->vert, this->numvert, this->tris, this->numtri);
	hgmin = 0;
	grading = 0.175;
	numnodes3s=numnodes3 = numvert;
}
BallPacking::~BallPacking()
{
	delete[] vert; // Release memory!!
	delete[] tris; // Release memory!!
	delete ptpoly; // Release memory!!
}
void BallPacking::SetGlobalH(double h)
{
	hglob = h/2.;
	global_size = hglob;
}
void BallPacking::SetMinH(double h)
{
	hgmin = h/2.;
	min_size = hgmin;
}
void  BallPacking::SetLocalH(const Poi<3> & pmin, const Poi<3> & pmax, double grading)
{
	Poi<3> c = Center(pmin, pmax);
	double d = max(pmax[0] - pmin[0],
		max(pmax[1] - pmin[1],
		pmax[2] - pmin[2]));
	d /= 2;
	Poi<3> pmin2 = c - Vec<3>(d, d, d);
	Poi<3> pmax2 = c + Vec<3>(d, d, d);

	if (lochfunc != NULL)
		delete lochfunc;
	lochfunc = new LocalH(pmin2, pmax2, grading);
}
void  BallPacking::RestrictLocalH(const Poi<3> & p, double hloc)
{
	if (hloc < hgmin)
		hloc = hgmin;

	//cout << "restrict h in " << p << " to " << hloc << endl;
	if (!lochfunc)
	{
		cout << ("RestrictLocalH called, creating mesh-size tree") << endl;
	}
	//在p点设置局部尺寸
	lochfunc->SetH(p, hloc);
}
double  BallPacking::GetH(const Poi<3> & p) const
{
	double hmi = hglob;
	if (lochfunc)
	{
		double hl = lochfunc->GetH(p);
		if (hl < hglob)
			hmi = hl;
	}
	return hmi;
}
void *callocate(unsigned nump, unsigned int  size){

	void *pointer;

	if (size <= 0) cout<<("allocate errc", 1);
	pointer = calloc(nump, size);
	if (pointer == NULL){
		cout<<("Memory allocation error --- heap overflow. \n", 1);
		return(NULL);
	}
	else return(pointer);
}
void BallPacking::creatLocalHFromSurface()
{
	Bo<3> boundingbox;
	Poi<3> point;
	point[0] = vert[0][0];
	point[1] = vert[0][1];
	point[2] = vert[0][2];
	boundingbox.Set(point);
	for (int i = 1; i < numvert; i++)
	{
		point[0] = vert[i][0];
		point[1] = vert[i][1];
		point[2] = vert[i][2];
		boundingbox.Add(point);
	}
	boundingbox.Increase(0.01*boundingbox.Diam());
	//开始创建尺寸场
	SetLocalH(boundingbox.PMin(), boundingbox.PMax(), grading);
	Poi<3> point1, point2, point3;
	for (int i = 0; i < numtri; i++)
	{
		point1[0] = vert[tris[i][0]][0];
		point1[1] = vert[tris[i][0]][1];
		point1[2] = vert[tris[i][0]][2];
		point2[0] = vert[tris[i][1]][0];
		point2[1] = vert[tris[i][1]][1];
		point2[2] = vert[tris[i][1]][2];
		point3[0] = vert[tris[i][2]][0];
		point3[1] = vert[tris[i][2]][1];
		point3[2] = vert[tris[i][2]][2];
		Poi<3> c = Center(point1, point2, point3);
		double d = (1.0 / 6.0)*((point1 - point2).Length() + (point3 - point2).Length() + (point1 - point3).Length());
		RestrictLocalH(c, d);
	}
	nodes3_r = (double *)callocate(numvert, sizeof(double));
	//创建表面节点的尺寸
	for (int i = 0; i < allpoints.size(); i++)
	{
		radis.push_back(GetH(allpoints[i]));
		nodes3_r[i] = radis[i];
	}

}

void BallPacking::getNearNodes(int p, Poi3dTree * pointtree, vector<int>& pis)
{
	double a = radis[p];
	double  dist = a + 1.02*a*(1 + grading);
	//double  dist = a + 2.01*a*(1 + grading) + a*(1 + grading)*(1 + grading);
	Poi<3> pmin = allpoints[p] - Vec<3>(dist, dist, dist);
	Poi<3> pmax = allpoints[p] + Vec<3>(dist, dist, dist);

	pointtree->GetIntersecting(pmin, pmax, pis);
}

#define sq3(x,y,z)    ((x)*(x)+(y)*(y)+(z)*(z))
#define sq2(x,y)       ((x)*(x)+(y)*(y))

int select_max_team(double(*rank)[4], int k, int n)
{
	double d, agd;   int l, i, j;

	for (l = k, d = rank[k][k], i = k + 1; i<n; i++)
		if (fabs(rank[i][k])>fabs(d)) { d = rank[i][k]; l = i; }

	if (fabs(d) <= 0.0001) return(0);      /* eps_1: bound of 0 team  :min_siz*3.14/60 */
	if (l == k) return(1);
	for (j = k; j<n + 1; j++)
	{
		agd = rank[l][j]; rank[l][j] = rank[k][j]; rank[k][j] = agd;
	}
	return(2);
}

int solve_linear_equation(double(*rank)[4], int n, double *x)
{
	int i, j, k;   double s;

	for (k = 0; k<n - 1; k++)
	{
		if (select_max_team(rank, k, n) == 0) return(0);
		for (j = k + 1; j<n + 1; j++)  rank[k][j] /= rank[k][k];
		for (j = k + 1; j<n + 1; j++)
			for (i = k + 1; i<n; i++)
				rank[i][j] -= rank[i][k] * rank[k][j];
	}
	if (fabs(rank[n - 1][n - 1]) <= 0.0001)  return(0);
	x[n - 1] = rank[n - 1][n] / rank[n - 1][n - 1];
	for (i = n - 2; i>-1; i--)
	{
		s = 0;
		for (j = i + 1; j<n; j++) s += rank[i][j] * x[j];
		x[i] = rank[i][n] - s;
	}
	return(1);
}

int BallPacking::set_node_3d(int ni, int nj, int nk, jf_point *p)
{
	jf_point norm;
	double rp, ris, rjs, rks, high, rank[3][4], pxc[3];

	rp = (nodes3_r[ni] + nodes3_r[nj] + nodes3_r[nk]) / 3.;/*2/3*/
	ris = (nodes3_r[ni] + rp)*(nodes3_r[ni] + rp);
	rjs = (nodes3_r[nj] + rp)*(nodes3_r[nj] + rp);
	rks = (nodes3_r[nk] + rp)*(nodes3_r[nk] + rp);
	rank[0][0] = coord_nodes3[ni][0] - coord_nodes3[nk][0];
	rank[0][1] = coord_nodes3[ni][1] - coord_nodes3[nk][1];
	rank[0][2] = coord_nodes3[ni][2] - coord_nodes3[nk][2];
	rank[1][0] = coord_nodes3[nj][0] - coord_nodes3[nk][0];
	rank[1][1] = coord_nodes3[nj][1] - coord_nodes3[nk][1];
	rank[1][2] = coord_nodes3[nj][2] - coord_nodes3[nk][2];
	vec_crop(rank[0], rank[1], rank[2]);
	rank[0][3] = sq3(rank[0][0], rank[0][1], rank[0][2]) / 2.;
	rank[1][3] = sq3(rank[1][0], rank[1][1], rank[1][2]) / 2.;
	rank[2][3] = 0;
	if (fabs(vec_dotp(rank[0], rank[1]))>1.98*sqrt(rank[0][3])*sqrt(rank[1][3]))
		return 0;                         /*2*cos6*/
	if (vec_uni(rank[2]) == 0) return 0;
	rank[0][3] += (rks - ris) / 2.; rank[1][3] += (rks - rjs) / 2.;
	norm[0] = rank[2][0];  norm[1] = rank[2][1];   norm[2] = rank[2][2];
	if (solve_linear_equation(rank, 3, pxc) == 0){
		/*  printf(" no node set " ); */     return(0);
	}
	high = rks - sq3(pxc[0], pxc[1], pxc[2]);
	if (high<-0.05*rks) return(0);
	if (high >= -0.05*rks&&high <= 0.05*rks)
		high = 0.;
	else
		high = sqrt(high);
	norm[0] *= high;   norm[1] *= high;   norm[2] *= high;
	pxc[0] += coord_nodes3[nk][0];
	pxc[1] += coord_nodes3[nk][1];
	pxc[2] += coord_nodes3[nk][2];
	p[0][0] = pxc[0] + norm[0];
	p[0][1] = pxc[1] + norm[1];
	p[0][2] = pxc[2] + norm[2];
	if (fabs(high)<0.0001) return 1;
	p[1][0] = pxc[0] - norm[0];
	p[1][1] = pxc[1] - norm[1];
	p[1][2] = pxc[2] - norm[2];
	return 2;
}

int BallPacking::is_near_3d(jf_point p, double r, vector<int>& pis, int nj, int nk)
{
	int i;    double wp;    Poi<3> point;
	point[0] = p[0];
	point[1] = p[1];
	point[2] = p[2];
	for (i = 0; i<pis.size(); i++){
		if (pis[i] == nj || pis[i] == nk) continue;
		wp = (point-allpoints[pis[i]]).Length2();
		if (1.20*wp < (radis[pis[i]] + r)*(radis[pis[i]] + r)) return(1);//2.35*2.35/4.
	}
	return(0);
}

vector<Poi<3>>& BallPacking::generateInnerPoints()
{
	//首先要创建一个节点树
	Poi3dTree * pointtree;

	Bo<3> boundingbox;
	Poi<3> point;

	boundingbox.Set(allpoints[0]);
	for (int i = 1; i < allpoints.size(); i++)
	{
		boundingbox.Add(allpoints[i]);
	}
	boundingbox.Increase(0.01*boundingbox.Diam());

	pointtree = new Poi3dTree(boundingbox.PMin(), boundingbox.PMax());

	for (int i = 0; i < allpoints.size(); i++)
	{
		pointtree->Insert(allpoints[i], i);
	}

	//其次遍历所有的节点，根据这个节点来寻找周围的节点，然后选取三个进行布置节点
	int itx;
	double p[2][3];
	for (int i = 0; i < allpoints.size(); i++)
	{
		int num = 0;
		cout << i <<" ";
		if (i == 494)
		{
			cout << "hello" << endl;
		}
		vector<int> pis;
		getNearNodes(i, pointtree, pis);
		cout << pis.size() << " ";
		for (int j = 0; j < pis.size() - 1; j++)
		{
			if (pis[j] <= i) continue;
			for (int k = j+1; k < pis.size(); k++)
			{
				if (pis[k] <= i) continue;
				//开始布置节点
				if ((itx = set_node_3d(i, pis[j], pis[k], p)) == 0) { continue; }
				for (int itw = 0; itw<itx; itw++){
					//					extern double GetDefinedEdgeSizeForPt(double p[3]);
					//					double rdefine=GetDefinedEdgeSizeForPt(p[itw]);
					if (ptpoly->isPinPolyhedron(p[itw])<0)
						continue;

					Poi<3> newp;
					newp[0] = p[itw][0];
					newp[1] = p[itw][1];
					newp[2] = p[itw][2];

					double r = GetH(newp);
					if (r>hglob) r = hglob;
					if (is_near_3d(p[itw], r, pis, pis[j], pis[k])) { continue; }

					allpoints.push_back(newp);
					radis.push_back(r);

					pointtree->Insert(newp,allpoints.size()-1);
					num++;
				}
			}
		}
		cout << num << endl;
	}
	return allpoints;
}

void BallPacking::outInnerPoints(char *strs)
{
	ofstream file(strs);
	for (int i = 0; i < allpoints.size(); i++)
	{
		file << allpoints[i][0] << " " << allpoints[i][1] << " " << allpoints[i][2] << " " << radis[i] << endl;
	}
}

void BallPacking::outInnerPoints2(char *strs)
{
	ofstream file(strs);
	cout << numnodes3 << endl;
	for (int i = 0; i <numnodes3; i++)
	{
		file << coord_nodes3[i][0] << " " << coord_nodes3[i][1] << " " << coord_nodes3[i][2] << " " << nodes3_r[i] << endl;
	}
}

void BallPacking::outInnernodes(char *strs)
{
	cout << numnodes3 << endl;
	//输出.node文件
	string filename1(strs);
	ofstream file1(filename1.c_str());
	double a, b, c;
	file1 << numnodes3-numnodes3s << "  " << 3 << "  " << 0 << "  " << 0 << endl;;
	for (int i = 1; i <= numnodes3 - numnodes3s; i++)
	{
		file1 << i << "  " << coord_nodes3[numnodes3s + i - 1][0] << "  " << coord_nodes3[numnodes3s + i - 1][1] << "  " << coord_nodes3[numnodes3s + i - 1][2] << endl;
	}
}


void wrapPointsUpasVert(void  ** &vti){

	vti = new void *[numnodes3];
	for (int i = 0; i<numnodes3; i++)
		vti[i] = startaddress + i;
}
int vertIndexFromPt(int *pta){
	return pta - startaddress;
}

void pofvforcoordnodes3(double p[3], void *pv){

	int n = vertIndexFromPt((int *)pv);
	p[0] = coord_nodes3[n][0];
	p[1] = coord_nodes3[n][1];
	p[2] = coord_nodes3[n][2];
}
int Record1NodeInRatioClass(std::list<int> **quevtclass, int numqvtclass, double r, int nd){

	int index = (int)(log(r / min_size) / log(ratio_g + 1.));
	if (index >= numqvtclass)
		cout<<("err recordvertex", 1);
	quevtclass[index]->push_back(nd);
	return index;
}
void RecordVertexInRatioClass(std::list<int> **&quevtclass, int &numqvtclass, double *rofnd, int nm){

	numqvtclass = (int)(log(global_size / min_size) / log(ratio_g + 1.) + 1);
	quevtclass = new std::list<int> *[numqvtclass];
	for (int i = 0; i<numqvtclass; i++)
		quevtclass[i] = new std::list<int>;
	for (int i = 0; i<nm; i++) Record1NodeInRatioClass(quevtclass, numqvtclass, rofnd[i], i);
}


void compCollisionScope(int node, Box &bd){

	double a = nodes3_r[node];
	double  dist0 = a + 2.01*a*g_maxincresize + a*g_maxincresize*g_maxincresize;

	if (dist0>4.*global_size) dist0 = 4.*global_size; //revise from 4.=>2.01, 06-09-23
	bd[0] = coord_nodes3[node][0] - dist0;
	bd[1] = coord_nodes3[node][1] - dist0;
	bd[2] = coord_nodes3[node][2] - dist0;
	bd[3] = coord_nodes3[node][0] + dist0;
	bd[4] = coord_nodes3[node][1] + dist0;
	bd[5] = coord_nodes3[node][2] + dist0;
}
double getMinDistOfNodeSet(int *nearnodes, int num, int &na, int &nb){

	double dist, dist0 = global_size;
	int ni, nj;
	for (int i = 0; i<num - 1; i++){
		ni = nearnodes[i];
		if (ni<numnodes3s) continue;
		for (int j = i + 1; j<num; j++){
			nj = nearnodes[j];
			if (nj<numnodes3s) continue;
			if ((dist = Distance3D(coord_nodes3[ni], coord_nodes3[nj]))<dist0){
				dist0 = dist;
				na = ni; nb = nj;
			}
		}
	}
	return dist0;
}
int sortnearnodes(const void *ai, const void *bi){
	int a, b;
	a = *((int *)ai); b = *(int *)bi;
	if (a<b) return -1;
	if (a>b) return 1;
	cout<<("errror sortnearnode", 1);
	return 0;
}
int maxctnear = 0;
void get_near_nodes3d(int node, int nearnodes[3200], int *pnum)
{
	int i;
	// double dist0 , dlx , drx , dly , dry , dlz , drz ;
	extern int ct;
	list<void *> lvert;

	*pnum = 0;  /* changing length to constant or search from point class ?*/
	double bd[6];
	compCollisionScope(node, bd);
	double a = nodes3_r[node];
	double  dist0 = a + 2.01*a*g_maxincresize + a*g_maxincresize*g_maxincresize;
	if (dist0>4.*global_size) dist0 = 4.*global_size;
	ktreeofnodes->collectVertsWithBox(bd, lvert);
	for (list<void *>::iterator pos = lvert.begin(); pos != lvert.end(); pos++){
		i = vertIndexFromPt((int *)*pos);

		/*	for( i=0; i<numnodes3; i++){
		if( coord_nodes3[i][0]<bd[0]||coord_nodes3[i][0]>bd[3]||
		coord_nodes3[i][1]<bd[1]||coord_nodes3[i][1]>bd[4]||
		coord_nodes3[i][2]<bd[2]||coord_nodes3[i][2]>bd[5] )continue ;
		*/
		double ab2sqrtab = nodes3_r[i] + nodes3_r[node] + 2.*sqrt(nodes3_r[i] * nodes3_r[node]);
		//	if( ( SqDistance3D( coord_nodes3[i],coord_nodes3[node] )<=dist0*dist0)&&i!=node){
		if ((0.9025*SqDistance3D(coord_nodes3[i], coord_nodes3[node])<ab2sqrtab*ab2sqrtab) && i != node){
			//	   nodes3_r[i]+ nodes3_r[node]+ 2*sqrt(nodes3_r[i]*nodes3_r[node]) )&&i!= node ){
			nearnodes[(*pnum)++] = i;
			if (*pnum >= 3200){                  /* test every time or here only 1 */
				cout<<(" error-near_node\n ", 1);
				int na, nb;
				double dmin = getMinDistOfNodeSet(nearnodes, *pnum, na, nb);
				qsort(nearnodes, *pnum, sizeof(int), sortnearnodes);
				dmin = dmin;
			}
		}
	}
	if (*pnum>maxctnear) maxctnear = *pnum;
}

int BallPacking::is_near_3d(jf_point p, double r, int nearnodes[3200], int numnear, int nj, int nk)
{
	int i;    double wp;
	for (i = 0; i<numnear; i++){
		if (nearnodes[i] == nj || nearnodes[i] == nk) continue;
		wp = SqDistance3D(p, coord_nodes3[nearnodes[i]]);
		if (1.11*wp < (nodes3_r[nearnodes[i]] + r)*(nodes3_r[nearnodes[i]] + r)) return(1);//2.35*2.35/4.
	}
	return(0);
}


bool *nodeProcessed;
void BallPacking::solid_fill()
{
	void **wvti;
	wrapPointsUpasVert(wvti);
	ktreeofnodes = new Kodtree(wvti, numnodes3, pofvforcoordnodes3, 2, 0.00001);
	delete wvti;
	numnodes3s = numnodes3;
	maxnodes3 = numnodes3s;
	nodeProcessed = (bool *)mallocate(sizeof(bool)*maxnodes3);
	for (int i = 0; i<maxnodes3; i++) nodeProcessed[i] = false;
	int nj, nk, itw, itx, numnear, nearnodes[3200];
	double p[2][3];
	std::list<int> **quevtclass;
	int numqvtclass;
	RecordVertexInRatioClass(quevtclass, numqvtclass, nodes3_r, numnodes3);
	for (int ilist = 0; ilist<numqvtclass; ilist++){
		std::list<int> *tlist = quevtclass[ilist];
		int backclass = ilist;
		while (!(*tlist).empty()){
			int i = (*tlist).front();
			(*tlist).pop_front();
			nodeProcessed[i] = true;
			get_near_nodes3d(i, nearnodes, &numnear);
			for (int j = 0; j<numnear - 1; j++){
				if (nodeProcessed[(nj = nearnodes[j])] == true) continue;
				for (int k = j + 1; k<numnear; k++){
					if (nodeProcessed[(nk = nearnodes[k])] == true) continue;
					if ((itx = set_node_3d(i, nj, nk, p)) == 0) { continue; }
					//				double r=nodes3_r[i]*(1.+ratio_g)<global_size?nodes3_r[i]*(1.+ratio_g):global_size;
					for (itw = 0; itw<itx; itw++){
						//					extern double GetDefinedEdgeSizeForPt(double p[3]);
						//					double rdefine=GetDefinedEdgeSizeForPt(p[itw]);
						double r = nodes3_r[i] + ratio_g*Distance3D(coord_nodes3[i], p[itw]);//(1-ratio_g);
						if (r>global_size) r = global_size;
						if (is_near_3d(p[itw], r, nearnodes, numnear, nj, nk)) { continue; }
						if (ptpoly->isPinPolyhedron(p[itw])<0)
							continue;
						if (numnodes3 == maxnodes3){
							maxnodes3 = 2 * maxnodes3;
							coord_nodes3 = (jf_point *)reallocat(coord_nodes3, maxnodes3 * sizeof(jf_point));
							nodes3_r = (double *)reallocat(nodes3_r, maxnodes3 * sizeof(double));
							nodeProcessed = (bool *)reallocat(nodeProcessed, maxnodes3 * sizeof(bool));
							for (int mi = maxnodes3 / 2; mi<maxnodes3; mi++)	nodeProcessed[mi] = false;
						}
						nearnodes[numnear++] = numnodes3;
						nodes3_r[numnodes3] = r;
						memcpy(coord_nodes3[numnodes3++], p[itw], 3 * sizeof(double));
						ktreeofnodes->insertVert(startaddress + numnodes3 - 1);
						if (numnodes3 == 6606)
							numnodes3 = 6606;
						int idclass = Record1NodeInRatioClass(quevtclass, numqvtclass, r, numnodes3 - 1);
						if (idclass<backclass) backclass = idclass;
						if (numnear >= 3200)
							std::cout<<("error nearnodes_3d\n ", 1);
					}
				}
			}
		}
		if (backclass<ilist)ilist = backclass - 1;
		//		delete tlist;
	}
	for (int i = 0; i<numqvtclass; i++)
		delete quevtclass[i];
	delete[] quevtclass;
	delete ktreeofnodes;
	delete[] nodeProcessed;
	cout << "maxctnear:" << maxctnear << " ";
}