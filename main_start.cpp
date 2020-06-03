#include "pinpolyhedron.h"
#include "stopwatch.h"
#include <float.h>
#include "ballpacking.h"
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

void main(){

	MyMesh mesh1;
	bool result1 = OpenMesh::IO::read_mesh(mesh1, "test.stl");
	//将网格的节点和三角网格给取出来

	cout << mesh1.n_vertices() << endl;

	int numvert = mesh1.n_vertices();
	double(*vert)[3];
	vert = (double(*)[3]) new double[3 * numvert];

	int numtri = mesh1.n_faces();
	int(*tris)[3];
	tris = (int(*)[3]) new int[3 * numtri];

	int i = 0;
	for (auto v_it = mesh1.vertices_begin(); v_it != mesh1.vertices_end(); v_it++, i++)
	{
		vert[i][0] = mesh1.point(*v_it)[0];
		vert[i][1] = mesh1.point(*v_it)[1];
		vert[i][2] = mesh1.point(*v_it)[2];
	}

	i = 0;
	for (auto f_it = mesh1.faces_begin(); f_it != mesh1.faces_end(); f_it++,i++)
	{
		int j = 0;
		for (auto fv_it = mesh1.fv_begin(*f_it); fv_it != mesh1.fv_end(*f_it); fv_it++,j++)
		{
			tris[i][j] = (*fv_it).idx();
		}
	}

	//创建ballpacking,生成球填充的节点
	BallPacking bp(vert, numvert, tris, numtri);
	bp.SetGlobalH(50.);
	bp.SetMinH(1.);
	bp.creatLocalHFromSurface();
	//bp.generateInnerPoints();
	bp.solid_fill();
	//bp.outInnerPoints2("packing.txt");
	bp.outInnernodes("add.node");
	delete[] vert; // Release memory!!
	delete[] tris; // Release memory!!

	cout << "complete!!!" << endl;
}