# include <iostream>
# include "FEMSolid.h"
# include "tet.h"

using namespace Eigen;
/*struct tetCorners
{
	std::array<Vector3d, 4> positions;
};*/

struct tetTetrahedras
{
	std::vector<Vector4i> indexes;
};

struct tetCorners
{
	std::vector<Vector3d> positions;
	std::vector<Vector3d> velocities;
	std::vector<Vector3d> elasticForces;
};

int main(int argc, char** argv)
{
	std::cout << "Just seeing if this works" << std::endl;
	Eigen::Matrix4f mat = Eigen::Matrix4f();
	std::cout << mat << std::endl;
	tetgenio ioOp = tetgenio();
	ioOp.initialize();
	//=========================================== test tetGen
	
	tetGen t1;

	tetTetrahedras tetlist;
	tetCorners cornerlist;
	int t;
	std::cout << t1.out.numberoftetrahedra << std::endl;
	std::cout << t1.out.numberofcorners << std::endl;

	for (int i = 0; i < t1.out.numberoftetrahedra; i++)
	{
		tetlist.indexes.push_back(Vector4i(t1.out.tetrahedronlist[4 * i],
										   t1.out.tetrahedronlist[4 * i + 1],
										   t1.out.tetrahedronlist[4 * i + 2],
										   t1.out.tetrahedronlist[4 * i + 3]));
	}

	for (int i = 0; i < t1.out.numberofpoints; i++)
	{
		Vector3d zero = Vector3d();
		cornerlist.positions.push_back(Vector3d(t1.out.pointlist[3 * i],
												t1.out.pointlist[3 * i + 1],
												t1.out.pointlist[3 * i + 2]));
		cornerlist.velocities.push_back(zero);
		cornerlist.elasticForces.push_back(zero);
	}

	//=======================================test output===========

	std::fstream file;
	file.open("test.poly", std::ios::out);
	objwriter::ObjWriter w(file);
	w.point();
	for (int i = 0; i < cornerlist.positions.size(); i++)
	{
		w.vertex(float(cornerlist.positions[i][0]), float(cornerlist.positions[i][1]), float(cornerlist.positions[i][2]), i);
	}
	w.polys();
	for (int i = 0; i < tetlist.indexes.size();i++)
	{
		w.tet(tetlist.indexes[i][0], tetlist.indexes[i][1], tetlist.indexes[i][2], tetlist.indexes[i][3], i);
	}
	w.end();
	//=============================================================

	//=================================================================
	//Main Algorithm
	std::vector<Matrix3d> Dm;
	std::vector<Matrix3d> Bm;
	std::vector<double> We;

	//For all tet, precompute the volume
	//for (int i = 0; i < tetlist.indexes.size(); i++)
	//{
	//	Matrix3d dmt;
	//	Vector3d& lastCorner = cornerlist.positions[tetlist.indexes[i][3]];

	//	for (int cornerInd = 0; cornerInd < 3; cornerInd++)
	//	{
	//		Vector3d& corner = cornerlist.positions[tetlist.indexes[i][cornerInd]];

	//		dmt.row(cornerInd) << corner - lastCorner;
	//	}
	//	Dm.push_back(dmt);
	//	Bm.push_back(dmt.inverse());
	//	We.push_back(dmt.determinant());
	//}
	//// Compute the elastic forces
	//// current positions
	//// todo: need to populate it
	//std::vector<tetCorners> currentTetList;
	//
	//// force on each corner
	//std::vector<tetCorners> tetForces;
	////stress tensor
	//std::vector<Matrix3d> F;
	////std::vector<Vector3d> H;
	////std::vector<Matrix3d> P;
	//for (int tetInd = 0; tetInd < tetlist.indexes.size(); tetInd++)
	//{
	//	Matrix3d dst;
	//	Vector3d& lastCorner = cornerlist.positions[tetlist.indexes[tetInd][3]];
	//	for (int cornerInd = 0; cornerInd < 3; cornerInd++)
	//	{
	//		Vector3d& corner = cornerlist.positions[tetlist.indexes[tetInd][cornerInd]];

	//		dst.row(cornerInd) << corner - lastCorner;
	//	}
	//	dst *= Bm[tetInd];
	//	F.push_back(dst);

	//	// todo: compute P
	//	//JacobiSVD<Matrix3d> svd(dst, ComputeFullU | ComputeFullV);
	//	Matrix3d p;

	//	// populate H
	//	Matrix3d H = - We[tetInd] * p * Bm[tetInd].transpose();

	//	Vector3d f3 = Vector3d();
	//	for (int cornerInd = 0; cornerInd < 3; cornerInd++)
	//	{
	//		cornerlist.elasticForces[tetlist.indexes[tetInd][cornerInd]] = H.col(cornerInd);
	//		f3 -= H.col(cornerInd);
	//	}
	//	
	//	cornerlist.elasticForces[tetlist.indexes[tetInd][3]] = f3;
	//}
}